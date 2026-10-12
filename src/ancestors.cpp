#include "common.h"
#include "genotypes.h"
#include "grouping.h"
#include "lowlevel.h"
#include "threads.h"
#include "zarr_r.h"

#include <algorithm>
#include <numeric>

using namespace Rcpp;

namespace {

struct AncestorResult {
    int64_t start_local = 0;
    int64_t end_local = 0;
};

/* Accumulates ancestor haplotypes for one sample chunk of the output
 * call_genotype array (sites x samples_chunk_size, row-major) and writes the
 * chunk to disk once it is complete. */
class ChunkWriter {
  public:
    ChunkWriter(const tsir::ZarrArray &array, size_t num_sites, int num_threads)
        : array_(array), num_sites_(num_sites), width_(static_cast<size_t>(array.chunks[1])),
          num_threads_(num_threads), buffer_(num_sites * width_, -1)
    {
    }

    void set(size_t site, int64_t ancestor, int8_t value)
    {
        buffer_[site * width_ + static_cast<size_t>(ancestor % static_cast<int64_t>(width_))]
            = value;
    }

    size_t width() const { return width_; }

    void flush(int64_t chunk)
    {
        const int64_t vchunk = array_.chunks[0];
        const int64_t num_vchunks
            = num_sites_ == 0 ? 0 : (static_cast<int64_t>(num_sites_) + vchunk - 1) / vchunk;
        tsir::parallel_for(static_cast<size_t>(num_vchunks), num_threads_, [&](size_t v) {
            std::vector<int8_t> data(static_cast<size_t>(vchunk) * width_, -1);
            for (int64_t r = 0; r < vchunk; r++) {
                size_t site = v * static_cast<size_t>(vchunk) + static_cast<size_t>(r);
                if (site >= num_sites_) {
                    break;
                }
                std::copy(buffer_.begin() + site * width_, buffer_.begin() + (site + 1) * width_,
                    data.begin() + r * static_cast<int64_t>(width_));
            }
            tsir::write_chunk(array_, { static_cast<int64_t>(v), chunk, 0 },
                reinterpret_cast<const uint8_t *>(data.data()));
        });
        std::fill(buffer_.begin(), buffer_.end(), static_cast<int8_t>(-1));
    }

  private:
    const tsir::ZarrArray &array_;
    size_t num_sites_;
    size_t width_;
    int num_threads_;
    std::vector<int8_t> buffer_;
};

} // namespace

// [[Rcpp::export(.build_ancestors)]]
List
build_ancestors(List view_spec, NumericVector times, IntegerMatrix intervals,
    NumericVector positions, List gt_meta, int num_threads, int genotype_encoding,
    bool verbose)
{
    tsir::GenotypeView view = tsir::view_from_list(view_spec);
    tsir::ZarrArray out = tsir::zarr_array_from_list(gt_meta);
    const size_t num_sites = view.num_positions;
    const size_t num_haps = view.num_haplotypes();
    if (static_cast<size_t>(times.size()) != num_sites
        || static_cast<size_t>(positions.size()) != num_sites) {
        stop("times and positions must have one entry per site");
    }
    const int threads = tsir::resolve_num_threads(num_threads);
    std::vector<double> anc_time;
    std::vector<int> anc_start, anc_end;
    std::vector<std::vector<int>> anc_focal;
    ChunkWriter writer(out, num_sites, threads);
    const int64_t width = static_cast<int64_t>(writer.width());
    int64_t written = 0;

    for (int iv = 0; iv < intervals.nrow(); iv++) {
        const size_t lo = static_cast<size_t>(intervals(iv, 0));
        const size_t hi = static_cast<size_t>(intervals(iv, 1));
        if (hi <= lo) {
            continue;
        }
        const size_t n_local = hi - lo;
        tsir::AncestorBuilderHandle builder(num_haps, n_local, genotype_encoding);
        ancestor_builder_t *ab = &builder.builder;
        std::vector<allele_t> derived(num_haps);
        int err = 0;
        try {
            tsir::iterate_genotypes(view, lo, hi, [&](size_t k, const int8_t *g) {
                if (err != 0) {
                    return;
                }
                for (size_t h = 0; h < num_haps; h++) {
                    derived[h] = g[h] > 0 ? 1 : g[h];
                }
                err = ancestor_builder_add_site(ab, times[k], derived.data());
            });
        } catch (std::exception &e) {
            stop(e.what());
        }
        tsir::check_tsi(err, "add_site");
        builder.finalise();
        const size_t n_anc = ab->num_ancestors;
        if (verbose) {
            Rcpp::message(wrap(tfm::format("Interval %d: %d sites (%d-%d), %d ancestors",
                iv, static_cast<int>(n_local), static_cast<int>(positions[lo]),
                static_cast<int>(positions[hi - 1]), static_cast<int>(n_anc))));
        }
        if (n_anc == 0) {
            continue;
        }
        const ancestor_descriptor_t *descr = ab->descriptors;

        /* Bounds pass: needed to group ancestors before writing them. */
        std::vector<AncestorResult> bounds(n_anc);
        try {
            tsir::parallel_for(n_anc, threads, [&](size_t j) {
                std::vector<allele_t> hap(n_local);
                tsk_id_t s, e;
                int ret = ancestor_builder_make_ancestor(
                    ab, descr[j].num_focal_sites, descr[j].focal_sites, &s, &e, hap.data());
                if (ret != 0) {
                    throw std::runtime_error(tsi_strerror(ret));
                }
                bounds[j].start_local = s;
                bounds[j].end_local = e;
            });
        } catch (std::exception &e) {
            stop("make_ancestor: %s", e.what());
        }
        Rcpp::checkUserInterrupt();

        std::vector<size_t> write_order(n_anc);
        std::iota(write_order.begin(), write_order.end(), 0);
        if (n_anc > 1) {
            std::vector<int64_t> starts(n_anc), ends(n_anc);
            std::vector<double> t(n_anc);
            for (size_t j = 0; j < n_anc; j++) {
                starts[j] = static_cast<int64_t>(positions[lo + bounds[j].start_local]);
                ends[j] = static_cast<int64_t>(positions[lo + bounds[j].end_local - 1]);
                t[j] = descr[j].time;
            }
            std::vector<int> group;
            try {
                group = tsir::group_by_linesweep(starts, ends, t);
            } catch (std::exception &e) {
                stop(e.what());
            }
            std::stable_sort(write_order.begin(), write_order.end(),
                [&](size_t a, size_t b) { return group[a] < group[b]; });
        }

        /* Write pass, one output sample chunk at a time. */
        size_t w = 0;
        while (w < n_anc) {
            const int64_t chunk = written / width;
            const size_t batch_end = std::min(
                n_anc, w + static_cast<size_t>((chunk + 1) * width - written));
            const size_t batch = batch_end - w;
            try {
                tsir::parallel_for(batch, threads, [&](size_t b) {
                    const size_t j = write_order[w + b];
                    std::vector<allele_t> hap(n_local);
                    tsk_id_t s, e;
                    int ret = ancestor_builder_make_ancestor(
                        ab, descr[j].num_focal_sites, descr[j].focal_sites, &s, &e, hap.data());
                    if (ret != 0) {
                        throw std::runtime_error(tsi_strerror(ret));
                    }
                    const int64_t ancestor = written + static_cast<int64_t>(b);
                    for (size_t l = 0; l < n_local; l++) {
                        writer.set(lo + l, ancestor, hap[l]);
                    }
                });
            } catch (std::exception &e) {
                stop("make_ancestor: %s", e.what());
            }
            for (size_t b = 0; b < batch; b++) {
                const size_t j = write_order[w + b];
                const AncestorResult &r = bounds[j];
                anc_time.push_back(descr[j].time);
                anc_start.push_back(static_cast<int>(positions[lo + r.start_local]));
                anc_end.push_back(static_cast<int>(positions[lo + r.end_local - 1]) + 1);
                std::vector<int> focal(descr[j].num_focal_sites);
                for (size_t f = 0; f < focal.size(); f++) {
                    focal[f] = static_cast<int>(positions[lo + descr[j].focal_sites[f]]);
                }
                anc_focal.push_back(std::move(focal));
            }
            written += static_cast<int64_t>(batch);
            w = batch_end;
            if (written % width == 0) {
                try {
                    writer.flush(chunk);
                } catch (std::exception &e) {
                    stop(e.what());
                }
                Rcpp::checkUserInterrupt();
            }
        }
    }
    if (written % width != 0) {
        try {
            writer.flush(written / width);
        } catch (std::exception &e) {
            stop(e.what());
        }
    }
    List focal(anc_focal.size());
    for (size_t j = 0; j < anc_focal.size(); j++) {
        focal[j] = IntegerVector(anc_focal[j].begin(), anc_focal[j].end());
    }
    return List::create(_["time"] = NumericVector(anc_time.begin(), anc_time.end()),
        _["start"] = IntegerVector(anc_start.begin(), anc_start.end()),
        _["end"] = IntegerVector(anc_end.begin(), anc_end.end()), _["focal_positions"] = focal);
}
