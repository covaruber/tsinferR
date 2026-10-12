#include "genotypes.h"
#include "zarr_r.h"

using namespace Rcpp;

namespace tsir {

size_t
GenotypeView::num_haplotypes() const
{
    size_t n = 0;
    for (const ViewSource &s : sources) {
        n += s.num_haplotypes();
    }
    return n;
}

GenotypeView
view_from_list(List spec)
{
    GenotypeView view;
    view.num_positions = static_cast<size_t>(as<double>(spec["num_positions"]));
    List sources = spec["sources"];
    for (R_xlen_t j = 0; j < sources.size(); j++) {
        List s = sources[j];
        ViewSource vs;
        vs.gt = get_zarr_array(s["gt"]);
        if (vs.gt->ndim() != 3) {
            stop("call_genotype must be 3-dimensional (variants, samples, ploidy)");
        }
        NumericVector rows = s["rows"];
        IntegerMatrix remap = s["remap"];
        NumericVector cols = s["cols"];
        vs.ploidy = as<int>(s["ploidy"]);
        if (static_cast<size_t>(rows.size()) != view.num_positions
            || static_cast<size_t>(remap.nrow()) != view.num_positions) {
            stop("genotype view: rows/remap must have one entry per position");
        }
        vs.rows.assign(rows.begin(), rows.end());
        vs.max_alleles = remap.ncol();
        vs.remap.resize(static_cast<size_t>(remap.nrow()) * vs.max_alleles);
        for (int k = 0; k < remap.nrow(); k++) {
            for (int a = 0; a < remap.ncol(); a++) {
                vs.remap[static_cast<size_t>(k) * vs.max_alleles + a]
                    = static_cast<int8_t>(remap(k, a));
            }
        }
        vs.cols.assign(cols.begin(), cols.end());
        for (size_t c = 1; c < vs.cols.size(); c++) {
            if (vs.cols[c] <= vs.cols[c - 1]) {
                stop("genotype view: sample columns must be strictly increasing");
            }
        }
        view.sources.push_back(std::move(vs));
    }
    return view;
}

namespace {

/* Raw genotypes of the selected samples for one variant chunk of a source. */
struct ChunkBuffer {
    int64_t chunk = -1;
    int64_t first_row = 0;
    std::vector<int32_t> data; /* rows x (cols * ploidy) */
};

void
load_variant_chunk(const ViewSource &s, int64_t vc, ChunkBuffer &buf)
{
    const ZarrArray &a = *s.gt;
    const int64_t vchunk = a.chunks[0];
    const int64_t schunk = a.chunks[1];
    const int64_t pchunk = a.chunks[2];
    const int64_t nrows = std::min(vchunk, a.shape[0] - vc * vchunk);
    const size_t width = s.num_haplotypes();
    buf.chunk = vc;
    buf.first_row = vc * vchunk;
    buf.data.assign(static_cast<size_t>(nrows) * width, -1);
    std::vector<int32_t> raw;
    size_t c = 0;
    while (c < s.cols.size()) {
        int64_t sc = s.cols[c] / schunk;
        size_t c_end = c;
        while (c_end < s.cols.size() && s.cols[c_end] / schunk == sc) {
            c_end++;
        }
        for (int64_t pc = 0; pc * pchunk < s.ploidy; pc++) {
            read_chunk_int32(a, { vc, sc, pc }, raw);
            for (int64_t r = 0; r < nrows; r++) {
                for (size_t k = c; k < c_end; k++) {
                    int64_t lc = s.cols[k] - sc * schunk;
                    for (int64_t lp = 0; lp < pchunk; lp++) {
                        int64_t p = pc * pchunk + lp;
                        if (p >= s.ploidy) {
                            break;
                        }
                        size_t src = static_cast<size_t>((r * schunk + lc) * pchunk + lp);
                        size_t dst = static_cast<size_t>(r) * width
                                     + k * static_cast<size_t>(s.ploidy)
                                     + static_cast<size_t>(p);
                        buf.data[dst] = raw[src];
                    }
                }
            }
        }
        c = c_end;
    }
}

} // namespace

void
iterate_genotypes(const GenotypeView &view, size_t begin, size_t end,
    const std::function<void(size_t, const int8_t *)> &f)
{
    const size_t total = view.num_haplotypes();
    std::vector<int8_t> out(total);
    std::vector<ChunkBuffer> buffers(view.sources.size());
    for (size_t k = begin; k < end; k++) {
        size_t offset = 0;
        for (size_t j = 0; j < view.sources.size(); j++) {
            const ViewSource &s = view.sources[j];
            const size_t width = s.num_haplotypes();
            int8_t *dest = out.data() + offset;
            offset += width;
            const int64_t row = s.rows[k];
            if (row < 0) {
                std::fill(dest, dest + width, static_cast<int8_t>(-1));
                continue;
            }
            ChunkBuffer &buf = buffers[j];
            int64_t vc = row / s.gt->chunks[0];
            if (buf.chunk != vc) {
                load_variant_chunk(s, vc, buf);
            }
            const int32_t *raw = buf.data.data() + static_cast<size_t>(row - buf.first_row) * width;
            const int8_t *remap = s.remap.data() + k * static_cast<size_t>(s.max_alleles);
            const int max_a = s.max_alleles;
            for (size_t h = 0; h < width; h++) {
                int32_t v = raw[h];
                if (v < 0) {
                    dest[h] = -1;
                } else {
                    /* Out-of-range allele indexes are clipped, as in tsinfer. */
                    int idx = v >= max_a ? max_a - 1 : v;
                    dest[h] = max_a > 0 ? remap[idx] : static_cast<int8_t>(-1);
                }
            }
        }
        f(k, out.data());
    }
}

} // namespace tsir

// [[Rcpp::export(.view_site_stats)]]
List
view_site_stats(List spec)
{
    tsir::GenotypeView view = tsir::view_from_list(spec);
    const size_t n = view.num_positions;
    const size_t haps = view.num_haplotypes();
    IntegerVector derived(n), non_missing(n);
    try {
        tsir::iterate_genotypes(view, 0, n, [&](size_t k, const int8_t *g) {
            int d = 0, m = 0;
            for (size_t h = 0; h < haps; h++) {
                if (g[h] > 0) {
                    d++;
                }
                if (g[h] >= 0) {
                    m++;
                }
            }
            derived[k] = d;
            non_missing[k] = m;
        });
    } catch (std::exception &e) {
        stop(e.what());
    }
    return List::create(_["derived_count"] = derived, _["non_missing"] = non_missing);
}

// [[Rcpp::export(.view_genotypes)]]
IntegerMatrix
view_genotypes(List spec, double begin, double end)
{
    tsir::GenotypeView view = tsir::view_from_list(spec);
    size_t b = static_cast<size_t>(begin), e = static_cast<size_t>(end);
    if (e > view.num_positions || b > e) {
        stop("position range out of bounds");
    }
    const size_t haps = view.num_haplotypes();
    IntegerMatrix out(static_cast<int>(e - b), static_cast<int>(haps));
    try {
        tsir::iterate_genotypes(view, b, e, [&](size_t k, const int8_t *g) {
            for (size_t h = 0; h < haps; h++) {
                out(static_cast<int>(k - b), static_cast<int>(h)) = g[h];
            }
        });
    } catch (std::exception &ex) {
        stop(ex.what());
    }
    return out;
}
