#include "common.h"
#include "lowlevel.h"
#include "tables.h"
#include "threads.h"
#include "zarr_r.h"

#include <algorithm>
#include <array>
#include <atomic>
#include <cfloat>
#include <cmath>
#include <set>
#include <thread>
#include <condition_variable>
#include <map>
#include <memory>
#include <mutex>
#include <unordered_map>

using namespace Rcpp;

namespace {

const tsk_flags_t NODE_IS_PC_ANCESTOR = 1u << 16;
const double PC_ANCESTOR_INCREMENT = 1.0 / 4294967296.0;

struct HapSource {
    tsir::ZarrArray gt;
    int ploidy = 1;
    int64_t sample_chunk = 1;
    std::vector<int64_t> sel_rows; /* selected source variant rows (ascending) */
    std::vector<int> sel_ref;      /* reference site of each selected row */
    std::vector<int8_t> remap;     /* sel_rows.size() x max_alleles */
    int max_alleles = 0;
    double recombination = 1e-2;
    double mismatch = 1e-20;
    /* For each variant chunk containing selected rows: (chunk, first selected
     * row index, end selected row index). */
    std::vector<std::array<int64_t, 3>> needed;
};

struct Job {
    int source;
    int64_t store_col;
    int ploidy_index;
    double time;
    tsk_flags_t flags;
    tsk_id_t individual;
    tsk_id_t population;
    std::string metadata;
};

using ChunkData = std::shared_ptr<const std::vector<int8_t>>;

/* Chunk cache shared by the matching threads. Entries are reference counted
 * using the read schedule and evicted once no further reads need them, or when
 * the memory budget is exceeded. */
class ChunkCache {
  public:
    explicit ChunkCache(size_t max_bytes) : max_bytes_(max_bytes) {}

    void add_reads(int64_t key, int64_t count)
    {
        std::lock_guard<std::mutex> lock(mutex_);
        refcount_[key] += count;
    }

    template <typename Loader> ChunkData get(int64_t key, Loader load)
    {
        std::unique_lock<std::mutex> lock(mutex_);
        while (true) {
            auto it = entries_.find(key);
            if (it != entries_.end()) {
                return it->second;
            }
            if (loading_.count(key) == 0) {
                break;
            }
            cond_.wait(lock);
        }
        loading_.insert(key);
        lock.unlock();
        ChunkData data;
        try {
            data = std::make_shared<const std::vector<int8_t>>(load());
        } catch (...) {
            lock.lock();
            loading_.erase(key);
            cond_.notify_all();
            throw;
        }
        lock.lock();
        loading_.erase(key);
        bytes_ += data->size();
        entries_[key] = data;
        evict_over_budget(key);
        cond_.notify_all();
        return data;
    }

    void record_read(int64_t key)
    {
        std::lock_guard<std::mutex> lock(mutex_);
        auto rc = refcount_.find(key);
        if (rc == refcount_.end()) {
            return;
        }
        if (--rc->second <= 0) {
            refcount_.erase(rc);
            auto it = entries_.find(key);
            if (it != entries_.end()) {
                bytes_ -= it->second->size();
                entries_.erase(it);
            }
        }
    }

    size_t bytes()
    {
        std::lock_guard<std::mutex> lock(mutex_);
        return bytes_;
    }

  private:
    void evict_over_budget(int64_t keep)
    {
        while (bytes_ > max_bytes_ && entries_.size() > 1) {
            int64_t victim = -1;
            int64_t best = INT64_MAX;
            for (const auto &e : entries_) {
                if (e.first == keep) {
                    continue;
                }
                auto rc = refcount_.find(e.first);
                int64_t refs = rc == refcount_.end() ? 0 : rc->second;
                if (refs < best) {
                    best = refs;
                    victim = e.first;
                }
            }
            if (victim < 0) {
                break;
            }
            bytes_ -= entries_[victim]->size();
            entries_.erase(victim);
        }
    }

    size_t max_bytes_;
    size_t bytes_ = 0;
    std::mutex mutex_;
    std::condition_variable cond_;
    std::unordered_map<int64_t, ChunkData> entries_;
    std::unordered_map<int64_t, int64_t> refcount_;
    std::set<int64_t> loading_;
};

struct MatchEngine {
    std::vector<double> positions;
    std::vector<tsk_size_t> num_alleles;
    std::vector<std::vector<std::string>> alleles;
    std::vector<HapSource> sources;
    std::vector<Job> jobs;
    std::unique_ptr<ChunkCache> cache;

    int64_t cache_key(const Job &job) const
    {
        const HapSource &s = sources[job.source];
        return static_cast<int64_t>(job.source) * (INT64_C(1) << 40)
               + job.store_col / s.sample_chunk;
    }

    std::vector<int8_t> load_chunk(int source, int64_t sc) const
    {
        const HapSource &s = sources[source];
        const tsir::ZarrArray &a = s.gt;
        const int64_t chunk_samples = std::min(s.sample_chunk, a.shape[1] - sc * s.sample_chunk);
        const size_t nsel = s.sel_rows.size();
        const size_t stride = static_cast<size_t>(chunk_samples) * s.ploidy;
        std::vector<int8_t> out(nsel * stride, -1);
        std::vector<int32_t> raw;
        const int64_t vchunk = a.chunks[0];
        const int64_t pchunk = a.chunks[2];
        for (const auto &need : s.needed) {
            const int64_t vc = need[0];
            for (int64_t pc = 0; pc * pchunk < s.ploidy; pc++) {
                tsir::read_chunk_int32(a, { vc, sc, pc }, raw);
                for (int64_t i = need[1]; i < need[2]; i++) {
                    const int64_t local = s.sel_rows[i] - vc * vchunk;
                    const int8_t *remap = s.remap.data() + i * s.max_alleles;
                    for (int64_t c = 0; c < chunk_samples; c++) {
                        for (int64_t lp = 0; lp < pchunk; lp++) {
                            const int64_t p = pc * pchunk + lp;
                            if (p >= s.ploidy) {
                                break;
                            }
                            int32_t v = raw[static_cast<size_t>(
                                (local * a.chunks[1] + c) * pchunk + lp)];
                            int8_t code = -1;
                            if (v >= 0 && s.max_alleles > 0) {
                                code = remap[v >= s.max_alleles ? s.max_alleles - 1 : v];
                            }
                            out[static_cast<size_t>(i) * stride
                                + static_cast<size_t>(c * s.ploidy + p)]
                                = code;
                        }
                    }
                }
            }
        }
        return out;
    }

    void read_haplotype(const Job &job, std::vector<allele_t> &hap)
    {
        const HapSource &s = sources[job.source];
        const int64_t sc = job.store_col / s.sample_chunk;
        const int64_t offset = job.store_col % s.sample_chunk;
        const int64_t key = cache_key(job);
        ChunkData data = cache->get(key, [&]() { return load_chunk(job.source, sc); });
        const int64_t chunk_samples
            = std::min(s.sample_chunk, s.gt.shape[1] - sc * s.sample_chunk);
        const size_t stride = static_cast<size_t>(chunk_samples) * s.ploidy;
        std::fill(hap.begin(), hap.end(), static_cast<allele_t>(-1));
        for (size_t i = 0; i < s.sel_rows.size(); i++) {
            hap[s.sel_ref[i]]
                = (*data)[i * stride + static_cast<size_t>(offset * s.ploidy + job.ploidy_index)];
        }
        cache->record_read(key);
    }
};

struct MatchResult {
    std::vector<tsk_id_t> left, right, parent;
    std::vector<std::pair<int, int8_t>> mutations;
    double mean_traceback = 0;
};

MatchEngine *
get_engine(SEXP xp)
{
    XPtr<MatchEngine> p(xp);
    return p.get();
}

} // namespace

// [[Rcpp::export(.match_engine_new)]]
SEXP
match_engine_new(NumericVector positions, List alleles, List sources, List jobs,
    double cache_bytes)
{
    auto engine = std::make_unique<MatchEngine>();
    engine->positions.assign(positions.begin(), positions.end());
    const size_t num_sites = engine->positions.size();
    if (static_cast<size_t>(alleles.size()) != num_sites) {
        stop("alleles must have one entry per site");
    }
    for (size_t j = 0; j < num_sites; j++) {
        CharacterVector a = alleles[j];
        std::vector<std::string> v;
        for (R_xlen_t k = 0; k < a.size(); k++) {
            v.push_back(Rf_translateCharUTF8(a[k]));
        }
        if (v.size() < 2) {
            v.resize(2);
        }
        if (v.size() > 127) {
            stop("Site %d has more than 127 alleles", static_cast<int>(j));
        }
        engine->num_alleles.push_back(static_cast<tsk_size_t>(a.size() < 2 ? 2 : a.size()));
        engine->alleles.push_back(std::move(v));
    }
    for (R_xlen_t j = 0; j < sources.size(); j++) {
        List s = sources[j];
        HapSource hs;
        hs.gt = tsir::zarr_array_from_list(s["gt_meta"]);
        if (hs.gt.ndim() != 3) {
            stop("call_genotype must be 3-dimensional");
        }
        hs.ploidy = as<int>(s["ploidy"]);
        hs.sample_chunk = hs.gt.chunks[1];
        NumericVector rows = s["sel_rows"];
        IntegerVector ref = s["sel_ref"];
        IntegerMatrix remap = s["remap"];
        hs.sel_rows.assign(rows.begin(), rows.end());
        hs.sel_ref.assign(ref.begin(), ref.end());
        hs.max_alleles = remap.ncol();
        hs.remap.resize(hs.sel_rows.size() * hs.max_alleles);
        for (size_t r = 0; r < hs.sel_rows.size(); r++) {
            if (hs.sel_ref[r] < 0 || static_cast<size_t>(hs.sel_ref[r]) >= num_sites) {
                stop("selected reference site out of range");
            }
            for (int a = 0; a < hs.max_alleles; a++) {
                hs.remap[r * hs.max_alleles + a]
                    = static_cast<int8_t>(remap(static_cast<int>(r), a));
            }
        }
        hs.recombination = as<double>(s["recombination"]);
        hs.mismatch = as<double>(s["mismatch"]);
        const int64_t vchunk = hs.gt.chunks[0];
        size_t i = 0;
        while (i < hs.sel_rows.size()) {
            int64_t vc = hs.sel_rows[i] / vchunk;
            size_t e = i;
            while (e < hs.sel_rows.size() && hs.sel_rows[e] / vchunk == vc) {
                e++;
            }
            hs.needed.push_back({ vc, static_cast<int64_t>(i), static_cast<int64_t>(e) });
            i = e;
        }
        engine->sources.push_back(std::move(hs));
    }
    IntegerVector src = jobs["source"];
    NumericVector col = jobs["store_col"];
    IntegerVector ploidy_index = jobs["ploidy_index"];
    NumericVector time = jobs["time"];
    NumericVector flags = jobs["node_flags"];
    IntegerVector individual = jobs["individual_id"];
    IntegerVector population = jobs["population_id"];
    CharacterVector metadata = jobs["metadata"];
    for (R_xlen_t j = 0; j < src.size(); j++) {
        Job job;
        job.source = src[j];
        job.store_col = static_cast<int64_t>(col[j]);
        job.ploidy_index = ploidy_index[j];
        job.time = time[j];
        job.flags = static_cast<tsk_flags_t>(flags[j]);
        job.individual = individual[j] == NA_INTEGER ? TSK_NULL : individual[j];
        job.population = population[j] == NA_INTEGER ? TSK_NULL : population[j];
        job.metadata = Rf_translateCharUTF8(metadata[j]);
        if (job.source < 0 || static_cast<size_t>(job.source) >= engine->sources.size()) {
            stop("job source out of range");
        }
        const HapSource &hs = engine->sources[job.source];
        if (job.store_col < 0 || job.store_col >= hs.gt.shape[1] || job.ploidy_index < 0
            || job.ploidy_index >= hs.ploidy) {
            stop("job sample column or ploidy index out of range");
        }
        engine->jobs.push_back(std::move(job));
    }
    engine->cache = std::make_unique<ChunkCache>(static_cast<size_t>(cache_bytes));
    XPtr<MatchEngine> xp(engine.release(), true);
    return xp;
}

// [[Rcpp::export(.match_engine_schedule)]]
void
match_engine_schedule(SEXP engine_xp, IntegerVector job_ids)
{
    MatchEngine *engine = get_engine(engine_xp);
    std::map<int64_t, int64_t> counts;
    for (int j : job_ids) {
        counts[engine->cache_key(engine->jobs.at(static_cast<size_t>(j)))]++;
    }
    for (const auto &c : counts) {
        engine->cache->add_reads(c.first, c.second);
    }
}

static void
compress_paths(tsk_table_collection_t *tables)
{
    const tsk_edge_table_t &edges = tables->edges;
    const double *node_time = tables->nodes.time;
    /* Path of each child: its (left, right, parent) edges, in edge table order
     * of first appearance of the child. */
    std::vector<tsk_id_t> child_order;
    std::unordered_map<tsk_id_t, std::vector<std::tuple<double, double, tsk_id_t>>> paths;
    for (tsk_size_t j = 0; j < edges.num_rows; j++) {
        tsk_id_t c = edges.child[j];
        auto it = paths.find(c);
        if (it == paths.end()) {
            child_order.push_back(c);
            it = paths.emplace(c, std::vector<std::tuple<double, double, tsk_id_t>>()).first;
        }
        it->second.emplace_back(edges.left[j], edges.right[j], edges.parent[j]);
    }
    std::map<std::vector<std::tuple<double, double, tsk_id_t>>, size_t> key_index;
    std::vector<std::vector<std::tuple<double, double, tsk_id_t>>> keys;
    std::vector<std::vector<tsk_id_t>> groups;
    for (tsk_id_t c : child_order) {
        auto path = paths[c];
        std::sort(path.begin(), path.end());
        auto it = key_index.find(path);
        if (it == key_index.end()) {
            key_index.emplace(path, keys.size());
            keys.push_back(path);
            groups.push_back({ c });
        } else {
            groups[it->second].push_back(c);
        }
    }
    std::vector<char> compressed(tables->nodes.num_rows, 0);
    std::vector<std::array<double, 2>> new_lr;
    std::vector<std::array<tsk_id_t, 2>> new_pc;
    for (size_t g = 0; g < keys.size(); g++) {
        const auto &path = keys[g];
        const auto &children = groups[g];
        if (children.size() < 2) {
            continue;
        }
        /* Children whose whole path is to a single path-compression ancestor
         * are already compressed; compressing again would only add unary
         * nodes on every call. */
        tsk_id_t first_parent = std::get<2>(path[0]);
        bool single_parent = true;
        for (const auto &e : path) {
            if (std::get<2>(e) != first_parent) {
                single_parent = false;
            }
        }
        if (single_parent && (tables->nodes.flags[first_parent] & NODE_IS_PC_ANCESTOR)) {
            continue;
        }
        double min_parent_time = INFINITY;
        for (const auto &e : path) {
            min_parent_time = std::min(min_parent_time, node_time[std::get<2>(e)]);
        }
        double pc_time = min_parent_time - PC_ANCESTOR_INCREMENT;
        double max_child_time = -INFINITY;
        for (tsk_id_t c : children) {
            max_child_time = std::max(max_child_time, node_time[c]);
        }
        if (pc_time <= max_child_time) {
            continue;
        }
        tsk_id_t pc = tsk_node_table_add_row(&tables->nodes, NODE_IS_PC_ANCESTOR, pc_time,
            TSK_NULL, TSK_NULL, "{}", 2);
        tsir::check_tsk(pc, "add path compression node");
        node_time = tables->nodes.time;
        compressed.resize(tables->nodes.num_rows, 0);
        for (const auto &e : path) {
            new_lr.push_back({ std::get<0>(e), std::get<1>(e) });
            new_pc.push_back({ std::get<2>(e), pc });
        }
        for (tsk_id_t c : children) {
            for (const auto &e : path) {
                new_lr.push_back({ std::get<0>(e), std::get<1>(e) });
                new_pc.push_back({ pc, c });
            }
            compressed[c] = 1;
        }
    }
    if (new_lr.empty()) {
        return;
    }
    tsk_edge_table_t kept;
    tsir::check_tsk(tsk_edge_table_init(&kept, 0), "edge table");
    for (tsk_size_t j = 0; j < edges.num_rows; j++) {
        if (!compressed[edges.child[j]]) {
            int ret = tsk_edge_table_add_row(&kept, edges.left[j], edges.right[j],
                edges.parent[j], edges.child[j], nullptr, 0);
            if (ret < 0) {
                tsk_edge_table_free(&kept);
                tsir::stop_tsk(ret, "add edge");
            }
        }
    }
    for (size_t j = 0; j < new_lr.size(); j++) {
        int ret = tsk_edge_table_add_row(
            &kept, new_lr[j][0], new_lr[j][1], new_pc[j][0], new_pc[j][1], nullptr, 0);
        if (ret < 0) {
            tsk_edge_table_free(&kept);
            tsir::stop_tsk(ret, "add edge");
        }
    }
    int ret = tsk_edge_table_copy(&kept, &tables->edges, TSK_NO_INIT);
    tsk_edge_table_free(&kept);
    tsir::check_tsk(ret, "copy edges");
    tsir::finalise_tables(tables, false);
}

// [[Rcpp::export(.compress_paths)]]
SEXP
compress_paths_r(SEXP tables_xp)
{
    tsk_table_collection_t *src = tsir::get_tables(tables_xp);
    tsk_table_collection_t *tables;
    SEXP out = PROTECT(tsir::new_tables_xptr(&tables));
    tsir::check_tsk(tsk_table_collection_copy(src, tables, 0), "copy");
    compress_paths(tables);
    UNPROTECT(1);
    return out;
}

// [[Rcpp::export(.match_group)]]
List
match_group(SEXP engine_xp, SEXP tables_xp, IntegerVector job_ids, int num_threads,
    bool path_compression, std::string provenance, bool return_results)
{
    MatchEngine *engine = get_engine(engine_xp);
    tsk_table_collection_t *tables = tsir::get_tables(tables_xp);
    const size_t num_sites = engine->positions.size();
    if (num_sites == 0) {
        stop("Cannot match haplotypes: there are no inference sites");
    }
    if (tables->sites.num_rows != num_sites) {
        stop("The tree sequence does not have one site per inference site");
    }
    if (!tsk_table_collection_has_index(tables, 0)) {
        tsir::check_tsk(tsk_table_collection_build_index(tables, 0), "build_index");
    }
    tsir::MatcherIndexesHandle mi(tables, engine->num_alleles.data());
    const size_t n = static_cast<size_t>(job_ids.size());
    std::vector<int> ids(job_ids.begin(), job_ids.end());
    for (int id : ids) {
        if (id < 0 || static_cast<size_t>(id) >= engine->jobs.size()) {
            stop("job id out of range");
        }
    }
    std::vector<MatchResult> results(n);
    const int threads = std::max(1, std::min(tsir::resolve_num_threads(num_threads),
                                        static_cast<int>(n == 0 ? 1 : n)));
    std::atomic<size_t> next(0);
    std::exception_ptr error = nullptr;
    std::mutex error_mutex;
    auto worker = [&]() {
        std::map<int, std::unique_ptr<tsir::AncestorMatcherHandle>> matchers;
        std::vector<allele_t> hap(num_sites), match(num_sites);
        std::vector<tsk_id_t> left(num_sites + 1), right(num_sites + 1), parent(num_sites + 1);
        while (true) {
            size_t i = next.fetch_add(1);
            if (i >= n) {
                break;
            }
            try {
                const Job &job = engine->jobs[ids[i]];
                auto &m = matchers[job.source];
                if (!m) {
                    const HapSource &s = engine->sources[job.source];
                    m = std::make_unique<tsir::AncestorMatcherHandle>(&mi.indexes,
                        std::vector<double>(num_sites, s.recombination),
                        std::vector<double>(num_sites, s.mismatch), DBL_MIN, 0);
                }
                engine->read_haplotype(job, hap);
                tsk_id_t start = 0, end = static_cast<tsk_id_t>(num_sites);
                size_t first = 0;
                while (first < num_sites && hap[first] < 0) {
                    first++;
                }
                if (first < num_sites) {
                    size_t last = num_sites;
                    while (hap[last - 1] < 0) {
                        last--;
                    }
                    start = static_cast<tsk_id_t>(first);
                    end = static_cast<tsk_id_t>(last);
                }
                size_t path_length = 0;
                int ret = ancestor_matcher_find_path(&m->matcher, start, end, hap.data(),
                    match.data(), &path_length, left.data(), right.data(), parent.data());
                if (ret != 0) {
                    throw std::runtime_error(tsi_strerror(ret));
                }
                MatchResult &r = results[i];
                r.left.assign(left.begin(), left.begin() + path_length);
                r.right.assign(right.begin(), right.begin() + path_length);
                r.parent.assign(parent.begin(), parent.begin() + path_length);
                for (tsk_id_t s = start; s < end; s++) {
                    if (hap[s] >= 0 && hap[s] != match[s]) {
                        r.mutations.emplace_back(s, hap[s]);
                    }
                }
                r.mean_traceback = ancestor_matcher_get_mean_traceback_size(&m->matcher);
            } catch (...) {
                std::lock_guard<std::mutex> lock(error_mutex);
                if (error == nullptr) {
                    error = std::current_exception();
                }
                next.store(n);
                break;
            }
        }
    };
    {
        std::vector<std::thread> pool;
        for (int t = 0; t < threads; t++) {
            pool.emplace_back(worker);
        }
        for (std::thread &t : pool) {
            t.join();
        }
    }
    if (error != nullptr) {
        try {
            std::rethrow_exception(error);
        } catch (std::exception &e) {
            stop("matching failed: %s", e.what());
        }
    }

    /* Extend the tables with the new nodes, edges and mutations. */
    std::vector<std::tuple<int, int8_t, tsk_id_t>> muts;
    double total_traceback = 0;
    for (size_t i = 0; i < n; i++) {
        const Job &job = engine->jobs[ids[i]];
        tsk_id_t node = tsk_node_table_add_row(&tables->nodes, job.flags, job.time,
            job.population, job.individual, job.metadata.data(), job.metadata.size());
        tsir::check_tsk(node, "add node");
        const MatchResult &r = results[i];
        for (size_t k = 0; k < r.left.size(); k++) {
            tsir::check_tsk(tsk_edge_table_add_row(&tables->edges, r.left[k], r.right[k],
                                r.parent[k], node, nullptr, 0),
                "add edge");
        }
        for (const auto &m : r.mutations) {
            muts.emplace_back(m.first, m.second, node);
        }
        total_traceback += r.mean_traceback;
    }
    for (const auto &m : muts) {
        const std::string &state = engine->alleles[std::get<0>(m)][std::get<1>(m)];
        tsir::check_tsk(tsk_mutation_table_add_row(&tables->mutations, std::get<0>(m),
                            std::get<2>(m), TSK_NULL, TSK_UNKNOWN_TIME, state.data(),
                            state.size(), nullptr, 0),
            "add mutation");
    }
    tsir::check_tsk(tsk_table_collection_sort(tables, nullptr, 0), "sort");
    tsir::check_tsk(tsk_table_collection_build_index(tables, 0), "build_index");
    tsir::check_tsk(tsk_table_collection_compute_mutation_parents(tables, 0),
        "compute_mutation_parents");
    tsir::add_provenance(tables, provenance);
    tsk_id_t check = tsk_table_collection_check_integrity(tables, TSK_CHECK_TREES);
    if (check < 0) {
        tsir::stop_tsk(static_cast<int>(check), "tree sequence integrity check");
    }
    if (path_compression) {
        compress_paths(tables);
    }
    List details;
    if (return_results) {
        details = List(n);
        for (size_t i = 0; i < n; i++) {
            const MatchResult &r = results[i];
            IntegerVector site(r.mutations.size());
            IntegerVector derived(r.mutations.size());
            for (size_t k = 0; k < r.mutations.size(); k++) {
                site[k] = r.mutations[k].first;
                derived[k] = r.mutations[k].second;
            }
            details[i] = List::create(_["left"] = IntegerVector(r.left.begin(), r.left.end()),
                _["right"] = IntegerVector(r.right.begin(), r.right.end()),
                _["parent"] = IntegerVector(r.parent.begin(), r.parent.end()),
                _["site"] = site, _["derived_state"] = derived);
        }
    }
    return List::create(_["num_matched"] = static_cast<double>(n),
        _["num_mutations"] = static_cast<double>(muts.size()),
        _["mean_traceback_size"] = n > 0 ? total_traceback / n : 0.0,
        _["cache_bytes"] = static_cast<double>(engine->cache->bytes()),
        _["results"] = details);
}

// [[Rcpp::export(.make_root_tables)]]
SEXP
make_root_tables(double sequence_length, NumericVector positions,
    CharacterVector ancestral, std::string metadata, double max_time,
    CharacterVector individuals, Nullable<CharacterVector> populations)
{
    tsk_table_collection_t *tables;
    SEXP xp = PROTECT(tsir::new_tables_xptr(&tables));
    tsir::check_tsk(tsk_table_collection_init(tables, 0), "init");
    tables->sequence_length = sequence_length;
    const std::string schema = "{\"codec\":\"json\"}";
    tsir::check_tsk(
        tsk_table_collection_set_metadata_schema(tables, schema.data(), schema.size()),
        "schema");
    tsir::check_tsk(tsk_table_collection_set_metadata(tables, metadata.data(), metadata.size()),
        "metadata");
    for (R_xlen_t j = 0; j < positions.size(); j++) {
        std::string a = Rf_translateCharUTF8(ancestral[j]);
        tsir::check_tsk(tsk_site_table_add_row(
                            &tables->sites, positions[j], a.data(), a.size(), nullptr, 0),
            "add site");
    }
    tsir::check_tsk(
        tsk_node_table_set_metadata_schema(&tables->nodes, schema.data(), schema.size()),
        "schema");
    tsir::check_tsk(tsk_individual_table_set_metadata_schema(
                        &tables->individuals, schema.data(), schema.size()),
        "schema");
    if (populations.isNotNull()) {
        CharacterVector pops(populations);
        tsir::check_tsk(tsk_population_table_set_metadata_schema(
                            &tables->populations, schema.data(), schema.size()),
            "schema");
        for (R_xlen_t j = 0; j < pops.size(); j++) {
            std::string md = Rf_translateCharUTF8(pops[j]);
            tsir::check_tsk(
                tsk_population_table_add_row(&tables->populations, md.data(), md.size()),
                "add population");
        }
    }
    for (R_xlen_t j = 0; j < individuals.size(); j++) {
        std::string md = Rf_translateCharUTF8(individuals[j]);
        tsir::check_tsk(tsk_individual_table_add_row(&tables->individuals, 0, nullptr, 0,
                            nullptr, 0, md.data(), md.size()),
            "add individual");
    }
    tsir::check_tsk(tsk_node_table_add_row(
                        &tables->nodes, 0, max_time + 2, TSK_NULL, TSK_NULL, "{}", 2),
        "add node");
    tsir::check_tsk(tsk_node_table_add_row(
                        &tables->nodes, 0, max_time + 1, TSK_NULL, TSK_NULL, "{}", 2),
        "add node");
    if (positions.size() > 0) {
        tsir::check_tsk(tsk_edge_table_add_row(
                            &tables->edges, positions[0], sequence_length, 0, 1, nullptr, 0),
            "add edge");
    }
    tsir::check_tsk(tsk_table_collection_sort(tables, nullptr, 0), "sort");
    tsir::check_tsk(tsk_table_collection_build_index(tables, 0), "build_index");
    UNPROTECT(1);
    return xp;
}

// [[Rcpp::export(.erase_flanks)]]
SEXP
erase_flanks(SEXP tables_xp, NumericMatrix intervals)
{
    tsk_table_collection_t *src = tsir::get_tables(tables_xp);
    tsk_table_collection_t *tables;
    SEXP out = PROTECT(tsir::new_tables_xptr(&tables));
    tsir::check_tsk(tsk_table_collection_copy(src, tables, 0), "copy");
    tsir::check_tsk(tsk_edge_table_clear(&tables->edges), "clear edges");
    const tsk_edge_table_t &e = src->edges;
    for (tsk_size_t j = 0; j < e.num_rows; j++) {
        for (int k = 0; k < intervals.nrow(); k++) {
            double l = std::max(e.left[j], intervals(k, 0));
            double r = std::min(e.right[j], intervals(k, 1));
            if (l < r) {
                tsir::check_tsk(tsk_edge_table_add_row(
                                    &tables->edges, l, r, e.parent[j], e.child[j], nullptr, 0),
                    "add edge");
            }
        }
    }
    tsir::check_tsk(tsk_table_collection_sort(tables, nullptr, 0), "sort");
    tsir::check_tsk(tsk_table_collection_build_index(tables, 0), "build_index");
    UNPROTECT(1);
    return out;
}

// [[Rcpp::export(.json_quote)]]
CharacterVector
json_quote_r(CharacterVector x)
{
    CharacterVector out(x.size());
    for (R_xlen_t j = 0; j < x.size(); j++) {
        out[j] = tsir::json_quote(Rf_translateCharUTF8(x[j]));
    }
    return out;
}

// [[Rcpp::export(.json_double)]]
CharacterVector
json_double_r(NumericVector x)
{
    CharacterVector out(x.size());
    for (R_xlen_t j = 0; j < x.size(); j++) {
        out[j] = tsir::json_double(x[j]);
    }
    return out;
}
