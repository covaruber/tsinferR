#include "common.h"
#include "genotypes.h"

#include <algorithm>
#include <cstdlib>

using namespace Rcpp;

namespace {

/* Owns a tree sequence and a tree over it for the duration of a call. */
struct TreeHolder {
    tsk_treeseq_t ts;
    tsk_tree_t tree;
    bool ts_ok = false;
    bool tree_ok = false;

    explicit TreeHolder(tsk_table_collection_t *tables)
    {
        int err = tsk_treeseq_init(&ts, tables, TSK_TS_INIT_BUILD_INDEXES);
        ts_ok = true;
        if (err == 0) {
            err = tsk_tree_init(&tree, &ts, 0);
            tree_ok = true;
        }
        if (err != 0) {
            release();
            tsir::stop_tsk(err, "tree sequence");
        }
    }

    void release()
    {
        if (tree_ok) {
            tsk_tree_free(&tree);
            tree_ok = false;
        }
        if (ts_ok) {
            tsk_treeseq_free(&ts);
            ts_ok = false;
        }
    }

    ~TreeHolder() { release(); }
};

} // namespace

// [[Rcpp::export(.augment_sites)]]
List
augment_sites_cpp(SEXP tables_xp, List view_spec, NumericVector positions, List alleles,
    IntegerVector permutation)
{
    tsk_table_collection_t *src = tsir::get_tables(tables_xp);
    tsir::GenotypeView view = tsir::view_from_list(view_spec);
    const size_t n_sites = view.num_positions;
    if (static_cast<size_t>(positions.size()) != n_sites
        || static_cast<size_t>(alleles.size()) != n_sites) {
        stop("positions and alleles must have one entry per site");
    }
    std::vector<std::vector<std::string>> site_alleles(n_sites);
    for (size_t k = 0; k < n_sites; k++) {
        CharacterVector a = alleles[k];
        for (R_xlen_t j = 0; j < a.size(); j++) {
            site_alleles[k].push_back(Rf_translateCharUTF8(a[j]));
        }
    }
    tsk_table_collection_t *tables;
    SEXP out = PROTECT(tsir::new_tables_xptr(&tables));
    tsir::check_tsk(tsk_table_collection_copy(src, tables, 0), "copy");
    TreeHolder holder(src);
    const tsk_size_t num_samples = tsk_treeseq_get_num_samples(&holder.ts);
    if (static_cast<tsk_size_t>(permutation.size()) != num_samples) {
        UNPROTECT(1);
        stop("permutation must have one entry per sample");
    }
    const size_t haps = view.num_haplotypes();
    std::vector<int32_t> genotypes(num_samples);
    int added = 0;
    try {
        tsir::iterate_genotypes(view, 0, n_sites, [&](size_t k, const int8_t *g) {
            bool any = false;
            for (tsk_size_t j = 0; j < num_samples; j++) {
                int p = permutation[j];
                genotypes[j] = (p >= 0 && static_cast<size_t>(p) < haps) ? g[p] : -1;
                any = any || genotypes[j] >= 0;
            }
            if (!any) {
                return;
            }
            const std::vector<std::string> &al = site_alleles[k];
            std::vector<std::string> states = al;
            tsk_flags_t options = 0;
            int32_t ancestral = 0;
            if (!al.empty() && !al[0].empty()) {
                options = TSK_MM_FIXED_ANCESTRAL_STATE;
            } else {
                /* No ancestral allele: drop empty alleles and let parsimony
                 * choose the root state. */
                states.clear();
                std::vector<int32_t> code_to_real(al.size(), -1);
                for (size_t c = 0; c < al.size(); c++) {
                    if (!al[c].empty()) {
                        code_to_real[c] = static_cast<int32_t>(states.size());
                        states.push_back(al[c]);
                    }
                }
                for (tsk_size_t j = 0; j < num_samples; j++) {
                    int32_t v = genotypes[j];
                    genotypes[j] = (v >= 0 && static_cast<size_t>(v) < al.size())
                                       ? code_to_real[v]
                                       : -1;
                }
            }
            int32_t max_state = 0;
            for (int32_t v : genotypes) {
                max_state = std::max(max_state, v);
            }
            if (max_state >= 64) {
                throw std::runtime_error("A maximum of 64 states is supported");
            }
            if (static_cast<size_t>(max_state) >= states.size()) {
                throw std::runtime_error("genotype refers to an unknown allele");
            }
            tsir::check_tsk(tsk_tree_seek(&holder.tree, positions[k], 0), "seek");
            tsk_size_t num_transitions = 0;
            tsk_state_transition_t *transitions = nullptr;
            int err = tsk_tree_map_mutations(&holder.tree, genotypes.data(), nullptr, options,
                &ancestral, &num_transitions, &transitions);
            if (err != 0) {
                std::free(transitions);
                throw std::runtime_error(std::string("map_mutations: ") + tsk_strerror(err));
            }
            const std::string &anc = states[ancestral];
            tsk_id_t site = tsk_site_table_add_row(
                &tables->sites, positions[k], anc.data(), anc.size(), nullptr, 0);
            std::vector<tsk_id_t> ids(num_transitions);
            for (tsk_size_t t = 0; t < num_transitions && site >= 0; t++) {
                const tsk_state_transition_t &tr = transitions[t];
                const std::string &d = states[tr.state];
                tsk_id_t parent = tr.parent == TSK_NULL ? TSK_NULL : ids[tr.parent];
                ids[t] = tsk_mutation_table_add_row(&tables->mutations, site, tr.node, parent,
                    TSK_UNKNOWN_TIME, d.data(), d.size(), nullptr, 0);
                if (ids[t] < 0) {
                    site = ids[t];
                }
            }
            std::free(transitions);
            if (site < 0) {
                throw std::runtime_error(tsk_strerror(site));
            }
            added++;
        });
    } catch (std::exception &e) {
        UNPROTECT(1);
        stop(e.what());
    }
    tsir::check_tsk(tsk_table_collection_sort(tables, nullptr, 0), "sort");
    tsir::check_tsk(tsk_table_collection_build_index(tables, 0), "build_index");
    tsk_id_t check = tsk_table_collection_check_integrity(tables, TSK_CHECK_TREES);
    if (check < 0) {
        UNPROTECT(1);
        tsir::stop_tsk(static_cast<int>(check), "tree sequence integrity check");
    }
    UNPROTECT(1);
    return List::create(_["tables"] = out, _["num_added"] = added);
}
