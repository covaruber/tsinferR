#include "common.h"
#include "lowlevel.h"

#include <cfloat>

using namespace Rcpp;

namespace tsir {

AncestorBuilderHandle::AncestorBuilderHandle(
    size_t num_samples, size_t max_sites, int genotype_encoding)
{
    int err = ancestor_builder_alloc(&builder, num_samples, max_sites, -1, genotype_encoding);
    if (err != 0) {
        ancestor_builder_free(&builder);
        stop_tsi(err, "AncestorBuilder");
    }
}

AncestorBuilderHandle::~AncestorBuilderHandle() { ancestor_builder_free(&builder); }

void
AncestorBuilderHandle::finalise()
{
    if (finalised) {
        return;
    }
    check_tsi(ancestor_builder_finalise(&builder), "AncestorBuilder");
    finalised = true;
}

MatcherIndexesHandle::MatcherIndexesHandle(
    const tsk_table_collection_t *tables, const tsk_size_t *num_alleles)
{
    std::memset(&indexes, 0, sizeof(indexes));
    int err = matcher_indexes_alloc(&indexes, tables, num_alleles, 0);
    if (err != 0) {
        matcher_indexes_free(&indexes);
        stop_tsi(err, "MatcherIndexes");
    }
}

MatcherIndexesHandle::~MatcherIndexesHandle() { matcher_indexes_free(&indexes); }

AncestorMatcherHandle::AncestorMatcherHandle(const matcher_indexes_t *indexes,
    const std::vector<double> &recombination, const std::vector<double> &mismatch,
    double likelihood_threshold, int flags)
{
    std::vector<double> rho(recombination);
    std::vector<double> mu(mismatch);
    int err = ancestor_matcher_alloc(
        &matcher, indexes, rho.data(), mu.data(), likelihood_threshold, flags);
    if (err != 0) {
        ancestor_matcher_free(&matcher);
        throw std::runtime_error(tsi_strerror(err));
    }
}

AncestorMatcherHandle::~AncestorMatcherHandle() { ancestor_matcher_free(&matcher); }

} // namespace tsir

using tsir::AncestorBuilderHandle;
using tsir::AncestorMatcherHandle;
using tsir::MatcherIndexesHandle;

// [[Rcpp::export(.ab_new)]]
SEXP
ab_new(int num_samples, int max_sites, int genotype_encoding)
{
    if (num_samples < 0 || max_sites < 0) {
        stop("num_samples and max_sites must be non-negative");
    }
    if (genotype_encoding != 0 && genotype_encoding != 1) {
        stop("genotype_encoding must be 0 (eight-bit) or 1 (one-bit)");
    }
    auto *handle = new AncestorBuilderHandle(
        static_cast<size_t>(num_samples), static_cast<size_t>(max_sites), genotype_encoding);
    XPtr<AncestorBuilderHandle> xp(handle, true);
    return xp;
}

// [[Rcpp::export(.ab_add_site)]]
void
ab_add_site(SEXP xp, double time, IntegerVector genotypes)
{
    XPtr<AncestorBuilderHandle> handle(xp);
    if (handle->finalised) {
        stop("Cannot add sites after ancestor descriptors have been computed");
    }
    ancestor_builder_t *ab = &handle->builder;
    if (static_cast<size_t>(genotypes.size()) != ab->num_samples) {
        stop("genotypes array wrong size.");
    }
    std::vector<allele_t> g(ab->num_samples);
    for (size_t j = 0; j < ab->num_samples; j++) {
        int v = genotypes[j];
        if (v == NA_INTEGER) {
            v = -1;
        }
        if (v < -1 || v > 1) {
            stop("genotypes must be 0 (ancestral), 1 (derived) or -1/NA (missing)");
        }
        g[j] = static_cast<allele_t>(v);
    }
    tsir::check_tsi(ancestor_builder_add_site(ab, time, g.data()), "add_site");
}

// [[Rcpp::export(.ab_descriptors)]]
List
ab_descriptors(SEXP xp)
{
    XPtr<AncestorBuilderHandle> handle(xp);
    handle->finalise();
    const ancestor_builder_t *ab = &handle->builder;
    NumericVector time(ab->num_ancestors);
    List focal(ab->num_ancestors);
    for (size_t j = 0; j < ab->num_ancestors; j++) {
        const ancestor_descriptor_t *d = &ab->descriptors[j];
        time[j] = d->time;
        IntegerVector sites(d->focal_sites, d->focal_sites + d->num_focal_sites);
        focal[j] = sites;
    }
    return List::create(_["time"] = time, _["focal_sites"] = focal);
}

// [[Rcpp::export(.ab_make_ancestor)]]
List
ab_make_ancestor(SEXP xp, IntegerVector focal_sites)
{
    XPtr<AncestorBuilderHandle> handle(xp);
    handle->finalise();
    const ancestor_builder_t *ab = &handle->builder;
    size_t num_focal = static_cast<size_t>(focal_sites.size());
    if (num_focal == 0 || num_focal > ab->num_sites) {
        stop("num_focal_sites must > 0 and <= num_sites");
    }
    std::vector<tsk_id_t> focal(num_focal);
    for (size_t j = 0; j < num_focal; j++) {
        int s = focal_sites[j];
        if (s == NA_INTEGER || s < 0 || static_cast<size_t>(s) >= ab->num_sites) {
            stop("focal site out of bounds");
        }
        if (j > 0 && s <= focal[j - 1]) {
            stop("focal sites must be strictly increasing");
        }
        focal[j] = s;
    }
    std::vector<allele_t> hap(ab->num_sites);
    tsk_id_t start = 0, end = 0;
    tsir::check_tsi(ancestor_builder_make_ancestor(
                        ab, num_focal, focal.data(), &start, &end, hap.data()),
        "make_ancestor");
    IntegerVector haplotype(hap.begin(), hap.end());
    return List::create(_["start"] = start, _["end"] = end, _["haplotype"] = haplotype);
}

// [[Rcpp::export(.ab_info)]]
List
ab_info(SEXP xp)
{
    XPtr<AncestorBuilderHandle> handle(xp);
    const ancestor_builder_t *ab = &handle->builder;
    return List::create(_["num_sites"] = static_cast<double>(ab->num_sites),
        _["num_samples"] = static_cast<double>(ab->num_samples),
        _["num_ancestors"] = static_cast<double>(ab->num_ancestors),
        _["mem_size"] = static_cast<double>(ancestor_builder_get_memsize(ab)));
}

// [[Rcpp::export(.mi_new)]]
SEXP
mi_new(SEXP tables_xp, Nullable<NumericVector> num_alleles)
{
    tsk_table_collection_t *tables = tsir::get_tables(tables_xp);
    std::vector<tsk_size_t> na;
    const tsk_size_t *na_ptr = nullptr;
    if (num_alleles.isNotNull()) {
        NumericVector v(num_alleles);
        if (static_cast<tsk_size_t>(v.size()) != tables->sites.num_rows) {
            stop("num_alleles must have one entry per site");
        }
        na.resize(v.size());
        for (R_xlen_t j = 0; j < v.size(); j++) {
            if (!(v[j] >= 2 && v[j] <= 127)) {
                stop("The number of alleles must be between 2 and 127");
            }
            na[j] = static_cast<tsk_size_t>(v[j]);
        }
        na_ptr = na.data();
    }
    if (!tsk_table_collection_has_index(tables, 0)) {
        tsir::check_tsk(tsk_table_collection_build_index(tables, 0), "build_index");
    }
    auto *handle = new MatcherIndexesHandle(tables, na_ptr);
    XPtr<MatcherIndexesHandle> xp(handle, true);
    return xp;
}

// [[Rcpp::export(.am_new)]]
SEXP
am_new(SEXP mi_xp, NumericVector recombination, NumericVector mismatch,
    double likelihood_threshold, bool weight_by_n, bool extended_checks)
{
    XPtr<MatcherIndexesHandle> mi(mi_xp);
    size_t num_sites = mi->indexes.num_sites;
    if (static_cast<size_t>(recombination.size()) != num_sites) {
        stop("Size of recombination array must be num_sites");
    }
    if (static_cast<size_t>(mismatch.size()) != num_sites) {
        stop("Size of mismatch array must be num_sites");
    }
    int flags = 0;
    if (extended_checks) {
        flags |= TSI_EXTENDED_CHECKS;
    }
    if (!weight_by_n) {
        flags |= TSI_DISABLE_WEIGHT_BY_N;
    }
    if (ISNAN(likelihood_threshold)) {
        likelihood_threshold = DBL_MIN;
    }
    std::vector<double> rho(recombination.begin(), recombination.end());
    std::vector<double> mu(mismatch.begin(), mismatch.end());
    AncestorMatcherHandle *handle = nullptr;
    try {
        handle = new AncestorMatcherHandle(&mi->indexes, rho, mu, likelihood_threshold, flags);
    } catch (std::exception &e) {
        stop("AncestorMatcher: %s", e.what());
    }
    /* Keep the matcher indexes alive for as long as the matcher exists. */
    XPtr<AncestorMatcherHandle> xp(handle, true, R_NilValue, mi_xp);
    return xp;
}

// [[Rcpp::export(.am_find_path)]]
List
am_find_path(SEXP am_xp, IntegerVector haplotype, int start, int end)
{
    XPtr<AncestorMatcherHandle> handle(am_xp);
    ancestor_matcher_t *am = &handle->matcher;
    size_t num_sites = am->num_sites;
    if (static_cast<size_t>(haplotype.size()) != num_sites) {
        stop("Incorrect size for input haplotype.");
    }
    std::vector<allele_t> h(num_sites);
    for (size_t j = 0; j < num_sites; j++) {
        int v = haplotype[j];
        if (v == NA_INTEGER) {
            v = -1;
        }
        if (v < -1 || v > 127) {
            stop("haplotype values must be in [-1, 127]");
        }
        h[j] = static_cast<allele_t>(v);
    }
    std::vector<allele_t> match(num_sites, -1);
    std::vector<tsk_id_t> left(num_sites + 1), right(num_sites + 1), parent(num_sites + 1);
    size_t path_length = 0;
    tsir::check_tsi(ancestor_matcher_find_path(am, start, end, h.data(), match.data(),
                        &path_length, left.data(), right.data(), parent.data()),
        "find_path");
    IntegerVector l(left.begin(), left.begin() + path_length);
    IntegerVector r(right.begin(), right.begin() + path_length);
    IntegerVector p(parent.begin(), parent.begin() + path_length);
    IntegerVector m(match.begin(), match.end());
    return List::create(_["left"] = l, _["right"] = r, _["parent"] = p, _["match"] = m,
        _["mean_traceback_size"] = ancestor_matcher_get_mean_traceback_size(am));
}
