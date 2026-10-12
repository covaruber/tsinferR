#ifndef TSINFERR_LOWLEVEL_H
#define TSINFERR_LOWLEVEL_H

#include "common.h"

namespace tsir {

struct AncestorBuilderHandle {
    ancestor_builder_t builder;
    bool finalised = false;

    AncestorBuilderHandle(size_t num_samples, size_t max_sites, int genotype_encoding);
    ~AncestorBuilderHandle();
    AncestorBuilderHandle(const AncestorBuilderHandle &) = delete;
    AncestorBuilderHandle &operator=(const AncestorBuilderHandle &) = delete;
    /* ancestor_builder_finalise() must run exactly once. */
    void finalise();
};

struct MatcherIndexesHandle {
    matcher_indexes_t indexes;

    MatcherIndexesHandle(const tsk_table_collection_t *tables, const tsk_size_t *num_alleles);
    ~MatcherIndexesHandle();
    MatcherIndexesHandle(const MatcherIndexesHandle &) = delete;
    MatcherIndexesHandle &operator=(const MatcherIndexesHandle &) = delete;
};

/* Safe to construct from worker threads: throws std::runtime_error. */
struct AncestorMatcherHandle {
    ancestor_matcher_t matcher;

    AncestorMatcherHandle(const matcher_indexes_t *indexes,
        const std::vector<double> &recombination, const std::vector<double> &mismatch,
        double likelihood_threshold, int flags);
    ~AncestorMatcherHandle();
    AncestorMatcherHandle(const AncestorMatcherHandle &) = delete;
    AncestorMatcherHandle &operator=(const AncestorMatcherHandle &) = delete;
};

} // namespace tsir

#endif
