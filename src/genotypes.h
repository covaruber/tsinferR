#ifndef TSINFERR_GENOTYPES_H
#define TSINFERR_GENOTYPES_H

#include "common.h"
#include "zarr.h"

#include <functional>

namespace tsir {

/* One source of a multi-source genotype view. For each requested position k,
 * rows[k] is the source variant row (or -1 if the source has no data there)
 * and remap[k * max_alleles + a] maps source allele a to the unified allele
 * code. cols are the selected store sample columns (sorted). */
struct ViewSource {
    const ZarrArray *gt = nullptr;
    std::vector<int64_t> rows;
    std::vector<int8_t> remap;
    int max_alleles = 0;
    std::vector<int64_t> cols;
    int ploidy = 1;

    size_t num_haplotypes() const { return cols.size() * static_cast<size_t>(ploidy); }
};

struct GenotypeView {
    std::vector<ViewSource> sources;
    size_t num_positions = 0;

    size_t num_haplotypes() const;
};

GenotypeView view_from_list(Rcpp::List spec);

/* Iterate positions [begin, end) in order, calling f(k, genotypes) with the
 * concatenated unified allele codes (-1 = missing) of all haplotypes. Reads
 * each Zarr chunk at most once. Does not call R. */
void iterate_genotypes(const GenotypeView &view, size_t begin, size_t end,
    const std::function<void(size_t, const int8_t *)> &f);

} // namespace tsir

#endif
