#ifndef TSINFERR_GROUPING_H
#define TSINFERR_GROUPING_H

#include <cstdint>
#include <vector>

namespace tsir {

/* Group haplotypes for parallel matching (tsinfer's linesweep + topological
 * sort). Returns the 0-based group of each haplotype; throws
 * std::runtime_error if the dependency graph has a cycle. */
std::vector<int> group_by_linesweep(const std::vector<int64_t> &start,
    const std::vector<int64_t> &end, const std::vector<double> &time);

} // namespace tsir

#endif
