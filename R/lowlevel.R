#' Low-level ancestor building and matching
#'
#' Direct access to the compiled 'tsinfer' engines. These functions mirror the
#' `_tsinfer` Python extension module and are mainly intended for testing and
#' for users who need fine control. Site and node indices are **0-based**, as
#' in the C library.
#'
#' * `ancestor_builder()` creates an ancestor builder for `num_samples`
#'   haplotypes and at most `max_sites` sites. Sites are added in genomic order
#'   with `ancestor_builder_add_site()`, using genotypes coded as 0
#'   (ancestral), 1 (derived) or -1 (missing).
#' * `ancestor_builder_descriptors()` returns the ancestors implied by the
#'   added sites as a list with elements `time` and `focal_sites`.
#' * `ancestor_builder_make_ancestor()` builds the haplotype of an ancestor
#'   with the given focal sites. It returns `start` and `end` (half-open site
#'   interval) and the `haplotype` (-1 outside the interval).
#' * `matcher_indexes()` prepares a tree sequence for matching. Node 0 must be
#'   the root of every tree.
#' * `ancestor_matcher()` creates a Li and Stephens matcher from per-site
#'   recombination and mismatch probabilities.
#' * `find_path()` finds the most likely copying path of a haplotype over the
#'   sites `[start, end)`. It returns the path as parallel vectors `left`,
#'   `right` (genomic coordinates) and `parent` (node IDs), plus the `match`
#'   haplotype.
#'
#' @param num_samples Number of sample haplotypes.
#' @param max_sites Maximum number of sites that will be added.
#' @param genotype_encoding 0 for eight-bit storage, 1 for one-bit storage
#'   (only valid without missing data).
#' @param builder An object returned by `ancestor_builder()`.
#' @param time The time (derived allele frequency) of the site.
#' @param genotypes Integer genotypes for one site.
#' @param focal_sites Integer vector of focal site indices.
#' @param ts A `tsinfer_ts` reference tree sequence.
#' @param num_alleles Optional number of alleles per site.
#' @param indexes An object returned by `matcher_indexes()`.
#' @param recombination,mismatch Per-site probabilities.
#' @param likelihood_threshold Likelihood floor (defaults to the smallest
#'   positive double).
#' @param weight_by_n Weight the recombination probability by the number of
#'   nodes, as in 'tsinfer'.
#' @param extended_checks Run expensive internal consistency checks.
#' @param matcher An object returned by `ancestor_matcher()`.
#' @param haplotype Integer haplotype (allele codes, -1 for missing).
#' @param start,end 0-based half-open site range to match over.
#' @return See Details.
#' @name lowlevel
#' @examples
#' ab <- ancestor_builder(num_samples = 4, max_sites = 3)
#' ancestor_builder_add_site(ab, 0.5, c(1, 1, 0, 0))
#' ancestor_builder_add_site(ab, 0.25, c(1, 0, 0, 0))
#' ancestor_builder_add_site(ab, 0.5, c(1, 1, 0, 0))
#' d <- ancestor_builder_descriptors(ab)
#' ancestor_builder_make_ancestor(ab, d$focal_sites[[1]])
NULL

#' @rdname lowlevel
#' @export
ancestor_builder <- function(num_samples, max_sites, genotype_encoding = 0L) {
  structure(
    list(ptr = .ab_new(as.integer(num_samples), as.integer(max_sites), as.integer(genotype_encoding))),
    class = "tsinfer_ancestor_builder"
  )
}

#' @rdname lowlevel
#' @export
ancestor_builder_add_site <- function(builder, time, genotypes) {
  .ab_add_site(builder$ptr, as.numeric(time), as.integer(genotypes))
  invisible(builder)
}

#' @rdname lowlevel
#' @export
ancestor_builder_descriptors <- function(builder) {
  .ab_descriptors(builder$ptr)
}

#' @rdname lowlevel
#' @export
ancestor_builder_make_ancestor <- function(builder, focal_sites) {
  .ab_make_ancestor(builder$ptr, as.integer(focal_sites))
}

#' @export
print.tsinfer_ancestor_builder <- function(x, ...) {
  info <- .ab_info(x$ptr)
  cat(sprintf(
    "<tsinfer_ancestor_builder> %d samples, %d sites, %d ancestors\n",
    as.integer(info$num_samples), as.integer(info$num_sites), as.integer(info$num_ancestors)
  ))
  invisible(x)
}

#' @rdname lowlevel
#' @export
matcher_indexes <- function(ts, num_alleles = NULL) {
  if (!is.null(num_alleles)) {
    num_alleles <- as.numeric(num_alleles)
  }
  structure(list(ptr = .mi_new(ts_ptr(ts), num_alleles)), class = "tsinfer_matcher_indexes")
}

#' @rdname lowlevel
#' @export
ancestor_matcher <- function(indexes, recombination, mismatch,
                             likelihood_threshold = NULL, weight_by_n = TRUE,
                             extended_checks = FALSE) {
  if (is.null(likelihood_threshold)) {
    likelihood_threshold <- NaN
  }
  structure(
    list(ptr = .am_new(
      indexes$ptr, as.numeric(recombination), as.numeric(mismatch),
      as.numeric(likelihood_threshold), isTRUE(weight_by_n), isTRUE(extended_checks)
    )),
    class = "tsinfer_ancestor_matcher"
  )
}

#' @rdname lowlevel
#' @export
find_path <- function(matcher, haplotype, start, end) {
  .am_find_path(matcher$ptr, as.integer(haplotype), as.integer(start), as.integer(end))
}
