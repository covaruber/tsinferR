#' tsinferR: infer tree sequences from genetic variation data
#'
#' A native R implementation of the 'tsinfer' pipeline. The package bundles
#' the 'tsinfer' and 'tskit' C libraries, so no 'Python' installation is
#' needed. It is an R mirror of the canonical tsinfer method; results produced
#' with this package should cite the tsinfer paper returned by
#' \code{citation("tsinferR")}.
#'
#' The inference pipeline has three stages:
#' \enumerate{
#'   \item [infer_ancestors()] builds ancestral haplotypes from sample
#'     genotypes and writes them to an ancestor VCZ store.
#'   \item [match_haplotypes()] matches ancestors against each other and then samples
#'     against the ancestors, building a tree sequence.
#'   \item [post_process()] cleans up the raw inferred tree sequence.
#' }
#' [run_tsinfer()] chains all stages, and [augment_sites()] places extra sites on an
#' inferred tree sequence by parsimony.
#'
#' Input data are VCF Zarr ('VCZ') stores. These can be read and written
#' with the built-in Zarr engine ([zarr_open()], [zarr_read()]), and VCF files
#' can be converted with [vcf_to_vcz()].
#'
#' @keywords internal
#' @useDynLib tsinferR, .registration = TRUE
#' @importFrom Rcpp sourceCpp
"_PACKAGE"
