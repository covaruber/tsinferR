#' Exchange tree sequences with companion R packages
#'
#' These helpers bridge through the standard `.trees` file format, so the
#' companion packages retain ownership of their own C pointers and table APIs.
#' `tskitR` and `tsdateR` are optional dependencies.
#'
#' * `as_tskit()` converts a `tsinfer_ts` to a `tskitR::TreeSequence`.
#' * `as_tsinfer()` converts a `tskitR::TreeSequence` to `tsinfer_ts`.
#' * `date_tsinfer()` sends a `tsinfer_ts` through `tsdateR::date_ts()` and
#'   returns the dated result as a new `tsinfer_ts`.
#'
#' @param ts A `tsinfer_ts` object.
#' @param tree_sequence A `tskitR::TreeSequence` object.
#' @param mutation_rate Mutation rate passed to `tsdateR::date_ts()`.
#' @param method Dating method passed to `tsdateR::date_ts()`.
#' @param progress Show the tsdateR progress indicator.
#' @param ... Additional arguments passed to `tsdateR::date_ts()`.
#' @return `as_tskit()` returns a `tskitR::TreeSequence`; `as_tsinfer()` and
#'   `date_tsinfer()` return `tsinfer_ts` objects.
#' @name r-package-interop
#' @examples
#' if (requireNamespace("tskitR", quietly = TRUE)) {
#'   ts <- ts_load(system.file("extdata", "example.trees", package = "tskitR"))
#'   tskit_ts <- as_tskit(ts)
#'   ts2 <- as_tsinfer(tskit_ts)
#'   ts_equals(ts, ts2)
#' }
NULL

#' @rdname r-package-interop
#' @export
as_tskit <- function(ts) {
  if (!requireNamespace("tskitR", quietly = TRUE)) {
    stop("Install the optional 'tskitR' package to convert tree sequences",
      call. = FALSE
    )
  }
  input <- tempfile(fileext = ".trees")
  on.exit(unlink(input), add = TRUE)
  ts_dump(ts, input)
  getExportedValue("tskitR", "tskit_load")(input)
}

#' @rdname r-package-interop
#' @export
as_tsinfer <- function(tree_sequence) {
  if (!inherits(tree_sequence, "TreeSequence")) {
    stop("'tree_sequence' must be a tskitR TreeSequence", call. = FALSE)
  }
  input <- tempfile(fileext = ".trees")
  on.exit(unlink(input), add = TRUE)
  tree_sequence$dump(input)
  ts_load(input)
}

#' @rdname r-package-interop
#' @export
date_tsinfer <- function(ts, mutation_rate, method = NULL, progress = FALSE, ...) {
  if (!requireNamespace("tsdateR", quietly = TRUE)) {
    stop("Install the optional 'tsdateR' package to date tree sequences",
      call. = FALSE
    )
  }
  input <- tempfile(fileext = ".trees")
  output <- tempfile(fileext = ".trees")
  on.exit(unlink(c(input, output)), add = TRUE)
  ts_dump(ts, input)
  dated <- tsdateR::ts_load(input)
  args <- list(ts = dated, mutation_rate = mutation_rate, progress = progress, ...)
  if (!is.null(method)) {
    args$method <- method
  }
  dated <- do.call(getExportedValue("tsdateR", "date_ts"), args)
  if (!inherits(dated, "tsdate_treeseq")) {
    stop("tsdateR::date_ts() returned a non-tree-sequence result; omit fit/likelihood outputs",
      call. = FALSE
    )
  }
  getExportedValue("tsdateR", "ts_dump")(dated, output)
  ts_load(output)
}
