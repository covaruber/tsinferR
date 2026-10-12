#' Tree sequences
#'
#' Inferred tree sequences are returned as `tsinfer_ts` objects, which wrap a
#' 'tskit' table collection held in compiled code. They can be saved to and
#' loaded from the standard `.trees` file format, so they interoperate with
#' 'tskit' (Python) and with R packages that read `.trees` files.
#'
#' Pointers do not survive saving an R session: use [ts_dump()] and
#' [ts_load()] to persist tree sequences.
#'
#' @param ts A `tsinfer_ts` object.
#' @param path Path to a `.trees` file.
#' @return `ts_load()` returns a `tsinfer_ts`; `ts_dump()` returns `path`
#'   invisibly.
#' @name tree-sequences
#' @examples
#' ts <- ts_from_tables(
#'   nodes = data.frame(flags = c(0, 1, 1), time = c(1, 0, 0)),
#'   edges = data.frame(left = 0, right = 10, parent = 0, child = 1:2),
#'   sequence_length = 10
#' )
#' path <- tempfile(fileext = ".trees")
#' ts_dump(ts, path)
#' ts2 <- ts_load(path)
#' ts_equals(ts, ts2)
NULL

new_ts <- function(ptr) {
  structure(list(ptr = ptr), class = "tsinfer_ts")
}

ts_ptr <- function(ts) {
  if (!inherits(ts, "tsinfer_ts")) {
    stop("expected a tsinfer_ts object", call. = FALSE)
  }
  ts$ptr
}

#' @rdname tree-sequences
#' @export
ts_load <- function(path) {
  path <- normalizePath(path, mustWork = TRUE)
  new_ts(.tc_load(path))
}

#' @rdname tree-sequences
#' @export
ts_dump <- function(ts, path) {
  .tc_dump(ts_ptr(ts), path.expand(path))
  invisible(path)
}

#' @rdname tree-sequences
#' @export
ts_copy <- function(ts) {
  new_ts(.tc_copy(ts_ptr(ts)))
}

#' @rdname tree-sequences
#' @param other Another `tsinfer_ts`.
#' @param ignore_provenance If `TRUE`, provenance tables are not compared.
#' @export
ts_equals <- function(ts, other, ignore_provenance = TRUE) {
  .tc_equals(ts_ptr(ts), ts_ptr(other), ignore_provenance, TRUE)
}

#' Summary information about a tree sequence
#'
#' @param ts A `tsinfer_ts` object.
#' @return `ts_info()` returns a named list with the sequence length and the
#'   number of rows in each table. `ts_num_trees()` returns the number of
#'   trees.
#' @export
ts_info <- function(ts) {
  info <- .tc_info(ts_ptr(ts))
  info[c(
    "sequence_length", "num_nodes", "num_edges", "num_sites",
    "num_mutations", "num_individuals", "num_populations", "num_provenances"
  )]
}

#' @rdname ts_info
#' @export
ts_num_trees <- function(ts) {
  .tc_num_trees(ts_ptr(ts))
}

#' @export
print.tsinfer_ts <- function(x, ...) {
  info <- ts_info(x)
  cat("<tsinfer_ts>\n")
  cat(sprintf("  sequence_length: %s\n", format(info$sequence_length)))
  for (name in setdiff(names(info), "sequence_length")) {
    cat(sprintf("  %s: %s\n", name, format(info[[name]])))
  }
  invisible(x)
}

decode_json_column <- function(x) {
  lapply(x, function(s) if (nzchar(s)) jsonlite::fromJSON(s, simplifyVector = FALSE) else NULL)
}

#' Tables of a tree sequence as data frames
#'
#' @param ts A `tsinfer_ts` object.
#' @param decode_metadata If `TRUE`, JSON metadata are decoded into a list
#'   column; otherwise the raw JSON strings are returned.
#' @return A data frame with one row per table row. Node, individual, site
#'   and population IDs are 0-based, as in 'tskit'.
#' @name ts-tables
NULL

table_df <- function(cols, decode_metadata) {
  cols$metadata_schema <- NULL
  if (!is.null(cols$metadata) && decode_metadata) {
    md <- decode_json_column(cols$metadata)
    cols$metadata <- NULL
    df <- if (length(cols) > 0) {
      as.data.frame(cols, stringsAsFactors = FALSE)
    } else {
      data.frame(row.names = seq_along(md))
    }
    df$metadata <- md
    return(df)
  }
  as.data.frame(cols, stringsAsFactors = FALSE)
}

#' @rdname ts-tables
#' @export
ts_nodes <- function(ts, decode_metadata = TRUE) {
  table_df(.tc_table(ts_ptr(ts), "nodes"), decode_metadata)
}

#' @rdname ts-tables
#' @export
ts_edges <- function(ts) {
  table_df(.tc_table(ts_ptr(ts), "edges"), FALSE)
}

#' @rdname ts-tables
#' @export
ts_sites <- function(ts, decode_metadata = TRUE) {
  table_df(.tc_table(ts_ptr(ts), "sites"), decode_metadata)
}

#' @rdname ts-tables
#' @export
ts_mutations <- function(ts, decode_metadata = TRUE) {
  table_df(.tc_table(ts_ptr(ts), "mutations"), decode_metadata)
}

#' @rdname ts-tables
#' @export
ts_individuals <- function(ts, decode_metadata = TRUE) {
  table_df(.tc_table(ts_ptr(ts), "individuals"), decode_metadata)
}

#' @rdname ts-tables
#' @export
ts_populations <- function(ts, decode_metadata = TRUE) {
  table_df(.tc_table(ts_ptr(ts), "populations"), decode_metadata)
}

#' @rdname ts-tables
#' @export
ts_provenances <- function(ts) {
  table_df(.tc_table(ts_ptr(ts), "provenances"), FALSE)
}

#' Top-level metadata of a tree sequence
#'
#' @param ts A `tsinfer_ts` object.
#' @return The decoded JSON metadata (a list), or `NULL` if there is none.
#' @export
ts_metadata <- function(ts) {
  md <- .tc_info(ts_ptr(ts))$metadata
  if (!nzchar(md)) {
    return(NULL)
  }
  jsonlite::fromJSON(md, simplifyVector = TRUE)
}

#' Genotype matrix of a tree sequence
#'
#' @param ts A `tsinfer_ts` object.
#' @return An integer matrix with one row per site and one column per sample
#'   node, holding allele indices (0 = ancestral state, -1 = missing). The
#'   allele strings for each site are attached as the `alleles` attribute and
#'   the sample node IDs as the `samples` attribute.
#' @export
ts_genotype_matrix <- function(ts) {
  res <- .tc_genotype_matrix(ts_ptr(ts))
  G <- res$genotypes
  attr(G, "alleles") <- res$alleles
  attr(G, "samples") <- res$samples
  G
}

#' Simplify a tree sequence
#'
#' @param ts A `tsinfer_ts` object.
#' @param samples Sample node IDs (0-based) to keep. Defaults to all nodes
#'   flagged as samples.
#' @param filter_sites,filter_populations,filter_individuals Remove unreferenced
#'   rows from the corresponding tables.
#' @param keep_unary Keep unary nodes.
#' @return A new, simplified `tsinfer_ts`.
#' @export
ts_simplify <- function(ts, samples = NULL, filter_sites = TRUE,
                        filter_populations = TRUE, filter_individuals = TRUE,
                        keep_unary = FALSE) {
  if (!is.null(samples)) {
    samples <- as.integer(samples)
  }
  new_ts(.tc_simplify(
    ts_ptr(ts), samples, filter_sites, filter_populations,
    filter_individuals, keep_unary
  ))
}

#' Build a tree sequence from data frames
#'
#' Mainly useful for testing and for constructing reference panels by hand.
#'
#' @param nodes Data frame with columns `flags` and `time`, and optionally
#'   `population`, `individual` and `metadata` (JSON strings).
#' @param edges Data frame with columns `left`, `right`, `parent`, `child`.
#' @param sites Optional data frame with columns `position` and
#'   `ancestral_state`.
#' @param mutations Optional data frame with columns `site`, `node`,
#'   `derived_state` and optionally `parent` and `time`.
#' @param sequence_length The sequence length.
#' @param sort If `TRUE`, the tables are sorted and indexed.
#' @return A `tsinfer_ts`.
#' @export
ts_from_tables <- function(nodes, edges, sites = NULL, mutations = NULL,
                           sequence_length, sort = TRUE) {
  ptr <- .tc_new(as.numeric(sequence_length))
  n <- nrow(nodes)
  .tc_set_table(ptr, "nodes", list(
    flags = as.numeric(nodes$flags),
    time = as.numeric(nodes$time),
    population = as.integer(if (is.null(nodes$population)) rep(-1L, n) else nodes$population),
    individual = as.integer(if (is.null(nodes$individual)) rep(-1L, n) else nodes$individual),
    metadata = as.character(if (is.null(nodes$metadata)) rep("", n) else nodes$metadata)
  ))
  .tc_set_table(ptr, "edges", list(
    left = as.numeric(edges$left), right = as.numeric(edges$right),
    parent = as.integer(edges$parent), child = as.integer(edges$child)
  ))
  if (!is.null(sites)) {
    .tc_set_table(ptr, "sites", list(
      position = as.numeric(sites$position),
      ancestral_state = as.character(sites$ancestral_state)
    ))
  }
  if (!is.null(mutations)) {
    m <- nrow(mutations)
    .tc_set_table(ptr, "mutations", list(
      site = as.integer(mutations$site),
      node = as.integer(mutations$node),
      parent = as.integer(if (is.null(mutations$parent)) rep(-1L, m) else mutations$parent),
      time = as.numeric(if (is.null(mutations$time)) rep(NaN, m) else mutations$time),
      derived_state = as.character(mutations$derived_state)
    ))
  }
  if (sort) {
    .tc_finalise(ptr, is.null(mutations$parent))
  }
  new_ts(ptr)
}
