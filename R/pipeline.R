#' Post-process an inferred tree sequence
#'
#' Cleans up the raw tree sequence produced by [match_haplotypes()]. Edges are
#' clipped to the inference intervals stored in the top-level metadata
#' (`erase_flanks`) and the tree sequence is simplified. Equivalent to
#' `tsinfer.post_process()`.
#'
#' @param ts A `tsinfer_ts` from [match_haplotypes()].
#' @param config A `tsinfer_config`. If it has no `post_process` section, `ts`
#'   is returned unchanged.
#' @param erase_flanks,split_ultimate Override the corresponding
#'   `post_process_config()` settings. Splitting the ultimate ancestor is not
#'   implemented yet (as in 'tsinfer'), so `split_ultimate` has no effect.
#' @return A `tsinfer_ts`.
#' @export
post_process <- function(ts, config, erase_flanks = NULL, split_ultimate = NULL) {
  pp <- config$post_process
  if (is.null(pp)) {
    return(ts)
  }
  if (is.null(erase_flanks)) {
    erase_flanks <- pp$erase_flanks
  }
  if (isTRUE(erase_flanks)) {
    md <- ts_metadata(ts)
    intervals <- md$sequence_intervals
    if (!is.null(intervals)) {
      intervals <- matrix(as.numeric(intervals), ncol = 2)
      ts <- new_ts(.erase_flanks(ts$ptr, intervals))
    }
  }
  simplified <- ts_simplify(ts)
  .tc_add_provenance(simplified$ptr, provenance_record("post_process", list()))
  simplified
}

#' Add sites to an inferred tree sequence by parsimony
#'
#' Places the variants of the sources listed in `config$augment_sites` that
#' are not already in `ts` (and that fall within its inference intervals)
#' onto the inferred trees, using Hartigan parsimony with the ancestral
#' allele fixed to that given by `config$ancestral_state`. Equivalent to
#' `tsinfer.augment_sites()`.
#'
#' Samples are matched to the sources through the `source`, `sample_id` and
#' `ploidy_index` keys of the sample node metadata written by
#' [match_haplotypes()].
#'
#' @param ts A `tsinfer_ts`.
#' @param config A `tsinfer_config`.
#' @param verbose Print progress messages.
#' @return A `tsinfer_ts` with the additional sites and mutations.
#' @export
augment_sites <- function(ts, config, verbose = FALSE) {
  aug <- config$augment_sites
  if (is.null(aug) || length(aug$sources) == 0) {
    return(ts)
  }
  md <- ts_metadata(ts)
  intervals <- md$sequence_intervals
  info <- ts_info(ts)
  intervals <- if (is.null(intervals)) {
    matrix(c(0, info$sequence_length), ncol = 2)
  } else {
    matrix(as.numeric(intervals), ncol = 2)
  }
  nodes <- ts_nodes(ts)
  samples <- which(bitwAnd(as.integer(nodes$flags %% 2^31), 1L) != 0)
  ident <- lapply(samples, function(u) {
    m <- nodes$metadata[[u]]
    if (is.null(m) || !all(c("source", "sample_id", "ploidy_index") %in% names(m))) {
      stop(sprintf(
        "Node %d is missing required metadata key(s) for augment_sites",
        u - 1L
      ), call. = FALSE)
    }
    m
  })
  sources <- lapply(aug$sources, function(name) {
    s <- config$sources[[name]]
    if (is.null(s)) stop("[augment_sites] references unknown source: '", name, "'", call. = FALSE)
    s
  })
  view <- multi_source_view(sources, config$ancestral_state)
  pos <- view$positions
  idx <- findInterval(pos, intervals[, 1])
  in_iv <- idx > 0 & pos < intervals[pmax(idx, 1), 2]
  pos <- pos[in_iv]
  pos <- pos[!pos %in% ts_sites(ts, FALSE)$position]
  if (length(pos) == 0) {
    if (verbose) message("augment_sites: no new sites to add")
    return(ts)
  }
  if (verbose) message(sprintf("augment_sites: %d new sites to place", length(pos)))
  vs <- view_spec(view, pos, apply_filters = TRUE)
  haps <- view_haplotypes(view)
  key_full <- paste(haps$source, haps$sample_id, haps$ploidy_index, sep = "\r")
  key_short <- paste(haps$sample_id, haps$ploidy_index, sep = "\r")
  src <- vapply(ident, function(m) as.character(m$source), character(1))
  sid <- vapply(ident, function(m) as.character(m$sample_id), character(1))
  pid <- vapply(ident, function(m) as.integer(m$ploidy_index), integer(1))
  perm <- match(paste(src, sid, pid, sep = "\r"), key_full)
  fallback <- is.na(perm)
  perm[fallback] <- match(paste(sid, pid, sep = "\r")[fallback], key_short)
  perm[is.na(perm)] <- 0L
  perm <- as.integer(perm - 1L)
  alleles <- lapply(seq_len(nrow(vs$table$forward)), function(i) {
    a <- vs$table$forward[i, ]
    n <- length(a)
    while (n > 0 && !nzchar(a[n])) n <- n - 1
    a[seq_len(n)]
  })
  res <- .augment_sites(ts$ptr, vs$spec, pos, alleles, perm)
  out <- new_ts(res$tables)
  if (verbose) {
    message(sprintf(
      "augment_sites complete: %d sites total (%d added)",
      ts_info(out)$num_sites, res$num_added
    ))
  }
  out
}

#' Run the full inference pipeline
#'
#' Runs [infer_ancestors()], [match_haplotypes()], [post_process()] and, if
#' configured, [augment_sites()]. Equivalent to `tsinfer.run()`.
#'
#' @param config A `tsinfer_config`, or the path of a TOML configuration file.
#' @param num_threads Number of threads.
#' @param cache_size Genotype chunk cache size in MiB used when matching.
#' @param genotype_encoding Ancestor builder genotype encoding (0 or 1).
#' @param match_file Optional JSON-lines file of copying paths.
#' @param output If not `NULL`, the final tree sequence is also written to
#'   this `.trees` file. Use `TRUE` to write to `config$match$output`.
#' @param verbose Print progress messages.
#' @param ... Passed to [match_haplotypes()] and [post_process()].
#' @return The inferred tree sequence (`tsinfer_ts`).
#' @export
run_tsinfer <- function(config, num_threads = 0L, cache_size = 256,
                        genotype_encoding = 0L, match_file = NULL, output = NULL,
                        verbose = FALSE, ...) {
  if (is.character(config)) {
    config <- read_config(config)
  }
  extra <- list(...)
  anc_cfg <- config$ancestors[[1]]
  if (is.null(anc_cfg$path)) {
    anc_cfg$path <- tempfile("ancestors-", fileext = ".vcz")
    on.exit(unlink(anc_cfg$path, recursive = TRUE), add = TRUE)
    config$ancestors[[1]] <- anc_cfg
    config$sources[[anc_cfg$name]]$path <- anc_cfg$path
  }
  sources <- lapply(anc_cfg$sources, function(name) config$sources[[name]])
  infer_ancestors(sources, anc_cfg, config$ancestral_state,
    num_threads = num_threads, genotype_encoding = genotype_encoding,
    verbose = verbose
  )
  ts <- match_haplotypes(config,
    num_threads = num_threads, cache_size = cache_size,
    match_file = match_file, path_compression = extra$path_compression,
    group_stop = extra$group_stop, verbose = verbose
  )
  ts <- post_process(ts, config,
    erase_flanks = extra$erase_flanks,
    split_ultimate = extra$split_ultimate
  )
  if (!is.null(config$augment_sites)) {
    ts <- augment_sites(ts, config, verbose = verbose)
  }
  if (isTRUE(output)) {
    output <- config$match$output
  }
  if (is.character(output)) {
    ts_dump(ts, output)
  }
  if (verbose) {
    info <- ts_info(ts)
    message(sprintf(
      "Pipeline complete: %d nodes, %d edges, %d sites",
      info$num_nodes, info$num_edges, info$num_sites
    ))
  }
  ts
}
