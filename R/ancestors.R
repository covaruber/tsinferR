#' Infer ancestral haplotypes
#'
#' Builds ancestral haplotypes from the sample genotypes of one or more VCZ
#' sources and writes them to an ancestor VCZ store. This is the first stage
#' of the inference pipeline and is equivalent to `tsinfer.infer_ancestors()`.
#'
#' Inference sites are the positions present in the ancestral state store
#' that are not duplicated within a source, polymorphic and not entirely
#' missing. The time of each site is the frequency of the derived alleles.
#' Sites are split into intervals wherever consecutive sites are more than
#' `max_gap_length` apart, and ancestors are built independently for each
#' interval using the compiled 'tsinfer' ancestor builder.
#'
#' @param sources A `tsinfer_source()` or a list of them.
#' @param config An `ancestors_config()`. Its `path` is the output store,
#'   which must not already exist.
#' @param ancestral_state An `ancestral_state()`.
#' @param num_threads Number of threads used to build ancestors (0 or 1 =
#'   single-threaded).
#' @param genotype_encoding 0 to store genotypes with eight bits per call,
#'   1 for one bit per call (less memory; no missing data allowed).
#' @param verbose Print progress messages.
#' @return The path of the ancestor store, invisibly.
#' @export
infer_ancestors <- function(sources, config, ancestral_state, num_threads = 0L,
                            genotype_encoding = 0L, verbose = FALSE) {
  if (inherits(sources, "tsinfer_source")) {
    sources <- list(sources)
  }
  out_path <- config$path
  if (is.null(out_path)) {
    out_path <- tempfile("ancestors-", fileext = ".vcz")
  }
  out_path <- path.expand(out_path)
  if (file.exists(out_path)) {
    stop("Ancestor output path already exists: ", out_path, call. = FALSE)
  }
  view <- multi_source_view(sources, ancestral_state)
  num_haps <- view_num_haplotypes(view)
  if (verbose) {
    message(sprintf("Inference sites identified: %d (%d total haplotypes)", length(view$positions), num_haps))
  }

  # Pass 1: site statistics.
  vs <- view_spec(view, view$positions, apply_filters = TRUE)
  stats <- .view_site_stats(vs$spec)
  keep <- stats$non_missing > 0 & stats$derived_count > 0 &
    stats$derived_count != stats$non_missing
  times <- stats$derived_count[keep] / num_haps
  positions <- view$positions[keep]
  alleles <- vs$table$forward[keep, , drop = FALSE]
  if (verbose) {
    message(sprintf("Pass 1: %d sites kept, %d dropped (fixed or all-missing)", sum(keep), sum(!keep)))
  }

  first <- view$sources[[1]]$group
  contig_id <- as.character(zarr_read(first, "contig_id"))[1]
  contig_length <- zarr_read(first, "contig_length")[1]
  intervals <- compute_sequence_intervals(positions, config$max_gap_length)
  compressor <- zarr_compressor("blosc",
    cname = config$compressor, clevel = config$compression_level,
    shuffle = "bitshuffle"
  )
  tmp_path <- paste0(out_path, ".", Sys.getpid(), ".tmp")
  unlink(tmp_path, recursive = TRUE)
  on.exit(unlink(tmp_path, recursive = TRUE), add = TRUE)
  g <- zarr_create_group(tmp_path)
  vchunk <- config$variants_chunk_size
  schunk <- config$samples_chunk_size

  num_sites <- length(positions)
  out_alleles <- matrix("", num_sites, 2)
  if (num_sites > 0) {
    out_alleles[, 1] <- alleles[, 1]
    rest <- alleles[, -1, drop = FALSE]
    has <- matrix(nzchar(rest), nrow = nrow(rest))
    first_derived <- max.col(has, ties.method = "first")
    any_derived <- rowSums(has) > 0
    out_alleles[any_derived, 2] <- rest[cbind(which(any_derived), first_derived[any_derived])]
  }

  if (num_sites == 0) {
    res <- list(time = numeric(0), start = integer(0), end = integer(0), focal_positions = list())
    zarr_write_array(g, "call_genotype", array(integer(0), c(0, 0, 1)),
      dtype = "int8",
      compressor = compressor, fill_value = -1, dims = c("variants", "samples", "ploidy")
    )
  } else {
    gt_meta <- list(
      path = file.path(tmp_path, "call_genotype"), zarr_format = 2L,
      shape = c(num_sites, 0, 1), chunks = c(vchunk, schunk, 1), dtype = "int8",
      itemsize = 1, big_endian = FALSE, serializer = "raw",
      codecs = list(compressor_json(compressor, "int8")), key_prefix = "",
      separator = ".", fill_value = -1, fill_string = ""
    )
    dir.create(gt_meta$path)
    spec <- view_spec(view, positions, apply_filters = TRUE)
    lo <- findInterval(intervals[, 1], positions, left.open = TRUE)
    hi <- findInterval(intervals[, 2] - 1, positions)
    site_ranges <- cbind(lo, hi)
    storage.mode(site_ranges) <- "integer"
    res <- .build_ancestors(
      spec$spec, times, site_ranges, positions, gt_meta,
      as.integer(num_threads), as.integer(genotype_encoding), isTRUE(verbose)
    )
    n_anc <- length(res$time)
    write_array_metadata(g, "call_genotype",
      shape = c(num_sites, n_anc, 1), chunks = c(vchunk, schunk, 1),
      dtype = "int8", compressor = compressor, fill_value = -1,
      dims = c("variants", "samples", "ploidy")
    )
  }
  n_anc <- length(res$time)
  zarr_write_array(g, "variant_position", as.integer(positions),
    chunks = vchunk,
    compressor = compressor, dims = "variants"
  )
  zarr_write_array(g, "variant_allele", out_alleles,
    chunks = c(vchunk, 2),
    compressor = compressor, dims = c("variants", "alleles")
  )
  storage.mode(intervals) <- "integer"
  zarr_write_array(g, "sequence_intervals", intervals,
    compressor = compressor,
    dims = c("intervals", "coords")
  )
  zarr_write_array(g, "contig_id", contig_id, compressor = compressor, dims = "contigs")
  zarr_write_array(g, "contig_length", contig_length,
    dtype = "int64", compressor = compressor,
    dims = "contigs"
  )
  zarr_write_array(g, "variant_contig", integer(num_sites),
    dtype = "int8",
    chunks = vchunk, compressor = compressor, dims = "variants"
  )
  zarr_write_array(g, "sample_time", as.numeric(res$time),
    chunks = schunk,
    compressor = compressor, dims = "samples"
  )
  zarr_write_array(g, "sample_start_position", as.integer(res$start),
    chunks = schunk,
    compressor = compressor, dims = "samples"
  )
  zarr_write_array(g, "sample_end_position", as.integer(res$end),
    chunks = schunk,
    compressor = compressor, dims = "samples"
  )
  zarr_write_array(g, "sample_id", if (n_anc > 0) paste0("a", seq_len(n_anc) - 1L) else character(0),
    compressor = compressor, dims = "samples"
  )
  max_focal <- max(1L, lengths(res$focal_positions))
  focal <- matrix(-2L, n_anc, max_focal)
  for (j in seq_len(n_anc)) {
    fp <- res$focal_positions[[j]]
    focal[j, seq_along(fp)] <- fp
  }
  zarr_write_array(g, "sample_focal_positions", focal,
    compressor = compressor,
    dims = c("samples", "focal_alleles")
  )
  if (!file.rename(tmp_path, out_path)) {
    stop("Could not move the ancestor store to ", out_path, call. = FALSE)
  }
  if (verbose) {
    message(sprintf("Ancestor inference complete: %d ancestors across %d sites", n_anc, num_sites))
  }
  invisible(out_path)
}

compute_sequence_intervals <- function(positions, max_gap_length) {
  n <- length(positions)
  if (n == 0) {
    return(matrix(0L, 0, 2))
  }
  breaks <- which(diff(positions) > max_gap_length)
  starts <- positions[c(1, breaks + 1)]
  ends <- positions[c(breaks, n)] + 1
  cbind(starts, ends, deparse.level = 0)
}

write_array_metadata <- function(group, name, shape, chunks, dtype, compressor,
                                 fill_value, dims) {
  dir <- file.path(group$path, name)
  meta <- list(
    chunks = as.list(chunks), compressor = compressor_json(compressor, dtype),
    dimension_separator = ".", dtype = dtype_v2_string(dtype),
    fill_value = fill_value, filters = NULL, order = "C",
    shape = as.list(shape), zarr_format = 2L
  )
  write_json_file(meta, file.path(dir, ".zarray"))
  write_json_file(list("_ARRAY_DIMENSIONS" = as.list(dims)), file.path(dir, ".zattrs"))
}
