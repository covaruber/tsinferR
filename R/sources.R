# Opening VCZ sources and building the unified multi-source views used by
# ancestor inference, matching and site augmentation.

# Read the per-source metadata needed by all pipeline steps.
open_source <- function(src) {
  g <- zarr_open(src$path)
  info <- zarr_info(g, "call_genotype")
  if (length(info$shape) != 3) {
    stop("Source '", src$name, "': call_genotype must have shape (variants, samples, ploidy)",
      call. = FALSE
    )
  }
  positions <- as.numeric(zarr_read(g, "variant_position"))
  if (is.unsorted(positions)) {
    stop("Source '", src$name, "': variant_position must be sorted", call. = FALSE)
  }
  alleles <- zarr_read(g, "variant_allele")
  if (is.null(dim(alleles))) {
    alleles <- matrix(alleles, ncol = 1)
  }
  alleles[is.na(alleles)] <- ""
  sample_ids <- as.character(zarr_read(g, "sample_id"))
  selection <- parse_samples_spec(src$samples, sample_ids)
  list(
    config = src, name = src$name, group = g, path = g$path,
    positions = positions, alleles = alleles,
    sample_ids = sample_ids, selection = selection,
    num_store_samples = info$shape[2], ploidy = as.integer(info$shape[3]),
    gt_meta = zarr_meta(g, "call_genotype")$meta,
    sample_chunk_size = info$chunks[2]
  )
}

# Resolve a sample-time (or other field) specification against a store.
resolve_field <- function(group, spec, join_key, n, fill_value = NULL) {
  if (is.null(spec)) {
    return(NULL)
  }
  if (is.character(spec) && length(spec) == 1) {
    return(zarr_read(group, spec))
  }
  if (is.list(spec)) {
    if (is.null(spec$path) || is.null(spec$field)) {
      stop("A field specification list must contain both 'path' and 'field'", call. = FALSE)
    }
    ann <- zarr_open(spec$path)
    keys <- zarr_read(group, join_key)
    ann_keys <- zarr_read(ann, join_key)
    ann_values <- zarr_read(ann, spec$field)
    idx <- match(keys, ann_keys)
    out <- ann_values[idx]
    if (!is.null(fill_value)) {
      out[is.na(idx)] <- fill_value
    }
    if (length(out) != n) {
      stop(sprintf("resolve_field: expected %d entries after join on '%s', got %d", n, join_key, length(out)),
        call. = FALSE
      )
    }
    return(out)
  }
  if (is.numeric(spec) && length(spec) == 1) {
    return(rep(spec, n))
  }
  stop("Invalid field specification", call. = FALSE)
}

duplicate_positions <- function(positions) {
  unique(positions[duplicated(positions)])
}

# Allele table: code 0 is the ancestral allele, then new alleles in order of
# first appearance across sources (source order, then variant row, then
# allele column), as in tsinfer's AlleleMapper.
build_allele_table <- function(positions, ancestral, per_source) {
  n <- length(positions)
  site <- seq_len(n)
  allele <- as.character(ancestral)
  for (ps in per_source) {
    rows <- which(ps$positions %in% positions)
    if (length(rows) == 0) {
      next
    }
    idx <- match(ps$positions[rows], positions)
    a <- ps$alleles[rows, , drop = FALSE]
    # Row-major order: all alleles of a row before the next row.
    site <- c(site, rep(idx, each = ncol(a)))
    allele <- c(allele, as.vector(t(a)))
  }
  keep <- !is.na(allele) & nzchar(allele)
  keep[seq_len(n)] <- TRUE
  site <- site[keep]
  allele <- allele[keep]
  first <- !duplicated(paste(site, allele, sep = "\r"))
  site <- site[first]
  allele <- allele[first]
  ord <- order(site, method = "radix")
  site <- site[ord]
  allele <- allele[ord]
  code <- sequence(tabulate(site, nbins = n)) - 1L
  max_a <- max(2L, if (length(code)) max(code) + 1L else 0L)
  forward <- matrix("", n, max_a)
  forward[cbind(site, code + 1L)] <- allele
  list(forward = forward, num_alleles = pmax(tabulate(site, nbins = n), 1L))
}

# Code of each (site, allele) pair in the allele table (-1 if unknown).
lookup_codes <- function(forward, site, alleles) {
  out <- rep(-1L, length(site))
  for (c in seq_len(ncol(forward))) {
    hit <- forward[cbind(site, c)] == alleles & nzchar(alleles)
    hit[is.na(hit)] <- FALSE
    out[hit] <- c - 1L
  }
  out
}

# Remap matrix for the given source rows (1-based) aligned with reference
# sites (1-based): entry [k, a] is the unified code of source allele a.
remap_matrix <- function(forward, src_alleles, rows, ref_sites) {
  max_a <- ncol(src_alleles)
  remap <- matrix(-1L, length(rows), max_a)
  if (length(rows) == 0) {
    return(remap)
  }
  for (j in seq_len(max_a)) {
    a <- src_alleles[rows, j]
    remap[, j] <- lookup_codes(forward, ref_sites, a)
  }
  remap
}

# The unified site set over several sources (tsinfer's MultiSourceView).
multi_source_view <- function(sources, anc_state) {
  srcs <- lapply(sources, open_source)
  ann <- zarr_open(anc_state$path)
  ann_pos <- as.numeric(zarr_read(ann, "variant_position"))
  ann_values <- if (isTRUE(anc_state$is_reference)) {
    a <- zarr_read(ann, "variant_allele")
    if (is.null(dim(a))) a else a[, 1]
  } else {
    zarr_read(ann, anc_state$field)
  }
  ann_values <- as.character(ann_values)
  valid <- lapply(srcs, function(s) {
    dups <- duplicate_positions(s$positions)
    if (length(dups) > 0) {
      message(sprintf(
        "Source '%s': %d duplicate position(s) excluded from inference",
        s$name, length(dups)
      ))
    }
    sort(unique(s$positions[!s$positions %in% dups & s$positions %in% ann_pos]))
  })
  positions <- sort(unique(unlist(valid)))
  has_site <- vapply(valid, function(v) positions %in% v, logical(length(positions)))
  if (!is.matrix(has_site)) {
    has_site <- matrix(has_site, nrow = length(positions))
  }
  anc_alleles <- ann_values[match(positions, ann_pos)]
  filtered <- lapply(srcs, function(s) {
    if (source_has_filters(s$config)) filtered_positions(s) else NULL
  })
  list(
    sources = srcs, positions = positions, has_site = has_site,
    ancestral = anc_alleles, filtered = filtered
  )
}

view_num_haplotypes <- function(view) {
  sum(vapply(view$sources, function(s) length(s$selection) * s$ploidy, numeric(1)))
}

# (source name, sample id, ploidy index) of each haplotype column of the view.
view_haplotypes <- function(view) {
  do.call(rbind, lapply(view$sources, function(s) {
    cols <- sort(s$selection)
    data.frame(
      source = s$name,
      sample_id = rep(s$sample_ids[cols + 1], each = s$ploidy),
      ploidy_index = rep(seq_len(s$ploidy) - 1L, times = length(cols)),
      stringsAsFactors = FALSE
    )
  }))
}

# Build the genotype view specification consumed by compiled code for a
# subset of the view's positions.
view_spec <- function(view, positions, apply_filters = FALSE) {
  idx <- match(positions, view$positions)
  if (anyNA(idx)) {
    stop("positions must be a subset of the view positions", call. = FALSE)
  }
  table <- build_allele_table(
    positions, view$ancestral[idx],
    lapply(view$sources, function(s) list(positions = s$positions, alleles = s$alleles))
  )
  specs <- vector("list", length(view$sources))
  for (j in seq_along(view$sources)) {
    s <- view$sources[[j]]
    has <- view$has_site[idx, j]
    if (apply_filters && !is.null(view$filtered[[j]])) {
      has <- has & positions %in% view$filtered[[j]]
    }
    rows <- rep(-1, length(positions))
    src_rows <- match(positions[has], s$positions)
    rows[has] <- src_rows - 1
    remap <- matrix(-1L, length(positions), ncol(s$alleles))
    remap[has, ] <- remap_matrix(table$forward, s$alleles, src_rows, which(has))
    specs[[j]] <- list(
      gt = .zarr_array(s$gt_meta), rows = rows, remap = remap,
      cols = as.numeric(sort(s$selection)), ploidy = s$ploidy
    )
  }
  list(
    spec = list(num_positions = length(positions), sources = specs),
    table = table
  )
}
