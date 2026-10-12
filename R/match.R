#' Match haplotypes to build a tree sequence
#'
#' Matches the ancestors against each other and then the samples against the
#' ancestors using the Li and Stephens model, building a tree sequence. This
#' is the second stage of the inference pipeline and is equivalent to
#' `tsinfer.match()`.
#'
#' All haplotypes (ancestors and samples from every source in
#' `config$match$sources`) are grouped so that each haplotype is matched only
#' against older haplotypes; groups are processed in order and the
#' haplotypes within a group are matched in parallel.
#'
#' If `config$match$workdir` is set, the list of match jobs and a checkpoint
#' tree sequence after each group are written there, and an interrupted run
#' resumes from the last checkpoint.
#'
#' @param config A `tsinfer_config`.
#' @param num_threads Number of matching threads (0 or 1 = single-threaded).
#' @param cache_size Size of the genotype chunk cache in MiB.
#' @param group_stop Stop before this group index (for partial runs).
#' @param match_file Optional path of a JSON-lines file recording the copying
#'   path and mutations of every haplotype.
#' @param path_compression Override `config$match$path_compression`.
#' @param verbose Print progress messages.
#' @return The inferred (raw) tree sequence as a `tsinfer_ts`.
#' @export
match_haplotypes <- function(config, num_threads = 0L, cache_size = 256,
                             group_stop = NULL, match_file = NULL,
                             path_compression = NULL, verbose = FALSE) {
  if (is.null(path_compression)) {
    path_compression <- config$match$path_compression
  }
  if (!is.null(group_stop) && group_stop < 0) {
    stop("group_stop must be non-negative", call. = FALSE)
  }
  if (length(config$ancestors) == 0) {
    stop("Matching against a reference tree sequence ([match].reference_ts) is not ",
      "supported yet; an [[ancestors]] section is required",
      call. = FALSE
    )
  }
  if (!is.null(group_stop) && is.null(config$match$workdir)) {
    warning("group_stop without workdir; partial result cannot be resumed", call. = FALSE)
  }
  anc_group <- zarr_open(config$ancestors[[1]]$path)
  positions <- as.numeric(zarr_read(anc_group, "variant_position"))
  site_alleles <- zarr_read(anc_group, "variant_allele")
  if (is.null(dim(site_alleles))) {
    site_alleles <- matrix(site_alleles, ncol = 2)
  }
  seq_intervals <- zarr_read(anc_group, "sequence_intervals")
  if (is.null(dim(seq_intervals))) {
    seq_intervals <- matrix(seq_intervals, ncol = 2)
  }
  if (length(positions) == 0) {
    stop("Cannot match haplotypes: the ancestor store has no inference sites", call. = FALSE)
  }

  workdir <- config$match$workdir
  last_group <- -1L
  start_ts <- NULL
  meta <- build_match_jobs(config, positions)
  if (!is.null(workdir)) {
    dir.create(workdir, recursive = TRUE, showWarnings = FALSE)
    jobs_path <- file.path(workdir, "match-jobs.json")
    if (file.exists(jobs_path)) {
      if (verbose) message("Loading existing groups from ", jobs_path)
      meta$jobs <- load_match_jobs(jobs_path, meta$jobs)
    } else {
      write_match_jobs(meta$jobs, jobs_path)
    }
    done <- list.files(workdir, pattern = "^group_[0-9]+\\.trees$")
    if (length(done) > 0) {
      last_group <- max(as.integer(sub("^group_([0-9]+)\\.trees$", "\\1", done)))
      start_ts <- ts_load(file.path(workdir, sprintf("group_%d.trees", last_group)))
      if (verbose) message("Resuming from group ", last_group)
    }
  }
  jobs <- meta$jobs

  anc_names <- vapply(config$ancestors, function(a) a$name, character(1))
  seq_len <- NULL
  for (name in names(config$match$sources)) {
    if (!name %in% anc_names) {
      seq_len <- as.numeric(zarr_read(zarr_open(config$sources[[name]]$path), "contig_length")[1])
      break
    }
  }
  if (is.null(seq_len)) {
    seq_len <- if (nrow(seq_intervals) > 0) max(seq_intervals) else 1
  }

  engine <- new_match_engine(config, positions, site_alleles[, 1], jobs, cache_size)

  groups <- sort(unique(jobs$group))
  todo <- groups[groups > last_group]
  if (!is.null(group_stop)) {
    todo <- todo[todo < group_stop]
  }
  scheduled <- which(jobs$group %in% todo) - 1L
  .match_engine_schedule(engine$ptr, scheduled)

  if (is.null(start_ts)) {
    anc_times <- zarr_read(anc_group, "sample_time")
    max_time <- if (length(anc_times) > 0) max(anc_times) else 0
    intervals_json <- canonical_json(list(sequence_intervals = lapply(
      seq_len(nrow(seq_intervals)),
      function(i) as.list(as.integer(seq_intervals[i, ]))
    )))
    ptr <- .make_root_tables(
      seq_len, positions, engine$forward[, 1], intervals_json, max_time,
      meta$individuals, if (length(meta$populations) > 0) meta$populations else NULL
    )
    ts <- new_ts(ptr)
  } else {
    ts <- start_ts
  }
  if (verbose) {
    message(sprintf("Match: %d haplotypes, %d sites, seq_len=%.0f", nrow(jobs), length(positions), seq_len))
  }

  fh <- NULL
  if (!is.null(match_file)) {
    fh <- file(match_file, open = "w")
    on.exit(close(fh), add = TRUE)
  }
  config_str <- format(config)
  prev_written <- NULL
  for (g in todo) {
    ids <- which(jobs$group == g)
    ids <- ids[order(jobs$haplotype_index[ids])] - 1L
    if (verbose) {
      message(sprintf("Group %d (%d of %d): %d haplotypes", g, match(g, groups), length(groups), length(ids)))
    }
    prov <- provenance_record("match", list(group = g, config = config_str))
    res <- .match_group(
      engine$ptr, ts$ptr, ids, as.integer(num_threads),
      isTRUE(path_compression), prov, !is.null(fh)
    )
    if (!is.null(fh)) {
      write_match_records(fh, g, jobs[ids + 1, , drop = FALSE], res$results, positions)
    }
    if (!is.null(workdir)) {
      path <- file.path(workdir, sprintf("group_%d.trees", g))
      ts_dump(ts, path)
      if (!config$match$keep_intermediates && !is.null(prev_written)) {
        unlink(file.path(workdir, sprintf("group_%d.trees", prev_written)))
      }
      prev_written <- g
    }
    Sys.sleep(0)
  }
  if (verbose) {
    info <- ts_info(ts)
    message(sprintf("Match complete: %d nodes, %d individuals", info$num_nodes, info$num_individuals))
  }
  ts
}

# Build one match job per haplotype, with individuals and populations.
build_match_jobs <- function(config, positions) {
  rows <- list()
  ind_keys <- character(0)
  ind_source <- character(0)
  ind_sample <- character(0)
  hap <- 0L
  for (name in names(config$match$sources)) {
    src <- config$sources[[name]]
    if (is.null(src)) {
      stop("[match.sources] references unknown source: '", name, "'", call. = FALSE)
    }
    mcfg <- config$match$sources[[name]]
    g <- zarr_open(src$path)
    info <- zarr_info(g, "call_genotype")
    all_ids <- as.character(zarr_read(g, "sample_id"))
    sel <- parse_samples_spec(src$samples, all_ids)
    ploidy <- as.integer(info$shape[3])
    st <- resolve_field(g, src$sample_time, "sample_id", info$shape[2], fill_value = 0)
    sample_time <- if (is.null(st)) rep(0, length(sel)) else as.numeric(st)[sel + 1]
    sample_ids <- all_ids[sel + 1]
    arrays <- zarr_list(g)
    if (all(c("sample_start_position", "sample_end_position") %in% arrays)) {
      start <- as.numeric(zarr_read(g, "sample_start_position"))[sel + 1]
      end <- as.numeric(zarr_read(g, "sample_end_position"))[sel + 1]
    } else {
      start <- rep(positions[1], length(sel))
      end <- rep(positions[length(positions)] + 1, length(sel))
    }
    individual <- rep(NA_integer_, length(sel))
    if (mcfg$create_individuals) {
      keys <- paste(name, sample_ids, sep = "\r")
      new <- unique(keys[!keys %in% ind_keys])
      new_idx <- match(new, keys)
      ind_keys <- c(ind_keys, new)
      ind_source <- c(ind_source, rep(name, length(new)))
      ind_sample <- c(ind_sample, sample_ids[new_idx])
      individual <- match(keys, ind_keys) - 1L
    }
    n <- length(sel)
    nh <- n * ploidy
    rows[[length(rows) + 1]] <- data.frame(
      haplotype_index = hap + seq_len(nh) - 1L,
      source = rep(name, nh),
      sample_id = rep(sample_ids, each = ploidy),
      ploidy_index = rep(seq_len(ploidy) - 1L, times = n),
      time = rep(sample_time, each = ploidy),
      start_position = rep(start, each = ploidy),
      end_position = rep(end, each = ploidy),
      group = 0L,
      node_flags = rep(mcfg$node_flags, nh),
      individual_id = rep(individual, each = ploidy),
      population_id = NA_integer_,
      sample_chunk = rep(sel %/% info$chunks[2], each = ploidy),
      store_col = rep(sel, each = ploidy),
      stringsAsFactors = FALSE
    )
    hap <- hap + nh
  }
  jobs <- do.call(rbind, rows)

  ind_cfg <- config$individual_metadata
  ind_md <- lapply(seq_along(ind_keys), function(i) {
    list(source = ind_source[i], sample_id = ind_sample[i])
  })
  pop_values <- rep(NA_character_, length(ind_keys))
  if (length(ind_keys) > 0 && !is.null(ind_cfg)) {
    for (name in unique(ind_source)) {
      g <- zarr_open(config$sources[[name]]$path)
      arrays <- zarr_list(g)
      all_ids <- as.character(zarr_read(g, "sample_id"))
      which_ind <- which(ind_source == name)
      sidx <- match(ind_sample[which_ind], all_ids)
      for (field in names(ind_cfg$fields)) {
        array_name <- ind_cfg$fields[[field]]
        if (array_name %in% arrays) {
          vals <- metadata_values(g, array_name, sidx)
          for (k in seq_along(which_ind)) {
            ind_md[[which_ind[k]]][[field]] <- vals[[k]]
          }
        }
      }
      if (!is.null(ind_cfg$population) && ind_cfg$population %in% arrays) {
        vals <- zarr_read(g, ind_cfg$population)
        pv <- vals[sidx]
        # Python's str() of the value.
        pop_values[which_ind] <- if (is.double(pv)) .json_double(pv) else as.character(pv)
      }
    }
  }
  populations <- character(0)
  if (!is.null(ind_cfg$population)) {
    pops <- unique(pop_values[!is.na(pop_values)])
    if (length(pops) > 0) {
      populations <- vapply(pops, function(p) canonical_json(list(name = p)), character(1), USE.NAMES = FALSE)
      ind_pop <- match(pop_values, pops) - 1L
      has_ind <- !is.na(jobs$individual_id)
      jobs$population_id[has_ind] <- ind_pop[jobs$individual_id[has_ind] + 1]
    }
  }
  groups <- .group_by_linesweep(jobs$start_position, jobs$end_position, jobs$time)
  jobs$group <- groups
  jobs <- jobs[order(jobs$group, jobs$haplotype_index), , drop = FALSE]
  rownames(jobs) <- NULL
  list(
    jobs = jobs,
    individuals = vapply(ind_md, canonical_json, character(1)),
    populations = populations
  )
}

# Values of a sample-dimensioned array, typed for JSON metadata.
metadata_values <- function(group, array_name, idx) {
  info <- zarr_info(group, array_name)
  vals <- zarr_read(group, array_name)
  if (!is.null(dim(vals))) {
    return(lapply(idx, function(i) as.list(vals[i, ])))
  }
  v <- vals[idx]
  if (is.numeric(v) && info$dtype %in% c("int64", "uint32", "uint64")) {
    return(lapply(v, json_int))
  }
  as.list(v)
}

match_job_columns <- c(
  "haplotype_index", "source", "sample_id", "ploidy_index", "time",
  "start_position", "end_position", "group", "node_flags", "individual_id",
  "population_id", "sample_chunk"
)

write_match_jobs <- function(jobs, path) {
  out <- jobs[, match_job_columns]
  json <- jsonlite::toJSON(out, dataframe = "rows", na = "null", digits = NA, pretty = TRUE)
  writeLines(json, path)
}

load_match_jobs <- function(path, fresh) {
  loaded <- jsonlite::fromJSON(path, simplifyVector = TRUE)
  key <- function(d) paste(d$source, d$sample_id, d$ploidy_index, sep = "\r")
  idx <- match(key(loaded), key(fresh))
  if (anyNA(idx)) {
    stop("match-jobs.json in the workdir does not match the configuration", call. = FALSE)
  }
  jobs <- fresh[idx, , drop = FALSE]
  for (col in intersect(match_job_columns, names(loaded))) {
    jobs[[col]] <- loaded[[col]]
  }
  jobs$individual_id <- as.integer(jobs$individual_id)
  jobs$population_id <- as.integer(jobs$population_id)
  rownames(jobs) <- NULL
  jobs
}

new_match_engine <- function(config, positions, ancestral, jobs, cache_size) {
  srcs <- lapply(names(config$match$sources), function(name) {
    src <- config$sources[[name]]
    g <- zarr_open(src$path)
    alleles <- zarr_read(g, "variant_allele")
    if (is.null(dim(alleles))) {
      alleles <- matrix(alleles, ncol = 1)
    }
    alleles[is.na(alleles)] <- ""
    list(
      name = name, group = g,
      positions = as.numeric(zarr_read(g, "variant_position")),
      alleles = alleles,
      gt_meta = zarr_meta(g, "call_genotype")$meta,
      ploidy = as.integer(zarr_info(g, "call_genotype")$shape[3])
    )
  })
  table <- build_allele_table(positions, ancestral, lapply(srcs, function(s) {
    list(positions = s$positions, alleles = s$alleles)
  }))
  source_specs <- lapply(seq_along(srcs), function(j) {
    s <- srcs[[j]]
    mcfg <- config$match$sources[[s$name]]
    # For duplicated positions the last source row is used, as in tsinfer.
    last <- !duplicated(s$positions, fromLast = TRUE)
    sel <- which(last & s$positions %in% positions)
    ref <- match(s$positions[sel], positions)
    list(
      gt_meta = s$gt_meta, ploidy = s$ploidy, sel_rows = as.numeric(sel - 1),
      sel_ref = as.integer(ref - 1), remap = remap_matrix(table$forward, s$alleles, sel, ref),
      recombination = if (is.null(mcfg$recombination)) 1e-2 else mcfg$recombination,
      mismatch = if (is.null(mcfg$mismatch)) 1e-20 else mcfg$mismatch
    )
  })
  src_index <- match(jobs$source, vapply(srcs, function(s) s$name, character(1))) - 1L
  node_md <- paste0(
    "{\"ploidy_index\":", jobs$ploidy_index,
    ",\"sample_id\":", .json_quote(jobs$sample_id),
    ",\"source\":", .json_quote(jobs$source), "}"
  )
  alleles <- lapply(seq_len(nrow(table$forward)), function(i) {
    a <- table$forward[i, ]
    a[seq_len(max(2L, table$num_alleles[i]))]
  })
  ptr <- .match_engine_new(
    positions, alleles, source_specs,
    list(
      source = src_index, store_col = as.numeric(jobs$store_col),
      ploidy_index = as.integer(jobs$ploidy_index), time = as.numeric(jobs$time),
      node_flags = as.numeric(jobs$node_flags), individual_id = as.integer(jobs$individual_id),
      population_id = as.integer(jobs$population_id), metadata = node_md
    ),
    cache_size * 1024 * 1024
  )
  list(ptr = ptr, forward = table$forward)
}

write_match_records <- function(fh, group, jobs, results, positions) {
  for (i in seq_len(nrow(jobs))) {
    r <- results[[i]]
    doc <- list(
      group = group, haplotype_index = jobs$haplotype_index[i],
      source = jobs$source[i], sample_id = jobs$sample_id[i],
      ploidy_index = jobs$ploidy_index[i], time = jobs$time[i],
      path = lapply(seq_along(r$left), function(k) {
        list(left = r$left[k], right = r$right[k], parent = r$parent[k])
      }),
      mutations = lapply(seq_along(r$site), function(k) {
        list(position = positions[r$site[k] + 1], derived_state = r$derived_state[k])
      })
    )
    writeLines(jsonlite::toJSON(doc, auto_unbox = TRUE, digits = NA), fh)
  }
}

provenance_record <- function(command, parameters) {
  record <- list(
    schema_version = "1.0.0",
    software = list(name = "tsinferR", version = as.character(utils::packageVersion("tsinferR"))),
    parameters = c(list(command = command), parameters),
    environment = list(
      os = list(system = Sys.info()[["sysname"]], release = Sys.info()[["release"]]),
      r = list(version = R.version.string)
    )
  )
  as.character(jsonlite::toJSON(record, auto_unbox = TRUE, digits = NA, null = "null"))
}
