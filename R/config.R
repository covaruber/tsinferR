#' Pipeline configuration
#'
#' Configuration objects for the inference pipeline, mirroring the 'tsinfer'
#' TOML configuration format. A configuration is usually read from a TOML
#' file with [read_config()], but can also be built directly from the
#' constructor functions below.
#'
#' * `tsinfer_source()` describes a named VCZ store. Variants can be filtered
#'   with bcftools-style `include`/`exclude` expressions (e.g.
#'   `"QUAL > 30"`, `'TYPE="snp"'`), `regions`/`targets` strings (e.g.
#'   `"1:1000-50000"`), and samples selected with `samples` (e.g.
#'   `"s1,s2"` or `"^s3"`). `sample_time` gives the time of each sample: the
#'   name of an array in the store, a `list(path=, field=)` pointing to an
#'   array in another store (joined on `sample_id`), or a number.
#' * `ancestral_state()` gives the store and array holding the ancestral
#'   allele of each variant (joined on `variant_position`). With
#'   `is_reference = TRUE` the reference allele is used.
#' * `ancestors_config()` configures [infer_ancestors()].
#' * `match_source_config()` and `match_config()` configure [match()].
#' * `post_process_config()`, `augment_sites_config()` and
#'   `individual_metadata_config()` configure the remaining steps.
#' * `tsinfer_config()` combines everything. Ancestor stores are added as
#'   sources automatically (with `sample_time = "sample_time"`).
#'
#' @param path Path to a VCZ store (or, for `ancestors_config()`, the output
#'   path of the ancestor store).
#' @param name Name of the source or ancestor set.
#' @param include,exclude Optional bcftools-style variant filter expression.
#' @param samples Optional bcftools-style samples selection string.
#' @param regions,targets Optional bcftools-style region strings.
#' @param sample_time Sample time specification (see Details).
#' @param field Array holding the ancestral alleles.
#' @param is_reference Use the reference allele as the ancestral allele.
#' @param sources For `ancestors_config()` and `augment_sites_config()`, a
#'   character vector of source names; for `match_config()`, a named list of
#'   `match_source_config()` objects (one per source to match, in matching
#'   order); for `tsinfer_config()`, a list of `tsinfer_source()` objects.
#' @param max_gap_length Inference sites further apart than this are split
#'   into separate intervals.
#' @param samples_chunk_size,variants_chunk_size Chunk sizes of the ancestor
#'   store.
#' @param compressor,compression_level Blosc compressor of the ancestor store.
#' @param node_flags Flags of the nodes created for this source (1 = sample).
#' @param create_individuals Create a 'tskit' individual for each sample.
#' @param recombination,mismatch Per-site recombination and mismatch
#'   probabilities used when matching this source (defaults 1e-2 and 1e-20).
#' @param output Output path of the inferred tree sequence (used by the
#'   command-line style workflow; not written by [match()]).
#' @param path_compression Apply path compression after each match group.
#' @param reference_ts Not supported yet; must be `NULL`.
#' @param workdir Optional directory for match checkpoints, allowing an
#'   interrupted [match()] to resume.
#' @param keep_intermediates Keep the checkpoint of every group.
#' @param split_ultimate,erase_flanks Post-processing options.
#' @param fields Named character vector/list mapping individual metadata keys
#'   to sample-dimensioned arrays.
#' @param population Name of a sample-dimensioned array whose values define
#'   populations.
#' @param ancestors A list of `ancestors_config()` objects.
#' @param ancestral_state For `tsinfer_config()`, an `ancestral_state()`.
#' @param match A `match_config()`.
#' @param individual_metadata,post_process,augment_sites Optional step
#'   configurations.
#' @return A configuration object (a classed list).
#' @name config
#' @examples
#' cfg <- tsinfer_config(
#'   sources = list(tsinfer_source("samples.vcz", "samples")),
#'   ancestral_state = ancestral_state("samples.vcz", "variant_ancestral_allele"),
#'   ancestors = list(ancestors_config("ancestors", "ancestors.vcz", "samples")),
#'   match = match_config(
#'     sources = list(
#'       ancestors = match_source_config(node_flags = 0, create_individuals = FALSE),
#'       samples = match_source_config()
#'     ),
#'     output = "out.trees"
#'   ),
#'   post_process = post_process_config()
#' )
#' cfg
NULL

DEFAULTS <- list(
  max_gap_length = 500000L,
  samples_chunk_size = 100L,
  variants_chunk_size = 50000L,
  compressor = "zstd",
  compression_level = 7L,
  node_flags = 1L,
  cache_size = 256
)

#' @rdname config
#' @export
tsinfer_source <- function(path, name, include = NULL, exclude = NULL,
                           samples = NULL, regions = NULL, targets = NULL,
                           sample_time = NULL) {
  if (missing(name) || is.null(name) || !nzchar(name)) {
    stop("Source must have a non-empty name", call. = FALSE)
  }
  if (!is.null(include) && !is.null(exclude)) {
    stop("Source '", name, "': include and exclude cannot both be set", call. = FALSE)
  }
  structure(list(
    path = path, name = name, include = include, exclude = exclude,
    samples = samples, regions = regions, targets = targets,
    sample_time = sample_time
  ), class = "tsinfer_source")
}

#' @rdname config
#' @export
ancestral_state <- function(path, field = NULL, is_reference = NULL) {
  if (is.null(is_reference) && is.null(field)) {
    stop("[ancestral_state] requires either 'field' or 'is_reference = true'", call. = FALSE)
  }
  if (isTRUE(is_reference) && !is.null(field)) {
    stop("[ancestral_state] field must not be set when is_reference is true", call. = FALSE)
  }
  structure(list(path = path, field = field, is_reference = is_reference),
    class = "tsinfer_ancestral_state"
  )
}

#' @rdname config
#' @export
ancestors_config <- function(name = "ancestors", path, sources,
                             max_gap_length = DEFAULTS$max_gap_length,
                             samples_chunk_size = DEFAULTS$samples_chunk_size,
                             variants_chunk_size = DEFAULTS$variants_chunk_size,
                             compressor = DEFAULTS$compressor,
                             compression_level = DEFAULTS$compression_level) {
  structure(list(
    name = as.character(name), path = path, sources = as.character(unlist(sources)),
    max_gap_length = as.integer(max_gap_length),
    samples_chunk_size = as.integer(samples_chunk_size),
    variants_chunk_size = as.integer(variants_chunk_size),
    compressor = as.character(compressor),
    compression_level = as.integer(compression_level)
  ), class = "tsinfer_ancestors_config")
}

#' @rdname config
#' @export
match_source_config <- function(node_flags = DEFAULTS$node_flags,
                                create_individuals = TRUE,
                                recombination = NULL, mismatch = NULL) {
  structure(list(
    node_flags = as.integer(node_flags),
    create_individuals = isTRUE(create_individuals),
    recombination = if (is.null(recombination)) NULL else as.numeric(recombination),
    mismatch = if (is.null(mismatch)) NULL else as.numeric(mismatch)
  ), class = "tsinfer_match_source_config")
}

#' @rdname config
#' @export
match_config <- function(sources, output = NULL, path_compression = FALSE,
                         reference_ts = NULL, workdir = NULL,
                         keep_intermediates = FALSE) {
  if (!is.list(sources) || is.null(names(sources)) || any(!nzchar(names(sources)))) {
    stop("[match] sources must be a named list of match_source_config()", call. = FALSE)
  }
  sources <- lapply(sources, function(s) {
    if (inherits(s, "tsinfer_match_source_config")) s else do.call(match_source_config, as.list(s))
  })
  structure(list(
    sources = sources, output = output,
    path_compression = isTRUE(path_compression),
    reference_ts = reference_ts, workdir = workdir,
    keep_intermediates = isTRUE(keep_intermediates)
  ), class = "tsinfer_match_config")
}

#' @rdname config
#' @export
post_process_config <- function(split_ultimate = TRUE, erase_flanks = TRUE) {
  structure(list(split_ultimate = isTRUE(split_ultimate), erase_flanks = isTRUE(erase_flanks)),
    class = "tsinfer_post_process_config"
  )
}

#' @rdname config
#' @export
augment_sites_config <- function(sources) {
  structure(list(sources = as.character(unlist(sources))), class = "tsinfer_augment_sites_config")
}

#' @rdname config
#' @export
individual_metadata_config <- function(fields = list(), population = NULL) {
  fields <- as.list(fields)
  if (length(fields) > 0 && (is.null(names(fields)) || any(!nzchar(names(fields))))) {
    stop("[individual_metadata] fields must be named", call. = FALSE)
  }
  structure(list(fields = fields, population = population),
    class = "tsinfer_individual_metadata_config"
  )
}

#' @rdname config
#' @export
tsinfer_config <- function(sources, ancestral_state, ancestors = list(), match,
                           individual_metadata = NULL, post_process = NULL,
                           augment_sites = NULL) {
  if (inherits(sources, "tsinfer_source")) {
    sources <- list(sources)
  }
  if (inherits(ancestors, "tsinfer_ancestors_config")) {
    ancestors <- list(ancestors)
  }
  src_names <- vapply(sources, function(s) s$name, character(1))
  if (anyDuplicated(src_names)) {
    stop("Duplicate source name: '", src_names[duplicated(src_names)][1], "'", call. = FALSE)
  }
  names(sources) <- src_names
  for (anc in ancestors) {
    if (is.null(sources[[anc$name]])) {
      sources[[anc$name]] <- tsinfer_source(anc$path, anc$name, sample_time = "sample_time")
    }
  }
  if (length(ancestors) == 0 && is.null(match$reference_ts)) {
    stop("Config must contain either an [[ancestors]] section or [match].reference_ts",
      call. = FALSE
    )
  }
  if (match$keep_intermediates && is.null(match$workdir)) {
    stop("keep_intermediates requires workdir to be set", call. = FALSE)
  }
  for (anc in ancestors) {
    if (!anc$name %in% names(match$sources)) {
      stop("Ancestor '", anc$name, "' must appear in [match.sources]", call. = FALSE)
    }
  }
  structure(list(
    sources = sources, ancestors = ancestors, match = match,
    ancestral_state = ancestral_state, individual_metadata = individual_metadata,
    post_process = post_process, augment_sites = augment_sites
  ), class = "tsinfer_config")
}

KNOWN_KEYS <- list(
  top = c(
    "source", "ancestral_state", "ancestors", "match", "individual_metadata",
    "post_process", "augment_sites"
  ),
  source = c("name", "path", "include", "exclude", "samples", "regions", "targets", "sample_time"),
  ancestral_state = c("path", "field", "is_reference"),
  ancestors = c(
    "name", "path", "sources", "max_gap_length", "samples_chunk_size",
    "variants_chunk_size", "compressor", "compression_level"
  ),
  match = c("sources", "output", "path_compression", "reference_ts", "workdir", "keep_intermediates"),
  match_source = c("node_flags", "create_individuals", "recombination", "mismatch"),
  individual_metadata = c("fields", "population"),
  post_process = c("split_ultimate", "erase_flanks"),
  augment_sites = c("sources")
)

check_keys <- function(section, entry, known) {
  unknown <- setdiff(names(entry), known)
  if (length(unknown) > 0) {
    stop(sprintf(
      "Unrecognised key(s) in [%s]: %s", section,
      paste(sort(unknown), collapse = ", ")
    ), call. = FALSE)
  }
}

#' Read a TOML configuration file
#'
#' Parses a 'tsinfer' TOML configuration file (the same format as used by the
#' 'tsinfer' Python package and command-line tool) into a `tsinfer_config`.
#' Relative paths are interpreted relative to the current working directory.
#'
#' @param path Path to the TOML file.
#' @param text Alternatively, the TOML content as a string.
#' @return A `tsinfer_config`.
#' @export
#' @examples
#' toml <- '
#' [[source]]
#' name = "samples"
#' path = "samples.vcz"
#'
#' [ancestral_state]
#' path = "samples.vcz"
#' field = "variant_ancestral_allele"
#'
#' [[ancestors]]
#' name = "ancestors"
#' path = "ancestors.vcz"
#' sources = ["samples"]
#'
#' [match]
#' output = "out.trees"
#' [match.sources.ancestors]
#' node_flags = 0
#' create_individuals = false
#' [match.sources.samples]
#'
#' [post_process]
#' '
#' cfg <- read_config(text = toml)
#' names(cfg$sources)
read_config <- function(path, text = NULL) {
  if (is.null(text)) {
    text <- paste(readLines(path, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  }
  raw <- unclass(RcppTOML::parseTOML(text, fromFile = FALSE, escape = FALSE))
  if (is.list(raw$match) && is.list(raw$match$sources)) {
    raw$match$sources <- raw$match$sources[toml_key_order(names(raw$match$sources), text)]
  }
  parse_config(raw)
}

# RcppTOML returns table keys sorted; recover the order in which the match
# sources appear in the file, which determines the haplotype order.
toml_key_order <- function(keys, text) {
  pos <- vapply(keys, function(k) {
    k <- gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", k)
    pattern <- sprintf(
      "sources\\s*\\.\\s*(\"%1$s\"|'%1$s'|%1$s)\\s*[]=.]|[{,]\\s*(\"%1$s\"|%1$s)\\s*=",
      k
    )
    m <- regexpr(pattern, text, perl = TRUE)
    if (m < 0) Inf else as.numeric(m)
  }, numeric(1))
  order(pos, seq_along(keys))
}

as_table_list <- function(x) {
  # [[section]] arrays parse to unnamed lists of tables; [section] to a table.
  if (is.null(x)) {
    return(list())
  }
  if (!is.null(names(x)) && all(nzchar(names(x)))) {
    return(list(x))
  }
  x
}

parse_config <- function(raw) {
  check_keys("top-level", raw, KNOWN_KEYS$top)
  sources <- list()
  for (entry in as_table_list(raw$source)) {
    check_keys("source", entry, KNOWN_KEYS$source)
    name <- entry$name
    if (is.null(name) || !nzchar(name)) {
      stop("Every [[source]] entry must have a 'name' field", call. = FALSE)
    }
    if (is.null(entry$path)) {
      stop("Source '", name, "' must have a 'path' field", call. = FALSE)
    }
    if (!is.null(sources[[name]])) {
      stop("Duplicate source name: '", name, "'", call. = FALSE)
    }
    sources[[name]] <- tsinfer_source(
      path = entry$path, name = name, include = entry$include,
      exclude = entry$exclude, samples = entry$samples,
      regions = entry$regions, targets = entry$targets,
      sample_time = entry$sample_time
    )
  }
  ancestors <- lapply(as_table_list(raw$ancestors), function(entry) {
    check_keys("ancestors", entry, KNOWN_KEYS$ancestors)
    for (key in c("path", "sources")) {
      if (is.null(entry[[key]])) {
        stop("[ancestors] missing required key: '", key, "'", call. = FALSE)
      }
    }
    args <- entry[intersect(names(entry), names(formals(ancestors_config)))]
    if (is.null(args$name)) {
      args$name <- "ancestors"
    }
    do.call(ancestors_config, args)
  })
  m <- raw$match
  if (is.null(m)) {
    stop("Config must contain a [match] section", call. = FALSE)
  }
  check_keys("match", m, KNOWN_KEYS$match)
  if (is.null(m$sources) || !is.list(m$sources)) {
    stop("[match] sources must be a table of tables, e.g. [match.sources.cohort]", call. = FALSE)
  }
  if (is.null(m$output)) {
    stop("[match] missing required key: 'output'", call. = FALSE)
  }
  msources <- list()
  for (name in names(m$sources)) {
    entry <- m$sources[[name]]
    if (is.list(entry)) {
      check_keys(paste0("match.sources.", name), entry, KNOWN_KEYS$match_source)
      msources[[name]] <- do.call(match_source_config, entry)
    } else {
      msources[[name]] <- match_source_config()
    }
  }
  match <- match_config(
    sources = msources, output = m$output,
    path_compression = isTRUE(m$path_compression), reference_ts = m$reference_ts,
    workdir = m$workdir, keep_intermediates = isTRUE(m$keep_intermediates)
  )
  for (anc in ancestors) {
    if (!is.null(sources[[anc$name]])) {
      stop("Ancestor name '", anc$name, "' conflicts with a [[source]] name", call. = FALSE)
    }
    sources[[anc$name]] <- tsinfer_source(anc$path, anc$name, sample_time = "sample_time")
  }
  as_entry <- raw$ancestral_state
  if (is.null(as_entry)) {
    stop("Config must contain an [ancestral_state] section", call. = FALSE)
  }
  check_keys("ancestral_state", as_entry, KNOWN_KEYS$ancestral_state)
  if (is.null(as_entry$path)) {
    stop("[ancestral_state] missing required key: 'path'", call. = FALSE)
  }
  anc_state <- ancestral_state(as_entry$path, as_entry$field, as_entry$is_reference)
  ind_md <- NULL
  if (!is.null(raw$individual_metadata)) {
    check_keys("individual_metadata", raw$individual_metadata, KNOWN_KEYS$individual_metadata)
    ind_md <- individual_metadata_config(
      fields = raw$individual_metadata$fields %||% list(),
      population = raw$individual_metadata$population
    )
  }
  pp <- NULL
  if (!is.null(raw$post_process)) {
    check_keys("post_process", raw$post_process, KNOWN_KEYS$post_process)
    pp <- post_process_config(
      split_ultimate = raw$post_process$split_ultimate %||% TRUE,
      erase_flanks = raw$post_process$erase_flanks %||% TRUE
    )
  }
  aug <- NULL
  if (!is.null(raw$augment_sites)) {
    check_keys("augment_sites", raw$augment_sites, KNOWN_KEYS$augment_sites)
    if (is.null(raw$augment_sites$sources)) {
      stop("[augment_sites] missing required key: 'sources'", call. = FALSE)
    }
    aug <- augment_sites_config(raw$augment_sites$sources)
  }
  tsinfer_config(
    sources = unname(sources), ancestral_state = anc_state, ancestors = ancestors,
    match = match, individual_metadata = ind_md, post_process = pp,
    augment_sites = aug
  )
}

#' Validate and describe a configuration
#'
#' `validate_config()` checks that all input paths exist and that every
#' referenced source is defined. `format()` and `print()` describe the
#' resolved configuration.
#'
#' @param config,x A `tsinfer_config`.
#' @param ... Unused.
#' @return `validate_config()` returns a character vector of problems (empty
#'   if the configuration is valid). `format()` returns a string.
#' @export
validate_config <- function(config) {
  errors <- character(0)
  anc_names <- vapply(config$ancestors, function(a) a$name, character(1))
  for (src in config$sources) {
    if (src$name %in% anc_names) {
      next
    }
    if (!is.null(src$path) && is.character(src$path) && !file.exists(src$path)) {
      errors <- c(errors, sprintf("Source '%s' path does not exist: %s", src$name, src$path))
    }
    if (is.list(src$sample_time) && !is.null(src$sample_time$path) &&
      !file.exists(src$sample_time$path)) {
      errors <- c(errors, sprintf(
        "Source '%s' sample_time path does not exist: %s",
        src$name, src$sample_time$path
      ))
    }
  }
  for (anc in config$ancestors) {
    for (s in anc$sources) {
      if (is.null(config$sources[[s]])) {
        errors <- c(errors, sprintf("Ancestors '%s' references unknown source: '%s'", anc$name, s))
      }
    }
  }
  for (s in names(config$match$sources)) {
    if (is.null(config$sources[[s]])) {
      errors <- c(errors, sprintf("[match.sources] references unknown source: '%s'", s))
    }
  }
  if (!is.null(config$augment_sites)) {
    for (s in config$augment_sites$sources) {
      if (is.null(config$sources[[s]])) {
        errors <- c(errors, sprintf("[augment_sites] references unknown source: '%s'", s))
      }
    }
  }
  if (!is.null(config$match$reference_ts) && !file.exists(config$match$reference_ts)) {
    errors <- c(errors, sprintf(
      "Match reference_ts path does not exist: %s",
      config$match$reference_ts
    ))
  }
  if (!file.exists(config$ancestral_state$path)) {
    errors <- c(errors, sprintf(
      "Ancestral state path does not exist: %s",
      config$ancestral_state$path
    ))
  }
  errors
}

fmt_value <- function(x) {
  if (is.null(x)) {
    return("None")
  }
  if (is.list(x)) {
    return(as.character(jsonlite::toJSON(x, auto_unbox = TRUE)))
  }
  if (is.logical(x)) {
    return(if (isTRUE(x)) "True" else "False")
  }
  if (length(x) > 1) {
    return(paste0("[", paste0("'", x, "'", collapse = ", "), "]"))
  }
  as.character(x)
}

#' @rdname validate_config
#' @export
format.tsinfer_config <- function(x, ...) {
  lines <- character(0)
  add <- function(...) lines <<- c(lines, paste0(...))
  for (src in x$sources) {
    add("[source.", src$name, "]")
    add("  path = ", src$path)
    for (key in c("include", "exclude", "samples", "regions", "targets", "sample_time")) {
      if (!is.null(src[[key]])) {
        add("  ", key, " = ", fmt_value(src[[key]]))
      }
    }
    add("")
  }
  add("[ancestral_state]")
  add("  path = ", x$ancestral_state$path)
  add("  field = ", fmt_value(x$ancestral_state$field))
  add("")
  for (anc in x$ancestors) {
    add("[[ancestors]] (name=", anc$name, ")")
    add("  path = ", anc$path)
    add("  sources = ", fmt_value(anc$sources))
    add("  max_gap_length = ", anc$max_gap_length)
    add("")
  }
  add("[match]")
  add("  output = ", fmt_value(x$match$output))
  for (name in names(x$match$sources)) {
    s <- x$match$sources[[name]]
    add("  [match.sources.", name, "]")
    add("    node_flags = ", s$node_flags)
    add("    create_individuals = ", fmt_value(s$create_individuals))
  }
  add("  path_compression = ", fmt_value(x$match$path_compression))
  if (!is.null(x$match$reference_ts)) add("  reference_ts = ", x$match$reference_ts)
  if (!is.null(x$match$workdir)) add("  workdir = ", x$match$workdir)
  if (x$match$keep_intermediates) add("  keep_intermediates = True")
  add("")
  if (!is.null(x$post_process)) {
    add("[post_process]")
    add("  split_ultimate = ", fmt_value(x$post_process$split_ultimate))
    add("  erase_flanks = ", fmt_value(x$post_process$erase_flanks))
    add("")
  }
  if (!is.null(x$individual_metadata)) {
    add("[individual_metadata]")
    add("  fields = ", fmt_value(x$individual_metadata$fields))
    if (!is.null(x$individual_metadata$population)) {
      add("  population = ", x$individual_metadata$population)
    }
    add("")
  }
  if (!is.null(x$augment_sites)) {
    add("[augment_sites]")
    add("  sources = ", fmt_value(x$augment_sites$sources))
    add("")
  }
  paste(lines, collapse = "\n")
}

#' @rdname validate_config
#' @export
print.tsinfer_config <- function(x, ...) {
  cat(format(x), "\n", sep = "")
  invisible(x)
}
