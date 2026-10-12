# Canonical JSON encoding compatible with 'tskit' (Python json.dumps with
# sort_keys=True, compact separators and ASCII escapes), so that metadata
# written by tsinferR is byte-identical to that written by 'tsinfer'.

canonical_json <- function(x) {
  if (is.null(x)) {
    return("null")
  }
  if (is.list(x)) {
    nm <- names(x)
    if (is.null(nm)) {
      return(paste0("[", paste(vapply(x, canonical_json, character(1)), collapse = ","), "]"))
    }
    nm_sorted <- sort(nm, method = "radix")
    parts <- vapply(nm_sorted, function(k) {
      paste0(.json_quote(k), ":", canonical_json(x[[k]]))
    }, character(1))
    return(paste0("{", paste(parts, collapse = ","), "}"))
  }
  if (length(x) != 1 || inherits(x, "AsIs")) {
    return(paste0("[", paste(vapply(as.list(unclass(x)), canonical_json, character(1)), collapse = ","), "]"))
  }
  if (is.character(x)) {
    return(.json_quote(x))
  }
  if (is.logical(x)) {
    return(if (is.na(x)) "null" else if (x) "true" else "false")
  }
  if (is.integer(x)) {
    return(if (is.na(x)) "null" else as.character(x))
  }
  if (inherits(x, "json_int")) {
    return(format(unclass(x), scientific = FALSE, trim = TRUE))
  }
  .json_double(x)
}

# Mark doubles that hold integer values (e.g. int64 Zarr data) so they are
# encoded without a decimal point.
json_int <- function(x) structure(x, class = "json_int")
