#' Zarr stores
#'
#' A self-contained engine for the 'Zarr' array storage format (versions 2 and
#' 3) on the local file system, used to read and write VCF Zarr ('VCZ')
#' stores. Chunk encoding and decoding is implemented in compiled code with
#' the bundled 'c-blosc', 'Zstandard', 'LZ4' and 'zlib' libraries, supporting
#' the `blosc`, `zstd`, `zlib`, `gzip`, `lz4` and `crc32c` codecs, the
#' `vlen-utf8` string encoding, fixed-width strings and `packbits` booleans.
#'
#' * `zarr_open()` opens an existing group (directory).
#' * `zarr_list()` lists the arrays in a group.
#' * `zarr_info()` returns the shape, chunks, data type and attributes of an
#'   array.
#' * `zarr_read()` reads an array, or a range of rows along its first
#'   dimension, into an R vector/matrix/array. Integer types up to 16 bits and
#'   32-bit signed integers become `integer`, other numeric types `double`,
#'   booleans `logical` and strings `character`.
#'
#' Arrays are stored in C (row-major) order on disk; `zarr_read()` returns R
#' arrays whose dimensions match the Zarr shape.
#'
#' @param path Path to a Zarr group directory.
#' @param group A `zarr_group` from `zarr_open()` or a path.
#' @param name Name of an array within the group.
#' @param rows Optional integer vector of (1-based) rows to read along the
#'   first dimension.
#' @return `zarr_open()` returns a `zarr_group`; `zarr_list()` a character
#'   vector; `zarr_info()` a list; `zarr_read()` an R vector or array.
#' @name zarr
#' @examples
#' path <- file.path(tempdir(), "example.zarr")
#' g <- zarr_create_group(path, overwrite = TRUE)
#' zarr_write_array(g, "x", matrix(1:6, 2), chunks = c(1, 2))
#' zarr_write_array(g, "s", c("a", "bb", "ccc"))
#' zarr_list(g)
#' zarr_read(g, "x")
#' zarr_read(g, "s", rows = 2:3)
NULL

#' @rdname zarr
#' @export
zarr_open <- function(path) {
  path <- normalizePath(path.expand(path), mustWork = FALSE)
  if (!dir.exists(path)) {
    stop("Zarr store not found: ", path, call. = FALSE)
  }
  format <- if (file.exists(file.path(path, "zarr.json"))) 3L else 2L
  structure(list(path = path, zarr_format = format), class = "zarr_group")
}

as_zarr_group <- function(group) {
  if (inherits(group, "zarr_group")) {
    return(group)
  }
  zarr_open(group)
}

#' @export
print.zarr_group <- function(x, ...) {
  arrays <- zarr_list(x)
  cat(sprintf("<zarr_group> %s (Zarr v%d, %d arrays)\n", x$path, x$zarr_format, length(arrays)))
  if (length(arrays) > 0) {
    cat(strwrap(paste(arrays, collapse = ", "), indent = 2, exdent = 2), sep = "\n")
  }
  invisible(x)
}

is_zarr_array_dir <- function(dir) {
  if (file.exists(file.path(dir, ".zarray"))) {
    return(TRUE)
  }
  meta_path <- file.path(dir, "zarr.json")
  if (file.exists(meta_path)) {
    meta <- jsonlite::fromJSON(meta_path, simplifyVector = FALSE)
    return(identical(meta$node_type, "array"))
  }
  FALSE
}

#' @rdname zarr
#' @export
zarr_list <- function(group) {
  group <- as_zarr_group(group)
  dirs <- list.dirs(group$path, full.names = FALSE, recursive = FALSE)
  dirs <- dirs[vapply(file.path(group$path, dirs), is_zarr_array_dir, logical(1))]
  sort(dirs)
}

zarr_has <- function(group, name) {
  group <- as_zarr_group(group)
  is_zarr_array_dir(file.path(group$path, name))
}

parse_v2_dtype <- function(dtype) {
  if (!is.character(dtype) || length(dtype) != 1) {
    stop("Structured Zarr data types are not supported", call. = FALSE)
  }
  order <- substr(dtype, 1, 1)
  kind <- substr(dtype, 2, 2)
  size <- suppressWarnings(as.integer(substring(dtype, 3)))
  big_endian <- order == ">"
  if (kind == "O") {
    return(list(dtype = "string", itemsize = 0, big_endian = FALSE))
  }
  if (kind == "U") {
    return(list(dtype = "fixed_unicode", itemsize = 4 * size, big_endian = big_endian))
  }
  if (kind == "S") {
    return(list(dtype = "fixed_bytes", itemsize = size, big_endian = FALSE))
  }
  if (kind == "b") {
    return(list(dtype = "bool", itemsize = 1, big_endian = FALSE))
  }
  name <- switch(kind,
    i = paste0("int", 8 * size),
    u = paste0("uint", 8 * size),
    f = paste0("float", 8 * size),
    stop("Unsupported Zarr dtype: ", dtype, call. = FALSE)
  )
  list(dtype = name, itemsize = size, big_endian = big_endian)
}

v3_itemsize <- function(data_type) {
  sizes <- c(
    int8 = 1, int16 = 2, int32 = 4, int64 = 8, uint8 = 1, uint16 = 2,
    uint32 = 4, uint64 = 8, float32 = 4, float64 = 8, bool = 1, string = 0
  )
  if (!data_type %in% names(sizes)) {
    stop("Unsupported Zarr v3 data type: ", data_type, call. = FALSE)
  }
  sizes[[data_type]]
}

fill_numeric <- function(x) {
  if (is.null(x)) {
    return(0)
  }
  if (is.character(x)) {
    return(switch(x,
      "NaN" = NaN,
      "Infinity" = Inf,
      "-Infinity" = -Inf,
      suppressWarnings(as.numeric(x))
    ))
  }
  if (is.logical(x)) {
    return(as.numeric(x))
  }
  as.numeric(x)
}

v2_codec <- function(spec) {
  id <- spec$id
  if (id == "blosc") {
    return(list(
      id = "blosc", cname = spec$cname %||% "lz4",
      clevel = as.integer(spec$clevel %||% 5), shuffle = as.integer(spec$shuffle %||% 1),
      blocksize = as.integer(spec$blocksize %||% 0)
    ))
  }
  if (id %in% c("zstd", "zlib", "gzip")) {
    return(list(id = id, level = as.integer(spec$level %||% 1)))
  }
  if (id == "lz4") {
    return(list(id = id, level = as.integer(spec$acceleration %||% 1)))
  }
  stop("Unsupported Zarr codec: ", id, call. = FALSE)
}

`%||%` <- function(a, b) if (is.null(a)) b else a

read_json_file <- function(path) {
  jsonlite::fromJSON(path, simplifyVector = FALSE)
}

zarr_meta <- function(group, name) {
  group <- as_zarr_group(group)
  dir <- file.path(group$path, name)
  if (file.exists(file.path(dir, ".zarray"))) {
    m <- read_json_file(file.path(dir, ".zarray"))
    attrs <- if (file.exists(file.path(dir, ".zattrs"))) read_json_file(file.path(dir, ".zattrs")) else list()
    if (!identical(m$order %||% "C", "C")) {
      stop("Only C-order Zarr arrays are supported (", name, ")", call. = FALSE)
    }
    dt <- parse_v2_dtype(m$dtype)
    serializer <- "raw"
    codecs <- list()
    for (f in m$filters %||% list()) {
      if (f$id %in% c("vlen-utf8", "vlen-bytes")) {
        serializer <- "vlen-utf8"
      } else if (f$id == "packbits") {
        serializer <- "packbits"
      } else {
        stop("Unsupported Zarr filter '", f$id, "' in array ", name, call. = FALSE)
      }
    }
    if (dt$dtype == "string") {
      serializer <- "vlen-utf8"
    }
    if (!is.null(m$compressor)) {
      codecs <- list(v2_codec(m$compressor))
    }
    fill <- m$fill_value
    fill_string <- ""
    if (dt$dtype %in% c("string", "fixed_unicode")) {
      fill_string <- if (is.character(fill)) fill else ""
      fill <- 0
    } else if (dt$dtype == "fixed_bytes") {
      fill_string <- if (is.character(fill) && nzchar(fill)) rawToChar(jsonlite::base64_dec(fill)) else ""
      fill <- 0
    }
    meta <- list(
      path = dir, zarr_format = 2L,
      shape = as.numeric(unlist(m$shape)), chunks = as.numeric(unlist(m$chunks)),
      dtype = dt$dtype, itemsize = dt$itemsize, big_endian = dt$big_endian,
      serializer = serializer, codecs = codecs, key_prefix = "",
      separator = m$dimension_separator %||% ".",
      fill_value = fill_numeric(fill), fill_string = fill_string
    )
    return(list(meta = meta, attrs = attrs))
  }
  if (file.exists(file.path(dir, "zarr.json"))) {
    m <- read_json_file(file.path(dir, "zarr.json"))
    if (!identical(m$node_type, "array")) {
      stop("'", name, "' is not a Zarr array", call. = FALSE)
    }
    data_type <- m$data_type
    if (!is.character(data_type)) {
      stop("Unsupported Zarr v3 data type in ", name, call. = FALSE)
    }
    itemsize <- v3_itemsize(data_type)
    grid <- m$chunk_grid
    if (!identical(grid$name, "regular")) {
      stop("Only regular chunk grids are supported (", name, ")", call. = FALSE)
    }
    key <- m$chunk_key_encoding %||% list(name = "default")
    sep <- key$configuration$separator %||% (if (identical(key$name, "v2")) "." else "/")
    prefix <- if (identical(key$name, "v2")) "" else "c"
    serializer <- "raw"
    big_endian <- FALSE
    codecs <- list()
    for (codec in m$codecs) {
      cname <- codec$name
      conf <- codec$configuration %||% list()
      if (cname == "bytes") {
        big_endian <- identical(conf$endian, "big")
      } else if (cname %in% c("vlen-utf8", "vlen-bytes")) {
        serializer <- "vlen-utf8"
      } else if (cname == "blosc") {
        shuffle <- switch(conf$shuffle %||% "noshuffle",
          noshuffle = 0L, shuffle = 1L, bitshuffle = 2L, 0L
        )
        codecs[[length(codecs) + 1]] <- list(
          id = "blosc", cname = conf$cname %||% "lz4",
          clevel = as.integer(conf$clevel %||% 5), shuffle = shuffle,
          blocksize = as.integer(conf$blocksize %||% 0),
          typesize = as.integer(conf$typesize %||% 0)
        )
      } else if (cname %in% c("zstd", "gzip")) {
        codecs[[length(codecs) + 1]] <- list(id = cname, level = as.integer(conf$level %||% 1))
      } else if (cname == "crc32c") {
        codecs[[length(codecs) + 1]] <- list(id = "crc32c")
      } else {
        stop("Unsupported Zarr v3 codec '", cname, "' in array ", name, call. = FALSE)
      }
    }
    fill <- m$fill_value
    fill_string <- ""
    if (data_type == "string") {
      fill_string <- if (is.character(fill)) fill else ""
      fill <- 0
    }
    meta <- list(
      path = dir, zarr_format = 3L,
      shape = as.numeric(unlist(m$shape)),
      chunks = as.numeric(unlist(grid$configuration$chunk_shape)),
      dtype = data_type, itemsize = itemsize, big_endian = big_endian,
      serializer = serializer, codecs = codecs, key_prefix = prefix,
      separator = sep, fill_value = fill_numeric(fill), fill_string = fill_string
    )
    attrs <- m$attributes %||% list()
    if (!is.null(m$dimension_names)) {
      attrs[["_ARRAY_DIMENSIONS"]] <- m$dimension_names
    }
    return(list(meta = meta, attrs = attrs))
  }
  stop("Array '", name, "' not found in ", group$path, call. = FALSE)
}

zarr_handle <- function(group, name) {
  m <- zarr_meta(group, name)
  list(ptr = .zarr_array(m$meta), meta = m$meta, attrs = m$attrs)
}

#' @rdname zarr
#' @export
zarr_info <- function(group, name) {
  m <- zarr_meta(group, name)
  list(
    shape = m$meta$shape, chunks = m$meta$chunks, dtype = m$meta$dtype,
    zarr_format = m$meta$zarr_format, attributes = m$attrs
  )
}

#' @rdname zarr
#' @export
zarr_read <- function(group, name, rows = NULL) {
  h <- zarr_handle(group, name)
  shape <- h$meta$shape
  if (length(shape) == 0) {
    return(.zarr_read(h$ptr, numeric(0), numeric(0)))
  }
  if (is.null(rows)) {
    return(.zarr_read(h$ptr, rep(0, length(shape)), shape))
  }
  rows <- as.integer(rows)
  if (length(rows) == 0) {
    count <- shape
    count[1] <- 0
    return(.zarr_read(h$ptr, rep(0, length(shape)), count))
  }
  if (anyNA(rows) || min(rows) < 1 || max(rows) > shape[1]) {
    stop("rows out of range", call. = FALSE)
  }
  lo <- min(rows)
  hi <- max(rows)
  start <- c(lo - 1, rep(0, length(shape) - 1))
  count <- c(hi - lo + 1, shape[-1])
  x <- .zarr_read(h$ptr, start, count)
  idx <- rows - lo + 1
  if (length(shape) == 1) {
    return(x[idx])
  }
  if (length(shape) == 2) {
    return(x[idx, , drop = FALSE])
  }
  dims <- dim(x)
  m <- matrix(x, nrow = dims[1])[idx, , drop = FALSE]
  array(m, dim = c(length(idx), dims[-1]))
}

# ---------------------------------------------------------------------------
# Writing

#' Write Zarr stores
#'
#' * `zarr_create_group()` creates a new (Zarr v2) group directory.
#' * `zarr_write_array()` creates an array in a group and writes `data` to it.
#' * `zarr_compressor()` describes a compressor for `zarr_write_array()`.
#'
#' @param path Directory of the new group.
#' @param overwrite If `TRUE`, an existing directory is deleted first.
#' @param attrs A named list of group or array attributes.
#' @param group A `zarr_group` or path.
#' @param name Array name.
#' @param data An R vector, matrix or array. Its dimensions give the array
#'   shape.
#' @param dtype Zarr data type: one of `"int8"`, `"int16"`, `"int32"`,
#'   `"int64"`, `"uint8"`, `"uint16"`, `"uint32"`, `"uint64"`, `"float32"`,
#'   `"float64"`, `"bool"` or `"string"`. Inferred from `data` by default.
#' @param chunks Chunk shape (defaults to the full shape).
#' @param compressor A compressor from `zarr_compressor()`, or `NULL` for no
#'   compression.
#' @param fill_value Fill value for missing chunks.
#' @param dims Dimension names, stored in the `_ARRAY_DIMENSIONS` attribute.
#' @param id Compressor id: `"blosc"`, `"zstd"`, `"zlib"`, `"gzip"` or
#'   `"lz4"`.
#' @param cname Blosc internal compressor (`"lz4"`, `"lz4hc"`, `"blosclz"`,
#'   `"zstd"` or `"zlib"`).
#' @param clevel,level Compression level.
#' @param shuffle Blosc shuffle: `"noshuffle"`, `"shuffle"`, `"bitshuffle"`
#'   or `"auto"`.
#' @return `zarr_create_group()` returns a `zarr_group`; `zarr_write_array()`
#'   returns the group invisibly; `zarr_compressor()` returns a compressor
#'   specification.
#' @name zarr-write
#' @examples
#' g <- zarr_create_group(file.path(tempdir(), "w.zarr"), overwrite = TRUE)
#' zarr_write_array(g, "pos", c(10L, 20L, 30L), dims = "variants",
#'   compressor = zarr_compressor("zstd", level = 3))
#' zarr_read(g, "pos")
NULL

#' @rdname zarr-write
#' @export
zarr_create_group <- function(path, overwrite = FALSE, attrs = NULL) {
  path <- path.expand(path)
  if (dir.exists(path) || file.exists(path)) {
    if (!overwrite) {
      stop("Path already exists: ", path, call. = FALSE)
    }
    unlink(path, recursive = TRUE)
  }
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  write_json_file(list(zarr_format = 2L), file.path(path, ".zgroup"))
  write_json_file(if (is.null(attrs)) structure(list(), names = character(0)) else attrs,
    file.path(path, ".zattrs"))
  zarr_open(path)
}

write_json_file <- function(x, path) {
  json <- jsonlite::toJSON(x, auto_unbox = TRUE, null = "null", digits = NA, pretty = TRUE)
  writeLines(json, path, useBytes = TRUE)
}

#' @rdname zarr-write
#' @export
zarr_compressor <- function(id = c("blosc", "zstd", "zlib", "gzip", "lz4"),
                            cname = "zstd", clevel = 7L, shuffle = "auto",
                            level = 5L) {
  id <- match.arg(id)
  if (id == "blosc") {
    shuffle_code <- switch(shuffle,
      noshuffle = 0L, shuffle = 1L, bitshuffle = 2L, auto = -1L,
      stop("invalid shuffle: ", shuffle, call. = FALSE)
    )
    if (!.blosc_has_compressor(cname)) {
      stop("Blosc compressor not available: ", cname, call. = FALSE)
    }
    return(structure(
      list(id = "blosc", cname = cname, clevel = as.integer(clevel), shuffle = shuffle_code, blocksize = 0L),
      class = "zarr_compressor"
    ))
  }
  if (id == "lz4") {
    return(structure(list(id = "lz4", acceleration = 1L), class = "zarr_compressor"))
  }
  structure(list(id = id, level = as.integer(level)), class = "zarr_compressor")
}

default_compressor <- function() {
  zarr_compressor("blosc", cname = "zstd", clevel = 7L, shuffle = "auto")
}

infer_dtype <- function(data) {
  if (is.character(data)) {
    return("string")
  }
  if (is.logical(data)) {
    return("bool")
  }
  if (is.integer(data)) {
    return("int32")
  }
  "float64"
}

dtype_v2_string <- function(dtype) {
  switch(dtype,
    int8 = "|i1", int16 = "<i2", int32 = "<i4", int64 = "<i8",
    uint8 = "|u1", uint16 = "<u2", uint32 = "<u4", uint64 = "<u8",
    float32 = "<f4", float64 = "<f8", bool = "|b1", string = "|O",
    stop("Unsupported dtype: ", dtype, call. = FALSE)
  )
}

compressor_json <- function(compressor, dtype) {
  if (is.null(compressor)) {
    return(NULL)
  }
  c <- unclass(compressor)
  if (identical(c$id, "blosc") && identical(c$shuffle, -1L)) {
    # Resolve numcodecs' AUTOSHUFFLE so the metadata is explicit.
    c$shuffle <- if (dtype %in% c("int8", "uint8", "bool", "string")) 2L else 1L
  }
  c
}

zarr_create_array <- function(group, name, shape, dtype, chunks = NULL,
                              compressor = default_compressor(), fill_value = NULL,
                              dims = NULL, attrs = NULL) {
  group <- as_zarr_group(group)
  shape <- as.numeric(shape)
  if (is.null(chunks)) {
    chunks <- pmax(shape, 1)
  }
  chunks <- pmax(as.numeric(chunks), 1)
  if (length(chunks) != length(shape)) {
    stop("chunks must have one entry per dimension", call. = FALSE)
  }
  dir <- file.path(group$path, name)
  if (dir.exists(dir)) {
    unlink(dir, recursive = TRUE)
  }
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  is_string <- dtype == "string"
  meta <- list(
    chunks = as.list(chunks),
    compressor = compressor_json(compressor, dtype),
    dimension_separator = ".",
    dtype = dtype_v2_string(dtype),
    fill_value = if (is_string) "" else fill_value,
    filters = if (is_string) list(list(id = "vlen-utf8")) else NULL,
    order = "C",
    shape = as.list(shape),
    zarr_format = 2L
  )
  write_json_file(meta, file.path(dir, ".zarray"))
  a <- attrs %||% list()
  if (!is.null(dims)) {
    a[["_ARRAY_DIMENSIONS"]] <- as.list(dims)
  }
  write_json_file(if (length(a) == 0) structure(list(), names = character(0)) else a,
    file.path(dir, ".zattrs"))
  invisible(zarr_handle(group, name))
}

#' @rdname zarr-write
#' @export
zarr_write_array <- function(group, name, data, dtype = NULL, chunks = NULL,
                             compressor = default_compressor(), fill_value = NULL,
                             dims = NULL, attrs = NULL) {
  group <- as_zarr_group(group)
  if (is.null(dtype)) {
    dtype <- infer_dtype(data)
  }
  shape <- if (is.null(dim(data))) length(data) else dim(data)
  if (is.null(fill_value) && dtype != "string") {
    fill_value <- 0
  }
  h <- zarr_create_array(group, name, shape, dtype,
    chunks = chunks, compressor = compressor,
    fill_value = fill_value, dims = dims, attrs = attrs
  )
  if (prod(shape) > 0) {
    if (is.factor(data)) {
      data <- as.character(data)
    }
    .zarr_write(h$ptr, data, rep(0, length(shape)), shape)
  }
  invisible(group)
}

#' @rdname zarr-write
#' @export
zarr_consolidate <- function(group) {
  group <- as_zarr_group(group)
  md <- list()
  root_files <- c(".zgroup", ".zattrs")
  for (f in root_files) {
    p <- file.path(group$path, f)
    if (file.exists(p)) {
      md[[f]] <- read_json_file(p)
    }
  }
  for (name in zarr_list(group)) {
    for (f in c(".zarray", ".zattrs")) {
      p <- file.path(group$path, name, f)
      if (file.exists(p)) {
        md[[paste0(name, "/", f)]] <- read_json_file(p)
      }
    }
  }
  write_json_file(list(metadata = md, zarr_consolidated_format = 1L), file.path(group$path, ".zmetadata"))
  invisible(group)
}
