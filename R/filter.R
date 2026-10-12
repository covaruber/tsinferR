# Variant filtering for VCZ sources: a bcftools-style expression language
# (the subset supported by 'vcztools'), regions/targets and sample selection.

parse_samples_spec <- function(samples, all_samples) {
  if (is.null(samples)) {
    return(seq_along(all_samples) - 1L)
  }
  exclude <- startsWith(samples, "^")
  ids <- strsplit(sub("^\\^", "", samples), ",", fixed = TRUE)[[1]]
  ids <- ids[nzchar(ids)]
  unknown <- setdiff(ids, all_samples)
  if (length(unknown) > 0) {
    warning("subset called for sample(s) not in header: ", paste(unknown, collapse = ","),
      call. = FALSE
    )
    ids <- setdiff(ids, unknown)
  }
  sel <- match(ids, all_samples) - 1L
  if (exclude) {
    sel <- setdiff(seq_along(all_samples) - 1L, sel)
  }
  sel
}

parse_region_strings <- function(regions) {
  parts <- strsplit(regions, ",", fixed = TRUE)[[1]]
  lapply(parts, function(r) {
    if (grepl(":[0-9]+-[0-9]*$", r)) {
      contig <- sub(":[0-9]+-[0-9]*$", "", r)
      se <- sub("^.*:", "", r)
      start <- as.numeric(sub("-.*$", "", se))
      end_s <- sub("^[0-9]+-", "", se)
      end <- if (nzchar(end_s)) as.numeric(end_s) else Inf
      list(contig = contig, start = start, end = end)
    } else if (grepl(":[0-9]+$", r)) {
      contig <- sub(":[0-9]+$", "", r)
      p <- as.numeric(sub("^.*:", "", r))
      list(contig = contig, start = p, end = p)
    } else {
      list(contig = r, start = 1, end = Inf)
    }
  })
}

# Logical mask of variants overlapping the regions (VCF 1-based, closed).
region_mask <- function(regions, contig_ids, variant_contig, pos, len, point = FALSE) {
  mask <- rep(FALSE, length(pos))
  for (r in parse_region_strings(regions)) {
    ci <- match(r$contig, contig_ids) - 1L
    if (is.na(ci)) {
      next
    }
    vstart <- pos
    vend <- if (point) pos else pos + len - 1
    mask <- mask | (variant_contig == ci & vend >= r$start & vstart <= r$end)
  }
  mask
}

# --- Expression language -----------------------------------------------------

tokenize_filter <- function(expr) {
  tokens <- list()
  i <- 1
  n <- nchar(expr)
  chars <- strsplit(expr, "")[[1]]
  ops <- c("&&", "||", "==", "!=", ">=", "<=", "!~", "=", "<", ">", "&", "|", "+", "-", "*", "/", "(", ")", "~", "!")
  while (i <= n) {
    ch <- chars[i]
    if (grepl("\\s", ch)) {
      i <- i + 1
      next
    }
    if (ch == "\"" || ch == "'") {
      j <- i + 1
      while (j <= n && chars[j] != ch) j <- j + 1
      if (j > n) {
        stop("Unterminated string in filter expression: ", expr, call. = FALSE)
      }
      tokens[[length(tokens) + 1]] <- list(type = "str", value = paste(chars[seq_len(j - i - 1) + i], collapse = ""))
      i <- j + 1
      next
    }
    if (grepl("[0-9.]", ch)) {
      m <- regmatches(substring(expr, i), regexpr("^[0-9]*\\.?[0-9]+([eE][-+]?[0-9]+)?", substring(expr, i)))
      tokens[[length(tokens) + 1]] <- list(type = "num", value = as.numeric(m))
      i <- i + nchar(m)
      next
    }
    if (grepl("[A-Za-z_@]", ch)) {
      m <- regmatches(substring(expr, i), regexpr("^[A-Za-z_@][A-Za-z0-9_./]*", substring(expr, i)))
      if (startsWith(m, "@")) {
        stop("Unsupported filtering feature: file references", call. = FALSE)
      }
      tokens[[length(tokens) + 1]] <- list(type = "id", value = m)
      i <- i + nchar(m)
      if (i <= n && chars[i] == "[") {
        stop("Unsupported filtering feature: array subscripts", call. = FALSE)
      }
      next
    }
    two <- substr(expr, i, i + 1)
    if (two %in% ops) {
      tokens[[length(tokens) + 1]] <- list(type = "op", value = two)
      i <- i + 2
      next
    }
    if (ch %in% ops) {
      tokens[[length(tokens) + 1]] <- list(type = "op", value = ch)
      i <- i + 1
      next
    }
    stop("Unexpected character '", ch, "' in filter expression: ", expr, call. = FALSE)
  }
  tokens
}

FILTER_PRECEDENCE <- list(
  "||" = 1, "|" = 2, "&&" = 3, "&" = 4,
  "=" = 5, "==" = 5, "!=" = 5, "<" = 5, ">" = 5, "<=" = 5, ">=" = 5, "~" = 5, "!~" = 5,
  "+" = 6, "-" = 6, "*" = 7, "/" = 7
)

parse_filter <- function(expr) {
  tokens <- tokenize_filter(expr)
  pos <- 1
  peek <- function() if (pos <= length(tokens)) tokens[[pos]] else NULL
  advance <- function() {
    t <- tokens[[pos]]
    pos <<- pos + 1
    t
  }
  parse_primary <- function() {
    t <- peek()
    if (is.null(t)) {
      stop("Unexpected end of filter expression: ", expr, call. = FALSE)
    }
    advance()
    if (t$type == "op" && t$value == "(") {
      node <- parse_expr(0)
      close <- peek()
      if (is.null(close) || close$value != ")") {
        stop("Missing ')' in filter expression: ", expr, call. = FALSE)
      }
      advance()
      return(node)
    }
    if (t$type == "op" && t$value == "-") {
      return(list(type = "neg", arg = parse_primary()))
    }
    if (t$type == "num") {
      return(list(type = "const", value = t$value))
    }
    if (t$type == "str") {
      if (t$value == ".") {
        stop("Unsupported filtering feature: missing data", call. = FALSE)
      }
      return(list(type = "const", value = t$value))
    }
    if (t$type == "id") {
      nxt <- peek()
      if (!is.null(nxt) && nxt$type == "op" && nxt$value == "(") {
        stop("Unsupported filtering feature: functions", call. = FALSE)
      }
      if (t$value == "GT") {
        stop("Unsupported filtering feature: genotype values", call. = FALSE)
      }
      return(list(type = "field", name = t$value))
    }
    stop("Unexpected token '", t$value, "' in filter expression: ", expr, call. = FALSE)
  }
  parse_expr <- function(min_prec) {
    lhs <- parse_primary()
    while (TRUE) {
      t <- peek()
      if (is.null(t) || t$type != "op" || is.null(FILTER_PRECEDENCE[[t$value]])) {
        break
      }
      prec <- FILTER_PRECEDENCE[[t$value]]
      if (prec < min_prec) {
        break
      }
      advance()
      rhs <- parse_expr(prec + 1)
      lhs <- list(type = "binary", op = t$value, lhs = lhs, rhs = rhs)
    }
    lhs
  }
  ast <- parse_expr(0)
  if (pos <= length(tokens)) {
    stop("Unexpected token '", tokens[[pos]]$value, "' in filter expression: ", expr, call. = FALSE)
  }
  ast
}

RESERVED_VCF_FIELDS <- c(
  CHROM = "variant_contig", POS = "variant_position", ID = "variant_id",
  REF = "variant_allele", ALT = "variant_allele", QUAL = "variant_quality",
  FILTER = "variant_filter"
)

vcz_field_for <- function(name, arrays) {
  if (name %in% names(RESERVED_VCF_FIELDS)) {
    return(RESERVED_VCF_FIELDS[[name]])
  }
  split <- strsplit(name, "/", fixed = TRUE)[[1]]
  base <- split[length(split)]
  candidates <- if (length(split) > 1) {
    if (split[1] %in% c("FORMAT", "FMT")) paste0("call_", base) else paste0("variant_", base)
  } else {
    c(paste0("call_", base), paste0("variant_", base))
  }
  found <- intersect(candidates, arrays)
  if (length(found) == 0) {
    stop("the tag \"", name, "\" is not defined", call. = FALSE)
  }
  if (length(found) > 1) {
    stop("ambiguous filtering expression: \"", name, "\", both INFO/", base,
      " and FORMAT/", base, " are defined",
      call. = FALSE
    )
  }
  found
}

variant_type_snp <- function(alleles, any) {
  ref <- alleles[, 1]
  alt <- alleles[, -1, drop = FALSE]
  is_snp <- nchar(alt) == 1 & nchar(ref) == 1 & alt != ref & grepl("^[ACGTN]$", toupper(alt))
  present <- nzchar(alt) & alt != "."
  is_snp[!present] <- NA
  if (any) {
    apply(is_snp, 1, function(r) any(r, na.rm = TRUE))
  } else {
    apply(is_snp, 1, function(r) all(r, na.rm = TRUE) && any(!is.na(r)))
  }
}

# Evaluate a parsed expression on a VCZ group, returning a logical mask over
# variants (2-D sample-level results are reduced with any()).
eval_filter <- function(ast, group, arrays) {
  cache <- new.env()
  field <- function(name) {
    vcz <- vcz_field_for(name, arrays)
    if (is.null(cache[[vcz]])) {
      x <- zarr_read(group, vcz)
      if (startsWith(vcz, "call_") && length(dim(x)) > 2) {
        stop("Unsupported filtering feature: higher dimensional FORMAT fields", call. = FALSE)
      }
      cache[[vcz]] <- x
    }
    cache[[vcz]]
  }
  ev <- function(node) {
    switch(node$type,
      const = node$value,
      field = {
        if (node$name == "TYPE") {
          stop("TYPE must be compared with a string", call. = FALSE)
        }
        x <- field(node$name)
        if (node$name %in% c("REF", "ALT")) {
          stop("Unsupported filtering feature: REF/ALT values", call. = FALSE)
        }
        x
      },
      neg = -ev(node$arg),
      binary = {
        op <- node$op
        lhs_field <- if (node$lhs$type == "field") node$lhs$name else ""
        if (lhs_field == "TYPE") {
          if (node$rhs$type != "const" || !node$rhs$value %in% c("snp", "ref")) {
            stop("Unsupported filtering feature: TYPE values other than \"snp\" and \"ref\"",
              call. = FALSE
            )
          }
          alleles <- field("REF")
          res <- if (node$rhs$value == "ref") {
            rowSums(nzchar(alleles[, -1, drop = FALSE])) == 0
          } else {
            variant_type_snp(alleles, any = op %in% c("~", "!~"))
          }
          if (op %in% c("!=", "!~")) res <- !res
          return(res)
        }
        if (lhs_field == "CHROM") {
          ids <- zarr_read(group, "contig_id")
          ci <- match(node$rhs$value, ids) - 1L
          x <- field("CHROM")
          res <- if (is.na(ci)) rep(FALSE, length(x)) else x == ci
          if (op == "!=") res <- !res
          return(res)
        }
        if (lhs_field == "FILTER") {
          f <- field("FILTER")
          ids <- zarr_read(group, "filter_id")
          wanted <- strsplit(node$rhs$value, "[;,]")[[1]]
          has <- f[, match(wanted, ids), drop = FALSE]
          n_set <- rowSums(f)
          res <- switch(op,
            "=" = ,
            "==" = rowSums(has) == length(wanted) & n_set == length(wanted),
            "!=" = !(rowSums(has) == length(wanted) & n_set == length(wanted)),
            "~" = rowSums(has) == length(wanted),
            "!~" = rowSums(has) != length(wanted),
            stop("Unsupported FILTER operator: ", op, call. = FALSE)
          )
          return(res)
        }
        if (op %in% c("~", "!~")) {
          stop("Unsupported filtering feature: regular expressions", call. = FALSE)
        }
        a <- ev(node$lhs)
        b <- ev(node$rhs)
        switch(op,
          "+" = a + b, "-" = a - b, "*" = a * b, "/" = a / b,
          "=" = , "==" = a == b, "!=" = a != b,
          "<" = a < b, ">" = a > b, "<=" = a <= b, ">=" = a >= b,
          "&" = , "&&" = combine_and(a, b),
          "|" = , "||" = combine_or(a, b)
        )
      }
    )
  }
  res <- ev(ast)
  if (is.matrix(res)) {
    res <- apply(res, 1, any)
  }
  res[is.na(res)] <- FALSE
  as.logical(res)
}

combine_and <- function(a, b) {
  if (is.matrix(a) || is.matrix(b)) {
    va <- if (is.matrix(a)) apply(a, 1, any) else a
    vb <- if (is.matrix(b)) apply(b, 1, any) else b
    return(va & vb)
  }
  a & b
}

combine_or <- function(a, b) {
  if (is.matrix(a) || is.matrix(b)) {
    va <- if (is.matrix(a)) apply(a, 1, any) else a
    vb <- if (is.matrix(b)) apply(b, 1, any) else b
    return(va | vb)
  }
  a | b
}

# Positions of the variants of a source passing its filters.
filtered_positions <- function(src) {
  g <- src$group
  n <- length(src$positions)
  keep <- rep(TRUE, n)
  s <- src$config
  arrays <- zarr_list(g)
  if (!is.null(s$include)) {
    keep <- keep & eval_filter(parse_filter(s$include), g, arrays)
  }
  if (!is.null(s$exclude)) {
    keep <- keep & !eval_filter(parse_filter(s$exclude), g, arrays)
  }
  if (!is.null(s$regions) || !is.null(s$targets)) {
    contig_ids <- zarr_read(g, "contig_id")
    vcontig <- zarr_read(g, "variant_contig")
    if (!is.null(s$regions)) {
      len <- if ("variant_length" %in% arrays) {
        zarr_read(g, "variant_length")
      } else {
        nchar(src$alleles[, 1])
      }
      keep <- keep & region_mask(s$regions, contig_ids, vcontig, src$positions, len)
    }
    if (!is.null(s$targets)) {
      complement <- startsWith(s$targets, "^")
      m <- region_mask(sub("^\\^", "", s$targets), contig_ids, vcontig, src$positions, 1, point = TRUE)
      keep <- keep & (if (complement) !m else m)
    }
  }
  src$positions[keep]
}

source_has_filters <- function(s) {
  !is.null(s$include) || !is.null(s$exclude) || !is.null(s$regions) || !is.null(s$targets)
}
