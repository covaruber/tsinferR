#' Convert a VCF file to a VCF Zarr (VCZ) store
#'
#' Converts a (optionally bgzip/gzip-compressed) VCF file into a VCZ store
#' following the VCF Zarr specification used by 'bio2zarr', ready for use as
#' a [tsinfer_source()]. Parsing and chunk compression are done in compiled
#' code, so large files can be converted without loading them into memory:
#' genotypes are streamed in chunks of `variants_chunk_size` variants.
#'
#' The store contains `call_genotype`, `call_genotype_phased`,
#' `call_genotype_mask`, `variant_position`, `variant_contig`,
#' `variant_allele`, `variant_id`, `variant_id_mask`, `variant_quality`,
#' `variant_filter`, `variant_length`, `filter_id`, `filter_description`,
#' `contig_id`, `contig_length` and `sample_id`, plus one `variant_<NAME>`
#' array per requested INFO field.
#'
#' Since 'tsinfer' infers one contig at a time, use `contig` to select a
#' single chromosome from a multi-contig VCF. Contigs without a length in the
#' header get the largest position seen plus one.
#'
#' @param vcf Path to the VCF file (`.vcf` or `.vcf.gz`).
#' @param vcz Path of the output store.
#' @param contig Optional name of the single contig to convert.
#' @param info_fields Character vector of INFO fields to store. Fields with
#'   `Number=1` are stored with their declared type (Integer, Float, String),
#'   flags as booleans, and other fields as strings.
#' @param ancestral_info Optional INFO field holding the ancestral allele
#'   (e.g. `"AA"`), stored upper-cased as `variant_ancestral_allele` with any
#'   `|`-separated suffix removed. Use it with
#'   `ancestral_state(vcz, "variant_ancestral_allele")`.
#' @param variants_chunk_size,samples_chunk_size Chunk sizes.
#' @param compressor A [zarr_compressor()].
#' @param num_threads Threads used to compress chunks.
#' @param overwrite Replace an existing store.
#' @return The path of the store, invisibly.
#' @export
#' @examples
#' vcf <- tempfile(fileext = ".vcf")
#' writeLines(c(
#'   "##fileformat=VCFv4.2",
#'   "##contig=<ID=1,length=1000>",
#'   "##INFO=<ID=AA,Number=1,Type=String,Description=\"Ancestral allele\">",
#'   "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB",
#'   "1\t10\t.\tA\tG\t.\tPASS\tAA=A\tGT\t0|1\t1|1",
#'   "1\t20\t.\tC\tT\t.\tPASS\tAA=t\tGT\t0|0\t1|0"
#' ), vcf)
#' vcz <- vcf_to_vcz(vcf, tempfile(fileext = ".vcz"), ancestral_info = "AA")
#' zarr_read(vcz, "call_genotype")
#' zarr_read(vcz, "variant_ancestral_allele")
vcf_to_vcz <- function(vcf, vcz, contig = NULL, info_fields = NULL,
                       ancestral_info = NULL, variants_chunk_size = 10000L,
                       samples_chunk_size = 1000L, compressor = default_compressor(),
                       num_threads = 1L, overwrite = FALSE) {
  vcf <- normalizePath(vcf, mustWork = TRUE)
  hdr <- .vcf_header(vcf)
  num_samples <- length(hdr$samples)
  if (num_samples == 0) {
    stop("The VCF file has no samples", call. = FALSE)
  }
  ploidy <- max(1L, hdr$ploidy)
  fields <- unique(c(info_fields, ancestral_info))
  unknown <- setdiff(fields, hdr$info_id)
  if (length(unknown) > 0) {
    warning("INFO field(s) not declared in the header: ", paste(unknown, collapse = ", "),
      call. = FALSE
    )
  }
  filter_ids <- unique(c("PASS", hdr$filter_id))
  filter_desc <- hdr$filter_description[match(filter_ids, hdr$filter_id)]
  filter_desc[is.na(filter_desc)] <- "All filters passed"

  g <- zarr_create_group(vcz,
    overwrite = overwrite,
    attrs = list(vcf_zarr_version = "0.4", source = "tsinferR")
  )
  vchunk <- as.integer(variants_chunk_size)
  schunk <- as.integer(min(samples_chunk_size, num_samples))
  meta_for <- function(name, dtype, chunks, fill) {
    dir.create(file.path(g$path, name), showWarnings = FALSE)
    list(
      path = file.path(g$path, name), zarr_format = 2L, shape = c(0, chunks[-1]),
      chunks = chunks, dtype = dtype, itemsize = 1, big_endian = FALSE,
      serializer = "raw", codecs = if (is.null(compressor)) list() else list(compressor_json(compressor, dtype)),
      key_prefix = "", separator = ".", fill_value = fill, fill_string = ""
    )
  }
  gt_meta <- meta_for("call_genotype", "int8", c(vchunk, schunk, ploidy), -1)
  ph_meta <- meta_for("call_genotype_phased", "bool", c(vchunk, schunk), 0)
  mk_meta <- meta_for("call_genotype_mask", "bool", c(vchunk, schunk, ploidy), 1)
  res <- .vcf_convert(
    vcf, gt_meta, ph_meta, mk_meta, num_samples, ploidy,
    if (is.null(contig)) "" else contig, as.character(fields), filter_ids,
    as.integer(num_threads)
  )
  n <- res$num_variants
  if (!is.null(contig) && n == 0) {
    stop("No variants found for contig '", contig, "'", call. = FALSE)
  }
  write_array_metadata(g, "call_genotype", c(n, num_samples, ploidy), c(vchunk, schunk, ploidy),
    "int8", compressor, -1, c("variants", "samples", "ploidy")
  )
  write_array_metadata(g, "call_genotype_phased", c(n, num_samples), c(vchunk, schunk),
    "bool", compressor, FALSE, c("variants", "samples")
  )
  write_array_metadata(g, "call_genotype_mask", c(n, num_samples, ploidy), c(vchunk, schunk, ploidy),
    "bool", compressor, TRUE, c("variants", "samples", "ploidy")
  )

  contigs <- hdr$contig_id
  lengths <- hdr$contig_length
  extra <- setdiff(res$contig_names, contigs)
  contigs <- c(contigs, extra)
  lengths <- c(lengths, rep(NA_real_, length(extra)))
  if (!is.null(contig)) {
    keep <- contigs == contig
    contigs <- contigs[keep]
    lengths <- lengths[keep]
  }
  vcontig <- match(res$contig_names, contigs)[res$contig + 1] - 1L
  for (j in which(is.na(lengths))) {
    p <- res$position[vcontig == j - 1]
    lengths[j] <- if (length(p) > 0) max(p) + 1 else 0
  }
  vc <- c(vchunk)
  zarr_write_array(g, "variant_position", as.integer(res$position),
    chunks = vc, compressor = compressor, dims = "variants"
  )
  zarr_write_array(g, "variant_contig", vcontig,
    dtype = if (length(contigs) > 127) "int16" else "int8",
    chunks = vc, compressor = compressor, dims = "variants"
  )
  zarr_write_array(g, "variant_allele", res$alleles,
    chunks = c(vchunk, ncol(res$alleles)), compressor = compressor,
    dims = c("variants", "alleles")
  )
  zarr_write_array(g, "variant_length", nchar(res$alleles[, 1]),
    chunks = vc, compressor = compressor, dims = "variants"
  )
  zarr_write_array(g, "variant_id", res$id, chunks = vc, compressor = compressor, dims = "variants")
  zarr_write_array(g, "variant_id_mask", !nzchar(res$id),
    chunks = vc, compressor = compressor, dims = "variants"
  )
  zarr_write_array(g, "variant_quality", res$quality,
    dtype = "float32", chunks = vc,
    compressor = compressor, fill_value = "NaN", dims = "variants"
  )
  zarr_write_array(g, "variant_filter", res$filters,
    chunks = c(vchunk, length(filter_ids)),
    compressor = compressor, dims = c("variants", "filters")
  )
  zarr_write_array(g, "filter_id", filter_ids, compressor = compressor, dims = "filters")
  zarr_write_array(g, "filter_description", filter_desc, compressor = compressor, dims = "filters")
  zarr_write_array(g, "contig_id", contigs, compressor = compressor, dims = "contigs")
  zarr_write_array(g, "contig_length", lengths,
    dtype = "int64", compressor = compressor,
    dims = "contigs"
  )
  zarr_write_array(g, "sample_id", hdr$samples, compressor = compressor, dims = "samples")
  if (n > 0) {
    # Per-chunk genomic extent, used by 'vcztools' for region queries.
    chunk <- (seq_len(n) - 1) %/% vchunk
    vend <- res$position + nchar(res$alleles[, 1]) - 1
    key <- paste(chunk, vcontig)
    first <- !duplicated(key)
    region_index <- cbind(
      chunk[first], vcontig[first],
      tapply(res$position, key, min)[key[first]],
      tapply(res$position, key, max)[key[first]],
      tapply(vend, key, max)[key[first]],
      tabulate(match(key, key[first]), sum(first))
    )
    storage.mode(region_index) <- "integer"
    dimnames(region_index) <- NULL
    zarr_write_array(g, "region_index", region_index,
      compressor = compressor,
      dims = c("region_index_values", "region_index_fields")
    )
  }
  for (field in info_fields) {
    write_info_field(g, field, res$info[[field]], hdr, vchunk, compressor)
  }
  if (!is.null(ancestral_info)) {
    aa <- res$info[[ancestral_info]]
    aa[aa %in% c("\x01", "\x02", ".")] <- ""
    aa <- toupper(sub("\\|.*$", "", aa))
    zarr_write_array(g, "variant_ancestral_allele", aa,
      chunks = vc, compressor = compressor, dims = "variants"
    )
  }
  zarr_consolidate(g)
  invisible(g$path)
}

write_info_field <- function(g, field, values, hdr, vchunk, compressor) {
  k <- match(field, hdr$info_id)
  type <- if (is.na(k)) "String" else hdr$info_type[k]
  number <- if (is.na(k)) "." else hdr$info_number[k]
  missing <- values == "\x02" | values == "."
  name <- paste0("variant_", field)
  if (type == "Flag") {
    zarr_write_array(g, name, values == "\x01", chunks = vchunk, compressor = compressor, dims = "variants")
  } else if (number == "1" && type == "Integer") {
    v <- suppressWarnings(as.integer(values))
    v[missing | is.na(v)] <- -1L
    zarr_write_array(g, name, v, chunks = vchunk, compressor = compressor, fill_value = -1, dims = "variants")
  } else if (number == "1" && type == "Float") {
    v <- suppressWarnings(as.numeric(values))
    v[missing] <- NaN
    zarr_write_array(g, name, v,
      dtype = "float32", chunks = vchunk, compressor = compressor,
      fill_value = "NaN", dims = "variants"
    )
  } else {
    values[missing | values == "\x01"] <- ""
    zarr_write_array(g, name, values, chunks = vchunk, compressor = compressor, dims = "variants")
  }
}
