test_that("Zarr v2 stores written by zarr-python are read correctly", {
  g <- zarr_open(fixture_path("case1.vcz"))
  expect_true(all(c("call_genotype", "variant_position", "sample_id") %in% zarr_list(g)))
  gt <- zarr_read(g, "call_genotype")
  pos <- zarr_read(g, "variant_position")
  ids <- zarr_read(g, "sample_id")
  expect_equal(dim(gt), c(length(pos), 16L, 1L))
  expect_equal(ids, paste0("s", 0:15))
  # The genotypes reproduce the simulated tree sequence's sites.
  final <- ts_load(fixture_path("case1_final.trees"))
  sites <- ts_sites(final)
  expect_true(all(sites$position %in% pos))
  alleles <- zarr_read(g, "variant_allele")
  expect_type(alleles, "character")
  expect_equal(nrow(alleles), length(pos))
  expect_equal(zarr_read(g, "variant_position", rows = c(3, 1)), pos[c(3, 1)])
  sub <- zarr_read(g, "call_genotype", rows = 2:4)
  expect_equal(sub, gt[2:4, , , drop = FALSE])
})

test_that("chunked ancestor stores are read correctly", {
  g <- zarr_open(fixture_path("case1_anc.vcz"))
  info <- zarr_info(g, "call_genotype")
  expect_equal(info$chunks, c(13, 7, 1))
  gt <- zarr_read(g, "call_genotype")
  expect_equal(dim(gt), as.integer(info$shape))
  expect_true(all(gt %in% c(-1L, 0L, 1L)))
  times <- zarr_read(g, "sample_time")
  expect_type(times, "double")
  expect_equal(length(times), info$shape[2])
})

test_that("bio2zarr stores (blosc zstd + bitshuffle, missing chunks) are read", {
  g <- zarr_open(fixture_path("bio2zarr.vcz"))
  pos <- zarr_read(g, "variant_position")
  expect_true(all(diff(pos) > 0))
  gt <- zarr_read(g, "call_genotype")
  expect_equal(dim(gt)[2:3], c(5L, 2L))
  mask <- zarr_read(g, "call_genotype_mask")
  expect_type(mask, "logical")
  expect_false(any(mask))
  phased <- zarr_read(g, "call_genotype_phased")
  expect_true(all(phased))
  expect_equal(zarr_read(g, "contig_id"), "1")
  vcf <- utils::read.table(fixture_path("bio2zarr.vcf"), comment.char = "#", sep = "\t")
  expect_equal(pos, vcf$V2)
})

test_that("legacy fixed-width string stores are read", {
  g <- zarr_open(fixture_path("legacy.vcz"))
  expect_equal(zarr_read(g, "sample_id"), c("S0", "S1", "S2"))
  alleles <- zarr_read(g, "variant_allele")
  expect_equal(dim(alleles), c(8L, 2L))
  expect_type(zarr_read(g, "call_genotype_phased"), "integer")
})

test_that("arrays round-trip for all data types and codecs", {
  dir <- file.path(tempdir(), "roundtrip.zarr")
  g <- zarr_create_group(dir, overwrite = TRUE)
  set.seed(1)
  m <- matrix(sample(-1:3, 60, replace = TRUE), 6, 10)
  comps <- list(
    NULL,
    zarr_compressor("blosc", cname = "lz4", shuffle = "shuffle"),
    zarr_compressor("blosc", cname = "zstd", shuffle = "bitshuffle"),
    zarr_compressor("blosc", cname = "zlib"),
    zarr_compressor("blosc", cname = "blosclz"),
    zarr_compressor("blosc", cname = "lz4hc"),
    zarr_compressor("zstd"), zarr_compressor("zlib"),
    zarr_compressor("gzip"), zarr_compressor("lz4")
  )
  for (k in seq_along(comps)) {
    for (dtype in c("int8", "int16", "int32", "int64", "float32", "float64")) {
      name <- paste0("a", k, dtype)
      zarr_write_array(g, name, m, dtype = dtype, chunks = c(4, 3), compressor = comps[[k]])
      x <- zarr_read(g, name)
      expect_equal(as.numeric(x), as.numeric(m), info = name)
    }
  }
  zarr_write_array(g, "u8", matrix(0:11, 3), dtype = "uint8")
  expect_equal(zarr_read(g, "u8"), matrix(0:11, 3))
  b <- c(TRUE, FALSE, TRUE, TRUE, FALSE)
  zarr_write_array(g, "b", b, chunks = 2)
  expect_equal(zarr_read(g, "b"), b)
  s <- matrix(c("A", "", "TTT", "\u00e9x", "C", "G"), 3)
  zarr_write_array(g, "s", s, chunks = c(2, 1))
  expect_equal(zarr_read(g, "s"), s)
  arr <- array(1:24, c(2, 3, 4))
  zarr_write_array(g, "arr", arr, chunks = c(1, 2, 3))
  expect_equal(zarr_read(g, "arr"), arr)
  expect_equal(zarr_read(g, "arr", rows = 2), arr[2, , , drop = FALSE])
  zarr_write_array(g, "empty", integer(0))
  expect_equal(zarr_read(g, "empty"), integer(0))
  zarr_consolidate(g)
  expect_true(file.exists(file.path(dir, ".zmetadata")))
})

test_that("Zarr v3 arrays are read", {
  dir <- file.path(tempdir(), "v3.zarr")
  unlink(dir, recursive = TRUE)
  dir.create(file.path(dir, "x", "c", "0"), recursive = TRUE)
  writeLines('{"zarr_format":3,"node_type":"group","attributes":{}}', file.path(dir, "zarr.json"))
  writeLines(paste0(
    '{"zarr_format":3,"node_type":"array","shape":[2,3],"data_type":"int16",',
    '"chunk_grid":{"name":"regular","configuration":{"chunk_shape":[2,3]}},',
    '"chunk_key_encoding":{"name":"default","configuration":{"separator":"/"}},',
    '"fill_value":0,"codecs":[{"name":"bytes","configuration":{"endian":"little"}}],',
    '"dimension_names":["a","b"]}'
  ), file.path(dir, "x", "zarr.json"))
  con <- file(file.path(dir, "x", "c", "0", "0"), "wb")
  writeBin(as.integer(0:5), con, size = 2, endian = "little")
  close(con)
  g <- zarr_open(dir)
  expect_equal(g$zarr_format, 3L)
  expect_equal(zarr_list(g), "x")
  expect_equal(zarr_read(g, "x"), matrix(0:5, 2, byrow = TRUE))
  expect_equal(zarr_info(g, "x")$attributes[["_ARRAY_DIMENSIONS"]], list("a", "b"))
})

test_that("Zarr errors are informative", {
  expect_error(zarr_open(file.path(tempdir(), "does-not-exist")), "not found")
  g <- zarr_create_group(file.path(tempdir(), "err.zarr"), overwrite = TRUE)
  expect_error(zarr_read(g, "nope"), "not found")
  expect_error(zarr_compressor("blosc", cname = "nope"), "not available")
})
