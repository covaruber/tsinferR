test_that("VCF conversion reproduces bio2zarr output", {
  b <- zarr_open(fixture_path("bio2zarr.vcz"))
  out <- file.path(tempdir(), "vcf-b2z-same.vcz")
  vcf_to_vcz(fixture_path("bio2zarr.vcf"), out, overwrite = TRUE)
  a <- zarr_open(out)
  for (name in zarr_list(b)) {
    expect_equal(zarr_read(a, name), zarr_read(b, name), info = name)
  }
  out <- file.path(tempdir(), "vcf-b2z.vcz")
  vcf_to_vcz(fixture_path("bio2zarr.vcf"), out,
    variants_chunk_size = 7, samples_chunk_size = 2,
    overwrite = TRUE
  )
  a <- zarr_open(out)
  for (name in setdiff(zarr_list(b), "region_index")) {
    expect_equal(zarr_read(a, name), zarr_read(b, name), info = name)
  }
  expect_equal(nrow(zarr_read(a, "region_index")), 14)
  expect_equal(zarr_info(a, "call_genotype")$chunks, c(7, 2, 2))
  expect_true(file.exists(file.path(out, ".zmetadata")))
})

test_that("gzip-compressed VCFs, INFO fields and missing data are handled", {
  vcf <- tempfile(fileext = ".vcf.gz")
  con <- gzfile(vcf, "w")
  writeLines(c(
    "##fileformat=VCFv4.2",
    "##contig=<ID=chr1,length=500>",
    "##contig=<ID=chr2>",
    "##FILTER=<ID=q10,Description=\"Quality below 10\">",
    "##INFO=<ID=AA,Number=1,Type=String,Description=\"Ancestral allele\">",
    "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">",
    "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Freq\">",
    "##INFO=<ID=DB,Number=0,Type=Flag,Description=\"dbSNP\">",
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\ts3",
    "chr1\t5\trs1\tA\tC,G\t30\tPASS\tAA=c|||;DP=10;AF=0.1,0.2;DB\tGT:DP\t0|1:3\t2/.:4\t.:5",
    "chr1\t9\t.\tT\t.\t.\tq10\tDP=7\tGT\t0|0\t0|0\t0",
    "chr2\t3\t.\tG\tA\t12.5\t.\t.\tDP:GT\t1:1|1\t2:0|1\t3:./."
  ), con)
  close(con)
  out <- vcf_to_vcz(vcf, file.path(tempdir(), "gz.vcz"),
    info_fields = c("DP", "AF", "DB"),
    ancestral_info = "AA", overwrite = TRUE, variants_chunk_size = 2
  )
  g <- zarr_open(out)
  gt <- zarr_read(g, "call_genotype")
  expect_equal(gt[1, , ], rbind(c(0L, 1L), c(2L, -1L), c(-1L, -2L)))
  expect_equal(gt[3, , ], rbind(c(1L, 1L), c(0L, 1L), c(-1L, -1L)))
  expect_equal(zarr_read(g, "call_genotype_phased")[1, ], c(TRUE, FALSE, FALSE))
  expect_equal(zarr_read(g, "variant_allele"), rbind(c("A", "C", "G"), c("T", "", ""), c("G", "A", "")))
  expect_equal(zarr_read(g, "variant_DP"), c(10L, 7L, -1L))
  expect_equal(zarr_read(g, "variant_AF"), c("0.1,0.2", "", ""))
  expect_equal(zarr_read(g, "variant_DB"), c(TRUE, FALSE, FALSE))
  expect_equal(zarr_read(g, "variant_ancestral_allele"), c("C", "", ""))
  expect_equal(zarr_read(g, "contig_id"), c("chr1", "chr2"))
  expect_equal(zarr_read(g, "contig_length"), c(500, 4))
  expect_equal(zarr_read(g, "variant_contig"), c(0L, 0L, 1L))
  expect_equal(zarr_read(g, "variant_quality"), c(30, NaN, 12.5))
  expect_equal(zarr_read(g, "filter_id"), c("PASS", "q10"))
  expect_equal(zarr_read(g, "variant_filter"), rbind(c(TRUE, FALSE), c(FALSE, TRUE), c(FALSE, FALSE)))

  one <- vcf_to_vcz(vcf, file.path(tempdir(), "chr2.vcz"), contig = "chr2", overwrite = TRUE)
  expect_equal(zarr_read(one, "variant_position"), 3L)
  expect_equal(zarr_read(one, "contig_id"), "chr2")
  expect_error(vcf_to_vcz(vcf, file.path(tempdir(), "none.vcz"), contig = "chrX", overwrite = TRUE),
    "No variants")
  expect_error(vcf_to_vcz(vcf, out), "already exists")
})

test_that("malformed VCFs give informative errors", {
  vcf <- tempfile(fileext = ".vcf")
  writeLines(c(
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1",
    "1\t5\t.\tA\tC\t.\t.\t.\tGT\t0|x"
  ), vcf)
  expect_error(vcf_to_vcz(vcf, tempfile()), "bad genotype")
  writeLines(c("##fileformat=VCFv4.2", "1\t5\t.\tA\tC\t.\t.\t."), vcf)
  expect_error(vcf_to_vcz(vcf, tempfile()), "#CHROM")
})

test_that("inference from a converted VCF matches inference from bio2zarr", {
  dir <- tempfile("vcfrun-")
  dir.create(dir)
  vcz <- vcf_to_vcz(fixture_path("bio2zarr.vcf"), file.path(dir, "s.vcz"), variants_chunk_size = 9)
  make_cfg <- function(store, anc) {
    tsinfer_config(
      sources = list(tsinfer_source(store, "s")),
      ancestral_state = ancestral_state(store, is_reference = TRUE),
      ancestors = list(ancestors_config("ancestors", anc, "s")),
      match = match_config(sources = list(
        ancestors = match_source_config(node_flags = 0, create_individuals = FALSE),
        s = match_source_config()
      )),
      post_process = post_process_config()
    )
  }
  a <- run_tsinfer(make_cfg(vcz, file.path(dir, "a1.vcz")))
  b <- run_tsinfer(make_cfg(fixture_path("bio2zarr.vcz"), file.path(dir, "a2.vcz")))
  expect_true(ts_equals(a, b))
  G <- ts_genotype_matrix(a)
  expect_equal(ncol(G), 10)
})
