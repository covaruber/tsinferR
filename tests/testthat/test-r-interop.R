test_that("tsinferR and tskitR exchange standard tree sequences", {
  skip_if_not_installed("tskitR")
  ts <- ts_load(fixture_path("case1_final.trees"))
  other <- as_tskit(ts)
  expect_s3_class(other, "TreeSequence")
  expect_equal(other$num_sites, ts_info(ts)$num_sites)
  expect_equal(other$num_nodes, ts_info(ts)$num_nodes)
  back <- as_tsinfer(other)
  expect_true(ts_equals(ts, back))
})

test_that("tsinferR and tsdateR exchange and date tree sequences", {
  skip_if_not_installed("tsdateR")
  ts <- ts_load(fixture_path("case1_final.trees"))
  dated <- date_tsinfer(ts, mutation_rate = 1e-8)
  expect_s3_class(dated, "tsinfer_ts")
  expect_equal(ts_info(dated)$num_sites, ts_info(ts)$num_sites)
  expect_equal(ncol(ts_genotype_matrix(dated)), ncol(ts_genotype_matrix(ts)))
  expect_true(all(is.finite(ts_nodes(dated)$time)))
})

test_that("citation points to the canonical tsinfer paper", {
  cited <- citation("tsinferR")
  expect_true(any(vapply(cited, function(entry) {
    inherits(entry, "bibentry") && grepl("Kelleher2019", entry$key, fixed = TRUE)
  }, logical(1))))
  expect_true(any(grepl(
    "10.1038/s41588-019-0483-y",
    capture.output(print(cited)), fixed = TRUE
  )))
})
