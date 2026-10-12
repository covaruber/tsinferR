test_that("tree sequences round-trip through .trees files", {
  ts <- ts_load(fixture_path("case1_final.trees"))
  path <- tempfile(fileext = ".trees")
  ts_dump(ts, path)
  ts2 <- ts_load(path)
  expect_true(ts_equals(ts, ts2, ignore_provenance = FALSE))
  expect_s3_class(ts, "tsinfer_ts")
  expect_output(print(ts), "num_nodes")
  expect_gt(ts_num_trees(ts), 0)
})

test_that("table accessors decode metadata", {
  ts <- ts_load(fixture_path("case1_raw.trees"))
  md <- ts_metadata(ts)
  expect_true("sequence_intervals" %in% names(md))
  nodes <- ts_nodes(ts)
  expect_equal(nrow(nodes), ts_info(ts)$num_nodes)
  sample_md <- nodes$metadata[[which(nodes$flags == 1)[1]]]
  expect_setequal(names(sample_md), c("source", "sample_id", "ploidy_index"))
  G <- ts_genotype_matrix(ts)
  expect_equal(nrow(G), ts_info(ts)$num_sites)
})

test_that("ts_from_tables builds valid tree sequences", {
  ts <- ts_from_tables(
    nodes = data.frame(flags = c(0, 1, 1), time = c(1, 0, 0)),
    edges = data.frame(left = 0, right = 10, parent = 0, child = 1:2),
    sites = data.frame(position = 5, ancestral_state = "A"),
    mutations = data.frame(site = 0, node = 1, derived_state = "T"),
    sequence_length = 10
  )
  G <- ts_genotype_matrix(ts)
  expect_equal(as.vector(G), c(1L, 0L))
  expect_equal(attr(G, "alleles")[[1]], c("A", "T"))
  s <- ts_simplify(ts, samples = 1L)
  expect_equal(ts_info(s)$num_nodes, 1)
  expect_error(ts_from_tables(
    nodes = data.frame(flags = c(0, 1), time = c(0, 1)),
    edges = data.frame(left = 0, right = 10, parent = 0, child = 1),
    sequence_length = 10
  ), "tskit error")
})
