fixture_env <- new.env()

# Zarr stores contain dot-files, so they are shipped as a tarball and
# extracted once per test session.
fixture_path <- function(...) {
  name <- file.path(...)
  if (grepl("\\.vcz$", name)) {
    if (is.null(fixture_env$dir)) {
      dir <- tempfile("tsinferR-fixtures-")
      dir.create(dir)
      utils::untar(testthat::test_path("fixtures", "stores.tar.gz"), exdir = dir)
      fixture_env$dir <- dir
    }
    return(file.path(fixture_env$dir, name))
  }
  testthat::test_path("fixtures", name)
}

read_fixture_json <- function(name) {
  jsonlite::fromJSON(fixture_path(name), simplifyVector = FALSE)
}

# Compare two tree sequences table-by-table, ignoring provenance.
expect_ts_equal <- function(actual, expected, metadata = TRUE) {
  testthat::expect_equal(ts_info(actual)[-8], ts_info(expected)[-8])
  testthat::expect_equal(ts_nodes(actual, FALSE)[c("flags", "time", "population", "individual")],
    ts_nodes(expected, FALSE)[c("flags", "time", "population", "individual")])
  testthat::expect_equal(ts_edges(actual), ts_edges(expected))
  testthat::expect_equal(ts_sites(actual, FALSE), ts_sites(expected, FALSE))
  testthat::expect_equal(ts_mutations(actual, FALSE)[c("site", "node", "parent", "derived_state")],
    ts_mutations(expected, FALSE)[c("site", "node", "parent", "derived_state")])
  if (metadata) {
    testthat::expect_equal(ts_nodes(actual, FALSE)$metadata, ts_nodes(expected, FALSE)$metadata)
    testthat::expect_equal(ts_individuals(actual, FALSE), ts_individuals(expected, FALSE))
    testthat::expect_equal(ts_populations(actual, FALSE), ts_populations(expected, FALSE))
    testthat::expect_equal(ts_metadata(actual), ts_metadata(expected))
  }
  testthat::expect_true(ts_equals(actual, expected, ignore_provenance = TRUE))
}

copy_fixture_dir <- function(name) {
  dest <- file.path(tempfile("tsinferR-"), name)
  dir.create(dirname(dest), recursive = TRUE)
  file.copy(fixture_path(name), dirname(dest), recursive = TRUE)
  dest
}
