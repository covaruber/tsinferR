case_config <- function(case, workdir, sample_store = NULL, anc_state_store = NULL,
                        path_compression = FALSE) {
  sample_store <- sample_store %||% fixture_path(paste0(case, ".vcz"))
  anc_state_store <- anc_state_store %||% sample_store
  tsinfer_config(
    sources = list(tsinfer_source(sample_store, "test")),
    ancestral_state = ancestral_state(anc_state_store, "variant_ancestral_allele"),
    ancestors = list(ancestors_config("ancestors", file.path(workdir, paste0(case, "_anc.vcz")),
      "test",
      samples_chunk_size = 7, variants_chunk_size = 13
    )),
    match = match_config(
      sources = list(
        ancestors = match_source_config(node_flags = 0, create_individuals = FALSE),
        test = match_source_config()
      ),
      output = file.path(workdir, "out.trees"), path_compression = path_compression
    ),
    post_process = post_process_config()
  )
}

`%||%` <- function(a, b) if (is.null(a)) b else a

expect_same_ancestors <- function(actual_path, expected_path) {
  a <- zarr_open(actual_path)
  e <- zarr_open(expected_path)
  for (name in c(
    "call_genotype", "variant_position", "variant_allele", "sample_time",
    "sample_start_position", "sample_end_position", "sample_focal_positions",
    "sequence_intervals", "sample_id", "contig_id", "contig_length"
  )) {
    expect_equal(zarr_read(a, name), zarr_read(e, name), info = name)
  }
}

test_that("case 1 (haploid) reproduces the Python pipeline exactly", {
  dir <- tempfile("case1-")
  dir.create(dir)
  cfg <- case_config("case1", dir)
  infer_ancestors(cfg$sources$test, cfg$ancestors[[1]], cfg$ancestral_state)
  expect_same_ancestors(cfg$ancestors[[1]]$path, fixture_path("case1_anc.vcz"))
  raw <- match_haplotypes(cfg)
  expect_ts_equal(raw, ts_load(fixture_path("case1_raw.trees")))
  final <- post_process(raw, cfg)
  expect_ts_equal(final, ts_load(fixture_path("case1_final.trees")))
})

test_that("case 2 (diploid, missing data, populations, threads) matches Python", {
  dir <- tempfile("case2-")
  dir.create(dir)
  cfg <- case_config("case2", dir)
  cfg$individual_metadata <- individual_metadata_config(
    fields = list(population_name = "sample_population", weight = "sample_weight"),
    population = "sample_population"
  )
  infer_ancestors(cfg$sources$test, cfg$ancestors[[1]], cfg$ancestral_state, num_threads = 3)
  expect_same_ancestors(cfg$ancestors[[1]]$path, fixture_path("case2_anc.vcz"))
  raw <- match_haplotypes(cfg, num_threads = 3)
  expect_ts_equal(raw, ts_load(fixture_path("case2_raw.trees")))
  final <- post_process(raw, cfg)
  expect_ts_equal(final, ts_load(fixture_path("case2_final.trees")))
  pops <- ts_populations(final)
  expect_equal(vapply(pops$metadata, function(m) m$name, ""), c("popA", "popB"))
})

test_that("case 3 (subset sites, path compression, augment_sites) matches Python", {
  dir <- tempfile("case3-")
  dir.create(dir)
  cfg <- case_config("case3", dir,
    sample_store = fixture_path("case3_subset.vcz"),
    anc_state_store = fixture_path("case3_full.vcz"), path_compression = TRUE
  )
  cfg <- tsinfer_config(
    sources = c(unname(cfg$sources), list(tsinfer_source(fixture_path("case3_full.vcz"), "full"))),
    ancestral_state = cfg$ancestral_state, ancestors = cfg$ancestors,
    match = cfg$match, post_process = cfg$post_process,
    augment_sites = augment_sites_config("full")
  )
  ts <- run_tsinfer(cfg, num_threads = 2)
  expect_ts_equal(ts, ts_load(fixture_path("case3_augmented.trees")))
  raw <- match_haplotypes(cfg)
  expect_ts_equal(raw, ts_load(fixture_path("case3_raw.trees")))
  expect_gt(sum(bitwAnd(as.integer(ts_nodes(raw)$flags), 65536L) != 0), 0)
})

test_that("checkpointing in a workdir resumes and gives the same result", {
  dir <- tempfile("wd-")
  dir.create(dir)
  cfg <- case_config("case1", dir)
  cfg$match$workdir <- file.path(dir, "work")
  infer_ancestors(cfg$sources$test, cfg$ancestors[[1]], cfg$ancestral_state)
  expect_warning(partial <- match_haplotypes(cfg, group_stop = 3), NA)
  expect_true(file.exists(file.path(cfg$match$workdir, "match-jobs.json")))
  expect_true(file.exists(file.path(cfg$match$workdir, "group_2.trees")))
  full <- match_haplotypes(cfg)
  expect_ts_equal(full, ts_load(fixture_path("case1_raw.trees")))
  jobs <- jsonlite::fromJSON(file.path(cfg$match$workdir, "match-jobs.json"))
  expect_setequal(names(jobs), c(
    "haplotype_index", "source", "sample_id", "ploidy_index", "time",
    "start_position", "end_position", "group", "node_flags", "individual_id",
    "population_id", "sample_chunk"
  ))
})

test_that("match_file records copying paths", {
  dir <- tempfile("mf-")
  dir.create(dir)
  cfg <- case_config("case1", dir)
  infer_ancestors(cfg$sources$test, cfg$ancestors[[1]], cfg$ancestral_state)
  mf <- file.path(dir, "match.jsonl")
  match_haplotypes(cfg, match_file = mf)
  lines <- readLines(mf)
  info <- zarr_info(cfg$ancestors[[1]]$path, "call_genotype")
  expect_equal(length(lines), info$shape[2] + 16)
  rec <- jsonlite::fromJSON(lines[length(lines)])
  expect_true(all(c("path", "mutations", "sample_id") %in% names(rec)))
})

test_that("a TOML configuration runs end to end", {
  dir <- tempfile("toml-")
  dir.create(dir)
  sample <- fixture_path("case1.vcz")
  toml <- sprintf('
[[source]]
name = "test"
path = "%s"

[ancestral_state]
path = "%s"
field = "variant_ancestral_allele"

[[ancestors]]
name = "ancestors"
path = "%s"
sources = ["test"]
samples_chunk_size = 7
variants_chunk_size = 13

[match]
output = "%s"
[match.sources.ancestors]
node_flags = 0
create_individuals = false
[match.sources.test]

[post_process]
', sample, sample, file.path(dir, "anc.vcz"), file.path(dir, "out.trees"))
  path <- file.path(dir, "config.toml")
  writeLines(toml, path)
  ts <- run_tsinfer(path, output = TRUE)
  expect_true(file.exists(file.path(dir, "out.trees")))
  expect_ts_equal(ts, ts_load(fixture_path("case1_final.trees")))
})
