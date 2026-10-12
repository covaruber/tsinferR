test_that("variant filters match vcztools", {
  cases <- read_fixture_json("filters.json")
  path <- fixture_path("bio2zarr.vcz")
  for (cs in cases) {
    args <- cs$args
    src <- do.call(tsinfer_source, c(list(path = path, name = "s"), args))
    info <- tsinferR:::open_source(src)
    got <- tsinferR:::filtered_positions(info)
    expect_equal(got, as.numeric(unlist(cs$positions)), info = jsonlite::toJSON(args, auto_unbox = TRUE))
  }
})

test_that("unsupported filter features give clear errors", {
  path <- fixture_path("bio2zarr.vcz")
  check <- function(expr) {
    src <- tsinfer_source(path, "s", include = expr)
    tsinferR:::filtered_positions(tsinferR:::open_source(src))
  }
  expect_error(check("GT = 1"), "genotype values")
  expect_error(check("strlen(REF) > 1"), "functions")
  expect_error(check("NOPE > 1"), "not defined")
  expect_error(check("POS > "), "Unexpected end")
  expect_error(check("(POS > 1"), "Missing")
  expect_error(check('TYPE="indel"'), "TYPE values")
  expect_error(check('ID ~ "rs"'), "regular expressions")
  expect_error(tsinfer_source(path, "s", include = "POS>1", exclude = "POS<1"), "cannot both")
})

test_that("sample selection follows bcftools conventions", {
  ids <- c("a", "b", "c", "d")
  expect_equal(tsinferR:::parse_samples_spec(NULL, ids), 0:3)
  expect_equal(tsinferR:::parse_samples_spec("c,a", ids), c(2L, 0L))
  expect_equal(tsinferR:::parse_samples_spec("^b,d", ids), c(0L, 2L))
  expect_warning(res <- tsinferR:::parse_samples_spec("a,zz", ids), "not in header")
  expect_equal(res, 0L)
})

test_that("source filters and sample subsets flow through the pipeline", {
  dir <- tempfile("filt-")
  dir.create(dir)
  store <- fixture_path("case1.vcz")
  pos <- zarr_read(store, "variant_position")
  cut <- pos[floor(length(pos) / 2)]
  cfg <- tsinfer_config(
    sources = list(tsinfer_source(store, "test",
      include = sprintf("POS <= %d", cut),
      samples = "^s0,s1"
    )),
    ancestral_state = ancestral_state(store, "variant_ancestral_allele"),
    ancestors = list(ancestors_config("ancestors", file.path(dir, "anc.vcz"), "test")),
    match = match_config(sources = list(
      ancestors = match_source_config(node_flags = 0, create_individuals = FALSE),
      test = match_source_config()
    )),
    post_process = post_process_config()
  )
  ts <- run_tsinfer(cfg)
  anc_pos <- zarr_read(file.path(dir, "anc.vcz"), "variant_position")
  expect_true(all(anc_pos <= cut))
  nodes <- ts_nodes(ts)
  sample_ids <- vapply(nodes$metadata[nodes$flags == 1], function(m) m$sample_id, "")
  expect_setequal(sample_ids, paste0("s", 2:15))
  expect_equal(ts_info(ts)$num_individuals, 14)
})

test_that("TOML configurations are parsed and validated", {
  toml <- '
[[source]]
name = "cohort"
path = "cohort.vcz"
include = "QUAL > 30"
sample_time = 0.5

[[source]]
name = "old"
path = "old.vcz"
sample_time = {path = "times.vcz", field = "t"}

[ancestral_state]
path = "anc.vcz"
is_reference = true

[[ancestors]]
name = "anc"
path = "ancestors.vcz"
sources = ["cohort", "old"]
max_gap_length = 1000

[match]
output = "out.trees"
path_compression = true
workdir = "work"
keep_intermediates = true
[match.sources.anc]
node_flags = 0
create_individuals = false
[match.sources.cohort]
recombination = 0.001
mismatch = 1e-10
[match.sources.old]

[individual_metadata]
fields = {pop = "sample_pop"}
population = "sample_pop"

[post_process]
erase_flanks = false

[augment_sites]
sources = ["cohort"]
'
  cfg <- read_config(text = toml)
  expect_s3_class(cfg, "tsinfer_config")
  expect_equal(names(cfg$sources), c("cohort", "old", "anc"))
  expect_equal(cfg$sources$cohort$sample_time, 0.5)
  expect_mapequal(cfg$sources$old$sample_time, list(path = "times.vcz", field = "t"))
  expect_equal(cfg$sources$anc$sample_time, "sample_time")
  expect_equal(cfg$ancestors[[1]]$max_gap_length, 1000L)
  expect_equal(names(cfg$match$sources), c("anc", "cohort", "old"))
  expect_equal(cfg$match$sources$cohort$recombination, 0.001)
  expect_true(cfg$match$path_compression)
  expect_false(cfg$post_process$erase_flanks)
  expect_true(cfg$post_process$split_ultimate)
  expect_equal(cfg$individual_metadata$fields, list(pop = "sample_pop"))
  expect_equal(cfg$augment_sites$sources, "cohort")
  expect_true(isTRUE(cfg$ancestral_state$is_reference))
  errors <- validate_config(cfg)
  expect_true(any(grepl("cohort.vcz", errors)))
  expect_true(any(grepl("times.vcz", errors)))
  expect_output(print(cfg), "\\[match\\]")

  bad <- function(text) read_config(text = text)
  ordered <- read_config(text = paste0(
    '[[source]]\nname="zeta"\npath="z"\n[[ancestors]]\nname="b.anc"\npath="a"\n',
    'sources=["zeta"]\n[ancestral_state]\npath="x"\nfield="f"\n[match]\noutput="o"\n',
    '[match.sources."b.anc"]\nnode_flags=0\n[match.sources.zeta]\n'
  ))
  expect_equal(names(ordered$match$sources), c("b.anc", "zeta"))
  ordered <- read_config(text = paste0(
    '[[source]]\nname="zeta"\npath="z"\n[[ancestors]]\nname="anc"\npath="a"\n',
    'sources=["zeta"]\n[ancestral_state]\npath="x"\nfield="f"\n[match]\noutput="o"\n',
    'sources = {zeta = {}, anc = {node_flags = 0}}\n'
  ))
  expect_equal(names(ordered$match$sources), c("zeta", "anc"))
  expect_error(bad('[foo]\nx=1\n'), "Unrecognised key")
  expect_error(bad('[[source]]\nname="a"\npath="x"\nbad=1\n'), "Unrecognised key")
  expect_error(bad('[[source]]\npath="x"\n'), "name")
  expect_error(
    bad('[ancestral_state]\npath="x"\nfield="f"\n[match]\noutput="o"\n[match.sources.a]\n'),
    "\\[\\[ancestors\\]\\]"
  )
  expect_error(
    bad('[[ancestors]]\npath="a"\nsources=["s"]\n[ancestral_state]\npath="x"\nfield="f"\n[match]\noutput="o"\n[match.sources.s]\n'),
    "must appear in"
  )
  expect_error(ancestral_state("x"), "requires either")
  expect_error(ancestral_state("x", "f", TRUE), "must not be set")
  expect_error(tsinfer_source("x", ""), "non-empty name")
})

test_that("canonical JSON matches Python's encoding", {
  cj <- tsinferR:::canonical_json
  expect_equal(cj(list(b = 1L, a = "x\u00e9\"")), '{"a":"x\\u00e9\\"","b":1}')
  expect_equal(cj(list(x = 1)), '{"x":1.0}')
  expect_equal(cj(list(x = 1e-20, y = 0.1, z = 1e16, w = 123456.5)),
    '{"w":123456.5,"x":1e-20,"y":0.1,"z":1e+16}')
  expect_equal(cj(list(a = TRUE, b = NULL, c = list(1L, 2L))), '{"a":true,"b":null,"c":[1,2]}')
  expect_equal(cj(list(s = "\U0001F600\n")), '{"s":"\\ud83d\\ude00\\n"}')
  expect_equal(tsinferR:::.json_double(c(1.5e-5, -2, NaN, Inf)), c("1.5e-05", "-2.0", "NaN", "Infinity"))
})

test_that("grouping matches tsinfer's linesweep", {
  g <- tsinferR:::.group_by_linesweep(c(0, 0, 5, 0), c(10, 4, 10, 10), c(3, 2, 2, 1))
  expect_equal(g, c(0L, 1L, 1L, 2L))
  expect_equal(tsinferR:::.group_by_linesweep(numeric(0), numeric(0), numeric(0)), integer(0))
  # Same-time overlapping haplotypes share a group.
  g <- tsinferR:::.group_by_linesweep(c(0, 5, 20), c(10, 15, 30), c(1, 1, 1))
  expect_equal(g, c(0L, 0L, 0L))
})

test_that("edge cases are reported clearly", {
  dir <- tempfile("edge-")
  dir.create(dir)
  store <- fixture_path("case1.vcz")
  cfg <- tsinfer_config(
    sources = list(tsinfer_source(store, "test", include = "POS < 0")),
    ancestral_state = ancestral_state(store, "variant_ancestral_allele"),
    ancestors = list(ancestors_config("ancestors", file.path(dir, "anc.vcz"), "test")),
    match = match_config(sources = list(
      ancestors = match_source_config(node_flags = 0, create_individuals = FALSE),
      test = match_source_config()
    ))
  )
  infer_ancestors(cfg$sources$test, cfg$ancestors[[1]], cfg$ancestral_state)
  anc <- zarr_open(file.path(dir, "anc.vcz"))
  expect_equal(zarr_info(anc, "call_genotype")$shape, c(0, 0, 1))
  expect_error(match_haplotypes(cfg), "no inference sites")
  expect_error(
    infer_ancestors(cfg$sources$test, cfg$ancestors[[1]], cfg$ancestral_state),
    "already exists"
  )
  cfg$ancestors <- list()
  expect_error(match_haplotypes(cfg), "not supported yet")
})
