test_that("AncestorBuilder matches the Python reference", {
  fx <- read_fixture_json("ancestor_builder.json")
  H <- do.call(rbind, lapply(fx$genotypes, unlist))
  times <- unlist(fx$times)
  for (enc in 0:1) {
    ab <- ancestor_builder(ncol(H), nrow(H), genotype_encoding = enc)
    for (j in seq_len(nrow(H))) {
      ancestor_builder_add_site(ab, times[j], H[j, ])
    }
    d <- ancestor_builder_descriptors(ab)
    expected <- fx[[paste0("encoding_", enc)]]
    expect_equal(length(d$time), length(expected))
    for (k in seq_along(expected)) {
      e <- expected[[k]]
      expect_equal(d$time[k], e$time)
      expect_equal(d$focal_sites[[k]], as.integer(unlist(e$focal)))
      a <- ancestor_builder_make_ancestor(ab, d$focal_sites[[k]])
      expect_equal(a$start, e$start)
      expect_equal(a$end, e$end)
      expect_equal(a$haplotype, as.integer(unlist(e$haplotype)))
    }
    # Descriptors are stable when requested twice.
    expect_equal(ancestor_builder_descriptors(ab), d)
  }
})

test_that("AncestorBuilder validates input", {
  expect_error(ancestor_builder(1, 10), "at least 2 samples")
  ab <- ancestor_builder(4, 1)
  expect_error(ancestor_builder_add_site(ab, 0.5, c(1, 0)), "wrong size")
  expect_error(ancestor_builder_add_site(ab, 0.5, c(2, 0, 0, 0)), "genotypes must be")
  ancestor_builder_add_site(ab, 0.5, c(1, 1, 0, 0))
  expect_error(ancestor_builder_add_site(ab, 0.5, c(1, 1, 0, 0)), "more sites")
  expect_error(ancestor_builder_make_ancestor(ab, integer(0)), "num_focal_sites")
  expect_error(ancestor_builder_make_ancestor(ab, 5L), "out of bounds")
  ab1 <- ancestor_builder(4, 2, genotype_encoding = 1)
  expect_error(ancestor_builder_add_site(ab1, 0.5, c(1, -1, 0, 0)), "binary")
})

test_that("AncestorMatcher matches the Python reference", {
  ts <- ts_load(fixture_path("case1_raw.trees"))
  mi <- matcher_indexes(ts)
  ns <- ts_info(ts)$num_sites
  cases <- read_fixture_json("ancestor_matcher.json")
  for (cs in cases) {
    am <- ancestor_matcher(mi, rep(cs$rho, ns), rep(cs$mu, ns))
    res <- find_path(am, unlist(cs$haplotype), cs$start, cs$end)
    expect_equal(res$left, as.integer(unlist(cs$left)))
    expect_equal(res$right, as.integer(unlist(cs$right)))
    expect_equal(res$parent, as.integer(unlist(cs$parent)))
    expect_equal(res$match[(cs$start + 1):cs$end], as.integer(unlist(cs$match)))
  }
})

test_that("AncestorMatcher validates input", {
  ts <- ts_load(fixture_path("case1_raw.trees"))
  mi <- matcher_indexes(ts)
  ns <- ts_info(ts)$num_sites
  expect_error(ancestor_matcher(mi, rep(0.1, ns - 1), rep(0.1, ns)), "recombination")
  am <- ancestor_matcher(mi, rep(1e-2, ns), rep(1e-20, ns))
  h <- integer(ns)
  expect_error(find_path(am, h[-1], 0, ns), "Incorrect size")
  expect_error(find_path(am, h, 0, 0), "Bad match range")
  expect_error(find_path(am, h, 0, ns + 1), "Bad match range")
  expect_error(find_path(am, h + 5L, 0, ns), "bad allele")
  expect_error(matcher_indexes(ts, num_alleles = rep(1, ns)), "between 2 and 127")
})
