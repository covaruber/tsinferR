# Run tsinferR on the dataset produced by tools/compare_large.py and compare
# with the Python results. Usage: Rscript tools/compare_large.R OUTDIR
library(tsinferR)
args <- commandArgs(trailingOnly = TRUE)
out <- normalizePath(args[1])
threads <- 4L
anc <- file.path(out, "r_anc.vcz")
unlink(anc, recursive = TRUE)
samples <- file.path(out, "samples.vcz")
cfg <- tsinfer_config(
  sources = list(tsinfer_source(samples, "samples")),
  ancestral_state = ancestral_state(samples, is_reference = TRUE),
  ancestors = list(ancestors_config("ancestors", anc, "samples",
    samples_chunk_size = 50, variants_chunk_size = 1000
  )),
  match = match_config(
    sources = list(
      ancestors = match_source_config(node_flags = 0, create_individuals = FALSE),
      samples = match_source_config()
    ),
    output = "x.trees"
  ),
  post_process = post_process_config()
)
t <- list()
t$infer_ancestors <- system.time(infer_ancestors(cfg$sources$samples, cfg$ancestors[[1]],
  cfg$ancestral_state, num_threads = threads))[["elapsed"]]
for (name in c("call_genotype", "sample_time", "sample_start_position", "sample_end_position",
               "sample_focal_positions", "variant_position", "variant_allele")) {
  same <- identical(zarr_read(anc, name), zarr_read(file.path(out, "python_anc.vcz"), name))
  cat(sprintf("ancestors %-24s identical: %s\n", name, same))
}
t$match <- system.time(raw <- match_haplotypes(cfg, num_threads = threads))[["elapsed"]]
t$post_process <- system.time(final <- post_process(raw, cfg))[["elapsed"]]
py <- jsonlite::fromJSON(file.path(out, "python_timings.json"))
cat("raw identical to Python:  ", ts_equals(raw, ts_load(file.path(out, "python_raw.trees"))), "\n")
cat("final identical to Python:", ts_equals(final, ts_load(file.path(out, "python_final.trees"))), "\n")
print(rbind(tsinferR = unlist(t), python = unlist(py)[names(t)]))
print(ts_info(final))
