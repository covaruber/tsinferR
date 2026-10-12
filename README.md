# tsinferR

`tsinferR` is a native R implementation of
[tsinfer](https://tskit.dev/tsinfer), which infers tree sequences (ancestral
recombination graphs) from genetic variation data. The `tsinfer` and `tskit`
C libraries are bundled with the package, and the Python orchestration code
is reimplemented in R and C++. No Python installation is required.

`tsinferR` is an R mirror of the canonical `tsinfer` method, not a new
inference method or a replacement for the original research software. If you
use this R mirror, please cite the canonical tsinfer paper by Kelleher et al.
(2019); `citation("tsinferR")` returns the recommended citation.

## Features

* The full `tsinfer` pipeline:
  * ancestor generation with `infer_ancestors()`;
  * Li & Stephens haplotype matching with `match_haplotypes()`, multi-threaded
    in C++;
  * post-processing with `post_process()`;
  * parsimony-based site placement with `augment_sites()`;
  * `run_tsinfer()` to run all of these in sequence.
* Output identical to the Python `tsinfer` package for the same input and
  configuration. The test suite checks this against reference outputs.
* TOML configuration files compatible with the `tsinfer` command-line tool,
  read with `read_config()`.
* A self-contained Zarr v2/v3 engine written in C/C++, used to read and write
  VCF Zarr (VCZ) stores: `zarr_open()`, `zarr_read()` and
  `zarr_write_array()`. It bundles c-blosc, Zstandard and LZ4, and links zlib.
* Fast conversion of VCF files (plain or gzip/bgzip-compressed) to VCZ with
  `vcf_to_vcz()`. Its output matches `bio2zarr`.
* bcftools-style variant filters (`include`, `exclude`, `regions`,
  `targets`) and sample selection, as supported by `vcztools`.
* Tree sequences are saved as standard `.trees` files (`ts_dump()`) and can
  be inspected as data frames (`ts_nodes()`, `ts_edges()`, ...).

## Installation

```r
# install.packages("remotes")
remotes::install_local("path/to/tsinferR")
```

You need a C/C++17 toolchain and zlib.

## Example

```r
library(tsinferR)
vcf <- system.file("extdata", "example.vcf.gz", package = "tsinferR")
samples <- vcf_to_vcz(vcf, "samples.vcz")

cfg <- tsinfer_config(
  sources = list(tsinfer_source(samples, "samples")),
  ancestral_state = ancestral_state(samples, is_reference = TRUE),
  ancestors = list(ancestors_config("ancestors", "ancestors.vcz", "samples")),
  match = match_config(sources = list(
    ancestors = match_source_config(node_flags = 0, create_individuals = FALSE),
    samples = match_source_config()
  )),
  post_process = post_process_config()
)
ts <- run_tsinfer(cfg, num_threads = 4)
ts_dump(ts, "inferred.trees")
```

## Python to R mapping

| Python (`tsinfer`)             | R (`tsinferR`)                   |
|--------------------------------|----------------------------------|
| `tsinfer.infer_ancestors()`    | `infer_ancestors()`              |
| `tsinfer.match()`              | `match_haplotypes()`             |
| `tsinfer.post_process()`       | `post_process()`                 |
| `tsinfer.augment_sites()`      | `augment_sites()`                |
| `tsinfer.run()`                | `run_tsinfer()`                  |
| `tsinfer.Config.from_toml()`   | `read_config()`                  |
| `tsinfer.Source`, `AncestralState`, `AncestorsConfig`, `MatchConfig`, `MatchSourceConfig`, `PostProcessConfig`, `AugmentSitesConfig`, `IndividualMetadataConfig`, `Config` | `tsinfer_source()`, `ancestral_state()`, `ancestors_config()`, `match_config()`, `match_source_config()`, `post_process_config()`, `augment_sites_config()`, `individual_metadata_config()`, `tsinfer_config()` |
| `_tsinfer.AncestorBuilder`     | `ancestor_builder()` and related |
| `_tsinfer.AncestorMatcher`     | `ancestor_matcher()`, `find_path()` |

The R function is called `match_haplotypes()` rather than `match()` because
`match()` would mask `base::match()`.

## R ecosystem interoperability

Tree sequences use the standard `.trees` format. The optional `tskitR`
package can load tsinferR output directly, or `as_tskit()` can create a
`tskitR` `TreeSequence` object:

```r
library(tsinferR)
tree_sequence <- ts_load("inferred.trees")
tskit_sequence <- as_tskit(tree_sequence)
tskit_sequence$num_trees
```

With `tsdateR` installed, `date_tsinfer()` dates the sequence and returns the
dated result as a `tsinfer_ts`:

```r
dated <- date_tsinfer(tree_sequence, mutation_rate = 1e-8)
ts_nodes(dated)$time
```

Simulation is outside the scope of tsinferR; supply VCF or prebuilt VCZ data
from an R simulation workflow or other source.

## Differences from Python tsinfer

* `compress_paths()` is idempotent. In `tsinfer` it adds a chain of unary
  path-compression nodes at every match group.
* Matcher indexes are freed after each group. In `tsinfer` they leak.
* Ancestor start and end positions are taken from the selected store
  columns when a sample subset is used.
* The haplotype labels and genotype columns of a multi-source view use the
  same sample order.
* When a `match()` run is resumed from a workdir before its first
  checkpoint, the individuals are rebuilt.
* `match.reference_ts` is rejected with a clear error. `tsinfer` ignores it
  and then crashes.
* The C library has extra input validation, including empty match ranges and
  unindexed tables.
* Only local file-system stores are supported. Remote stores (`s3://`,
  `https://`) are not.
