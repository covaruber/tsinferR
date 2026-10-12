# tsinferR 0.1.0

* First release. A native R implementation of the tsinfer pipeline that does
  not depend on Python:
  * `infer_ancestors()`, `match_haplotypes()`, `post_process()`,
    `augment_sites()` and `run_tsinfer()`;
  * TOML configuration files;
  * a Zarr v2/v3 engine;
  * VCF to VCZ conversion.
* Results are identical to tsinfer (Python) for the same input.
* Fixes for problems found in tsinfer:
  * the matcher position buffer was allocated with the wrong element size;
  * matcher indexes leaked memory after every match group;
  * `find_path()` with an empty site range read out of bounds;
  * calling `ancestor_descriptors()` twice returned no ancestors;
  * path compression was not idempotent;
  * sample start and end positions were misaligned when a sample subset was
    used;
  * individuals were lost when a match resumed from a workdir without a
    checkpoint;
  * ancestor-builder arena allocations could misalign 32-bit focal-site IDs;
  * zero-length tskit memory comparisons could pass null pointers to libc;
  * an error message was missing a space.
