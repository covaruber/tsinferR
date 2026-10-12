"""
Larger-scale parity/performance benchmark between tsinfer (Python) and
tsinferR. Development tool only.

    cd ../tsinfer
    uv run python ../tsinferR/tools/compare_large.py OUTDIR [num_individuals] [seq_len]

Writes OUTDIR/samples.vcz (via bio2zarr from a simulated VCF), runs the Python
pipeline and stores OUTDIR/python_raw.trees, python_final.trees and timings.
Then run tools/compare_large.R OUTDIR to run tsinferR and compare.
"""

import json
import pathlib
import shutil
import sys
import time

import msprime
import numpy as np
from bio2zarr import vcf as b2z

import tsinfer
from tsinfer import config, pipeline

out = pathlib.Path(sys.argv[1]).resolve()
n = int(sys.argv[2]) if len(sys.argv) > 2 else 100
L = int(float(sys.argv[3])) if len(sys.argv) > 3 else 1_000_000
threads = 4
out.mkdir(parents=True, exist_ok=True)

ts = msprime.sim_ancestry(n, ploidy=2, sequence_length=L, recombination_rate=1e-8,
                          population_size=10_000, random_seed=1234)
ts = msprime.sim_mutations(ts, rate=1e-8, random_seed=1234)
vcf = out / "samples.vcf"
with open(vcf, "w") as f:
    ts.write_vcf(f, contig_id="1")
shutil.rmtree(out / "samples.vcz", ignore_errors=True)
b2z.convert([str(vcf)], str(out / "samples.vcz"), worker_processes=0,
            variants_chunk_size=500, samples_chunk_size=40)

anc_path = out / "python_anc.vcz"
shutil.rmtree(anc_path, ignore_errors=True)
cfg = config.Config(
    sources={
        "samples": config.Source(path=str(out / "samples.vcz"), name="samples"),
        "ancestors": config.Source(path=str(anc_path), name="ancestors",
                                   sample_time="sample_time"),
    },
    ancestors=[config.AncestorsConfig(name="ancestors", path=str(anc_path),
                                      sources=["samples"], samples_chunk_size=50,
                                      variants_chunk_size=1000)],
    match=config.MatchConfig(
        sources={"ancestors": config.MatchSourceConfig(node_flags=0,
                                                       create_individuals=False),
                 "samples": config.MatchSourceConfig()},
        output="x.trees"),
    ancestral_state=config.AncestralState(path=str(out / "samples.vcz"),
                                          is_reference=True),
    post_process=config.PostProcessConfig(),
)
timings = {}
t0 = time.time()
tsinfer.infer_ancestors([cfg.sources["samples"]], cfg.ancestors[0], cfg.ancestral_state,
                        num_threads=threads)
timings["infer_ancestors"] = time.time() - t0
t0 = time.time()
raw = pipeline.match(cfg, num_threads=threads)
timings["match"] = time.time() - t0
t0 = time.time()
final = pipeline.post_process(raw, cfg)
timings["post_process"] = time.time() - t0
raw.dump(out / "python_raw.trees")
final.dump(out / "python_final.trees")
(out / "python_timings.json").write_text(json.dumps(timings))
print("sites", ts.num_sites, "samples", ts.num_samples, timings)
print("raw", raw.num_nodes, raw.num_edges, "final", final.num_nodes, final.num_edges)
