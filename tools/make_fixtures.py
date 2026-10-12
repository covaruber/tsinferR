"""
Generate reference fixtures for tsinferR from the Python tsinfer implementation.

Development tool only (not part of the R package build). Run from the tsinfer
checkout so that the Python reference implementation is importable:

    cd ../tsinfer
    uv run python ../tsinferR/tools/make_fixtures.py ../tsinferR/tests/testthat/fixtures
"""

import json
import pathlib
import shutil
import sys
import tarfile

import msprime
import numpy as np
import tskit
import zarr
from zarr.core.dtype.npy.string import VariableLengthUTF8

import _tsinfer
import tsinfer
from tsinfer import arg_ops, config, pipeline

OUT = pathlib.Path(sys.argv[1]).resolve()
OUT.mkdir(parents=True, exist_ok=True)
# Zarr stores are shipped packed (they contain dot-files, which CRAN flags).
if (OUT / "stores.tar.gz").exists():
    with tarfile.open(OUT / "stores.tar.gz") as tar:
        tar.extractall(OUT, filter="data")

_VLEN = VariableLengthUTF8()

# tsinferR fixes a bug in tsinfer's compress_paths: it is not idempotent and
# adds a chain of unary path-compression nodes on every match group. The
# reference is generated with the same fix: groups whose whole path is to a
# single path-compression ancestor are already compressed and are skipped.
_original_compress_paths = arg_ops.compress_paths


def _fixed_compress_paths(ts):
    import collections

    tables = ts.dump_tables()
    node_time = tables.nodes.time
    flags = tables.nodes.flags
    edges_by_child = collections.defaultdict(list)
    for edge in tables.edges:
        edges_by_child[edge.child].append((edge.left, edge.right, edge.parent))
    path_map = collections.defaultdict(list)
    for child, edges in edges_by_child.items():
        path_map[tuple(sorted(edges))].append(child)
    compressed_children = set()
    new_edges = []
    for path, children in path_map.items():
        if len(children) < 2:
            continue
        parents = {p for _, _, p in path}
        if len(parents) == 1 and arg_ops.is_pc_ancestor(flags[next(iter(parents))]):
            continue
        min_parent_time = min(node_time[p] for _, _, p in path)
        pc_time = min_parent_time - arg_ops.PC_ANCESTOR_INCREMENT
        if pc_time <= max(node_time[c] for c in children):
            continue
        pc_node = tables.nodes.add_row(flags=arg_ops.NODE_IS_PC_ANCESTOR, time=pc_time)
        for left, right, parent in path:
            new_edges.append((left, right, parent, pc_node))
        for child in children:
            for left, right, _ in path:
                new_edges.append((left, right, pc_node, child))
        compressed_children.update(children)
    tables.edges.clear()
    for edge in ts.edges():
        if edge.child not in compressed_children:
            tables.edges.add_row(left=edge.left, right=edge.right, parent=edge.parent,
                                 child=edge.child)
    for left, right, parent, child in new_edges:
        tables.edges.add_row(left=left, right=right, parent=parent, child=child)
    tables.sort()
    return tables.tree_sequence()


arg_ops.compress_paths = _fixed_compress_paths


def write_sample_vcz(path, genotypes, positions, alleles, ancestral, seq_len,
                     sample_ids, extra=None, chunks=None):
    shutil.rmtree(path, ignore_errors=True)
    root = zarr.open_group(str(path), mode="w", zarr_format=2)

    def arr(name, data, dims, chunk=None):
        kw = {} if chunk is None else {"chunks": chunk}
        a = root.create_array(name, data=np.asarray(data), **kw)
        a.attrs["_ARRAY_DIMENSIONS"] = dims

    def sarr(name, data, dims, chunk=None):
        data = np.asarray(data)
        kw = {} if chunk is None else {"chunks": chunk}
        a = root.create_array(name, shape=data.shape, dtype=_VLEN, **kw)
        a[:] = data
        a.attrs["_ARRAY_DIMENSIONS"] = dims

    nsites, nsamp, _ = genotypes.shape
    gchunks = None if chunks is None else (chunks[0], chunks[1], genotypes.shape[2])
    arr("call_genotype", genotypes.astype(np.int8), ["variants", "samples", "ploidy"], gchunks)
    vch = None if chunks is None else (chunks[0],)
    arr("variant_position", positions.astype(np.int32), ["variants"], vch)
    arr("variant_contig", np.zeros(nsites, dtype=np.int8), ["variants"], vch)
    sarr("variant_allele", alleles, ["variants", "alleles"],
         None if chunks is None else (chunks[0], alleles.shape[1]))
    sarr("variant_ancestral_allele", ancestral, ["variants"], vch)
    sarr("contig_id", np.array(["1"]), ["contigs"])
    arr("contig_length", np.array([seq_len], dtype=np.int64), ["contigs"])
    sarr("sample_id", np.asarray(sample_ids), ["samples"],
         None if chunks is None else (chunks[1],))
    for name, (data, dims) in (extra or {}).items():
        data = np.asarray(data)
        if data.dtype.kind in "OUS":
            sarr(name, data, dims)
        else:
            arr(name, data, dims)


def ts_to_arrays(ts, ploidy):
    n_ind = ts.num_samples // ploidy
    samples = ts.samples()
    max_alleles = max(len(v.alleles) for v in ts.variants())
    G = np.zeros((ts.num_sites, n_ind, ploidy), dtype=np.int8)
    alleles = np.empty((ts.num_sites, max_alleles), dtype=object)
    alleles[:] = ""
    anc = np.empty(ts.num_sites, dtype=object)
    pos = np.zeros(ts.num_sites, dtype=np.int32)
    for j, v in enumerate(ts.variants(samples=samples)):
        pos[j] = int(v.site.position)
        for k, a in enumerate(v.alleles):
            alleles[j, k] = a if a is not None else ""
        anc[j] = v.site.ancestral_state
        G[j] = v.genotypes.reshape(n_ind, ploidy)
    return G, pos, alleles, anc


def simulate(n, ploidy, seq_len, seed, rate=1e-8):
    ts = msprime.sim_ancestry(n, ploidy=ploidy, sequence_length=seq_len,
                              recombination_rate=rate, population_size=10_000,
                              random_seed=seed)
    return msprime.sim_mutations(ts, rate=rate, random_seed=seed)


def base_config(sample_path, anc_path, out, anc_state=None, **match_kw):
    if anc_state is None:
        anc_state = config.AncestralState(path=str(sample_path),
                                          field="variant_ancestral_allele")
    return config.Config(
        sources={
            "test": config.Source(path=str(sample_path), name="test"),
            "ancestors": config.Source(path=str(anc_path), name="ancestors",
                                       sample_time="sample_time"),
        },
        ancestors=[config.AncestorsConfig(name="ancestors", path=str(anc_path),
                                          sources=["test"],
                                          samples_chunk_size=7,
                                          variants_chunk_size=13)],
        match=config.MatchConfig(
            sources={
                "ancestors": config.MatchSourceConfig(node_flags=0,
                                                      create_individuals=False),
                "test": config.MatchSourceConfig(),
            },
            output=str(out), **match_kw),
        ancestral_state=anc_state,
        post_process=config.PostProcessConfig(),
    )


def run_case(name, cfg, num_threads=0):
    anc_cfg = cfg.ancestors[0]
    shutil.rmtree(anc_cfg.path, ignore_errors=True)
    tsinfer.infer_ancestors([cfg.sources[s] for s in anc_cfg.sources], anc_cfg,
                            cfg.ancestral_state, num_threads=num_threads)
    raw = pipeline.match(cfg, num_threads=num_threads)
    raw.dump(OUT / f"{name}_raw.trees")
    final = pipeline.post_process(raw, cfg)
    final.dump(OUT / f"{name}_final.trees")
    if cfg.augment_sites is not None:
        aug = pipeline.augment_sites(final, cfg)
        aug.dump(OUT / f"{name}_augmented.trees")
    return raw


# ---------------------------------------------------------------------------
# Case 1: haploid, simple
ts1 = simulate(16, 1, 100_000, 42)
G, pos, alleles, anc = ts_to_arrays(ts1, 1)
write_sample_vcz(OUT / "case1.vcz", G, pos, alleles, anc, 100_000,
                 [f"s{i}" for i in range(G.shape[1])])
cfg1 = base_config(OUT / "case1.vcz", OUT / "case1_anc.vcz", OUT / "case1.trees")
run_case("case1", cfg1)

# ---------------------------------------------------------------------------
# Case 2: diploid with missing data, chunked, populations, sample time field
ts2 = simulate(6, 2, 150_000, 7)
G, pos, alleles, anc = ts_to_arrays(ts2, 2)
rng = np.random.default_rng(3)
G[rng.random(G.shape) < 0.04] = -1
pops = np.array(["popA", "popB", "popA", "popB", "popA", "popB"])
write_sample_vcz(OUT / "case2.vcz", G, pos, alleles, anc, 150_000,
                 [f"ind{i}" for i in range(G.shape[1])],
                 extra={"sample_population": (pops, ["samples"]),
                        "sample_weight": (np.arange(6, dtype=np.float64) / 4,
                                          ["samples"])},
                 chunks=(17, 4))
cfg2 = base_config(OUT / "case2.vcz", OUT / "case2_anc.vcz", OUT / "case2.trees")
cfg2.individual_metadata = config.IndividualMetadataConfig(
    fields={"population_name": "sample_population", "weight": "sample_weight"},
    population="sample_population")
run_case("case2", cfg2, num_threads=2)

# ---------------------------------------------------------------------------
# Case 3: inference on a subset of sites, ancestral state from a separate
# store, path compression, augment_sites from the full store.
ts3 = simulate(10, 1, 120_000, 11)
G, pos, alleles, anc = ts_to_arrays(ts3, 1)
ids = [f"h{i}" for i in range(G.shape[1])]
write_sample_vcz(OUT / "case3_full.vcz", G, pos, alleles, anc, 120_000, ids)
keep = np.arange(len(pos)) % 3 != 0
write_sample_vcz(OUT / "case3_subset.vcz", G[keep], pos[keep], alleles[keep],
                 anc[keep], 120_000, ids)
cfg3 = base_config(OUT / "case3_subset.vcz", OUT / "case3_anc.vcz", OUT / "case3.trees",
                   anc_state=config.AncestralState(path=str(OUT / "case3_full.vcz"),
                                                   field="variant_ancestral_allele"),
                   path_compression=True)
cfg3.sources["full"] = config.Source(path=str(OUT / "case3_full.vcz"), name="full")
cfg3.augment_sites = config.AugmentSitesConfig(sources=["full"])
run_case("case3", cfg3)

# ---------------------------------------------------------------------------
# Low-level AncestorBuilder fixture
ab_ts = simulate(12, 1, 50_000, 5)
H = ab_ts.genotype_matrix()
keep = (H.sum(axis=1) > 1) & (H.sum(axis=1) < H.shape[1]) & (H.max(axis=1) == 1)
H = H[keep]
times = H.sum(axis=1) / H.shape[1]
ab_fixture = {"genotypes": H.tolist(), "times": times.tolist()}
for enc in (0, 1):
    ab = _tsinfer.AncestorBuilder(num_samples=H.shape[1], max_sites=H.shape[0],
                                  genotype_encoding=enc)
    for t, g in zip(times, H):
        ab.add_site(time=float(t), genotypes=g.astype(np.int8))
    descr = ab.ancestor_descriptors()
    ancs = []
    for t, focal in descr:
        a = np.zeros(H.shape[0], dtype=np.int8)
        s, e = ab.make_ancestor(focal, a)
        ancs.append({"time": t, "focal": focal.tolist(), "start": s, "end": e,
                     "haplotype": a.tolist()})
    ab_fixture[f"encoding_{enc}"] = ancs
(OUT / "ancestor_builder.json").write_text(json.dumps(ab_fixture))

# ---------------------------------------------------------------------------
# Low-level AncestorMatcher fixture against the case1 raw match output
raw = tskit.load(OUT / "case1_raw.trees")
tables = raw.dump_tables()
ll = _tsinfer.LightweightTableCollection(tables.sequence_length)
ll.fromdict(tables.asdict())
mi = _tsinfer.MatcherIndexes(ll)
ns = raw.num_sites
rng = np.random.default_rng(1)
haps = []
for k in range(8):
    h = rng.integers(0, 2, size=ns).astype(np.int8)
    if k % 2 == 1:
        h[rng.random(ns) < 0.1] = -1
    start = int(rng.integers(0, ns // 3))
    end = int(rng.integers(2 * ns // 3, ns + 1))
    for rho, mu in ((1e-2, 1e-20), (1e-3, 1e-2)):
        am = _tsinfer.AncestorMatcher(mi, recombination=np.full(ns, rho),
                                      mismatch=np.full(ns, mu))
        n, left, right, parent, match = am.find_path(h, start, end)
        haps.append({"haplotype": h.tolist(), "start": start, "end": end,
                     "rho": rho, "mu": mu, "left": left[:n].tolist(),
                     "right": right[:n].tolist(), "parent": parent[:n].tolist(),
                     "match": match[start:end].tolist()})
(OUT / "ancestor_matcher.json").write_text(json.dumps(haps))

with tarfile.open(OUT / "stores.tar.gz", "w:gz") as tar:
    for store in sorted(OUT.glob("*.vcz")):
        tar.add(store, arcname=store.name)
for store in OUT.glob("*.vcz"):
    shutil.rmtree(store)
print("fixtures written to", OUT)
