"""Reference results of vcztools variant filtering (development tool)."""
import json, sys
import zarr
from tsinfer import vcz

store = zarr.open(sys.argv[1], mode="r")
cases = [
    {"include": "POS > 3000"},
    {"include": "POS > 3000 && POS < 8000"},
    {"include": "POS < 1000 || POS >= 9000"},
    {"include": "POS*2 >= 10000 | POS < 500"},
    {"include": 'TYPE="snp"'},
    {"include": 'FILTER="PASS"'},
    {"include": 'CHROM="1"'},
    {"exclude": "POS < 5000"},
    {"exclude": "(POS > 2000 && POS < 4000) || POS > 9000"},
    {"include": "-POS > -5000"},
    {"regions": "1:1000-5000"},
    {"regions": "1:3000-"},
    {"targets": "1:1000-5000"},
    {"targets": "^1:1000-5000"},
    {"regions": "1"},
]
out = []
for c in cases:
    pos = [int(v["variant_position"]) for v in vcz.iter_variants(store, fields=["variant_position"], **c)]
    out.append({"args": c, "positions": pos})
json.dump(out, open(sys.argv[2], "w"))
print(len(out), "cases")
