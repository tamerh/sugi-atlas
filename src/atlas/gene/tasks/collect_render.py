#!/usr/bin/env python3
"""Atlas §1-12 collect + render. Reads typed params from $ENJU_RUN_DIR/context.json."""
import json, os
from atlas.gene import collect as C
from atlas.gene import render as R
from atlas.biobtree import CALLS
from atlas.pipeline import datasets_from_calls

ctx = json.load(open(os.path.join(os.environ["ENJU_RUN_DIR"], "context.json")))
symbol = ctx["iteration"]["symbol"]

os.makedirs(f"build/{symbol}", exist_ok=True)
CALLS.clear()  # capture the datasets actually queried for this page
bundle = {s: C.SECTIONS[s](symbol) for s in C.SECTIONS}
with open(f"build/{symbol}/bundle.json", "w") as f:
    json.dump(bundle, f, indent=2, sort_keys=True)
json.dump(datasets_from_calls(CALLS), open(f"build/{symbol}/datasets.json", "w"))
body = "\n\n".join(R.RENDER[s](bundle[s]) for s in R.RENDER)
open(f"build/{symbol}/body.md", "w").write(body)
print(f"collect_render {symbol}: bundle={len(json.dumps(bundle))}c body={len(body)}c")
