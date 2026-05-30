#!/usr/bin/env bash
set -euo pipefail
cd "$ENJU_PROJECT_DIR"
mkdir -p "build/$SYMBOL"
python - <<PY
import json
from atlas.gene import collect as C
from atlas.gene import render as R
bundle = {s: C.SECTIONS[s]("$SYMBOL") for s in C.SECTIONS}
with open("build/$SYMBOL/bundle.json", "w") as f:
    json.dump(bundle, f, indent=2, sort_keys=True)
body = "\n\n".join(R.RENDER[s](bundle[s]) for s in R.RENDER)
open("build/$SYMBOL/body.md", "w").write(body)
print(f"collect_render $SYMBOL: bundle={len(json.dumps(bundle))}c body={len(body)}c")
PY
