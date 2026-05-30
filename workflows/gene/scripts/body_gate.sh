#!/usr/bin/env bash
set -euo pipefail
cd "$ENJU_PROJECT_DIR"
python - <<PY
import json, sys
from atlas.validation import body_gate
bundle = json.load(open("build/$SYMBOL/bundle.json"))
r = body_gate.check("$SYMBOL", bundle)
with open("build/$SYMBOL/body_gate.json", "w") as f:
    json.dump(r, f, indent=2)
print(f"body_gate $SYMBOL: {r['verdict']} ({r['summary']})")
sys.exit(0 if r["verdict"] in {"clean", "drift", "first_run"} else 2)
PY
