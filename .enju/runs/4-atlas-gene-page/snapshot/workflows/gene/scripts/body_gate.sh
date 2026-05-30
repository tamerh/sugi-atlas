#!/usr/bin/env bash
set -euo pipefail
SYMBOL="${ENJU_PARAM_symbol:-${SYMBOL:?need SYMBOL}}"
SUMMARY_MODEL="${ENJU_PARAM_summary_model:-${SUMMARY_MODEL:-qwen/qwen3-235b-a22b-2507|Together}}"
SKIP_SUMMARY="${ENJU_PARAM_skip_summary:-${SKIP_SUMMARY:-false}}"
DIST_ROOT="${ENJU_PARAM_dist_root:-${DIST_ROOT:-/data/sugi-atlas-dist}}"
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
