#!/usr/bin/env bash
set -euo pipefail
SYMBOL="${SYMBOL:-${ENJU_PARAM_symbol:?need SYMBOL or ENJU_PARAM_symbol}}"
SUMMARY_MODEL="${SUMMARY_MODEL:-${ENJU_PARAM_summary_model:-qwen/qwen3-235b-a22b-2507|Together}}"
SKIP_SUMMARY="${SKIP_SUMMARY:-${ENJU_PARAM_skip_summary:-false}}"
DIST_ROOT="${DIST_ROOT:-${ENJU_PARAM_dist_root:-/data/sugi-atlas-dist}}"
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
