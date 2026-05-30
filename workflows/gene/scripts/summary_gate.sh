#!/usr/bin/env bash
set -euo pipefail
SYMBOL="${ENJU_PARAM_symbol:-${SYMBOL:?need SYMBOL}}"
SUMMARY_MODEL="${ENJU_PARAM_summary_model:-${SUMMARY_MODEL:-qwen/qwen3-235b-a22b-2507|Together}}"
SKIP_SUMMARY="${ENJU_PARAM_skip_summary:-${SKIP_SUMMARY:-false}}"
DIST_ROOT="${ENJU_PARAM_dist_root:-${DIST_ROOT:-/data/sugi-atlas-dist}}"
cd "$ENJU_PROJECT_DIR"
if [ "${SKIP_SUMMARY:-false}" = "true" ]; then
    echo "summary_gate $SYMBOL: SKIPPED"
    echo '{"verdict":"skipped"}' > "build/$SYMBOL/judge.json"
    exit 0
fi
python - <<PY
import json
from atlas.validation import summary_gate
body = open("build/$SYMBOL/body.md").read()
summary = open("build/$SYMBOL/summary.md").read()
r = summary_gate.check_summary(body, summary)
with open("build/$SYMBOL/judge.json", "w") as f:
    json.dump(r, f, indent=2)
print(f"summary_gate $SYMBOL: {r['verdict']} both={len(r.get('both', []))} single={len(r.get('single', []))}")
PY
