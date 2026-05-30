#!/usr/bin/env bash
set -euo pipefail
SYMBOL="${SYMBOL:-${ENJU_PARAM_symbol:?need SYMBOL or ENJU_PARAM_symbol}}"
SUMMARY_MODEL="${SUMMARY_MODEL:-${ENJU_PARAM_summary_model:-qwen/qwen3-235b-a22b-2507|Together}}"
SKIP_SUMMARY="${SKIP_SUMMARY:-${ENJU_PARAM_skip_summary:-false}}"
DIST_ROOT="${DIST_ROOT:-${ENJU_PARAM_dist_root:-/data/sugi-atlas-dist}}"
cd "$ENJU_PROJECT_DIR"
if [ "${SKIP_SUMMARY:-false}" = "true" ]; then
    echo "summary $SYMBOL: SKIPPED"
    : > "build/$SYMBOL/summary.md"
    exit 0
fi
python - <<PY
from atlas.bench import summary as B
body = open("build/$SYMBOL/body.md").read()
key = B.api_key()
prompt = B.INSTR.format(g="$SYMBOL") + "\n\nBODY:\n" + body
d, dt = B.call("$SUMMARY_MODEL", prompt, key, max_tokens=600)
if "choices" not in d:
    raise SystemExit("API: " + str(d)[:200])
txt = str(d["choices"][0]["message"].get("content") or "").strip()
open("build/$SYMBOL/summary.md", "w").write(txt + "\n")
print(f"summary $SYMBOL: {len(txt)}c in {dt:.1f}s ($SUMMARY_MODEL)")
PY
