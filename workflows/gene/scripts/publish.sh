#!/usr/bin/env bash
set -euo pipefail
SYMBOL="${SYMBOL:-${ENJU_PARAM_symbol:?need SYMBOL or ENJU_PARAM_symbol}}"
SUMMARY_MODEL="${SUMMARY_MODEL:-${ENJU_PARAM_summary_model:-qwen/qwen3-235b-a22b-2507|Together}}"
SKIP_SUMMARY="${SKIP_SUMMARY:-${ENJU_PARAM_skip_summary:-false}}"
DIST_ROOT="${DIST_ROOT:-${ENJU_PARAM_dist_root:-/data/sugi-atlas-dist}}"
cd "$ENJU_PROJECT_DIR"
OUT="$DIST_ROOT/atlas/gene/$SYMBOL"
mkdir -p "$OUT"
python - <<PY
import json, os, shutil
from datetime import datetime, timezone
from atlas import __version__ as V
from atlas.pipeline import assemble_page, biobtree_version
bundle  = json.load(open("build/$SYMBOL/bundle.json"))
body    = open("build/$SYMBOL/body.md").read()
summary = open("build/$SYMBOL/summary.md").read().strip()
meta = {
    "title": "$SYMBOL", "symbol": "$SYMBOL", "entity_type": "gene",
    "generated_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
    "atlas_version": V, "biobtree_version": biobtree_version(),
}
page = assemble_page("$SYMBOL", summary, body, meta)
out = "$OUT"
open(out + "/page.md", "w").write(page)
for f in ("bundle.json", "body.md", "summary.md", "judge.json", "body_gate.json"):
    src = f"build/$SYMBOL/{f}"
    if os.path.exists(src):
        shutil.copy(src, out + "/" + f)
print(f"publish $SYMBOL: page={len(page)}c -> {out}")
PY
