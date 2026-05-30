#!/usr/bin/env python3
"""Assemble Hugo-frontmatter page.md and write to $dist_root/atlas/gene/<symbol>/."""
import json, os, shutil
from datetime import datetime, timezone
from atlas import __version__ as V
from atlas.pipeline import assemble_page, biobtree_version

ctx = json.load(open(os.path.join(os.environ["ENJU_RUN_DIR"], "context.json")))
symbol = ctx["iteration"]["symbol"]
dist_root = ctx["params"].get("dist_root") or "/data/sugi-atlas-dist"

out = f"{dist_root}/atlas/gene/{symbol}"
os.makedirs(out, exist_ok=True)

bundle  = json.load(open(f"build/{symbol}/bundle.json"))
body    = open(f"build/{symbol}/body.md").read()
summary = open(f"build/{symbol}/summary.md").read().strip()
meta = {
    "title": symbol, "symbol": symbol, "entity_type": "gene",
    "generated_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
    "atlas_version": V, "biobtree_version": biobtree_version(),
}
page = assemble_page(symbol, summary, body, meta)
open(f"{out}/page.md", "w").write(page)
for f in ("bundle.json", "body.md", "summary.md", "judge.json", "body_gate.json"):
    src = f"build/{symbol}/{f}"
    if os.path.exists(src):
        shutil.copy(src, f"{out}/{f}")
print(f"publish {symbol}: page={len(page)}c -> {out}")
