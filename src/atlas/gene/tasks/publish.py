#!/usr/bin/env python3
"""Assemble Hugo-frontmatter page.md and write to $dist_root/atlas/gene/<symbol>/."""
import json, os, shutil
from datetime import datetime, timezone
from atlas import __version__ as V
from atlas.pipeline import assemble_page, biobtree_version
from atlas.page.jsonld import build_jsonld, as_jsonld_string

ctx = json.load(open(os.path.join(os.environ["ENJU_RUN_DIR"], "context.json")))
symbol = ctx["iteration"]["symbol"]
dist_root = ctx["params"].get("dist_root") or "/data/sugi-atlas-dist"

out = f"{dist_root}/atlas/gene/{symbol}"
os.makedirs(out, exist_ok=True)

bundle  = json.load(open(f"build/{symbol}/bundle.json"))
body    = open(f"build/{symbol}/body.md").read()
summary = open(f"build/{symbol}/summary.md").read().strip()
summary_model = ctx["params"].get("summary_model") or "qwen/qwen3-235b-a22b-2507|Together"
# strip provider suffix + path for the human-facing disclosure line
display_model = summary_model.split("|", 1)[0].split("/", 1)[-1]
meta = {
    "title": symbol, "symbol": symbol, "entity_type": "gene",
    "generated_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
    "atlas_version": V, "biobtree_version": biobtree_version(),
    "summary_model": display_model,
}
page = assemble_page(symbol, summary, body, meta, bundle=bundle)
open(f"{out}/page.md", "w").write(page)

# schema.org Gene JSON-LD sidecar — machine clients fetch this directly via
# the <link rel="alternate" type="application/ld+json"> hint (added by Hugo
# theme; the file itself lives here next to page.md).
open(f"{out}/entity.jsonld", "w").write(as_jsonld_string(build_jsonld(bundle)))

for f in ("bundle.json", "body.md", "summary.md", "judge.json", "body_gate.json"):
    src = f"build/{symbol}/{f}"
    if os.path.exists(src):
        shutil.copy(src, f"{out}/{f}")
print(f"publish {symbol}: page={len(page)}c -> {out}")
