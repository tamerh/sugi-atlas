#!/usr/bin/env python3
"""Assemble Hugo-frontmatter page.md and write to $dist_root/atlas/disease/<slug>/.

No JSON-LD or provenance sidecar yet — those are gene-specific (build_jsonld
expects HGNC anchors). Disease equivalents are a follow-up polish item."""
import json, os, shutil
from datetime import datetime, timezone
from atlas import __version__ as V
from atlas.pipeline import assemble_page, biobtree_version
from atlas.disease.slug import slugify

ctx = json.load(open(os.path.join(os.environ["ENJU_RUN_DIR"], "context.json")))
disease = ctx["iteration"]["disease"]
name = disease["name"] if isinstance(disease, dict) else str(disease)
slug = (disease.get("slug") if isinstance(disease, dict) else None) or slugify(name)
dist_root = ctx["params"].get("dist_root") or "/data/sugi-atlas-dist"

out = f"{dist_root}/atlas/disease/{slug}"
os.makedirs(out, exist_ok=True)

body    = open(f"build/{slug}/body.md").read()
summary = open(f"build/{slug}/summary.md").read().strip()
summary_model = ctx["params"].get("summary_model") or "qwen/qwen3-235b-a22b-2507|Together"
display_model = summary_model.split("|", 1)[0].split("/", 1)[-1]

am_path = f"build/{slug}/anchors_meta.json"
am = json.load(open(am_path)) if os.path.exists(am_path) else {}
title = am.get("canonical_name") or name

meta = {
    "title": title,
    "symbol": slug,
    "entity_type": "disease",
    "generated_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
    "atlas_version": V,
    "biobtree_version": biobtree_version(),
    "summary_model": display_model,
}
# bundle=None → assemble_page skips the gene-specific lead (jsonld + declarative
# sentence). Disease-equivalent lead is a follow-up.
page = assemble_page(slug, summary, body, meta, bundle=None)
open(f"{out}/page.md", "w").write(page)

for f in ("bundle.json", "body.md", "summary.md", "body_gate.json",
          "anchors_meta.json"):
    src = f"build/{slug}/{f}"
    if os.path.exists(src):
        shutil.copy(src, f"{out}/{f}")
print(f"publish {slug}: page={len(page)}c -> {out}")
