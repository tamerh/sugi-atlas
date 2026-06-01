#!/usr/bin/env python3
"""Assemble Hugo-frontmatter page.md and write to $dist_root/atlas/drug/<slug>/.

Sidecars (parity with gene/disease publish):
  entity.jsonld    — schema.org Drug (sameAs ChEMBL/PubChem/ChEBI/ATC + target
                     Gene nodes + treats MedicalCondition nodes)
  provenance.json  — schema.org Dataset (per-section datasets + chains +
                     upstream-source URLs)
"""
import json, os, shutil
from datetime import datetime, timezone
from atlas import __version__ as V
from atlas.pipeline import assemble_page, biobtree_version, GENERATED_BY, datasets_union
from atlas.drug.collect import REGISTRY
from atlas.drug.slug import slugify
from atlas.page.drug_jsonld import build_jsonld, as_jsonld_string

ctx = json.load(open(os.path.join(os.environ["ENJU_RUN_DIR"], "context.json")))
drug = ctx["iteration"]["drug"]
name = drug["name"] if isinstance(drug, dict) else str(drug)
slug = (drug.get("slug") if isinstance(drug, dict) else None) or slugify(name)
dist_root = ctx["params"].get("dist_root") or "/data/sugi-atlas-dist"

out = f"{dist_root}/atlas/drug/{slug}"
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
    "entity_type": "drug",
    "generated_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
    "atlas_version": V,
    "biobtree_version": biobtree_version(),
    "generated_by": GENERATED_BY,
    "datasets": datasets_union(REGISTRY),
    "summary_model": display_model,
}
bundle = json.load(open(f"build/{slug}/bundle.json"))
page = assemble_page(slug, summary, body, meta, bundle=bundle)
open(f"{out}/page.md", "w").write(page)
open(f"{out}/entity.jsonld", "w").write(as_jsonld_string(build_jsonld(bundle, slug)))
# bundle.json + provenance.json intentionally NOT published (data dump / api
# trail kept internal); transparency = frontmatter datasets + generated_by.
for f in ("summary.md",):
    src = f"build/{slug}/{f}"
    if os.path.exists(src):
        shutil.copy(src, f"{out}/{f}")
print(f"publish {slug}: page={len(page)}c -> {out}")
