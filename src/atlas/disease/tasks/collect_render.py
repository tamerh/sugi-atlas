#!/usr/bin/env python3
"""Atlas disease §1-§14 collect + render. Reads typed params from
$ENJU_RUN_DIR/context.json; iteration var is `disease` (a record with
`name` and `slug` fields)."""
import json, os
from atlas.disease import collect as DC
from atlas.disease import render as DR
from atlas.disease.anchors import resolve as resolve_anchors
from atlas.disease.slug import slugify
from atlas.biobtree import CALLS
from atlas.pipeline import datasets_from_calls

ctx = json.load(open(os.path.join(os.environ["ENJU_RUN_DIR"], "context.json")))
disease = ctx["iteration"]["disease"]
name = disease["name"] if isinstance(disease, dict) else str(disease)
# Trust caller's slug if provided, else compute. Either way, normalise.
slug = (disease.get("slug") if isinstance(disease, dict) else None) or slugify(name)

os.makedirs(f"build/{slug}", exist_ok=True)
CALLS.clear()  # capture datasets actually queried
a = resolve_anchors(name)
bundle = {sid: DC.REGISTRY[sid].collect_fn(a) for sid in DC.REGISTRY}
with open(f"build/{slug}/bundle.json", "w") as f:
    json.dump(bundle, f, indent=2, sort_keys=True, default=str)
json.dump(datasets_from_calls(CALLS), open(f"build/{slug}/datasets.json", "w"))
body = DR.render_all(bundle)
open(f"build/{slug}/body.md", "w").write(body)
# Persist resolved name + slug + mondo_id so downstream tasks don't re-resolve.
open(f"build/{slug}/anchors_meta.json", "w").write(json.dumps({
    "name": name, "slug": slug, "mondo_id": a.mondo_id,
    "canonical_name": a.canonical_name, "is_cancer": a.is_cancer,
}, indent=2))
print(f"collect_render {slug}: cohort={len(a.cohort)} "
      f"bundle={len(json.dumps(bundle, default=str))}c body={len(body)}c")
