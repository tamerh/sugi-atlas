#!/usr/bin/env python3
"""Atlas drug §1-§12 collect + render. Reads typed params from
$ENJU_RUN_DIR/context.json; iteration var is `drug` (a record with `name`
and `slug` fields; `name` may be a ChEMBL id for unambiguous resolution)."""
import json, os
from atlas.drug import collect as DRC
from atlas.drug import render as DRR
from atlas.drug.anchors import resolve as resolve_anchors
from atlas.drug.slug import slugify

ctx = json.load(open(os.path.join(os.environ["ENJU_RUN_DIR"], "context.json")))
drug = ctx["iteration"]["drug"]
name = drug["name"] if isinstance(drug, dict) else str(drug)
slug = (drug.get("slug") if isinstance(drug, dict) else None) or slugify(name)

os.makedirs(f"build/{slug}", exist_ok=True)
a = resolve_anchors(name)
bundle = {sid: DRC.REGISTRY[sid].collect_fn(a) for sid in DRC.REGISTRY}
with open(f"build/{slug}/bundle.json", "w") as f:
    json.dump(bundle, f, indent=2, sort_keys=True, default=str)
body = DRR.render_all(bundle)
open(f"build/{slug}/body.md", "w").write(body)
open(f"build/{slug}/anchors_meta.json", "w").write(json.dumps({
    "name": name, "slug": slug, "chembl_id": a.chembl_id,
    "canonical_name": a.canonical_name, "molecule_type": a.molecule_type,
    "max_phase": a.max_phase,
}, indent=2))
print(f"collect_render {slug}: targets={len(a.targets)} "
      f"bundle={len(json.dumps(bundle, default=str))}c body={len(body)}c")
