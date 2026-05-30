#!/usr/bin/env python3
"""Body-gate regression check. Snapshot dir = <dist>/snapshots/disease/.
Exits non-zero on regression so Enju marks the task failed_retryable."""
import json, os, sys
from atlas.validation import body_gate
from atlas.disease.slug import slugify

ctx = json.load(open(os.path.join(os.environ["ENJU_RUN_DIR"], "context.json")))
disease = ctx["iteration"]["disease"]
slug = (disease.get("slug") if isinstance(disease, dict) else None) \
       or slugify(disease["name"] if isinstance(disease, dict) else str(disease))
dist_root = ctx["params"].get("dist_root") or "/data/sugi-atlas-dist"
snap_dir = body_gate.snap_dir_for(dist_root, "disease")

bundle = json.load(open(f"build/{slug}/bundle.json"))
r = body_gate.check(slug, bundle, snap_dir)
with open(f"build/{slug}/body_gate.json", "w") as f:
    json.dump(r, f, indent=2)
print(f"body_gate {slug}: {r['verdict']} ({r['summary']})")
sys.exit(0 if r["verdict"] in {"clean", "drift", "first_run"} else 2)
