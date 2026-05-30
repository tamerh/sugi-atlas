#!/usr/bin/env python3
"""Disease-pipeline orchestrator. Thin façade over the per-section modules.

  - `resolve(name)` once -> DiseaseAnchors (atlas.disease.anchors.resolve)
  - each section's collect_fn(anchors) -> section bundle
  - `SECTIONS[<id>](name) -> bundle` for backward-compat / bench helpers

  python3 -m atlas.disease.collect "age-related macular degeneration" 1
  python3 -m atlas.disease.collect "endometrial cancer" all
"""
import json, sys
from atlas.biobtree import CALLS
from atlas.disease.anchors import DiseaseAnchors, resolve as resolve_anchors
from atlas.disease.sections import REGISTRY

def _run(sid, name):
    return REGISTRY[sid].collect_fn(resolve_anchors(name))

SECTIONS = {sid: (lambda s, sid=sid: _run(sid, s)) for sid in REGISTRY}

def collect_all(name):
    """Resolve anchors once, run all sections, return {section_id: bundle}."""
    a = resolve_anchors(name)
    return {sid: REGISTRY[sid].collect_fn(a) for sid in REGISTRY}

if __name__ == "__main__":
    name = sys.argv[1] if len(sys.argv) > 1 else "age-related macular degeneration"
    sec = sys.argv[2] if len(sys.argv) > 2 else "1"
    out = collect_all(name) if sec == "all" else SECTIONS[sec](name)
    print(json.dumps(out, indent=2, default=str))
    print(f"\n--- {len(CALLS)} api calls ---", file=sys.stderr)
