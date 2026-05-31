#!/usr/bin/env python3
"""Drug-pipeline orchestrator. Thin façade over the per-section modules.

  - `resolve(name)` once -> DrugAnchors (atlas.drug.anchors.resolve)
  - each section's collect_fn(anchors) -> section bundle
  - `SECTIONS[<id>](name)` is preserved for backward compatibility with
    body_gate / bench, mirroring the gene + disease orchestrators.

  python3 -m atlas.drug.collect Imatinib 1     # collect §1 for Imatinib
  python3 -m atlas.drug.collect Imatinib all    # collect every wired section
"""
import json, sys
from atlas.drug.anchors import resolve as resolve_anchors
from atlas.drug.sections import REGISTRY


def _run(sid, name):
    return REGISTRY[sid].collect_fn(resolve_anchors(name))


SECTIONS = {sid: (lambda n, sid=sid: _run(sid, n)) for sid in REGISTRY}


def collect_all(name):
    """Resolve anchors once, run every wired section. Returns
    {section_id: bundle}."""
    a = resolve_anchors(name)
    return {sid: REGISTRY[sid].collect_fn(a) for sid in REGISTRY}


def main(argv=None):
    a = argv or sys.argv[1:]
    if len(a) < 2:
        print("usage: python -m atlas.drug.collect <name|CHEMBLid> <section|all>")
        sys.exit(2)
    name, sec = a[0], a[1]
    out = collect_all(name) if sec == "all" else {sec: SECTIONS[sec](name)}
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
