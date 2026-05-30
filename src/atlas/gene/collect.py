#!/usr/bin/env python3
"""Gene-pipeline orchestrator. Thin façade over the per-section modules.

  - `resolve(symbol)` once -> Anchors (atlas.gene.anchors.resolve)
  - each section's collect_fn(anchors) -> section bundle
  - `SECTIONS[<id>](symbol) -> bundle` is preserved for backward compatibility
    with coverage/stress/bench/body_gate.

For full bundle assembly use `collect_all(symbol)` which resolves anchors once
and runs all 12 sections.

  python3 -m atlas.gene.collect TP53 7       # collect §7 for TP53
  python3 -m atlas.gene.collect TP53 all     # collect every section
"""
import json, sys
from atlas.biobtree import (
    search, entry, bbmap, map_all, rows, map_targets, CALLS, xref_counts,
)
from atlas.gene.anchors import Anchors, resolve as resolve_anchors, resolve_hgnc
from atlas.gene.sections import REGISTRY, Section

# Section ID -> callable(symbol) -> bundle  (backward-compatible API).
def _run(sid, symbol):
    return REGISTRY[sid].collect_fn(resolve_anchors(symbol))

SECTIONS = {sid: (lambda s, sid=sid: _run(sid, s)) for sid in REGISTRY}

# Individual collect_* names — used by atlas.validation.stress and others.
collect_gene_ids      = SECTIONS["1"]
collect_transcripts   = SECTIONS["2"]
collect_protein_ids   = SECTIONS["3"]
collect_structure     = SECTIONS["4"]
collect_orthologs     = SECTIONS["5"]
collect_variants      = SECTIONS["6"]
collect_pathways      = SECTIONS["7"]
collect_interactions  = SECTIONS["8"]
collect_tf_regulation = SECTIONS["9"]
collect_drugs         = SECTIONS["10"]
collect_expression    = SECTIONS["11"]
collect_diseases      = SECTIONS["12"]

def collect_all(symbol):
    """Resolve anchors once, run all 12 sections, return {section_id: bundle}."""
    a = resolve_anchors(symbol)
    return {sid: REGISTRY[sid].collect_fn(a) for sid in REGISTRY}

if __name__ == "__main__":
    sym = sys.argv[1] if len(sys.argv) > 1 else "TP53"
    sec = sys.argv[2] if len(sys.argv) > 2 else "1"
    if sec == "all":
        out = collect_all(sym)
    else:
        out = SECTIONS[sec](sym)
    print(json.dumps(out, indent=2))
    print(f"\n--- {len(CALLS)} api calls ---", file=sys.stderr)
    for c in CALLS:
        print(f"  {c['path']}({c['params']})", file=sys.stderr)
