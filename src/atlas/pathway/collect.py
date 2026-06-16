#!/usr/bin/env python3
"""Pathway-pipeline orchestrator. Resolves the anchor once and assembles the
single page bundle. Pathways are structurally simple (members + hierarchy), so
unlike gene/drug/disease there's one collect, not a section registry.

  python -m atlas.pathway.collect R-HSA-109581
  python -m atlas.pathway.collect Apoptosis
"""
import json
import sys

from atlas.pathway.anchors import resolve as resolve_anchors


def collect(id_or_name):
    a = resolve_anchors(id_or_name)
    return {
        "reactome_id": a.reactome_id,
        "name": a.name,
        "members": a.members,
        "member_count": len(a.members),
        "parent": a.parent,
        "children": a.children,
        "child_count": len(a.children),
        "go_id": a.go_id,
    }


def main(argv=None):
    a = argv or sys.argv[1:]
    if not a:
        print("usage: python -m atlas.pathway.collect <R-HSA-id|name>")
        sys.exit(2)
    print(json.dumps(collect(a[0]), indent=1))


if __name__ == "__main__":
    main()
