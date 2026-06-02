"""Cross-entity mesh integrity over the built corpus."""
import json
import os

import pytest

from ._harness import ATLAS, _INTERNAL_LINK, report

pytestmark = pytest.mark.integration


def test_no_dangling_internal_links(pages):
    """Every /atlas/<entity>/<slug>/ link points at a page that exists in the
    dist — the mesh must never dead-end."""
    bad = []
    for p in pages:
        for _full, entity, slug in _INTERNAL_LINK.findall(p.body):
            if not os.path.isdir(os.path.join(ATLAS, entity, slug)):
                bad.append(f"{p.entity}/{p.slug} → /atlas/{entity}/{slug}/ (missing)")
    assert not bad, report(bad)


def test_drug_mesh_labels_not_shouting(pages):
    """Drug links in the Related block render de-SHOUTed (Cisplatin, not
    CISPLATIN)."""
    bad = []
    for p in pages:
        if "{#related}" not in p.body:
            continue
        rel = p.body.split("{#related}", 1)[-1]
        for full, entity, slug in _INTERNAL_LINK.findall(rel):
            if entity != "drug":
                continue
            idx = rel.find(f"]({full})")          # within the Related block only
            label = rel[:idx].rsplit("[", 1)[-1] if idx >= 0 else ""
            if len(label) > 3 and label == label.upper():
                bad.append(f"{p.entity}/{p.slug}: SHOUTING drug label {label!r}")
    assert not bad, report(bad)


def test_manifest_and_reverse_index_parse(dist_root):
    for name in ("manifest.json", "reverse_edges.json"):
        path = os.path.join(dist_root, "atlas", name)
        assert os.path.exists(path), f"missing {name}"
        with open(path) as f:
            data = json.load(f)
        assert isinstance(data, dict) and data, f"{name} empty/invalid"
