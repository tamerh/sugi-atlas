"""JSON-LD sidecar + inline-script validity over the built corpus."""
import json

import pytest

from ._harness import report

pytestmark = pytest.mark.integration

_TYPE = {"gene": "Gene", "disease": "MedicalCondition", "drug": "Drug"}


def test_sidecar_parses_and_typed(pages):
    bad = []
    for p in pages:
        j = None
        try:
            j = p.jsonld()
        except (json.JSONDecodeError, ValueError) as e:
            bad.append(f"{p.entity}/{p.slug}: unparseable ({e})")
            continue
        if j is None:
            bad.append(f"{p.entity}/{p.slug}: no entity.jsonld")
            continue
        if j.get("@type") != _TYPE[p.entity]:
            bad.append(f"{p.entity}/{p.slug}: @type={j.get('@type')}")
        if not j.get("@id") or not j.get("identifier"):
            bad.append(f"{p.entity}/{p.slug}: missing @id/identifier")
    assert not bad, report(bad)


def test_gene_protein_nodes_have_id_and_reciprocal_edge(pages):
    """Protein-coding genes emit Protein node(s) with a #protein-<acc> @id and
    the reciprocal isEncodedByBioChemEntity (layer A)."""
    bad = []
    for p in pages:
        if p.entity != "gene":
            continue
        j = p.jsonld() or {}
        enc = j.get("encodesBioChemEntity")
        if not enc:
            continue                                    # ncRNA — no product
        for node in (enc if isinstance(enc, list) else [enc]):
            if "#protein-" not in (node.get("@id") or ""):
                bad.append(f"gene/{p.slug}: protein @id={node.get('@id')}")
            if not (node.get("isEncodedByBioChemEntity") or {}).get("@id"):
                bad.append(f"gene/{p.slug}: no isEncodedByBioChemEntity")
    assert not bad, report(bad)


def test_inline_script_jsonld_parses(pages):
    """The compacted inline <script type=application/ld+json> in the body
    parses (audit #6 capping must not produce invalid JSON)."""
    bad = []
    tag = '<script type="application/ld+json">'
    for p in pages:
        if tag not in p.body:
            bad.append(f"{p.entity}/{p.slug}: no inline JSON-LD")
            continue
        block = p.body.split(tag, 1)[1].split("</script>", 1)[0].strip()
        try:
            json.loads(block)
        except (json.JSONDecodeError, ValueError) as e:
            bad.append(f"{p.entity}/{p.slug}: inline JSON-LD invalid ({e})")
    assert not bad, report(bad)
