"""Cross-entity link mesh: a synonym/alias must NEVER shadow another entity's own
canonical name or ID (the bug that mislinked 1,214+ pages — Panobinostat →
/drug/vorinostat/, and generic cancers → their pediatric-only subtype pages).
Covers BOTH manifest builders: batch `_merge_manifest` (prod) and `upsert` (single-page)."""
import json
import os

import pytest

from atlas.page import links
from atlas.batch import _merge_manifest


def _rec(entity, slug, canonical, id_keys=(), synonyms=()):
    # name_keys mirrors the pipeline: canonical name first, then synonyms/aliases.
    return {"entity": entity, "slug": slug, "canonical": canonical,
            "id_keys": list(id_keys), "name_keys": [canonical, *synonyms]}


# The two colliding records: vorinostat lists "Panobinostat" as an alias (real bug).
_PANO = _rec("drug", "panobinostat", "Panobinostat", id_keys=["CHEMBL483254"])
_VORI = _rec("drug", "vorinostat", "Vorinostat", id_keys=["CHEMBL98"],
             synonyms=["Panobinostat", "suberoylanilide hydroxamic acid", "SID144206870"])
# Disease: pediatric subtype lists the generic cancer name as a synonym.
_LYMPH = _rec("disease", "lymphoma", "lymphoma", id_keys=["MONDO:0005062"])
_PED = _rec("disease", "pediatric-lymphoma", "pediatric lymphoma", id_keys=["MONDO:0006000"],
            synonyms=["lymphoma", "childhood lymphoma"])


@pytest.mark.parametrize("order", [
    [_PANO, _VORI, _LYMPH, _PED],      # authoritative owner built first
    [_VORI, _PANO, _PED, _LYMPH],      # alias-holder built first (the order that broke it)
])
def test_merge_manifest_alias_never_shadows_canonical(order, tmp_path):
    os.makedirs(tmp_path / "atlas")
    _merge_manifest(order, str(tmp_path))
    links.reset()
    links.load(str(tmp_path))
    # Each entity's own name resolves to ITS page, regardless of build order.
    assert links.drug_url(name="Panobinostat") == "/atlas/drug/panobinostat/"
    assert links.drug_url(name="Vorinostat") == "/atlas/drug/vorinostat/"
    assert links.disease_url(name="lymphoma") == "/atlas/disease/lymphoma/"
    assert links.disease_url(name="pediatric lymphoma") == "/atlas/disease/pediatric-lymphoma/"
    # A genuine alias with no authoritative owner still resolves (to its owner).
    assert links.drug_url(name="suberoylanilide hydroxamic acid") == "/atlas/drug/vorinostat/"
    # Registry-ID junk is NOT a resolvable alias.
    assert links.drug_url(name="SID144206870") is None
    # IDs resolve verbatim.
    assert links.drug_url(chembl_id="CHEMBL483254") == "/atlas/drug/panobinostat/"


def test_upsert_alias_never_shadows_canonical(tmp_path):
    links.reset()
    d = str(tmp_path)
    # vorinostat upserted AFTER panobinostat — the order that clobbered before.
    links.upsert(d, "drug", "panobinostat", id_keys=["CHEMBL483254"],
                 name_keys=["Panobinostat"], canonical="Panobinostat")
    links.upsert(d, "drug", "vorinostat", id_keys=["CHEMBL98"],
                 name_keys=["Vorinostat", "Panobinostat", "SID144206870"], canonical="Vorinostat")
    links.load(d)
    assert links.drug_url(name="Panobinostat") == "/atlas/drug/panobinostat/"
    assert links.drug_url(name="Vorinostat") == "/atlas/drug/vorinostat/"
    assert links.drug_url(name="SID144206870") is None


def test_is_name_junk():
    assert links._is_name_junk("sid144206870")
    assert links._is_name_junk("nsc125066")
    assert links._is_name_junk("12345")
    assert not links._is_name_junk("panobinostat")
    assert not links._is_name_junk("suberoylanilide hydroxamic acid")
