"""Sparse-subtype parent-drug fallback: only a parent that is a specific clinical
entity (Orphanet/GARD/OMIM xref) lends its drugs to a subtype. v1.11.7 rendered
"hereditary disease" / "syndromic disease" drugs on ~1,900 unrelated pages."""
import json

from atlas import batch
from atlas.disease.render import parent_is_clinical_entity


def test_grouping_classes_rejected():
    assert not parent_is_clinical_entity({"mondochild": 1921, "clinvar": 6, "mesh": 1})  # hereditary disease
    assert not parent_is_clinical_entity({"gencc": 1, "icd10cm": 1})                     # eye / nervous system disorder
    assert not parent_is_clinical_entity({})
    assert not parent_is_clinical_entity(None)


def test_clinical_entities_accepted():
    assert parent_is_clinical_entity({"orphanet": 1, "gard": 1, "mim": 1})   # Noonan syndrome
    assert parent_is_clinical_entity({"mim": 1, "gencc": 3})                 # myopia
    assert parent_is_clinical_entity({"orphanet": 1})                        # alopecia


def test_batch_gate_reads_parent_cache(tmp_path):
    d = tmp_path / "disease"; d.mkdir()
    (d / "noonan-syndrome.json").write_text(json.dumps({"bundle": {"1": {"xref_counts": {"orphanet": 1}}}}))
    (d / "hereditary-disease.json").write_text(json.dumps({"bundle": {"1": {"xref_counts": {"mesh": 1}}}}))
    assert batch._parent_is_clinical(str(tmp_path), "noonan-syndrome")
    assert not batch._parent_is_clinical(str(tmp_path), "hereditary-disease")
    assert not batch._parent_is_clinical(str(tmp_path), "missing-parent")      # no cache → no fallback
