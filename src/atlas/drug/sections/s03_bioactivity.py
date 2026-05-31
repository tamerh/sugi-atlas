"""§3 — bioactivity: the drug's own ChEMBL activity rows
(>>chembl_molecule>>chembl_activity). pchembl is the standardized potency
metric (-log10(M): 10 = 0.1 nM, 6 = 1 µM). Keep rows at pchembl ≥ 5 (real
binding), sorted by potency, top 30. Empty for biologics (no small-molecule
assay rows) — block elides."""
from atlas.biobtree import map_all
from atlas.section import Section


def _pchembl(r):
    try:
        return float(r.get("pchembl") or 0)
    except (TypeError, ValueError):
        return 0.0


def collect(a):
    rows = map_all(a.chembl_id, ">>chembl_molecule>>chembl_activity")
    potent = [r for r in rows if _pchembl(r) >= 5.0]
    potent.sort(key=_pchembl, reverse=True)
    return {
        "section": "03_bioactivity",
        "activity_total": len(rows),
        "potent_count": len(potent),
        "activities": [{
            "id": r.get("id"),
            "type": r.get("standard_type"),
            "value": r.get("standard_value"),
            "unit": r.get("standard_units"),
            "pchembl": r.get("pchembl"),
        } for r in potent[:30]],
    }


SECTION = Section(
    id="3", name="bioactivity",
    description=("ChEMBL bioactivity for the drug (chembl_molecule→chembl_activity); "
                 "pchembl-ranked, ≥5 (real binding), top 30 by potency"),
    needs=("chembl_id",),
    produces=("activity_total", "potent_count", "activities"),
    datasets=("chembl_molecule", "chembl_activity"),
    chains=(">>chembl_molecule>>chembl_activity",),
    collect_fn=collect,
)
