"""§3 — bioactivity: the drug's own ChEMBL activity rows
(>>chembl_molecule>>chembl_activity). pchembl is the standardized potency
metric (-log10(M): 10 = 0.1 nM, 6 = 1 µM). Keep rows at pchembl ≥ 5 (real
binding), sorted by potency, top 50. Empty for biologics (no small-molecule
assay rows) — block elides."""
from atlas.biobtree import map_all, entry
from atlas.section import Section


def _pchembl(r):
    try:
        return float(r.get("pchembl") or 0)
    except (TypeError, ValueError):
        return 0.0


_TGT_CACHE = {}


def _activity_target(act_id):
    """Resolve a ChEMBL activity → its target gene symbol. The chembl_activity
    map row omits the target, but the activity links to one UniProt — resolve
    that to an HGNC symbol (so the potency table says WHAT each row is against,
    instead of being target-agnostic). Cached; falls back to the accession."""
    if not act_id:
        return None
    try:
        u = map_all(act_id, ">>chembl_activity>>uniprot")
        uni = u[0].get("id") if u else None
        if not uni:
            return None
        if uni in _TGT_CACHE:
            return _TGT_CACHE[uni]
        sym = uni
        h = map_all(uni, ">>uniprot>>hgnc")
        if h:
            he = entry(h[0]["id"], "hgnc")
            syms = ((he.get("Attributes") or {}).get("Hgnc") or {}).get("symbols") or []
            sym = syms[0] if syms else uni
        _TGT_CACHE[uni] = sym
        return sym
    except Exception:
        return None


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
            "target": _activity_target(r.get("id")),
            "type": r.get("standard_type"),
            "value": r.get("standard_value"),
            "unit": r.get("standard_units"),
            "pchembl": r.get("pchembl"),
        } for r in potent[:50]],
    }


SECTION = Section(
    id="3", name="bioactivity",
    description=("ChEMBL bioactivity for the drug (chembl_molecule→chembl_activity); "
                 "pchembl-ranked, ≥5 (real binding), top 50 by potency"),
    needs=("chembl_id",),
    produces=("activity_total", "potent_count", "activities"),
    datasets=("chembl_molecule", "chembl_activity"),
    chains=(">>chembl_molecule>>chembl_activity",),
    collect_fn=collect,
)
