"""§3 — bioactivity: the drug's own ChEMBL activity rows
(>>chembl_molecule>>chembl_activity). pchembl is the standardized potency
metric (-log10(M): 10 = 0.1 nM, 6 = 1 µM). Keep rows at pchembl ≥ 5 (real
binding), sorted by potency, top 100. Empty for biologics (no small-molecule
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
    """Resolve a ChEMBL activity → (gene_symbol_or_None, display_label). The
    activity row omits the target but links to one UniProt; resolve that to an
    HGNC symbol when possible (so the symbol can link to its gene page), else
    fall back to the UniProt protein NAME — never a bare accession, so every
    potency row reads as a real target. Cached."""
    if not act_id:
        return (None, None)
    try:
        u = map_all(act_id, ">>chembl_activity>>uniprot")
        uni = u[0].get("id") if u else None
        if not uni:
            return (None, None)
        if uni in _TGT_CACHE:
            return _TGT_CACHE[uni]
        symbol, label = None, uni
        h = map_all(uni, ">>uniprot>>hgnc")
        if h:
            he = entry(h[0]["id"], "hgnc")
            syms = ((he.get("Attributes") or {}).get("Hgnc") or {}).get("symbols") or []
            if syms:
                symbol = label = syms[0]
        if symbol is None:                       # no gene symbol → protein name, not the bare accession
            ue = entry(uni, "uniprot")
            nm = ((ue.get("Attributes") or {}).get("Uniprot") or {}).get("name")
            if nm:
                label = nm
        _TGT_CACHE[uni] = (symbol, label)
        return (symbol, label)
    except Exception:
        return (None, None)


def collect(a):
    rows = map_all(a.chembl_id, ">>chembl_molecule>>chembl_activity")
    potent = [r for r in rows if _pchembl(r) >= 5.0]
    potent.sort(key=_pchembl, reverse=True)
    # Fill to 100 DISTINCT rows by potency. ChEMBL carries many activities for the
    # same drug↔target↔value (different assays) that rendered as identical repeats;
    # we dedup by (target, type, value, unit) AND resolve past the top 100 to
    # backfill what dedup drops — up to 100 distinct or a 200-row safety cap
    # (heavy-dup drugs simply show fewer).
    acts, seen = [], set()
    for r in potent[:200]:
        if len(acts) >= 100:
            break
        symbol, label = _activity_target(r.get("id"))
        key = (label, r.get("standard_type"), r.get("standard_value"), r.get("standard_units"))
        if key in seen:
            continue
        seen.add(key)
        acts.append({
            "id": r.get("id"),
            "target": label,
            "target_symbol": symbol,
            "type": r.get("standard_type"),
            "value": r.get("standard_value"),
            "unit": r.get("standard_units"),
            "pchembl": r.get("pchembl"),
        })
    return {
        "section": "03_bioactivity",
        "activity_total": len(rows),
        "potent_count": len(potent),
        "activities": acts,
    }


SECTION = Section(
    id="3", name="bioactivity",
    description=("ChEMBL bioactivity for the drug (chembl_molecule→chembl_activity); "
                 "pchembl-ranked, ≥5 (real binding), top 100 by potency"),
    needs=("chembl_id",),
    produces=("activity_total", "potent_count", "activities"),
    datasets=("chembl_molecule", "chembl_activity"),
    chains=(">>chembl_molecule>>chembl_activity",),
    collect_fn=collect,
)
