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
    """Resolve a ChEMBL activity → (gene_symbol, protein_name, uniprot_acc,
    organism), any of which may be None. The activity links to one UniProt —
    that's the assayed PROTEIN (accession + name); the gene symbol (when it maps)
    is the encoding gene, used to link the Atlas gene page. `organism` is the
    species common name and is set ONLY for NON-human targets (no HGNC symbol) —
    many assays run on animal orthologs (sheep/cattle COX-1), which correctly have
    no human gene; labelling the species explains the blank gene cell. Cached."""
    if not act_id:
        return (None, None, None, None)
    try:
        u = map_all(act_id, ">>chembl_activity>>uniprot")
        uni = u[0].get("id") if u else None
        if not uni:
            return (None, None, None, None)
        if uni in _TGT_CACHE:
            return _TGT_CACHE[uni]
        # the assayed protein's name, from the UniProt entry
        ue = entry(uni, "uniprot")
        pname = ((ue.get("Attributes") or {}).get("Uniprot") or {}).get("name") or None
        # the encoding gene's HGNC symbol (for the Atlas gene link), when it maps
        symbol = None
        h = map_all(uni, ">>uniprot>>hgnc")
        if h:
            he = entry(h[0]["id"], "hgnc")
            syms = ((he.get("Attributes") or {}).get("Hgnc") or {}).get("symbols") or []
            if syms:
                symbol = syms[0]
        organism = None
        if symbol is None:                       # non-human ortholog → label the species
            tx = map_all(uni, ">>uniprot>>taxonomy")
            if tx:
                organism = tx[0].get("common_name") or tx[0].get("name") or None
        result = (symbol, pname, uni, organism)
        _TGT_CACHE[uni] = result
        return result
    except Exception:
        return (None, None, None, None)


def collect(a):
    rows = map_all(a.chembl_id, ">>chembl_molecule>>chembl_activity")
    potent = [r for r in rows if _pchembl(r) >= 5.0]
    potent.sort(key=_pchembl, reverse=True)
    # Fill to 100 DISTINCT rows by potency. ChEMBL carries many activities for the
    # same drug↔protein↔value (different assays) that rendered as identical repeats;
    # dedup by (uniprot, type, value, unit) AND resolve past the top 100 to backfill
    # what dedup drops — up to 100 distinct or a 200-row safety cap (heavy-dup drugs
    # simply show fewer).
    acts, seen = [], set()
    for r in potent[:200]:
        if len(acts) >= 100:
            break
        symbol, pname, uni, organism = _activity_target(r.get("id"))
        key = (uni, r.get("standard_type"), r.get("standard_value"), r.get("standard_units"))
        if key in seen:
            continue
        seen.add(key)
        acts.append({
            "id": r.get("id"),
            "target_symbol": symbol,
            "protein_name": pname,
            "uniprot": uni,
            "organism": organism,
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
    datasets=("chembl_molecule", "chembl_activity", "uniprot", "hgnc", "taxonomy"),
    chains=(">>chembl_molecule>>chembl_activity",),
    collect_fn=collect,
)
