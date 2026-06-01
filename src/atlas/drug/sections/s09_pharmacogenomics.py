"""§9 — pharmacogenomics. Drug-level PGx = CPIC / DPWG genotype-guided dosing
guidelines for THIS drug (drug × metabolizing-gene, e.g. atorvastatin × SLCO1B1
statin myopathy). NOT the drug's *target* gene — PGx is about the patient's
pharmacogene genotype.

Route: the drug has a `pharmgkb` chemical node (PA id), reached by name (there's
no chembl_molecule>>pharmgkb edge); `>>pharmgkb>>pharmgkb_guideline` then yields
the dosing guidelines. Clinical-annotation + variant tables are gene-keyed (not
reachable from the chemical node), so they live on the relevant gene page.
Empty for drugs without a curated guideline (e.g. newer targeted agents)."""
from atlas.biobtree import search, rows, map_all
from atlas.section import Section


def _pharmgkb_chemical_id(name):
    """Resolve the drug's PharmGKB chemical id (PA…) by name — prefer exact
    case-insensitive match, else highest-xref hit. {None} if not in PharmGKB."""
    res = rows(search(name, source="pharmgkb"))
    if not res:
        return None
    q = name.lower()
    exact = [r for r in res if (r.get("name") or "").lower() == q]
    pick = (exact or sorted(res, key=lambda r: int(r.get("xref_count") or 0),
                            reverse=True))[0]
    return pick.get("id")


def collect(a):
    pa = _pharmgkb_chemical_id(a.canonical_name)
    guidelines = []
    if pa:
        for r in map_all(pa, ">>pharmgkb>>pharmgkb_guideline"):
            guidelines.append({
                "id": r.get("id"),
                "name": r.get("name"),
                "source": r.get("source"),          # CPIC / DPWG / ...
                "genes": r.get("gene_symbols"),
                "chemicals": r.get("chemical_names"),
                "has_dosing": r.get("has_dosing_info") == "true",
                "has_recommendation": r.get("has_recommendation") == "true",
            })
    return {
        "section": "09_pharmacogenomics",
        "pharmgkb_chemical_id": pa,
        "guidelines": guidelines,
        "guideline_count": len(guidelines),
    }


SECTION = Section(
    id="9", name="pharmacogenomics",
    description=("Drug-level pharmacogenomics: CPIC / DPWG genotype-guided dosing "
                 "guidelines (drug × pharmacogene) via pharmgkb→pharmgkb_guideline"),
    needs=("canonical_name",),
    produces=("pharmgkb_chemical_id", "guidelines", "guideline_count"),
    datasets=("pharmgkb", "pharmgkb_guideline"),
    chains=(">>pharmgkb>>pharmgkb_guideline",),
    collect_fn=collect,
)
