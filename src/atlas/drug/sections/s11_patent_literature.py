"""§11 — patent literature (SureChEMBL). Each chembl_molecule maps to 0-N
patent_compound records; each carries an xref count of patents mentioning the
compound. Sum across the molecule's patent_compounds for an IP-intensity
signal. Counts attach to the compound (not the gene–compound relationship), so
promiscuous compounds score high — surfaced as a drug-level total. Empty for
biologics (no SureChEMBL small-molecule extraction)."""
from atlas.biobtree import map_all, entry, xref_counts
from atlas.section import Section


def collect(a):
    pcids = [t.get("id") for t in map_all(a.chembl_id, ">>chembl_molecule>>patent_compound", cap=2)
             if t.get("id")]
    # Per-matched-structure patent counts. The headline total is often
    # dominated by one promiscuous/reference structure; surfacing the split
    # lets the page quantify that (the only honest enrichment available —
    # assignee/CPC/family-rollup/date-facets aren't exposed; see
    # BIOBTREE_ISSUES #25/#26/#27).
    breakdown = []
    for pcid in pcids:
        try:
            n = xref_counts(entry(pcid, "patent_compound")).get("patent", 0)
        except Exception:
            n = 0
        breakdown.append({"id": pcid, "patent_count": n})
    breakdown.sort(key=lambda x: -x["patent_count"])
    total = sum(b["patent_count"] for b in breakdown)
    return {
        "section": "11_patent_literature",
        "patent_compound_ids": pcids,
        "patent_compound_breakdown": breakdown,
        "patent_total": total,
    }


SECTION = Section(
    id="11", name="patent_literature",
    description=("SureChEMBL patent coverage: sum of patent mentions across the "
                 "molecule's patent_compound records (IP-intensity signal)"),
    needs=("chembl_id",),
    produces=("patent_compound_ids", "patent_compound_breakdown", "patent_total"),
    datasets=("chembl_molecule", "patent_compound"),
    chains=(">>chembl_molecule>>patent_compound",),
    collect_fn=collect,
)
