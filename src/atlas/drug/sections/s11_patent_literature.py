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
    # Per-matched-structure counts: raw patent mentions AND the distinct
    # patent_family count (BIOBTREE_ISSUES #26 — now a one-call xref). Families
    # are the honest dedup metric (one invention, many jurisdictions); the raw
    # mention total is inflated by promiscuous/reference structures, so surfacing
    # both — plus the per-structure split — lets the page quantify that.
    # (assignee/CPC attributes now exist too (#25), but a *representative*
    # landscape still needs date-sort / server-side facets (#27, still
    # id-ordered), so we don't sample them — see the issues log.)
    breakdown = []
    for pcid in pcids:
        try:
            xc = xref_counts(entry(pcid, "patent_compound"))
            n, fam = xc.get("patent", 0), xc.get("patent_family", 0)
        except Exception:
            n, fam = 0, 0
        breakdown.append({"id": pcid, "patent_count": n, "family_count": fam})
    breakdown.sort(key=lambda x: -x["patent_count"])
    total = sum(b["patent_count"] for b in breakdown)
    family_total = sum(b["family_count"] for b in breakdown)
    return {
        "section": "11_patent_literature",
        "patent_compound_ids": pcids,
        "patent_compound_breakdown": breakdown,
        "patent_total": total,
        "patent_family_total": family_total,
    }


SECTION = Section(
    id="11", name="patent_literature",
    description=("SureChEMBL patent coverage: sum of patent mentions across the "
                 "molecule's patent_compound records (IP-intensity signal)"),
    needs=("chembl_id",),
    produces=("patent_compound_ids", "patent_compound_breakdown", "patent_total",
              "patent_family_total"),
    datasets=("chembl_molecule", "patent_compound"),
    chains=(">>chembl_molecule>>patent_compound",),
    collect_fn=collect,
)
