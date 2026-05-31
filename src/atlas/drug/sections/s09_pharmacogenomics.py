"""§9 — pharmacogenomics. The direct drug→PGx route
(>>chembl_molecule>>pharmgkb_drug) is empty in biobtree (same root cause as
BIOBTREE_ISSUES #13 — pharmgkb clinical/guideline/variant edges declared but
unpopulated). Fallback: surface PharmGKB *gene-level* coverage for the drug's
target genes (VIP / CPIC-guideline flags via >>hgnc>>pharmgkb_gene). When #13
lands, swap to the drug-level guideline/variant content."""
from atlas.biobtree import map_all
from atlas.section import Section


def collect(a):
    rows, seen = [], set()
    for t in a.targets:
        if not t.hgnc_id or t.hgnc_id in seen:
            continue
        seen.add(t.hgnc_id)
        for r in map_all(t.hgnc_id, ">>hgnc>>pharmgkb_gene"):
            rows.append({"gene": t.gene_symbol, "pharmgkb_id": r.get("id"),
                         "vip": r.get("is_vip"),
                         "cpic_guideline": r.get("has_cpic_guideline")})
    return {
        "section": "09_pharmacogenomics",
        "source": "target_gene_fallback",
        "pgx_entries": rows,
        "pgx_count": len(rows),
    }


SECTION = Section(
    id="9", name="pharmacogenomics",
    description=("PharmGKB pharmacogenomics. Direct drug→PGx is blocked (biobtree "
                 "#13); fallback surfaces target-gene PharmGKB coverage "
                 "(VIP / CPIC flags)"),
    needs=("targets",),
    produces=("pgx_entries", "pgx_count", "source"),
    datasets=("hgnc", "pharmgkb_gene"),
    chains=(">>hgnc>>pharmgkb_gene",),
    collect_fn=collect,
)
