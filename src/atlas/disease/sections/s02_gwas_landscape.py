"""§2 — GWAS landscape: mondo → GWAS assocs + studies. Counts, top hits,
study-level metadata.

NEW collector. Skeleton — collect() to be filled by an agent in the next pass."""
from atlas.section import Section
from atlas.biobtree import map_all, entry

CHAINS   = (">>mondo>>gwas", ">>mondo>>gwas_study", ">>mondo>>gwas>>hgnc")
DATASETS = ("mondo", "gwas", "gwas_study", "hgnc")

def collect(a):
    bundle = {"section": "02_gwas_landscape", "mondo_id": a.mondo_id,
              "_todo": "Implement: top assocs with p-values + study attrs"}
    return bundle

SECTION = Section(
    id="2", name="gwas_landscape",
    description=("GWAS associations + studies for the disease — total counts, "
                 "top assocs by p-value, study-level meta (lead author, year, "
                 "case/control counts)."),
    needs=("mondo_id", "xref_counts"),
    produces=("assoc_total", "study_total", "top_assocs", "studies"),
    datasets=DATASETS, chains=CHAINS, collect_fn=collect,
)
