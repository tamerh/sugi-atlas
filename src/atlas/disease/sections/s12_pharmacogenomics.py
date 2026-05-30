"""§12 — pharmacogenomics: per-cohort-gene PharmGKB coverage (gene-drug
interactions, CPIC guidelines, VIP flags). REUSE wrapper over gene §10
(pharmgkb_gene field)."""
from atlas.section import Section

CHAINS   = (">>hgnc>>pharmgkb_gene",)
DATASETS = ("hgnc", "pharmgkb_gene", "pharmgkb")

def collect(a):
    bundle = {"section": "12_pharmacogenomics", "mondo_id": a.mondo_id,
              "_todo": "Implement: fan gene §10 (pharmgkb only) over cohort"}
    return bundle

SECTION = Section(
    id="12", name="pharmacogenomics",
    description=("Per-cohort-gene PharmGKB coverage: gene-drug interactions, "
                 "CPIC guidelines, VIP flags. Gene-level PGx signal."),
    needs=("cohort",),
    produces=("per_gene_pgx", "pgx_genes", "cpic_count", "vip_count"),
    datasets=DATASETS, chains=CHAINS, collect_fn=collect,
)
