"""§14 — pathway analysis: per-cohort-gene Reactome pathways + cohort-level
aggregation (which pathways are most over-represented?). REUSE wrapper over
gene §7."""
from atlas.section import Section

CHAINS   = (">>hgnc>>reactome", ">>hgnc>>ensembl>>reactome")
DATASETS = ("hgnc", "ensembl", "reactome")

def collect(a):
    bundle = {"section": "14_pathways", "mondo_id": a.mondo_id,
              "_todo": "Implement: fan gene §7 over cohort + reactome rollup"}
    return bundle

SECTION = Section(
    id="14", name="pathways",
    description=("Cohort-level Reactome pathway aggregation: top enriched "
                 "pathways, per-pathway gene counts, druggability per pathway."),
    needs=("cohort",),
    produces=("per_gene_pathways", "top_pathways", "pathway_gene_counts"),
    datasets=DATASETS, chains=CHAINS, collect_fn=collect,
)
