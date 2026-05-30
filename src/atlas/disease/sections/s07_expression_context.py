"""§7 — expression context: per-cohort-gene tissue + single-cell expression
(Bgee, SCXA). REUSE wrapper over gene §11.

Lower priority: surfacing the disease-relevant tissue is heuristic; will refine
after the §11 reuse is in place."""
from atlas.section import Section

CHAINS   = (">>hgnc>>ensembl>>bgee", ">>hgnc>>ensembl>>scxa")  # via gene §11
DATASETS = ("hgnc", "ensembl", "bgee", "scxa")

def collect(a):
    bundle = {"section": "07_expression_context", "mondo_id": a.mondo_id,
              "_todo": "Implement: fan gene §11 over cohort + tissue summary"}
    return bundle

SECTION = Section(
    id="7", name="expression_context",
    description=("Per-cohort-gene Bgee tissue expression + SCXA single-cell. "
                 "Tissue-specificity ranking; flag genes with no expression "
                 "in disease-relevant tissue."),
    needs=("cohort",),
    produces=("per_gene_expression", "tissue_overlap"),
    datasets=DATASETS, chains=CHAINS, collect_fn=collect,
)
