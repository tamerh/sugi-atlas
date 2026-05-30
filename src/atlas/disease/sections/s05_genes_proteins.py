"""§5 — genes → proteins: per-cohort-gene id mapping (HGNC, Ensembl, UniProt,
protein name). REUSE wrapper — fans gene §1 (gene_ids) + §3 (protein_ids) over
the cohort and aggregates.

Skeleton — implementation will use cohort.fan() against gene §1/§3 collectors."""
from atlas.section import Section

CHAINS   = (">>hgnc>>ensembl>>uniprot",)  # reused via gene collectors
DATASETS = ("hgnc", "ensembl", "uniprot")

def collect(a):
    bundle = {"section": "05_genes_proteins", "mondo_id": a.mondo_id,
              "_todo": "Implement: fan gene §1+§3 over a.cohort"}
    return bundle

SECTION = Section(
    id="5", name="genes_proteins",
    description=("Per-cohort-gene identifier mapping (HGNC → Ensembl → UniProt + "
                 "protein name + evidence tier + Mendelian overlap flag). "
                 "REUSE wrapper over gene §1/§3."),
    needs=("cohort", "cohort_evidence"),
    produces=("genes", "protein_count", "gene_count"),
    datasets=DATASETS, chains=CHAINS, collect_fn=collect,
)
