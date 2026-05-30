"""§4 — Mendelian disease overlap: which cohort genes also cause Mendelian
forms (via OMIM / Orphanet / GenCC). For cancers: which genes also carry
somatic-driver evidence (CIViC / intogen at gene level).

Genes with BOTH GWAS + Mendelian evidence = highest-confidence targets.

NEW collector. Skeleton."""
from atlas.section import Section

CHAINS   = (">>mondo>>orphanet", ">>orphanet>>hgnc", ">>orphanet>>mim",
            ">>mondo>>gencc>>hgnc", ">>hgnc>>intogen", ">>hgnc>>civic")
DATASETS = ("mondo", "orphanet", "mim", "gencc", "hgnc", "intogen", "civic")

def collect(a):
    bundle = {"section": "04_mendelian_overlap", "mondo_id": a.mondo_id,
              "_todo": "Implement: gencc/orphanet/omim genes + intersection"}
    return bundle

SECTION = Section(
    id="4", name="mendelian_overlap",
    description=("Cohort genes with Mendelian evidence (GenCC, Orphanet, OMIM); "
                 "for cancers: cohort genes with somatic-driver evidence "
                 "(intOGen, CIViC). High-confidence target subset."),
    needs=("mondo_id", "cohort", "cohort_evidence", "is_cancer", "orphanet_ids"),
    produces=("gencc_genes", "orphanet_genes", "omim_genes",
              "somatic_driver_genes", "dual_evidence_genes"),
    datasets=DATASETS, chains=CHAINS, collect_fn=collect,
)
