"""§5 — genes → proteins: per-cohort-gene id mapping (HGNC, Ensembl, UniProt,
protein name). REUSE wrapper — fans gene §1 (gene_ids) + §3 (protein_ids) over
the cohort and aggregates.

No new chains/datasets — everything comes from gene §1 + §3 bundles via the
cohort fan-out helper."""
from atlas.section import Section
from atlas.disease.cohort import fan
from atlas.gene.sections import s01_gene_ids, s03_protein_ids

CHAINS   = (">>hgnc>>ensembl>>uniprot",)  # reused via gene collectors
DATASETS = ("hgnc", "ensembl", "uniprot")

_EVIDENCE_KEYS = ("gwas", "gencc", "clinvar", "civic_evidence")


def _classify(ev: dict) -> str:
    """Partition a gene's evidence flags into a single bucket name."""
    g = bool(ev.get("gwas"))
    gc = bool(ev.get("gencc"))
    cv = bool(ev.get("clinvar"))
    ci = bool(ev.get("civic_evidence"))
    n = sum((g, gc, cv, ci))
    if n >= 3:
        return "multi_evidence"
    if g and gc and n == 2:
        return "gwas_and_gencc"
    if g and cv and n == 2:
        return "gwas_and_clinvar"
    if g and not (gc or cv or ci):
        return "gwas_only"
    if ci and not (g or gc or cv):
        return "civic_only"
    # Two-evidence combos not listed above (e.g. gencc+clinvar, gwas+civic)
    # collapse into multi_evidence; single-source non-GWAS/non-CIViC into
    # multi_evidence too so the partition stays exhaustive.
    if n >= 2:
        return "multi_evidence"
    return "multi_evidence"


def collect(a):
    g1_bundles = fan(s01_gene_ids.SECTION.collect_fn, a.cohort)
    g3_bundles = fan(s03_protein_ids.SECTION.collect_fn, a.cohort)
    g3_by = {b.get("symbol"): b for b in g3_bundles}

    summary = {"gwas_only": 0, "gwas_and_gencc": 0, "gwas_and_clinvar": 0,
               "civic_only": 0, "multi_evidence": 0}
    genes = []
    seen_canonical = set()
    for b1 in g1_bundles:
        sym = b1.get("symbol")
        b3 = g3_by.get(sym, {})
        ev_raw = a.cohort_evidence.get(b1.get("hgnc_id")) or {}
        ev = {k: bool(ev_raw.get(k)) for k in _EVIDENCE_KEYS}
        canonical = b3.get("canonical_uniprot")
        if canonical:
            seen_canonical.add(canonical)
        genes.append({
            "symbol": sym,
            "hgnc_id": b1.get("hgnc_id"),
            "ensembl_id": b1.get("ensembl_id"),
            "hgnc_name": b1.get("name"),
            "canonical_uniprot": canonical,
            "all_uniprots": b3.get("uniprot_all", []),
            "protein_name": None,  # not surfaced by §3; would need extra uniprot fetch
            "evidence": ev,
        })
        summary[_classify(ev)] += 1

    return {
        "section": "05_genes_proteins",
        "mondo_id": a.mondo_id,
        "genes": genes,
        "gene_count": len(genes),
        "protein_count": len(seen_canonical),
        "evidence_summary": summary,
    }


SECTION = Section(
    id="5", name="genes_proteins",
    description=("Per-cohort-gene identifier mapping (HGNC → Ensembl → UniProt + "
                 "protein name + evidence tier + Mendelian overlap flag). "
                 "REUSE wrapper over gene §1/§3."),
    needs=("cohort", "cohort_evidence"),
    produces=("genes", "protein_count", "gene_count"),
    datasets=DATASETS, chains=CHAINS, collect_fn=collect,
)
