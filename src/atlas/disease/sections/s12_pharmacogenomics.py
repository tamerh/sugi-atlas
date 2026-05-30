"""§12 — pharmacogenomics: per-cohort-gene PharmGKB coverage (gene-drug
interactions, CPIC guidelines, VIP flags).

Calls `map_all(hgnc_id, >>hgnc>>pharmgkb_gene)` directly per cohort gene
rather than fanning the full gene §10 (which is heavy: ChEMBL targets,
molecules, patents, BindingDB, ...). Cross-section sharing of the §10
pharmgkb slice isn't supported yet, and a direct map_all per gene is
~1 call/gene — cheap enough to repeat."""
from atlas.biobtree import map_all
from atlas.section import Section

CHAINS   = (">>hgnc>>pharmgkb_gene",)
DATASETS = ("hgnc", "pharmgkb_gene")

def collect(a):
    bundle = {"section": "12_pharmacogenomics", "mondo_id": a.mondo_id}

    per_gene_pgx = []
    for ga in a.cohort:
        try:
            rows = map_all(ga.hgnc_id, ">>hgnc>>pharmgkb_gene")
        except Exception as e:
            per_gene_pgx.append({"symbol": ga.symbol, "hgnc_id": ga.hgnc_id,
                                 "pharmgkb_entries": [], "_error": str(e)})
            continue
        # biobtree returns these flags as the literal strings 'true'/'false'
        # (not booleans) — normalize defensively.
        def _truthy(v):
            return v is True or (isinstance(v, str) and v.lower() == "true")
        entries = [{"id": r["id"], "vip": _truthy(r.get("is_vip")),
                    "cpic_guideline": _truthy(r.get("has_cpic_guideline"))}
                   for r in rows if r.get("id")]
        per_gene_pgx.append({"symbol": ga.symbol, "hgnc_id": ga.hgnc_id,
                             "pharmgkb_entries": entries})

    pgx_genes = []
    vip_count = 0
    cpic_count = 0
    for g in per_gene_pgx:
        entries = g["pharmgkb_entries"]
        if not entries:
            continue
        v = sum(1 for e in entries if e["vip"])
        c = sum(1 for e in entries if e["cpic_guideline"])
        pgx_genes.append({"symbol": g["symbol"], "vip_count": v, "cpic_count": c})
        if v: vip_count += 1
        if c: cpic_count += 1

    bundle["per_gene_pgx"]   = per_gene_pgx
    bundle["pgx_genes"]      = pgx_genes
    bundle["pgx_gene_count"] = len(pgx_genes)
    bundle["vip_count"]      = vip_count
    bundle["cpic_count"]     = cpic_count
    return bundle

SECTION = Section(
    id="12", name="pharmacogenomics",
    description=("Per-cohort-gene PharmGKB coverage: gene-drug interactions, "
                 "CPIC guidelines, VIP flags. Gene-level PGx signal."),
    needs=("cohort",),
    produces=("per_gene_pgx", "pgx_genes", "cpic_count", "vip_count"),
    datasets=DATASETS, chains=CHAINS, collect_fn=collect,
)
