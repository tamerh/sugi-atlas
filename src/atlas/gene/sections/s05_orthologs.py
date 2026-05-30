"""§5 — orthologs: cross-species orthologs (Ensembl Compara) + paralogs."""
from atlas.biobtree import map_all
from atlas.gene.sections.base import Section

CHAINS = (">>ensembl>>ortholog", ">>ensembl>>paralog")
DATASETS = ("ensembl", "ortholog", "paralog")

def collect(a):
    bundle = {"section": "05_orthologs", "symbol": a.symbol, "ensembl_id": a.ensembl_id}
    orths = map_all(a.ensembl_id, ">>ensembl>>ortholog") if a.ensembl_id else []
    bundle["orthologs"] = [{"id": t["id"], "symbol": t.get("name"),
                            "organism": t.get("genome")} for t in orths]
    bundle["ortholog_count"] = len(orths)
    paras = map_all(a.ensembl_id, ">>ensembl>>paralog") if a.ensembl_id else []
    bundle["paralogs"] = [{"id": t["id"], "symbol": t.get("name")} for t in paras]
    bundle["paralog_count"] = len(paras)
    return bundle

SECTION = Section(
    id="5", name="orthologs",
    description="Orthologous genes in model organisms + paralogs (Ensembl Compara)",
    needs=("ensembl_id",),
    produces=("orthologs", "paralogs"),
    datasets=DATASETS, chains=CHAINS, collect_fn=collect,
)
