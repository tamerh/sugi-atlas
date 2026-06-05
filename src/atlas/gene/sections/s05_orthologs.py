"""§5 — orthologs: cross-species orthologs (Ensembl Compara) + paralogs."""
from atlas.biobtree import map_all
from atlas.gene.sections.base import Section

CHAINS = (">>ensembl>>ortholog", ">>ensembl>>paralog")
DATASETS = ("ensembl", "ortholog", "paralog")

# biobtree's >>ensembl>>ortholog edge leaves name/genome EMPTY for some model-
# organism-database namespaces (WormBase C. elegans confirmed; FlyBase/zebrafish/
# mouse/rat come back populated). Recover the organism from the id prefix so the
# row isn't half-blank; the symbol isn't resolvable via biobtree, so it stays
# empty (the id remains in its own column). See docs/internal/BIOBTREE_ISSUES.md.
_NS_ORGANISM = (("WBGene", "caenorhabditis_elegans"),
                ("FBgn", "drosophila_melanogaster"))


def _organism_from_id(gid):
    g = (gid or "").upper()
    return next((org for pre, org in _NS_ORGANISM if g.startswith(pre.upper())), "")


def collect(a):
    bundle = {"section": "05_orthologs", "symbol": a.symbol, "ensembl_id": a.ensembl_id}
    orths = map_all(a.ensembl_id, ">>ensembl>>ortholog") if a.ensembl_id else []
    bundle["orthologs"] = [{"id": t["id"], "symbol": t.get("name"),
                            "organism": t.get("genome") or _organism_from_id(t["id"])}
                           for t in orths]
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
