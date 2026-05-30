"""§4 — structure: experimental PDB (method + resolution), AlphaFold predicted model (pLDDT)."""
from atlas.biobtree import map_all
from atlas.gene.sections.base import Section

CHAINS = (">>uniprot>>pdb", ">>uniprot>>alphafold")
DATASETS = ("pdb", "alphafold", "uniprot")

def collect(a):
    bundle = {"section": "04_structure", "symbol": a.symbol,
              "hgnc_id": a.hgnc_id, "reviewed_uniprot": list(a.reviewed_uniprots)}

    pdb, af = {}, []
    for u in a.reviewed_uniprots:
        for t in map_all(u, ">>uniprot>>pdb"):
            pdb[t["id"]] = {"id": t["id"], "method": t.get("method"),
                            "resolution": t.get("resolution")}
        # Every UniProt protein has an AlphaFold model (id AF-<acc>-F1). biobtree's
        # alphafold map is EMPTY for very large/fragmented proteins (>~2700 aa) —
        # see BIOBTREE_ISSUES.md #10. Construct AF-<acc>-F1 and attach pLDDT when present.
        m = map_all(u, ">>uniprot>>alphafold")
        af.append({"id": f"AF-{u}-F1", "uniprot": u,
                   "plddt": m[0].get("global_metric") if m else None,
                   "fraction_plddt_very_high": m[0].get("fraction_plddt_very_high") if m else None})
    bundle["pdb"] = list(pdb.values())
    bundle["pdb_count"] = len(pdb)
    bundle["alphafold"] = af
    return bundle

SECTION = Section(
    id="4", name="structure",
    description="Experimental PDB structures (method + resolution) and AlphaFold predicted model (pLDDT)",
    needs=("hgnc_id", "reviewed_uniprots"),
    produces=("pdb", "pdb_count", "alphafold"),
    datasets=DATASETS, chains=CHAINS, collect_fn=collect,
)
