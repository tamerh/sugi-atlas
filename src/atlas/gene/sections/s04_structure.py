"""§4 — structure: experimental PDB (method + resolution), AlphaFold predicted model (pLDDT)."""
from atlas.biobtree import map_all
from atlas.gene.sections.base import Section

CHAINS = (">>uniprot>>pdb", ">>uniprot>>alphafold", ">>uniprot>>pdb>>antibody")
DATASETS = ("pdb", "alphafold", "uniprot", "antibody")

def collect(a):
    bundle = {"section": "04_structure", "symbol": a.symbol,
              "hgnc_id": a.hgnc_id, "reviewed_uniprot": list(a.reviewed_uniprots)}

    pdb, af, ab_pdbs = {}, [], set()
    for u in a.reviewed_uniprots:
        for t in map_all(u, ">>uniprot>>pdb"):
            pdb[t["id"]] = {"id": t["id"], "method": t.get("method"),
                            "resolution": t.get("resolution"),
                            "title": (t.get("title") or "").strip()}
        # Antibody-complex structures (SAbDab) bound to this protein, reached via
        # the bidirectional PDB↔antibody edge (>>uniprot>>pdb>>antibody). The ids
        # are PDB_Hchain_Lchain — dedup to the distinct PDB. A high count flags a
        # validated antibody target (PD-1, EGFR, CD20). Map-only (no entry) so the
        # disease cohort-fan stays cheap; therapeutic INN names aren't edge-linked.
        for t in map_all(u, ">>uniprot>>pdb>>antibody"):
            pid = (t.get("id") or "").split("_")[0]
            if pid:
                ab_pdbs.add(pid)
        # biobtree 2026-05-31 refresh resolved BIOBTREE_ISSUES #10 — alphafold
        # coverage now extends to ~3000 aa (MTOR works at 2549 aa). Remaining
        # empties (ATM/BRCA2/DMD/TTN/MUC16 — all >~3000 aa) reflect AlphaFold-DB
        # upstream truly not having a model for these proteins, not a biobtree
        # gap. Trust biobtree: emit the entry only when it returns one, set a
        # `present` flag the renderer keys off to elide cleanly + add a footnote.
        m = map_all(u, ">>uniprot>>alphafold")
        if m:
            af.append({"id": f"AF-{u}-F1", "uniprot": u, "present": True,
                       "plddt": m[0].get("global_metric"),
                       "fraction_plddt_very_high": m[0].get("fraction_plddt_very_high")})
        else:
            af.append({"id": None, "uniprot": u, "present": False,
                       "plddt": None, "fraction_plddt_very_high": None,
                       "note": "no AlphaFold model — AlphaFold DB does not provide a model for proteins > ~3000 aa"})
    bundle["pdb"] = list(pdb.values())
    bundle["pdb_count"] = len(pdb)
    bundle["alphafold"] = af
    bundle["antibody_structures"] = sorted(ab_pdbs)
    return bundle

SECTION = Section(
    id="4", name="structure",
    description="Experimental PDB structures (method + resolution) and AlphaFold predicted model (pLDDT)",
    needs=("hgnc_id", "reviewed_uniprots"),
    produces=("pdb", "pdb_count", "alphafold", "antibody_structures"),
    datasets=DATASETS, chains=CHAINS, collect_fn=collect,
)
