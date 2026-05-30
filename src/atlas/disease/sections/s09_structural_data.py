"""§9 — structural data: per-cohort-gene PDB + AlphaFold availability.
REUSE wrapper over gene §4 + aggregate.

Structure availability drives druggability (a kinase with PDB + co-crystal
is more druggable than one with AlphaFold only)."""
from atlas.section import Section
from atlas.disease.cohort import fan
from atlas.gene.sections import s04_structure

CHAINS   = (">>uniprot>>pdb", ">>uniprot>>alphafold")  # via gene §4
DATASETS = ("uniprot", "pdb", "alphafold")

# Limit on PDB / AlphaFold-only lists surfaced in the bundle. Mirrors the
# "top 30" cap the migrated AMD/EC pages used so downstream renderers don't
# explode for cohorts where every gene has a structure.
_LIST_CAP = 30


def collect(a):
    g4_bundles = fan(s04_structure.SECTION.collect_fn, a.cohort)

    per_gene = []
    pdb_genes = []
    af_only_genes = []
    pdb_count = af_only_count = no_struct_count = 0

    # Index cohort gene anchors by hgnc_id to recover canonical_uniprot
    # (gene §4 surfaces reviewed_uniprot list but not the canonical pick).
    by_hgnc = {g.hgnc_id: g for g in a.cohort}

    for b in g4_bundles:
        if b.get("_error"):
            continue
        symbol = b.get("symbol")
        hgnc_id = b.get("hgnc_id")
        ga = by_hgnc.get(hgnc_id)
        canonical_u = ga.canonical_uniprot if ga else None

        pdb_n = int(b.get("pdb_count") or 0)
        af_list = b.get("alphafold") or []
        # Prefer the AF model for canonical uniprot; else first non-null pLDDT.
        plddt = None
        af_present = False
        canon_af = next((af for af in af_list if af.get("uniprot") == canonical_u), None)
        pick = canon_af or (af_list[0] if af_list else None)
        if pick is not None:
            af_present = True
            plddt = pick.get("plddt")

        if pdb_n > 0:
            tier = "PDB"
            pdb_count += 1
            pdb_genes.append({"symbol": symbol, "uniprot": canonical_u,
                              "pdb_count": pdb_n})
        elif af_present:
            tier = "AlphaFold only"
            af_only_count += 1
            af_only_genes.append({"symbol": symbol, "uniprot": canonical_u,
                                  "plddt": plddt})
        else:
            tier = "No structure"
            no_struct_count += 1

        per_gene.append({
            "symbol": symbol,
            "hgnc_id": hgnc_id,
            "canonical_uniprot": canonical_u,
            "pdb_count": pdb_n,
            "has_alphafold": af_present,
            "alphafold_plddt": plddt,
            "structure_tier": tier,
        })

    pdb_genes.sort(key=lambda r: r["pdb_count"], reverse=True)
    def _plddt_key(r):
        v = r.get("plddt")
        try:
            return float(v) if v is not None else -1.0
        except (TypeError, ValueError):
            return -1.0
    af_only_genes.sort(key=_plddt_key, reverse=True)

    return {
        "section": "09_structural_data",
        "mondo_id": a.mondo_id,
        "per_gene_structure": per_gene,
        "pdb_count": pdb_count,
        "alphafold_only_count": af_only_count,
        "no_structure_count": no_struct_count,
        "pdb_genes": pdb_genes[:_LIST_CAP],
        "alphafold_only_genes": af_only_genes[:_LIST_CAP],
    }


SECTION = Section(
    id="9", name="structural_data",
    description=("Per-cohort-gene PDB + AlphaFold structure availability. "
                 "PDB-present / AlphaFold-only / no-structure split feeds §16."),
    needs=("cohort",),
    produces=("per_gene_structure", "pdb_count", "alphafold_only_count",
              "no_structure_count"),
    datasets=DATASETS, chains=CHAINS, collect_fn=collect,
)
