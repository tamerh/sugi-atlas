"""§11 — bioactivity & enzyme data: per-cohort-gene ChEMBL assay depth +
BRENDA EC classification. REUSE wrapper over gene §3 (brenda_ec) + gene §10
(chembl_assay_total + type counts).

Bioactivity counts are a 'studied-ness' signal that helps rank undrugged
targets for §17 (high assay count = easier starting point). Drugged-vs-
undrugged resolution is deferred to render-time §17, which has access to the
gene §10 molecules list."""
from atlas.section import Section
from atlas.disease.cohort import fan
from atlas.gene.sections import s03_protein_ids, s10_drugs

CHAINS   = (">>uniprot>>chembl_target>>chembl_assay", ">>uniprot>>brenda",
            ">>uniprot>>pubchem_activity")
DATASETS = ("uniprot", "chembl_target", "chembl_assay", "brenda",
            "pubchem_activity")

# Threshold for "high screening signal" — anything ≥100 assays has enough
# chemical matter to seed a hit-to-lead campaign even if no drug is approved
# yet. Tuned from migrated EC page: ESR1/TP53/PIK3CA/PPARG/KRAS all clear it;
# the undrugged-but-screened bucket (BPTF, MAP2K5) lands here too.
HIGH_ASSAY_THRESHOLD = 100


def collect(a):
    # Fan both gene §3 (for brenda_ec) and gene §10 (for chembl_assay_*).
    # ~600 calls + ~30s on a 50-gene cohort; both fans are needed because
    # §3 and §10 produce disjoint fields and we don't want a third gene
    # collector just for the union.
    s03 = {b.get("hgnc_id"): b for b in fan(s03_protein_ids.SECTION.collect_fn, a.cohort)}
    s10 = {b.get("hgnc_id"): b for b in fan(s10_drugs.SECTION.collect_fn, a.cohort)}

    per_gene_bioactivity = []
    enzyme_genes = []
    undrugged_starting_points = []

    for ga in a.cohort:
        hg = ga.hgnc_id
        sym = ga.symbol
        b03 = s03.get(hg) or {}
        b10 = s10.get(hg) or {}

        assay_total = int(b10.get("chembl_assay_total") or 0)
        assay_types = dict(b10.get("chembl_assay_type_counts") or {})

        brenda = b03.get("brenda_ec") or []
        # Dedupe EC numbers across multiple uniprot accessions; keep the first
        # name we see for each EC.
        ec_map = {}
        for r in brenda:
            ec = r.get("ec")
            if not ec or ec in ec_map:
                continue
            ec_map[ec] = r.get("name")
        ec_numbers = list(ec_map.keys())
        ec_names = list(ec_map.values())

        per_gene_bioactivity.append({
            "symbol": sym,
            "hgnc_id": hg,
            "chembl_assay_total": assay_total,
            "chembl_assay_types": assay_types,
            "ec_numbers": ec_numbers,
            "ec_names": ec_names,
        })

        if ec_numbers:
            enzyme_genes.append({
                "symbol": sym,
                "ec_numbers": ec_numbers,
                "ec_names": ec_names,
            })

        if assay_total >= HIGH_ASSAY_THRESHOLD:
            # Drugged-vs-undrugged status resolved at render-time §17 via the
            # gene §10 molecules list (would require a third fan here).
            undrugged_starting_points.append({
                "symbol": sym,
                "chembl_assay_total": assay_total,
                "note": "high-screening signal, drugged status TBD",
            })

    top_studied_genes = sorted(
        ({"symbol": g["symbol"],
          "total": g["chembl_assay_total"],
          "type_summary": g["chembl_assay_types"]}
         for g in per_gene_bioactivity
         if g["chembl_assay_total"] > 0),
        key=lambda g: g["total"], reverse=True,
    )[:10]

    return {
        "section": "11_bioactivity_enzyme",
        "mondo_id": a.mondo_id,
        "per_gene_bioactivity": per_gene_bioactivity,
        "top_studied_genes": top_studied_genes,
        "enzyme_genes": enzyme_genes,
        "enzyme_count": len(enzyme_genes),
        "undrugged_starting_points": undrugged_starting_points,
    }


SECTION = Section(
    id="11", name="bioactivity_enzyme",
    description=("Per-cohort-gene ChEMBL assay depth + BRENDA enzyme "
                 "classification. Bioactivity signal for undrugged-target "
                 "prioritisation."),
    needs=("cohort",),
    produces=("per_gene_bioactivity", "top_studied_genes", "enzyme_genes",
              "enzyme_count", "undrugged_starting_points"),
    datasets=DATASETS, chains=CHAINS, collect_fn=collect,
)
