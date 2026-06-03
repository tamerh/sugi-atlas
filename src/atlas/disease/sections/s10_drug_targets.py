"""§10 — drug targets: per-cohort-gene ChEMBL targets + phased molecules.
REUSE wrapper over gene §10 (chembl_targets + molecules subset) + aggregate
into approved / phase≥3 / phase≥1 buckets and drugged-vs-undrugged split.

Key derived stat for §16: 'genes with ≥1 approved drug' vs 'undrugged'."""
from collections import Counter
from atlas.section import Section
from atlas.disease.cohort import fan, enrichment_fan
from atlas.gene.sections import s10_drugs

CHAINS   = (">>uniprot>>chembl_target",
            ">>chembl_target>>chembl_molecule[highestDevelopmentPhase>=1]")
DATASETS = ("uniprot", "chembl_target", "chembl_molecule")


def _phase(d):
    try:
        return int(d) if (str(d) or "").isdigit() else 0
    except (TypeError, ValueError):
        return 0


def collect(a):
    bundles = fan(s10_drugs.SECTION.collect_fn, a.cohort)

    per_gene_drugs = []
    approved_count = 0
    phase3_count = 0
    phased_count = 0
    undrugged_count = 0
    approved_genes = []
    drugs_agg = {}  # id -> {id, name, max_phase, gene_targets:set}

    for b in bundles:
        sym = b.get("symbol")
        molecules = b.get("molecules") or []
        sorted_mols = sorted(molecules, key=lambda m: _phase(m.get("phase")),
                             reverse=True)
        max_phase = _phase(sorted_mols[0].get("phase")) if sorted_mols else 0
        is_target = bool(b.get("is_drug_target"))
        mol_count = len(sorted_mols)

        top_mols = [{"id": m.get("id"), "name": m.get("name"),
                     "phase": _phase(m.get("phase"))}
                    for m in sorted_mols[:5]]

        per_gene_drugs.append({
            "symbol": sym,
            "hgnc_id": b.get("hgnc_id"),
            "is_drug_target": is_target,
            "chembl_target_count": len(b.get("chembl_targets") or []),
            "molecule_count": mol_count,
            "max_phase": max_phase,
            "top_molecules": top_mols,
        })

        if max_phase >= 4:
            approved_count += 1
            top_name = top_mols[0].get("name") if top_mols else None
            approved_genes.append({"symbol": sym, "drug": top_name})
        if max_phase >= 3:
            phase3_count += 1
        if max_phase >= 1:
            phased_count += 1
        if (not is_target) or mol_count == 0:
            undrugged_count += 1

        for m in sorted_mols:
            mid = m.get("id")
            if not mid:
                continue
            p = _phase(m.get("phase"))
            rec = drugs_agg.get(mid)
            if rec is None:
                drugs_agg[mid] = {"id": mid, "name": m.get("name"),
                                  "max_phase": p, "gene_targets": {sym} if sym else set()}
            else:
                if p > rec["max_phase"]:
                    rec["max_phase"] = p
                if sym:
                    rec["gene_targets"].add(sym)

    drugs = sorted(
        ({"id": r["id"], "name": r["name"], "max_phase": r["max_phase"],
          "gene_targets": sorted(r["gene_targets"])}
         for r in drugs_agg.values()),
        key=lambda d: d["max_phase"], reverse=True,
    )[:30]

    top_targets = sorted(
        ({"symbol": g["symbol"], "molecule_count": g["molecule_count"],
          "max_phase": g["max_phase"]} for g in per_gene_drugs),
        key=lambda g: g["molecule_count"], reverse=True,
    )[:10]

    # Druggability BREADTH over the wide enrichment cohort: of all the
    # evidence-associated genes (not just the displayed 75), how many have a
    # ChEMBL target at all? A cheap presence fan (one chain/gene over ~250),
    # complementing the deep per-gene molecule data above on the display cohort.
    enrichment = enrichment_fan(a.enrichment_cohort,
                                ">>hgnc>>ensembl>>uniprot>>chembl_target")
    enrichment_druggable = sum(1 for _h, _s, rows in enrichment if rows)

    return {
        "section": "10_drug_targets",
        "mondo_id": a.mondo_id,
        "per_gene_drugs": per_gene_drugs,
        "enrichment_size": len(a.enrichment_cohort),
        "enrichment_druggable": enrichment_druggable,
        "approved_count": approved_count,
        "phase3_count": phase3_count,
        "phased_count": phased_count,
        "undrugged_count": undrugged_count,
        "approved_genes": approved_genes,
        "top_targets": top_targets,
        "drugs": drugs,
    }


SECTION = Section(
    id="10", name="drug_targets",
    description=("Per-cohort-gene ChEMBL targets + phased molecules. "
                 "Approved / Phase ≥3 / Phase ≥1 buckets per gene. "
                 "Drugged-vs-undrugged split."),
    needs=("cohort", "enrichment_cohort"),
    produces=("per_gene_drugs", "approved_count", "phase3_count",
              "phased_count", "undrugged_count", "approved_genes",
              "top_targets", "drugs", "enrichment_size", "enrichment_druggable"),
    datasets=DATASETS, chains=CHAINS, collect_fn=collect,
)
