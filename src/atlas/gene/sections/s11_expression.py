"""§11 — expression: Bgee gene-level + bgee_evidence per-tissue scores,
FANTOM5 CAGE (gene + alternative promoters), single-cell datasets (scxa)."""
from atlas.biobtree import map_all
from atlas.gene.sections.base import Section

CHAINS = (
    ">>ensembl>>bgee", ">>ensembl>>bgee>>bgee_evidence",
    ">>ensembl>>fantom5_gene", ">>ensembl>>fantom5_promoter",
    ">>ensembl>>scxa",
)
DATASETS = ("bgee", "bgee_evidence", "fantom5_gene", "fantom5_promoter", "scxa", "ensembl")

def _num(x):
    try: return float(x)
    except (TypeError, ValueError): return 0.0

def collect(a):
    bundle = {"section": "11_expression", "symbol": a.symbol}
    ens = a.ensembl_id

    bgee = map_all(ens, ">>ensembl>>bgee") if ens else []
    if bgee:
        b = bgee[0]
        bundle["bgee"] = {"breadth": b.get("expression_breadth"),
                          "present_calls": b.get("total_present_calls"),
                          "max_expression_score": b.get("max_expression_score")}

    # Per-tissue expression (Bgee evidence) — top by expression score.
    tissues = map_all(ens, ">>ensembl>>bgee>>bgee_evidence") if ens else []
    tissues.sort(key=lambda t: _num(t.get("expression_score")), reverse=True)
    bundle["top_tissues"] = [{"tissue": t.get("anatomical_entity_name"),
                              "score": t.get("expression_score"),
                              "rank": t.get("expression_rank"),
                              "quality": t.get("call_quality")} for t in tissues[:30]]
    bundle["tissue_count"] = len(tissues)

    # FANTOM5 CAGE gene-level expression (TPM + breadth).
    f5 = map_all(ens, ">>ensembl>>fantom5_gene") if ens else []
    if f5:
        x = f5[0]
        bundle["fantom5"] = {"tpm_average": x.get("tpm_average"), "tpm_max": x.get("tpm_max"),
                             "samples_expressed": x.get("samples_expressed"),
                             "breadth": x.get("expression_breadth")}

    # FANTOM5 alternative promoters (TSS usage) — sorted by activity.
    proms = map_all(ens, ">>ensembl>>fantom5_promoter") if ens else []
    proms.sort(key=lambda t: _num(t.get("tpm_average")), reverse=True)
    bundle["fantom5_promoters"] = [{"id": t["id"], "tpm_average": t.get("tpm_average"),
                                    "samples_expressed": t.get("samples_expressed")} for t in proms]

    bundle["single_cell_datasets"] = [{"id": t["id"], "description": t.get("description"),
                                       "cells": t.get("number_of_cells")}
                                      for t in (map_all(ens, ">>ensembl>>scxa") if ens else [])]
    return bundle

SECTION = Section(
    id="11", name="expression",
    description="Bgee per-tissue expression + FANTOM5 CAGE (gene + alt promoters) + single-cell datasets",
    needs=("ensembl_id",),
    produces=("bgee", "top_tissues", "fantom5", "fantom5_promoters", "single_cell_datasets"),
    datasets=DATASETS, chains=CHAINS, collect_fn=collect,
)
