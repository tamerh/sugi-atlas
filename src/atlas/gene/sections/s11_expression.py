"""§11 — expression: Bgee gene-level + bgee_evidence per-tissue scores,
FANTOM5 CAGE (gene + alternative promoters), single-cell datasets (scxa)."""
from atlas.biobtree import map_all, entry
from atlas.gene.sections.base import Section

CHAINS = (
    ">>ensembl>>bgee", ">>ensembl>>bgee>>bgee_evidence",
    ">>ensembl>>fantom5_gene", ">>ensembl>>fantom5_promoter",
    ">>ensembl>>scxa_expression",
    ">>ensembl>>scxa_expression>>scxa_gene_experiment",
)
DATASETS = ("bgee", "bgee_evidence", "fantom5_gene", "fantom5_promoter",
            "scxa_expression", "scxa_gene_experiment", "ensembl")

def _num(x):
    try: return float(x)
    except (TypeError, ValueError): return 0.0

def _int(x):
    try: return int(x)
    except (TypeError, ValueError): return 0

def _gtex_egene(ens):
    """GTEx v10 eGene calls (biobtree v2.16 `gtex_egene`): tissues where this
    gene's expression (eQTL) or splicing (sQTL) is significantly associated with
    a nearby variant (q < 0.05). The entry carries every tested tissue; keep the
    significant ones, strongest first."""
    if not ens:
        return {}
    try:
        e = (entry(ens, "gtex_egene").get("Attributes") or {}).get("GtexEgene") or {}
    except Exception:
        return {}
    rows = e.get("tissues") or []
    sig = []
    for t in rows:
        try:
            q = float(t.get("qval"))
        except (TypeError, ValueError):
            continue
        if t.get("is_egene") is True or str(t.get("is_egene")).lower() == "true" or q < 0.05:
            sig.append({"tissue": (t.get("tissue") or "").replace("_", " "),
                        "type": (t.get("qtl_type") or "").lower(), "q": q,
                        "slope": t.get("slope"), "lead_variant": t.get("lead_variant")})
    sig.sort(key=lambda t: (t["q"], t["tissue"]))
    return {"tested": len(rows), "egenes": sig}


def collect(a):
    bundle = {"section": "11_expression", "symbol": a.symbol}
    ens = a.ensembl_id

    bundle["gtex_egene"] = _gtex_egene(ens)
    bgee = map_all(ens, ">>ensembl>>bgee") if ens else []
    if bgee:
        b = bgee[0]
        bundle["bgee"] = {"breadth": b.get("expression_breadth"),
                          "present_calls": b.get("total_present_calls"),
                          "max_expression_score": b.get("max_expression_score")}

    # Per-tissue expression (Bgee evidence) — top by expression score.
    # bgee_evidence id is "<ensembl_id>|<UBERON:N or CL:N>"; the tail is the
    # canonical anatomy/cell-type identifier we expose for federated lookup.
    def _anatomy_id(rid):
        last = (rid or "").rsplit("|", 1)[-1]
        return last if last.startswith(("UBERON:", "CL:")) else None
    tissues = map_all(ens, ">>ensembl>>bgee>>bgee_evidence") if ens else []
    tissues.sort(key=lambda t: _num(t.get("expression_score")), reverse=True)
    bundle["top_tissues"] = [{"tissue": t.get("anatomical_entity_name"),
                              "anatomy_id": _anatomy_id(t.get("id")),
                              "score": t.get("expression_score"),
                              "rank": t.get("expression_rank"),
                              "quality": t.get("call_quality")} for t in tissues[:100]]
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

    # Single-cell (scxa) — per-gene marker + max expression across SC experiments
    # via the scxa_expression node (the old >>ensembl>>scxa was dataset-
    # degenerate — biobtree #31). scxa_gene_experiment carries the per-experiment
    # marker flag + max mean expression. No per-experiment description fetch
    # (keeps the disease cohort-fan cheap; the accession is a stable EBI ref).
    sc = {}
    if ens:
        se = map_all(ens, ">>ensembl>>scxa_expression")
        if se:
            exps = [{"experiment_id": e.get("experiment_id"),
                     "is_marker": str(e.get("is_marker_in_experiment")).lower() == "true",
                     "max_expression": e.get("max_mean_expression")}
                    for e in map_all(ens, ">>ensembl>>scxa_expression>>scxa_gene_experiment")]
            exps.sort(key=lambda r: (not r["is_marker"], -_num(r.get("max_expression"))))
            sc = {"total_experiments": _int(se[0].get("total_experiments")),
                  "marker_experiments": _int(se[0].get("marker_experiment_count")),
                  "experiments": exps[:20]}
    bundle["single_cell"] = sc
    return bundle

SECTION = Section(
    id="11", name="expression",
    description="Bgee per-tissue expression + FANTOM5 CAGE (gene + alt promoters) + single-cell datasets",
    needs=("ensembl_id",),
    produces=("bgee", "top_tissues", "fantom5", "fantom5_promoters", "single_cell", "gtex_egene"),
    datasets=DATASETS, chains=CHAINS, collect_fn=collect,
)
