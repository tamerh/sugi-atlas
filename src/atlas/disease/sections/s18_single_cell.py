"""§18 — single-cell datasets (CZ CELLxGENE): disease-level scRNA-seq data
availability. The disease anchors directly to CELLxGENE datasets annotated with
its Mondo term (>>mondo>>cellxgene) — a research-resource pointer, NOT a per-gene
or cohort signal (id 18 because 15-17 are render-only sections; see sections/__init__).

Only the DATASET list is surfaced (count + largest by cell number). The
`cellxgene_celltype` edge is deliberately NOT used: it is contaminated by
multi-disease atlases — e.g. type-2-diabetes returns cortical-neuron cell types
carried by brain atlases that merely co-annotate the disease — so a disease-level
cell-type list would mislead. Framed honestly: these are datasets *annotated with*
the disease, some of them broad multi-disease atlases.
"""
from atlas.biobtree import map_all
from atlas.section import Section

CHAINS   = (">>mondo>>cellxgene",)
DATASETS = ("mondo", "cellxgene")

_TOP_DATASETS = 8


def _cells(d):
    try:
        return int(d.get("cell_count") or 0)
    except (TypeError, ValueError):
        return 0


def collect(a):
    ds = map_all(a.mondo_id, ">>mondo>>cellxgene")
    seen, rows, total_cells = set(), [], 0
    for d in ds:
        did = d.get("id")
        if not did or did in seen:
            continue
        seen.add(did)
        n = _cells(d)
        total_cells += n
        rows.append({"id": did, "title": (d.get("title") or "").strip(),
                     "organism": d.get("organism"), "cells": n})
    rows.sort(key=lambda r: -r["cells"])
    return {
        "section": "18_single_cell",
        "mondo_id": a.mondo_id,
        "dataset_count": len(rows),
        "total_cells": total_cells,
        "top_datasets": rows[:_TOP_DATASETS],
    }


SECTION = Section(
    id="18", name="single_cell",
    description=("Disease-level single-cell RNA-seq data availability from CZ "
                 "CELLxGENE (mondo→cellxgene): dataset count + the largest datasets "
                 "by cell number. A research-resource pointer, not a cohort signal; "
                 "the cell-type edge is excluded as multi-disease-atlas contaminated."),
    needs=("mondo_id",),
    produces=("dataset_count", "total_cells", "top_datasets"),
    datasets=DATASETS, chains=CHAINS, collect_fn=collect,
)
