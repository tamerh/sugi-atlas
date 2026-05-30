"""Section contract — the data each gene section module declares.

Sections are first-class objects: their CHAINS / DATASETS / NEEDS metadata is
machine-readable, so graph.py can emit Mermaid diagrams of the data flow with
zero drift between code and figure."""
from dataclasses import dataclass, field
from typing import Callable, Tuple

@dataclass(frozen=True)
class Section:
    id: str            # "1".."12"
    name: str          # short identifier: gene_ids, transcripts, ...
    description: str   # one-line what this section covers
    needs: Tuple[str, ...]      # anchors this section reads (hgnc_id, ensembl_id, ...)
    produces: Tuple[str, ...]   # primary bundle keys it sets
    datasets: Tuple[str, ...]   # biobtree datasets touched (for the graph)
    chains: Tuple[str, ...]     # biobtree chain strings used (for the graph)
    collect_fn: Callable        # collect(anchors) -> bundle dict
