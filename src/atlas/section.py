"""Shared Section contract — used by both gene and disease entity pipelines.

A Section is a first-class metadata record so graph.py / provenance.py can
introspect "what this section reads" without importing the collect function.
The contract is identical across entity types; only the anchors type passed
to `collect_fn` differs (atlas.gene.anchors.Anchors vs
atlas.disease.anchors.DiseaseAnchors)."""
from dataclasses import dataclass
from typing import Callable, Tuple

@dataclass(frozen=True)
class Section:
    id: str                     # "1".."N"
    name: str                   # short identifier
    description: str            # one-line what this section covers
    needs: Tuple[str, ...]      # anchor fields this section reads
    produces: Tuple[str, ...]   # primary bundle keys it sets
    datasets: Tuple[str, ...]   # biobtree datasets touched (for the graph)
    chains: Tuple[str, ...]     # biobtree chain strings used (for the graph)
    collect_fn: Callable        # collect(anchors) -> bundle dict
    # Bundle keys whose count is allowed to shrink without flagging a
    # regression in body_gate. Use for fields known to legitimately
    # fluctuate when biobtree (or its upstream) recurates — e.g. the
    # 2026-05-30 RefSeq REVIEWED-only filter dropped TP53 mRNA 46→25.
    # Default empty: every key strictly monitored.
    shrinkable: Tuple[str, ...] = ()
