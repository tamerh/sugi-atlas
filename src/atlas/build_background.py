#!/usr/bin/env python3
"""Precompute the Reactome ORA background — genome-wide pathway sizes (K) and the
annotated universe size (N) — by fanning the SAME `>>hgnc>>ensembl>>reactome`
chain the disease cohort uses over every protein-coding HGNC gene. Counting K the
identical way the cohort counts its overlap k keeps the hypergeometric test
self-consistent. Refresh per biobtree DATA release (alongside seed regeneration),
since the biobtree version stamp does not move on a data-only re-index.

  python -m atlas.build_background

Output: data/background/reactome.json = {chain, biobtree_version, universe_n, sizes}
"""
import json
import os
from collections import Counter
from multiprocessing import Pool

from atlas.biobtree import map_all
from atlas import pipeline

SEEDS = "data/corpus/seeds/genes_hgnc_protein_coding.txt"
OUT = "data/background/reactome.json"
CHAIN = ">>hgnc>>ensembl>>reactome"


def _pathways(sym):
    """Distinct Reactome pathway ids for one gene (one count per gene/pathway)."""
    try:
        rows = map_all(sym, CHAIN, cap=10)
    except Exception:
        return (sym, [])
    return (sym, sorted({p["id"] for p in rows if p.get("id")}))


def main():
    syms = [l.strip() for l in open(SEEDS) if l.strip()]
    print(f"[background] fanning {CHAIN} over {len(syms)} protein-coding genes …", flush=True)
    sizes: Counter = Counter()
    annotated = 0
    with Pool(16) as pool:
        for i, (_sym, pids) in enumerate(pool.imap_unordered(_pathways, syms, chunksize=8), 1):
            if pids:
                annotated += 1
                for pid in pids:
                    sizes[pid] += 1
            if i % 2500 == 0:
                print(f"  {i}/{len(syms)} — {annotated} annotated, {len(sizes)} pathways",
                      flush=True)
    out = {"chain": CHAIN,
           "biobtree_version": pipeline.biobtree_version(),
           "universe_n": annotated,
           "sizes": dict(sizes)}
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w") as f:
        json.dump(out, f, sort_keys=True)
    print(f"[background] universe_n={annotated}, {len(sizes)} pathways → {OUT}", flush=True)


if __name__ == "__main__":
    main()
