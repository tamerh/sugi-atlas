#!/usr/bin/env python3
"""Precompute the ORA backgrounds — genome-wide category sizes (K) + the annotated
universe size (N) — for the disease cohort enrichment sections. Each background
fans the SAME chain/classifier the cohort uses over every protein-coding HGNC
gene, so K and the cohort's overlap k are counted identically (the hypergeometric
test stays self-consistent). Refresh per biobtree DATA release (alongside seed
regeneration), since the biobtree version stamp does not move on a data-only
re-index.

  python -m atlas.build_background           # all three
  python -m atlas.build_background reactome  # just one

Outputs: data/background/{reactome,go,family}.json = {chain, biobtree_version,
universe_n, sizes}.
"""
import gzip
import json
import os
import sys
from collections import Counter
from multiprocessing import Pool

from atlas.biobtree import map_all
from atlas import pipeline
from atlas.disease.sections.s06_protein_families import _classify

SEEDS = "data/corpus/seeds/genes_hgnc_protein_coding.txt"
OUTDIR = "data/background"


# Each worker returns (annotated_bool, [category_ids]) for one gene. annotated =
# the gene contributes to the universe N (has >=1 category); for the family
# background every gene is classified, so it always counts toward N.
def _reactome(sym):
    try:
        ids = sorted({p["id"] for p in map_all(sym, ">>hgnc>>ensembl>>reactome", cap=10) if p.get("id")})
    except Exception:
        ids = []
    return (bool(ids), ids)


def _go_bp(sym):
    try:
        ids = sorted({t["id"] for t in map_all(sym, ">>hgnc>>uniprot>>go", cap=10)
                      if t.get("id") and t.get("type") == "biological_process"})
    except Exception:
        ids = []
    return (bool(ids), ids)


def _family(sym):
    try:
        names = [t.get("short_name") for t in map_all(sym, ">>hgnc>>uniprot>>interpro", cap=20)]
        has_ec = bool(map_all(sym, ">>hgnc>>uniprot>>brenda", cap=1))
    except Exception:
        names, has_ec = [], False
    return (True, [_classify([n for n in names if n], has_ec)])   # every gene gets a family


_BUILDERS = {"reactome": _reactome, "go": _go_bp, "family": _family}


def build(name):
    fn = _BUILDERS[name]
    syms = [l.strip() for l in open(SEEDS) if l.strip()]
    print(f"[background:{name}] fanning over {len(syms)} protein-coding genes …", flush=True)
    sizes: Counter = Counter()
    annotated = 0
    with Pool(16) as pool:
        for i, (is_ann, ids) in enumerate(pool.imap_unordered(fn, syms, chunksize=8), 1):
            if is_ann:
                annotated += 1
            for cid in ids:
                sizes[cid] += 1
            if i % 5000 == 0:
                print(f"  {i}/{len(syms)} — {annotated} annotated, {len(sizes)} categories", flush=True)
    out = {"chain": name, "biobtree_version": pipeline.biobtree_version(),
           "universe_n": annotated, "sizes": dict(sizes)}
    os.makedirs(OUTDIR, exist_ok=True)
    path = os.path.join(OUTDIR, f"{name}.json")
    with open(path, "w") as f:
        json.dump(out, f, sort_keys=True)
    print(f"[background:{name}] universe_n={annotated}, {len(sizes)} categories → {path}", flush=True)


def _membership(sym):
    """Per-gene Reactome pathway + GO-BP (id, name) pairs (the same chains the
    backgrounds use). Powers interaction-partner ORA on gene pages: a gene's
    partner set is looked up against this table — no per-gene fan at corpus build."""
    try:
        rc = [(p["id"], p.get("name")) for p in map_all(sym, ">>hgnc>>ensembl>>reactome", cap=10) if p.get("id")]
        go = [(t["id"], t.get("name")) for t in map_all(sym, ">>hgnc>>uniprot>>go", cap=10)
              if t.get("id") and t.get("type") == "biological_process"]
    except Exception:
        rc, go = [], []
    return (sym, rc, go)


def build_membership():
    syms = [l.strip() for l in open(SEEDS) if l.strip()]
    print(f"[membership] reactome + GO-BP per gene over {len(syms)} genes …", flush=True)
    mem, names = {}, {}
    with Pool(16) as pool:
        for i, (sym, rc, go) in enumerate(pool.imap_unordered(_membership, syms, chunksize=8), 1):
            rc_ids = sorted({i for i, _ in rc})
            go_ids = sorted({i for i, _ in go})
            if rc_ids or go_ids:
                mem[sym] = {"reactome": rc_ids, "go": go_ids}
            for cid, nm in rc + go:
                if nm and cid not in names:
                    names[cid] = nm
            if i % 5000 == 0:
                print(f"  {i}/{len(syms)} — {len(mem)} with annotation", flush=True)
    os.makedirs(OUTDIR, exist_ok=True)
    path = os.path.join(OUTDIR, "membership.json.gz")
    blob = json.dumps({"biobtree_version": pipeline.biobtree_version(),
                       "genes": mem, "names": names}, sort_keys=True).encode()
    with open(path, "wb") as raw, gzip.GzipFile(fileobj=raw, mode="wb", mtime=0) as gz:
        gz.write(blob)                   # mtime=0 → byte-reproducible across rebuilds
    print(f"[membership] {len(mem)} genes, {len(names)} names → {path} "
          f"({os.path.getsize(path)//1024} KB)", flush=True)


def main(argv=None):
    argv = argv if argv is not None else sys.argv[1:]
    names = argv or (list(_BUILDERS) + ["membership"])
    for name in names:
        if name == "membership":
            build_membership()
        else:
            build(name)


if __name__ == "__main__":
    main()
