#!/usr/bin/env python3
"""Regenerate the corpus SEED lists from biobtree's own ingest sources.

Writes corpus/seeds/*.txt (gitignored — regenerable, large). These are the
canonical entity universes the Atlas can build, derived from the exact files
biobtree ingests (see /data/biobtree/conf/source1.dataset.json):

  genes   — HGNC complete set (hgnc dataset, id 10). Our gene slug == HGNC
            symbol and the anchor resolves via HGNC, so HGNC is the complete +
            correct gene corpus: a locus with no HGNC symbol can't be an Atlas
            page. Split protein-coding (full §-content) vs ncRNA (RNAcentral).
  drugs   — ChEMBL molecules (chembl_molecule dataset, id 22) local JSONL,
            filtered by max_phase: approved (4) and clinical (1-3), named only.
  disease — the ranked Mondo corpus already built by atlas.disease.corpus
            (build/mondo_corpus.json), ordered by signal_score.

Usage:  python scripts/build_corpus.py [--out corpus/seeds]
"""
import argparse
import json
import os
import urllib.request

HGNC_URL = "https://storage.googleapis.com/public-download-files/hgnc/json/json/hgnc_complete_set.json"
CHEMBL_JSONL = "/data/biobtree/raw_data/chembl/extracted/chembl_molecules.jsonl"
MONDO_CORPUS = os.path.join(os.path.dirname(__file__), "..", "build", "mondo_corpus.json")


def _write(path, items):
    with open(path, "w") as f:
        f.write("\n".join(items) + "\n")
    print(f"  {path}: {len(items)}")


def genes(out):
    print("HGNC complete set …")
    with urllib.request.urlopen(HGNC_URL, timeout=60) as r:
        docs = json.load(r).get("response", {}).get("docs", [])
    pc = sorted({d["symbol"] for d in docs
                 if d.get("locus_group") == "protein-coding gene"
                 and d.get("status") == "Approved" and d.get("symbol")})
    nc = sorted({d["symbol"] for d in docs
                 if d.get("locus_group") == "non-coding RNA" and d.get("symbol")})
    _write(os.path.join(out, "genes_hgnc_protein_coding.txt"), pc)
    _write(os.path.join(out, "genes_hgnc_ncrna.txt"), nc)


def drugs(out):
    print(f"ChEMBL molecules ({CHEMBL_JSONL}) …")
    approved, clinical = set(), set()
    with open(CHEMBL_JSONL) as f:
        for ln in f:
            try:
                r = json.loads(ln)
            except ValueError:
                continue
            nm = (r.get("name") or "").strip()
            if not nm:
                continue
            mp = r.get("max_phase")
            if mp == 4:
                approved.add(nm)
            elif mp in (1, 2, 3):
                clinical.add(nm)
    _write(os.path.join(out, "drugs_chembl_approved.txt"), sorted(approved))
    _write(os.path.join(out, "drugs_chembl_clinical.txt"), sorted(clinical))


def diseases(out):
    print(f"Mondo ranked corpus ({MONDO_CORPUS}) …")
    d = json.load(open(MONDO_CORPUS))["diseases"]  # already signal-ranked
    _write(os.path.join(out, "diseases_mondo_ranked.txt"),
           [e["canonical_name"] for e in d if e.get("canonical_name")])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=os.path.join(os.path.dirname(__file__), "..", "corpus", "seeds"))
    ap.add_argument("--only", choices=["genes", "drugs", "diseases"])
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    if a.only in (None, "genes"):
        genes(a.out)
    if a.only in (None, "drugs"):
        drugs(a.out)
    if a.only in (None, "diseases"):
        diseases(a.out)


if __name__ == "__main__":
    main()
