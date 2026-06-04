#!/usr/bin/env python3
"""Regenerate the corpus SEED lists from biobtree's own ingest sources.

Writes data/corpus/seeds/*.txt (gitignored — regenerable, large). These are the
canonical entity universes the Atlas can build, derived from the exact files
biobtree ingests (see /data/biobtree/conf/source1.dataset.json):

  genes   — HGNC complete set (hgnc dataset, id 10). Our gene slug == HGNC
            symbol and the anchor resolves via HGNC, so HGNC is the complete +
            correct gene corpus: a locus with no HGNC symbol can't be an Atlas
            page. Split protein-coding (full §-content) vs ncRNA (RNAcentral).
  drugs   — ChEMBL molecules (chembl_molecule dataset, id 22) local JSONL,
            filtered by max_phase: approved (4) and clinical (1-3), named only.
  disease — the ranked Mondo corpus already built by atlas.disease.corpus
            (dist/build/mondo_corpus.json), ordered by signal_score.

Usage:  python -m atlas.build_corpus [--out data/corpus/seeds]
"""
import argparse
import json
import os
import urllib.request
from multiprocessing import Pool

# Repo root = two levels up from this module (src/atlas/build_corpus.py).
_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
HGNC_URL = "https://storage.googleapis.com/public-download-files/hgnc/json/json/hgnc_complete_set.json"
CHEMBL_JSONL = "/data/biobtree/raw_data/chembl/extracted/chembl_molecules.jsonl"
MONDO_CORPUS = os.path.join(_ROOT, "dist", "build", "mondo_corpus.json")

# Gate 2: ChEBI roles that are NON-therapeutic — a molecule whose ONLY roles are
# these (and which has no ATC and no curated target) is a reagent/excipient/
# metabolite, not a drug (water's roles = solvent/greenhouse gas/metabolite).
# Substring match, lowercased. Conservative — anything NOT here counts as a
# (possibly) pharmacological role, so we keep generously.
_NONPHARMA_ROLE = (
    "solvent", "metabolite", "greenhouse gas", "fertilis", "fertiliz", "fuel",
    "reference compound", "nmr", "chemical shift", "food", "nutrient",
    "contaminant", "pollutant", "reagent", "dye", "stain", "buffer",
    "cosmetic", "flavour", "flavor", "fragrance",
)


def _is_pharma_role(role: str) -> bool:
    r = (role or "").lower()
    return r != "" and not any(k in r for k in _NONPHARMA_ROLE)


def _gate_one(cid: str) -> dict:
    """Resolve a ChEMBL id and compute the gate-2/3 signals. Top-level for Pool."""
    from atlas.drug.anchors import resolve as resolve_drug
    try:
        a = resolve_drug(cid)
    except (Exception, SystemExit) as e:   # SystemExit too (harden M1)
        return {"id": cid, "ok": False, "err": repr(e)[:80]}
    roles = tuple(getattr(a, "chebi_roles", None) or ())
    # Drop ONLY a clear reagent: HAS ChEBI roles and ALL are non-pharmacological
    # (water = solvent/metabolite/greenhouse-gas). KEEP if it has ATC, a curated
    # target, a pharma role, OR no roles at all — because many real (older/
    # obscure) drugs are simply un-annotated in ChEMBL (mivacurium/loracarbef:
    # atc=() targets=0 roles=()), and under-gate > over-gate.
    therapeutic = (bool(a.atc_codes) or bool(a.targets) or not roles
                   or any(_is_pharma_role(r) for r in roles))
    return {"id": cid, "ok": True, "name": a.canonical_name,
            "therapeutic": therapeutic, "parent": a.parent_chembl,
            "had_roles": bool(roles)}


def gate_drugs(ids, workers=16):
    """Gate 2 (therapeutic-value filter) + Gate 3 (salt→parent collapse) over a
    set of ChEMBL ids. Returns (kept_ids, dropped[(id,name,reason)])."""
    ids = sorted(set(ids))
    with Pool(workers) as p:
        res = p.map(_gate_one, ids)
    present = {r["id"] for r in res if r.get("ok")}
    kept, dropped = [], []
    for r in res:
        if not r.get("ok"):
            dropped.append((r["id"], "", f"resolve-failed: {r.get('err','')}"))
        elif r["parent"] and r["parent"] in present:   # salt check first → correct label
            dropped.append((r["id"], r["name"], f"salt/child of {r['parent']}"))
        elif not r["therapeutic"]:
            dropped.append((r["id"], r["name"], "non-therapeutic (no ATC / target / pharma ChEBI role)"))
        else:
            kept.append(r["id"])
    return kept, dropped


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
    # Combined seed (protein-coding + ncRNA) — the full buildable gene corpus.
    # ncRNA pages are thinner (no protein §) but supported (RNAcentral + ncRNA
    # biotypes). Pseudogene/other excluded as too thin.
    _write(os.path.join(out, "genes_hgnc.txt"), sorted(set(pc) | set(nc)))


def drugs(out):
    print(f"ChEMBL molecules ({CHEMBL_JSONL}) …")
    # Seed by molecule_id (CHEMBLxxx), NOT name — resolve_drug() then does a
    # direct entry fetch (reliable) instead of a flaky name search (weird
    # chemical names like "3,3',4',5-TETRACHLOROSALICYLANILIDE" don't round-trip
    # through search). Require a name so the page has a human label/slug.
    approved, clinical = set(), set()
    with open(CHEMBL_JSONL) as f:
        for ln in f:
            try:
                r = json.loads(ln)
            except ValueError:
                continue
            mid = (r.get("molecule_id") or "").strip()
            if not mid or not (r.get("name") or "").strip():
                continue
            mp = r.get("max_phase")
            if mp == 4:
                approved.add(mid)
            elif mp in (1, 2, 3):
                clinical.add(mid)
    _write(os.path.join(out, "drugs_chembl_clinical.txt"), sorted(clinical))
    # Gates 2+3: refine the approved set — resolve each, drop reagents/excipients
    # (no ATC/target/pharma-role) and salt-form children whose parent is also
    # approved. Emits an audit drop-list to eyeball before a full run.
    print(f"  gating {len(approved)} approved (resolve + filter) …", flush=True)
    kept, dropped = gate_drugs(approved)
    n_reagent = sum(1 for _, _, why in dropped if why.startswith("non-therapeutic"))
    n_salt = sum(1 for _, _, why in dropped if why.startswith("salt"))
    n_fail = sum(1 for _, _, why in dropped if why.startswith("resolve"))
    print(f"  gate2/3: kept {len(kept)}/{len(approved)} — dropped "
          f"{n_reagent} reagents, {n_salt} salts, {n_fail} resolve-fails")
    _write(os.path.join(out, "drugs_chembl_approved.txt"), sorted(kept))
    _write(os.path.join(out, "drugs_chembl_approved.dropped.txt"),
           [f"{i}\t{nm}\t{why}" for i, nm, why in sorted(dropped)])


def diseases(out):
    print(f"Mondo ranked corpus ({MONDO_CORPUS}) …")
    # Seed by MONDO id (signal-ranked), NOT canonical_name — resolve_disease()
    # then does a direct lookup instead of a name search (long subtype names
    # like "epidermolysis bullosa simplex 5C, with pyloric atresia" don't
    # round-trip through search). Page slug still derives from canonical_name.
    d = [e for e in json.load(open(MONDO_CORPUS))["diseases"] if e.get("id")]

    # Admission gates: (1) the disease-characteristic qualifier subtree
    # (inherited/acquired/sporadic/… — not diseases), and (2) non-human /
    # veterinary diseases (achondroplasia-in-cattle, canine rhabdomyosarcoma —
    # out of scope, render empty). Surgical subtrees, NOT the broad
    # disease_grouping subset (which holds real human hubs).
    from atlas.disease.corpus import (parse_obo, characteristic_ids,
                                      non_human_ids, OBO_PATH)
    obo = os.path.abspath(OBO_PATH)   # canonical cache path (corpus json records a stale bin/.. path)
    if not os.path.exists(obo):
        raise FileNotFoundError(f"Mondo OBO not found at {obo} — run atlas.disease.corpus ensure_obo")
    terms = parse_obo(obo)
    qual = characteristic_ids(terms)
    nonhuman = non_human_ids(terms)
    drop = qual | nonhuman
    kept = [e for e in d if e["id"] not in drop]
    dropped = [e for e in d if e["id"] in drop]
    print(f"  gate1: dropped {sum(1 for e in d if e['id'] in qual)} qualifier nodes")
    print(f"  gate2: dropped {sum(1 for e in d if e['id'] in nonhuman)} non-human / "
          f"veterinary nodes (e.g. "
          f"{[e['canonical_name'] for e in d if e['id'] in nonhuman][:4]})")
    _write(os.path.join(out, "diseases_mondo_ranked.txt"), [e["id"] for e in kept])
    # Audit trail — eyeball before a full run (id  name  signal).
    _write(os.path.join(out, "diseases_mondo_ranked.dropped.txt"),
           [f"{e['id']}\t{e.get('canonical_name')}\t{e.get('signal_score')}" for e in dropped])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=os.path.join(_ROOT, "data", "corpus", "seeds"))
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
