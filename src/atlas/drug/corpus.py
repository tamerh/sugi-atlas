#!/usr/bin/env python3
"""Drug corpus discovery + filtering + batch-build driver.

Mirrors atlas.disease.corpus, but the enumeration source differs. biobtree has
no list-all (retracted #16), so the canonical corpus comes from ChEMBL's own
release: the production `build` parses the ChEMBL SQLite/flat-file dump for
(chembl_id, pref_name, max_phase, atc) and tiers by development phase
(D1 = approved / max_phase 4, D2 ≥ 3, ...). That dump is ~GBs and not wired
yet — until it is, `build` seeds from the curated REFERENCE_DRUGS validation
set (the spec §10 list), resolving each through the live ChEMBL entry to fill
chembl_id / canonical_name / max_phase / slug.

Tiers (by max development phase):
  D1  max_phase == 4   approved
  D2  max_phase >= 3
  D3  max_phase >= 2
  D4  max_phase >= 1   any clinical
  D5  max_phase  < 1   preclinical

Examples:
  python -m atlas.drug.corpus build                      # seed → drug_corpus.json
  python -m atlas.drug.corpus filter --tier D1           # approved subset → stdout
  python -m atlas.drug.corpus run --tier D1 --limit 10   # dev-batch first 10
"""
import argparse, json, os, sys, time
from collections import Counter
from typing import List, Dict, Optional

from atlas.drug.slug import slugify
from atlas.drug.anchors import resolve_chembl, _mol_attrs, _phase

# Curated reference / smoke set (spec §10). Mix of kinase inhibitors, antibodies,
# small molecules, and non-cancer drugs — exercises every section's data shape.
REFERENCE_DRUGS = [
    "Imatinib", "Trastuzumab", "Sotorasib", "Osimertinib", "Pembrolizumab",
    "Olaparib", "Vemurafenib", "Crizotinib", "Alectinib", "Venetoclax",
    "Ibrutinib", "Palbociclib", "Bevacizumab", "Nivolumab", "Lenalidomide",
    "Erlotinib", "Gefitinib", "Sorafenib", "Metformin", "Atorvastatin",
]

DEFAULT_CORPUS_PATH = os.path.join(os.path.dirname(__file__), "..", "..", "..",
                                   "build", "drug_corpus.json")

_TIERS = [("D1", 4), ("D2", 3), ("D3", 2), ("D4", 1), ("D5", 0)]
_TIER_ORDER = ["D1", "D2", "D3", "D4", "D5"]


def tier_of(max_phase: int) -> str:
    for label, cut in _TIERS:
        if max_phase >= cut:
            return label
    return "D5"


def _resolve_record(name: str) -> Optional[Dict]:
    """Lightweight resolve for the corpus row — chembl entry only (no targets/
    chemistry). Returns None if the drug can't be resolved."""
    try:
        chembl_id, en = resolve_chembl(name)
    except Exception as e:
        print(f"  ! {name}: {type(e).__name__}: {e}", flush=True)
        return None
    mol = _mol_attrs(en)
    canonical = mol.get("name") or chembl_id
    max_phase = _phase(mol.get("highestDevelopmentPhase"))
    return {
        "id": chembl_id,
        "input_name": name,
        "canonical_name": canonical,
        "slug": slugify(canonical),
        "molecule_type": mol.get("type"),
        "max_phase": max_phase,
        "atc_codes": list(mol.get("atcClassification") or ()),
        "tier": tier_of(max_phase),
    }


def cmd_build(args):
    seed = REFERENCE_DRUGS
    if args.limit:
        seed = seed[: args.limit]
    print(f"Resolving {len(seed)} seed drugs through ChEMBL ...", flush=True)
    t0 = time.time()
    records = []
    for name in seed:
        rec = _resolve_record(name)
        if rec:
            records.append(rec)
    records.sort(key=lambda r: (-r["max_phase"], r["canonical_name"]))
    print(f"  resolved {len(records)}/{len(seed)} in {time.time()-t0:.1f}s")

    tcounts = Counter(r["tier"] for r in records)
    print("Tier distribution:")
    for label, cut in _TIERS:
        print(f"  {label} (max_phase ≥ {cut}): {tcounts.get(label, 0)}")

    out_path = args.out or os.path.abspath(DEFAULT_CORPUS_PATH)
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, "w") as f:
        json.dump({
            "source": "seed:REFERENCE_DRUGS (validation set); production: ChEMBL release D1",
            "tier_thresholds": {l: c for l, c in _TIERS},
            "total": len(records),
            "tier_counts": dict(tcounts),
            "drugs": records,
        }, f, indent=2)
    print(f"\nCorpus → {out_path} ({len(records)} drugs)")


def load_corpus(path: Optional[str] = None) -> dict:
    p = path or os.path.abspath(DEFAULT_CORPUS_PATH)
    with open(p) as f:
        return json.load(f)


def filter_corpus(corpus: dict, tier: Optional[str] = None,
                  limit: int = 0) -> List[Dict]:
    """Return [{id, name, slug}] records, optionally filtered to tier-or-above
    (D1 = D1 only; D3 = D1+D2+D3). `id` is the ChEMBL id — pass it to resolve()
    to avoid name-collision against the broader search index."""
    rows = corpus.get("drugs", [])
    if tier:
        keep = set(_TIER_ORDER[: _TIER_ORDER.index(tier) + 1])
        rows = [r for r in rows if r["tier"] in keep]
    if limit:
        rows = rows[:limit]
    return [{"id": r["id"], "name": r["canonical_name"], "slug": r["slug"]}
            for r in rows]


def cmd_filter(args):
    records = filter_corpus(load_corpus(args.corpus), tier=args.tier, limit=args.limit)
    print(json.dumps(records, indent=2))


def cmd_run(args):
    """Dev-time serial batch driver — for production use the Enju workflow."""
    from atlas.pipeline import run_drug

    records = filter_corpus(load_corpus(args.corpus), tier=args.tier, limit=args.limit)
    print(f"Driving {len(records)} drug pages ({args.tier or 'all'}) → {args.dist}"
          + ("  [with LLM summary]" if args.summary else "  [deterministic only]"))

    existing = set()
    if not args.force:
        snap_dir = os.path.join(args.dist, "snapshots", "drug")
        if os.path.isdir(snap_dir):
            existing = {f[:-5] for f in os.listdir(snap_dir) if f.endswith(".json")}
        print(f"Already-snapshotted: {len(existing)}")

    results = {"ok": [], "fail": [], "skipped": []}
    t_start = time.time()
    for i, rec in enumerate(records, 1):
        # Pass the ChEMBL id (not the name) to resolve — unambiguous.
        chembl_id, slug = rec["id"], rec["slug"]
        if slug in existing and not args.force:
            results["skipped"].append(slug)
            continue
        t = time.time()
        try:
            run_drug(chembl_id, args.dist, do_summary=args.summary, accept_first_run=True)
            dt = time.time() - t
            print(f"  [{i:>3d}/{len(records)}] {slug} OK ({dt:.1f}s)")
            results["ok"].append({"slug": slug, "seconds": round(dt, 1)})
        except Exception as e:
            print(f"  [{i:>3d}/{len(records)}] {slug} FAIL: {type(e).__name__}: {e}")
            results["fail"].append({"slug": slug, "error": str(e)[:200]})

    elapsed = (time.time() - t_start) / 60
    print(f"\nDone in {elapsed:.1f} min: ok={len(results['ok'])} "
          f"fail={len(results['fail'])} skipped={len(results['skipped'])}")
    if results["fail"]:
        print("\nFailures:")
        for f in results["fail"][:25]:
            print(f"  {f['slug']}: {f['error']}")
    with open("/tmp/drug_batch_run.json", "w") as f:
        json.dump(results, f, indent=2)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    b = sub.add_parser("build", help="Build the corpus (seed = REFERENCE_DRUGS)")
    b.add_argument("--out", default=None)
    b.add_argument("--limit", type=int, default=0, help="resolve first N seed drugs (testing)")
    b.set_defaults(func=cmd_build)

    f = sub.add_parser("filter", help="Filter corpus to a tier + emit [{id,name,slug}] to stdout")
    f.add_argument("--corpus", default=None)
    f.add_argument("--tier", choices=_TIER_ORDER, default=None)
    f.add_argument("--limit", type=int, default=0)
    f.set_defaults(func=cmd_filter)

    r = sub.add_parser("run", help="Dev-time serial batch builder (production: Enju workflow)")
    r.add_argument("--dist", default="/data/sugi-atlas-dist")
    r.add_argument("--corpus", default=None)
    r.add_argument("--tier", choices=_TIER_ORDER, default=None)
    r.add_argument("--limit", type=int, default=0)
    r.add_argument("--force", action="store_true")
    r.add_argument("--summary", action="store_true",
                   help="run the LLM executive summary step (default: deterministic only)")
    r.set_defaults(func=cmd_run)

    args = ap.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
