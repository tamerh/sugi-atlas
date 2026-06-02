#!/usr/bin/env python3
"""Disease corpus discovery + filtering + batch-build driver.

Phase 1 (`build`): parses Mondo's OBO release, probes biobtree per term for
xref counts, computes Atlas signal score, emits a ranked JSON corpus.

Phase 2 (`filter` / `run`): consumes the corpus to drive page builds.
  - `filter` emits a tiered subset as JSON (feedable to the Enju workflow's
    `diseases` param, which is `list<{name, slug}>`)
  - `run` is the dev-time serial driver — calls atlas.pipeline.run_disease
    per disease, no Enju. Production runs should use the Enju workflow at
    src/atlas/disease/enju.yaml for retries + parallelism + audit trail.

Examples:
  python -m atlas.disease.corpus build                       # 28k Mondo → corpus.json
  python -m atlas.disease.corpus build --refresh-obo         # re-download OBO
  python -m atlas.disease.corpus filter --tier T1            # 3.3k T1 records → stdout
  python -m atlas.disease.corpus run --tier T1 --limit 20    # dev-batch first 20

The corpus JSON shape — one record per Mondo class with material signal:
  {
    "id": "MONDO:0005233",
    "canonical_name": "non-small cell lung carcinoma",
    "slug": "non-small-cell-lung-carcinoma",
    "signal_score": 1247.3,
    "tier": "T1",
    "xref_counts": {"gwas": 120, "civic_evidence": 849, ...}
  }
"""
import argparse, json, os, re, sys, time, urllib.parse, urllib.request
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import List, Dict, Optional

from atlas.disease.slug import slugify

# ---------------------------------------------------------------------------
# OBO acquisition + parsing
# ---------------------------------------------------------------------------

OBO_URL = "http://purl.obolibrary.org/obo/mondo.obo"
CACHE_DIR = os.path.join(os.path.dirname(__file__), "..", "..", "..", "dist", "build", "cache")
OBO_PATH = os.path.join(CACHE_DIR, "mondo.obo")
DEFAULT_CORPUS_PATH = os.path.join(os.path.dirname(__file__), "..", "..", "..",
                                   "dist", "build", "mondo_corpus.json")
API = "http://127.0.0.1:8000/api"


def ensure_obo(refresh: bool = False) -> str:
    """Download Mondo OBO once, cache locally. ~50 MB; refresh quarterly."""
    os.makedirs(CACHE_DIR, exist_ok=True)
    path = os.path.abspath(OBO_PATH)
    if not refresh and os.path.exists(path) and os.path.getsize(path) > 1_000_000:
        return path
    print(f"Downloading Mondo OBO from {OBO_URL} ...", flush=True)
    t = time.time()
    with urllib.request.urlopen(OBO_URL, timeout=120) as r, open(path, "wb") as f:
        f.write(r.read())
    print(f"  cached {os.path.getsize(path)/1e6:.1f} MB in {time.time()-t:.1f}s -> {path}",
          flush=True)
    return path


# OBO term-block parsing.
_TERM_BLOCK = re.compile(r"^\[Term\]\s*$", re.M)
_ID_LINE    = re.compile(r"^id:\s*(MONDO:\d+)\s*$", re.M)
_NAME_LINE  = re.compile(r"^name:\s*(.+?)\s*$", re.M)
_OBSOLETE   = re.compile(r"^is_obsolete:\s*true\s*$", re.M)
_ISA_LINE   = re.compile(r"^is_a:\s*(MONDO:\d+)", re.M)

# Admission gate 1: the "disease characteristic" subtree — PATO-rooted
# qualities (inherited / acquired / sporadic / congenital / X-linked) that are
# NOT diseases but accrue word-matched evidence and float to the top of the
# ranked corpus (e.g. MONDO:0021152 "inherited" rendered a full junk page).
# We drop ONLY this subtree. Deliberately NOT gating on the broad
# `disease_grouping` subset — that also flags real hub diseases (cardiomyopathy,
# AML, renal cell carcinoma). See docs/research/06_admission_gates.md.
DISEASE_CHARACTERISTIC_ROOT = "MONDO:0021125"


def parse_obo(path: str) -> List[Dict]:
    """Parse all non-obsolete MONDO term blocks → [{id, name, parents}, ...]."""
    with open(path) as f:
        text = f.read()
    parts = _TERM_BLOCK.split(text)[1:]  # drop header before first [Term]
    out = []
    for block in parts:
        if _OBSOLETE.search(block):
            continue
        mid = _ID_LINE.search(block)
        if not mid:
            continue
        nm = _NAME_LINE.search(block)
        out.append({"id": mid.group(1), "name": (nm.group(1) if nm else None),
                    "parents": _ISA_LINE.findall(block)})
    return out


def characteristic_ids(terms: List[Dict], root: str = DISEASE_CHARACTERISTIC_ROOT) -> set:
    """MONDO ids in the `disease characteristic` qualifier subtree (root + all
    is_a descendants) — the admission-gate-1 exclusion set. Pure OBO graph walk,
    no network. `terms` is parse_obo() output."""
    children = {}
    for t in terms:
        for p in t.get("parents") or ():
            children.setdefault(p, []).append(t["id"])
    out, stack = set(), [root]
    while stack:
        n = stack.pop()
        for c in children.get(n, ()):
            if c not in out:
                out.add(c)
                stack.append(c)
    out.add(root)
    return out


# ---------------------------------------------------------------------------
# Signal scoring (weighted xref sum)
# ---------------------------------------------------------------------------

# Curated > literature-mined > raw counts. clinical_trials damped /100 below
# so high-trial diseases (breast cancer ≈11k trials) don't drown out
# mechanistically rich classes with fewer trials.
SIGNAL_WEIGHTS = {
    "gencc":           3.0,
    "civic_evidence":  2.0,
    "intogen":         2.0,
    "gwas":            1.0,
    "clinvar":         1.0,
    "hpo":             1.0,
    "civic_assertion": 1.0,
    "gwas_study":      0.5,
}

TIERS = [
    ("T1", 100.0),
    ("T2",  50.0),
    ("T3",  20.0),
    ("T4",   5.0),
    ("T5",   0.001),
]


def tier_of(score: float) -> str:
    for label, cut in TIERS:
        if score >= cut:
            return label
    return "T6"


def signal_score(xref_counts: dict) -> float:
    """Weighted sum of signal-carrying xref counts."""
    s = 0.0
    for key, w in SIGNAL_WEIGHTS.items():
        s += w * (xref_counts.get(key) or 0)
    s += 1.0 * ((xref_counts.get("clinical_trials") or 0) / 100.0)
    return round(s, 1)


# ---------------------------------------------------------------------------
# biobtree probing
# ---------------------------------------------------------------------------

def _entry_xref_counts(mondo_id: str) -> Dict[str, int]:
    """Single biobtree /entry call → xref count dict. Returns {} on error."""
    url = f"{API}/entry?i={urllib.parse.quote(mondo_id)}&s=mondo"
    try:
        with urllib.request.urlopen(url, timeout=15) as r:
            d = json.load(r)
    except Exception:
        return {}
    out = {}
    for line in (d.get("xrefs") or {}).get("data") or []:
        try:
            k, v = line.split("|", 1)
            out[k] = int(v)
        except ValueError:
            continue
    return out


def probe_parallel(terms: List[Dict], workers: int = 8) -> List[Dict]:
    """For each term, fetch xref counts and compute signal score. Returns
    enriched records (zero-signal terms included; caller filters)."""
    out = []
    n = len(terms)
    t0 = time.time()
    last_report = t0
    with ThreadPoolExecutor(max_workers=workers) as pool:
        futs = {pool.submit(_entry_xref_counts, t["id"]): t for t in terms}
        for i, fut in enumerate(as_completed(futs), 1):
            term = futs[fut]
            xrefs = fut.result()
            score = signal_score(xrefs)
            out.append({
                "id": term["id"],
                "canonical_name": term["name"],
                "slug": slugify(term["name"] or term["id"]),
                "signal_score": score,
                "tier": tier_of(score),
                "xref_counts": xrefs,
            })
            now = time.time()
            if now - last_report > 5:
                rate = i / max(now - t0, 1e-9)
                eta = (n - i) / rate
                print(f"  probed {i:>5d}/{n} ({i*100//n:>2d}%), "
                      f"rate={rate:.0f}/s, eta={eta/60:.1f}min", flush=True)
                last_report = now
    return out


# ---------------------------------------------------------------------------
# Top-level commands
# ---------------------------------------------------------------------------

def cmd_build(args):
    obo = args.obo or ensure_obo(refresh=args.refresh_obo)
    print(f"Parsing OBO {obo} ...", flush=True)
    terms = parse_obo(obo)
    print(f"  {len(terms)} non-obsolete MONDO terms")
    # Admission gate 1: drop the disease-characteristic qualifier subtree
    # (inherited/acquired/… — not diseases) before probing. Saves the probes
    # and keeps the corpus clean at the source.
    qual = characteristic_ids(terms)
    before = len(terms)
    terms = [t for t in terms if t["id"] not in qual]
    print(f"  gate1: dropped {before - len(terms)} disease-characteristic qualifier nodes")
    if args.limit:
        terms = terms[: args.limit]
        print(f"  (--limit {args.limit} → first {len(terms)})")
    print(f"Probing biobtree (parallel={args.parallel}) ...", flush=True)
    t0 = time.time()
    enriched = probe_parallel(terms, workers=args.parallel)
    print(f"  done in {(time.time()-t0)/60:.1f} min")

    nonzero = sorted([e for e in enriched if e["signal_score"] > 0],
                     key=lambda e: -e["signal_score"])
    print(f"\nKept {len(nonzero)} of {len(enriched)} terms (signal > 0)")
    tcounts = Counter(e["tier"] for e in nonzero)
    print("\nTier distribution:")
    for label, cut in TIERS:
        print(f"  {label} (score ≥ {cut:g}): {tcounts.get(label, 0):>5d}")
    print("\nTop 10 by signal:")
    for e in nonzero[:10]:
        print(f"  {e['signal_score']:>9.1f}  {e['id']:14s}  {e['tier']}  {e['canonical_name']}")

    out_path = args.out or os.path.abspath(DEFAULT_CORPUS_PATH)
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, "w") as f:
        json.dump({
            "obo_version_file": obo,
            "weights": SIGNAL_WEIGHTS,
            "tier_thresholds": dict(TIERS),
            "total_terms_in_obo": len(enriched),
            "kept_nonzero": len(nonzero),
            "tier_counts": dict(tcounts),
            "diseases": nonzero,
        }, f, indent=2)
    print(f"\nCorpus → {out_path} ({os.path.getsize(out_path)/1e6:.1f} MB)")


def load_corpus(path: Optional[str] = None) -> dict:
    """Read corpus JSON; default location is build/mondo_corpus.json."""
    p = path or os.path.abspath(DEFAULT_CORPUS_PATH)
    with open(p) as f:
        return json.load(f)


# Tier order for "≥ tier" filtering — keep T1 includes T1 only (smaller is
# stricter). To include "T1 and better" use cumulative across the list.
_TIER_ORDER = ["T1", "T2", "T3", "T4", "T5"]


def filter_corpus(corpus: dict, tier: Optional[str] = None,
                  limit: int = 0) -> List[Dict]:
    """Return [{id, name, slug}] records, optionally filtered to tier-or-above.

    `tier='T2'` means T1+T2 (stricter tier + everything above it gets included).
    Record shape:
      id    — MONDO:NNN (the canonical mondo id; pass this to resolve() to
              avoid name-collision against the broader search index — e.g.
              biobtree's name-search for 'cardiomyopathy' returns DCM-1G
              by xref_count rank, but the corpus knows the umbrella's id)
      name  — canonical_name from Mondo
      slug  — URL-safe slug derived from canonical_name
    Compatible with the Enju workflow's `diseases` list<record> (extra
    fields silently ignored by the workflow)."""
    rows = corpus.get("diseases", [])
    if tier:
        keep = set(_TIER_ORDER[: _TIER_ORDER.index(tier) + 1])
        rows = [r for r in rows if r["tier"] in keep]
    if limit:
        rows = rows[:limit]
    return [{"id": r["id"], "name": r["canonical_name"], "slug": r["slug"]}
            for r in rows]


def cmd_filter(args):
    corpus = load_corpus(args.corpus)
    records = filter_corpus(corpus, tier=args.tier, limit=args.limit)
    print(json.dumps(records, indent=2))


def cmd_run(args):
    """Dev-time serial batch driver — for production use the Enju workflow."""
    from atlas.pipeline import run_disease
    from atlas.disease.anchors import resolve as resolve_disease

    corpus = load_corpus(args.corpus)
    records = filter_corpus(corpus, tier=args.tier, limit=args.limit)
    print(f"Driving {len(records)} pages ({args.tier or 'all'}) → {args.dist}")

    existing = set()
    if not args.force:
        snap_dir = os.path.join(args.dist, "snapshots", "disease")
        if os.path.isdir(snap_dir):
            existing = {f[:-5] for f in os.listdir(snap_dir) if f.endswith(".json")}
        print(f"Already-snapshotted: {len(existing)}")

    results = {"ok": [], "fail": [], "skipped": []}
    t_start = time.time()
    for i, rec in enumerate(records, 1):
        # Pass the MONDO id (not the name) to resolve — avoids the
        # name-collision where biobtree's search picks a high-xref
        # subtype over an umbrella term ('cardiomyopathy' resolved to
        # 'dilated cardiomyopathy 1G' before this fix). The corpus has
        # the authoritative id.
        mondo_id, name, slug = rec["id"], rec["name"], rec["slug"]
        if slug in existing and not args.force:
            results["skipped"].append(slug)
            continue
        t = time.time()
        try:
            run_disease(mondo_id, args.dist, do_summary=False, accept_first_run=True)
            dt = time.time() - t
            print(f"  [{i:>5d}/{len(records)}] {slug} OK ({dt:.1f}s)")
            results["ok"].append({"slug": slug, "seconds": round(dt, 1)})
        except Exception as e:
            print(f"  [{i:>5d}/{len(records)}] {slug} FAIL: {type(e).__name__}: {e}")
            results["fail"].append({"slug": slug, "error": str(e)[:200]})

    elapsed = (time.time() - t_start) / 60
    print(f"\nDone in {elapsed:.1f} min: ok={len(results['ok'])} "
          f"fail={len(results['fail'])} skipped={len(results['skipped'])}")
    if results["fail"]:
        print("\nFailures:")
        for f in results["fail"][:25]:
            print(f"  {f['slug']}: {f['error']}")
    with open("/tmp/disease_batch_run.json", "w") as f:
        json.dump(results, f, indent=2)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    b = sub.add_parser("build", help="Build the corpus from Mondo OBO + biobtree probes")
    b.add_argument("--out", default=None, help="output JSON path (default: build/mondo_corpus.json)")
    b.add_argument("--obo", default=None, help="path to mondo.obo (auto-download if absent)")
    b.add_argument("--refresh-obo", action="store_true")
    b.add_argument("--parallel", type=int, default=16)
    b.add_argument("--limit", type=int, default=0, help="probe first N terms only (testing)")
    b.set_defaults(func=cmd_build)

    f = sub.add_parser("filter", help="Filter corpus to a tier + emit [{name, slug}] list to stdout")
    f.add_argument("--corpus", default=None)
    f.add_argument("--tier", choices=_TIER_ORDER, default=None,
                   help="include this tier and all stricter tiers (T1 = T1 only; T3 = T1+T2+T3)")
    f.add_argument("--limit", type=int, default=0)
    f.set_defaults(func=cmd_filter)

    r = sub.add_parser("run", help="Dev-time serial batch builder (production: use Enju workflow)")
    r.add_argument("--dist", default="/data/sugi-atlas-dist")
    r.add_argument("--corpus", default=None)
    r.add_argument("--tier", choices=_TIER_ORDER, default=None)
    r.add_argument("--limit", type=int, default=0)
    r.add_argument("--force", action="store_true",
                   help="rebuild diseases already snapshotted")
    r.set_defaults(func=cmd_run)

    args = ap.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
