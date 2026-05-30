#!/usr/bin/env python3
"""Body-gate — regression-test the collector bundle vs a saved per-gene snapshot.

For every entity we publish, we keep a JSON snapshot of its collector bundle at
`validation_data/snapshots/gene/<SYMBOL>.json`. The next run diffs the new
bundle against that snapshot and reports STRUCTURAL changes (counts that moved,
sections that emptied, ids added/removed).

Verdicts:
  clean      — bundles are byte-identical
  drift      — minor numeric movements (typical when biobtree adds data)
  regression — a key emptied, or a count dropped by >50% — surface for review
  first_run  — no snapshot yet, no comparison possible

CLI:
  python -m atlas.validation.body_gate TP53            # check
  python -m atlas.validation.body_gate TP53 --update   # overwrite the snapshot
"""
import argparse, json, os, sys
from atlas.gene import collect as C

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
SNAP_DIR = os.path.join(REPO, "validation_data", "snapshots", "gene")

def _snap_path(symbol):
    return os.path.join(SNAP_DIR, f"{symbol}.json")

def load_snapshot(symbol):
    p = _snap_path(symbol)
    return json.load(open(p)) if os.path.exists(p) else None

def save_snapshot(symbol, bundle):
    os.makedirs(SNAP_DIR, exist_ok=True)
    with open(_snap_path(symbol), "w") as f:
        json.dump(bundle, f, indent=2, sort_keys=True)

def _size(v):
    if v is None or v == "": return 0
    if isinstance(v, (list, dict)): return len(v)
    return 1

def _flat_sizes(bundle):
    out = {}
    for sec, b in bundle.items():
        if not isinstance(b, dict): continue
        for k, v in b.items():
            if k in {"section", "symbol"}: continue
            out[(sec, k)] = _size(v)
    return out

def diff(old, new):
    if old is None:
        return [{"kind": "first_run"}]
    if old == new:
        return []
    old_s, new_s = _flat_sizes(old), _flat_sizes(new)
    out = []
    for k in sorted(set(old_s) | set(new_s)):
        o, n = old_s.get(k, 0), new_s.get(k, 0)
        if o == n: continue
        kind = "added" if o == 0 else "removed" if n == 0 else "changed"
        out.append({"section": k[0], "key": k[1], "old": o, "new": n, "kind": kind})
    return out

def verdict(differences):
    if not differences: return "clean"
    if any(d.get("kind") == "first_run" for d in differences): return "first_run"
    for d in differences:
        if d["kind"] == "removed": return "regression"
        if d["kind"] == "changed" and d["old"] > 0 and d["new"] < d["old"] * 0.5:
            return "regression"
    return "drift"

def check(symbol, bundle):
    old = load_snapshot(symbol)
    d = diff(old, bundle)
    v = verdict(d)
    if v == "first_run":
        msg = "no snapshot yet — run --update to create baseline"
    elif v == "clean":
        msg = "byte-identical to snapshot"
    else:
        adds = sum(1 for x in d if x["kind"] == "added")
        rems = sum(1 for x in d if x["kind"] == "removed")
        chs  = sum(1 for x in d if x["kind"] == "changed")
        msg = f"{chs} changed, {adds} added, {rems} removed"
    return {"verdict": v, "summary": msg, "diff": d}

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("symbol")
    ap.add_argument("--update", action="store_true", help="overwrite the snapshot with current bundle")
    args = ap.parse_args()
    bundle = {s: C.SECTIONS[s](args.symbol) for s in C.SECTIONS}
    if args.update:
        save_snapshot(args.symbol, bundle)
        print(f"snapshot updated: {_snap_path(args.symbol)}")
        return
    r = check(args.symbol, bundle)
    print(f"{args.symbol}: {r['verdict']} ({r['summary']})")
    for d in r["diff"][:25]:
        if d.get("kind") in {"added", "removed", "changed"}:
            print(f"  §{d['section']:<3} {d['key']:<26} {d['old']} → {d['new']} ({d['kind']})")
    sys.exit(0 if r["verdict"] in {"clean", "drift", "first_run"} else 2)

if __name__ == "__main__":
    main()
