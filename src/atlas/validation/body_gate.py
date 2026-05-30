#!/usr/bin/env python3
"""Body-gate — regression-test the collector bundle vs a saved per-gene snapshot.

Snapshots live in the PRIVATE dist repo (alongside the published pages), not in
the public atlas code repo. Default snap_dir is `$ATLAS_DIST_ROOT/snapshots/gene/`
(default ATLAS_DIST_ROOT=/data/sugi-atlas-dist). The workflow body_gate.py script
passes ctx["params"]["dist_root"] explicitly; ad-hoc CLI use picks up the env.

Verdicts:
  clean      — bundles are byte-identical
  drift      — minor numeric movements (typical when biobtree adds data)
  regression — a key emptied, or a count dropped by >50% — surface for review
  first_run  — no snapshot yet, no comparison possible

CLI:
  python -m atlas.validation.body_gate TP53                # check
  python -m atlas.validation.body_gate TP53 --update       # overwrite the snapshot
  python -m atlas.validation.body_gate TP53 --dist /path   # override dist root
"""
import argparse, json, os, sys

def default_dist_root():
    return os.environ.get("ATLAS_DIST_ROOT", "/data/sugi-atlas-dist")

def snap_dir_for(dist_root, entity="gene"):
    return os.path.join(dist_root, "snapshots", entity)

def _snap_path(snap_dir, symbol):
    return os.path.join(snap_dir, f"{symbol}.json")

def load_snapshot(snap_dir, symbol):
    p = _snap_path(snap_dir, symbol)
    return json.load(open(p)) if os.path.exists(p) else None

def save_snapshot(snap_dir, symbol, bundle):
    os.makedirs(snap_dir, exist_ok=True)
    with open(_snap_path(snap_dir, symbol), "w") as f:
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

def check(symbol, bundle, snap_dir=None):
    if snap_dir is None:
        snap_dir = snap_dir_for(default_dist_root())
    old = load_snapshot(snap_dir, symbol)
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
    ap.add_argument("--entity", default="gene", choices=["gene", "disease"],
                    help="entity type (selects which collector + snapshots/<entity>/ dir)")
    ap.add_argument("--update", action="store_true",
                    help="overwrite the snapshot with current bundle")
    ap.add_argument("--dist", default=default_dist_root(),
                    help="dist repo root (snapshots live at <dist>/snapshots/<entity>/)")
    args = ap.parse_args()
    snap_dir = snap_dir_for(args.dist, args.entity)
    if args.entity == "gene":
        from atlas.gene import collect as C
        bundle = {s: C.SECTIONS[s](args.symbol) for s in C.SECTIONS}
        key = args.symbol
    else:
        from atlas.disease import collect as DC
        from atlas.disease.anchors import resolve as resolve_disease
        # For disease the snapshot key is the slug (deterministic, filename-safe).
        # Resolve once to use as collect input AND derive slug.
        a = resolve_disease(args.symbol)
        bundle = {sid: DC.REGISTRY[sid].collect_fn(a) for sid in DC.REGISTRY}
        from atlas.disease.slug import slugify
        key = slugify(a.canonical_name or args.symbol)
    if args.update:
        save_snapshot(snap_dir, key, bundle)
        print(f"snapshot updated: {_snap_path(snap_dir, key)}")
        return
    r = check(key, bundle, snap_dir)
    print(f"{key}: {r['verdict']} ({r['summary']})")
    for d in r["diff"][:25]:
        if d.get("kind") in {"added", "removed", "changed"}:
            print(f"  §{d['section']:<3} {d['key']:<26} {d['old']} → {d['new']} ({d['kind']})")
    sys.exit(0 if r["verdict"] in {"clean", "drift", "first_run"} else 2)

if __name__ == "__main__":
    main()
