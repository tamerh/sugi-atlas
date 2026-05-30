#!/usr/bin/env python3
"""Dev-time backlog driver for disease pages.

Two modes:
- Default: walk `biobtree-content/biobtree/disease/entities.yaml` (the parked
  curated 61-disease list).
- `--all`:  walk every Mondo entry with substantial cohort signal (gwas or
  civic_evidence or clinvar or gencc xref counts > 0). Used for the
  full-Mondo scale-out — same engine, different input list.

For each disease:
  resolve → collect §1-§14 → render → body_gate (first_run accepted) →
  publish to <dist>/atlas/disease/<slug>/. No LLM summary by default.

Skips diseases that already have a snapshot in <dist>/snapshots/disease/
(re-run with `--force` to rebuild).

This is the dev driver. Production runs should use the Enju workflow at
`src/atlas/disease/enju.yaml`, which gives per-disease retries, parallel
DAGs, and a full audit trail.

  python -m bin.run_disease_backlog               # parked 61
  python -m bin.run_disease_backlog --limit 5     # first 5 only
  python -m bin.run_disease_backlog --force       # rebuild already-snapshotted
  python -m bin.run_disease_backlog --all         # full Mondo enumeration

Logs are appended to /tmp/disease_backlog_run.json.
"""
import argparse, json, os, re, sys, time, traceback

# Make atlas importable even when called as a plain script.
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

DEFAULT_DIST = "/data/sugi-atlas-dist"
ENT_YAML = "/data/biobtree-content/biobtree/disease/entities.yaml"


def load_curated():
    """Parse the parked entities.yaml (`- {symbol: "...", status: ...}` lines)."""
    diseases = []
    for line in open(ENT_YAML):
        m = re.match(r'^- \{symbol: "([^"]+)"', line)
        if m:
            diseases.append(m.group(1))
    return diseases


def load_full_mondo():
    """Enumerate every Mondo node with material cohort signal. NOT IMPLEMENTED
    yet — biobtree doesn't expose a 'list all mondo ids' endpoint; we'd need
    to either page through search('mondo', s='mondo') with no input filter or
    drive off an external Mondo dump (mondo.obo). Placeholder for the
    all-diseases scale-out."""
    raise NotImplementedError(
        "Full-Mondo enumeration not yet wired. Needs either:\n"
        "  - a biobtree 'list ids in dataset N' endpoint (file as upstream req)\n"
        "  - external Mondo .obo download + parse for MONDO:NNN ids\n"
        "  - filter to nodes with >=1 of {gwas, civic_evidence, clinvar, gencc}\n"
        "    xref counts to skip leaf-most ontology nodes with no Atlas signal.")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dist", default=DEFAULT_DIST,
                    help="dist repo root (default %(default)s)")
    ap.add_argument("--limit", type=int, default=0,
                    help="stop after N successful builds (0 = no limit)")
    ap.add_argument("--force", action="store_true",
                    help="rebuild diseases that already have a snapshot")
    ap.add_argument("--all", action="store_true",
                    help="full Mondo enumeration (not yet implemented)")
    args = ap.parse_args()

    if args.all:
        diseases = load_full_mondo()
    else:
        diseases = load_curated()
    print(f"Loaded {len(diseases)} diseases")

    from atlas.pipeline import run_disease
    from atlas.disease.anchors import resolve as resolve_disease
    from atlas.disease.slug import slugify

    existing = set()
    if not args.force:
        snap_dir = os.path.join(args.dist, "snapshots", "disease")
        if os.path.isdir(snap_dir):
            existing = {f[:-5] for f in os.listdir(snap_dir) if f.endswith(".json")}
        print(f"Already-snapshotted: {len(existing)}")

    results = {"ok": [], "fail": [], "skipped": []}
    t_start = time.time()
    ok_count = 0
    for i, name in enumerate(diseases, 1):
        try:
            a = resolve_disease(name)
            slug = slugify(a.canonical_name or name)
        except Exception as e:
            print(f"  [{i:>2}/{len(diseases)}] {name!r}: RESOLVE FAILED — {e}")
            results["fail"].append({"name": name, "stage": "resolve", "error": str(e)})
            continue
        if slug in existing:
            print(f"  [{i:>2}/{len(diseases)}] {name!r} → {slug}: SKIP")
            results["skipped"].append(slug)
            continue

        t = time.time()
        try:
            run_disease(name, args.dist, do_summary=False, accept_first_run=True)
            results["ok"].append({"name": name, "slug": slug,
                                  "seconds": round(time.time() - t, 1)})
            ok_count += 1
            if args.limit and ok_count >= args.limit:
                print(f"  -- limit {args.limit} reached, stopping --")
                break
        except SystemExit as e:
            print(f"  [{i:>2}/{len(diseases)}] {name!r} → {slug}: SYS_EXIT {e}")
            results["fail"].append({"name": name, "slug": slug, "stage": "pipeline",
                                    "error": f"sys.exit {e}"})
        except Exception as e:
            print(f"  [{i:>2}/{len(diseases)}] {name!r} → {slug}: "
                  f"PIPELINE FAILED — {type(e).__name__}: {e}")
            traceback.print_exc(limit=2)
            results["fail"].append({"name": name, "slug": slug, "stage": "pipeline",
                                    "error": f"{type(e).__name__}: {e}"})

    elapsed = time.time() - t_start
    print(f"\n--- DONE in {elapsed/60:.1f}min ---")
    print(f"ok={len(results['ok'])}  fail={len(results['fail'])}  "
          f"skipped={len(results['skipped'])}")
    for f in results["fail"]:
        print(f"  FAIL — {f['name']!r}: [{f['stage']}] {f['error']}")

    with open("/tmp/disease_backlog_run.json", "w") as f:
        json.dump(results, f, indent=2)
    print(f"Run log: /tmp/disease_backlog_run.json")


if __name__ == "__main__":
    main()
