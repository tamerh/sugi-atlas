#!/usr/bin/env python3
"""Atlas one-entity end-to-end pipeline.

  collect §1..12  →  render  →  body_gate  →  summary (LLM)  →  summary_gate  →  publish

Writes to <dist>/atlas/<entity>/<symbol>/{page.md, bundle.json, summary.md, judge.json}.

CLI:
  python -m atlas.pipeline gene TP53 --dist /data/sugi-atlas-dist
  python -m atlas.pipeline gene TP53 --dist . --no-summary           # deterministic only
  python -m atlas.pipeline gene TP53 --dist . --first-run            # accept body_gate first_run
"""
import argparse, json, os, sys, time, urllib.request
from datetime import datetime, timezone

from atlas import __version__ as ATLAS_VERSION
from atlas.gene import collect as C
from atlas.gene import render as R
from atlas.disease import collect as DC
from atlas.disease import render as DR
from atlas.disease.anchors import resolve as resolve_disease
from atlas.disease.slug import slugify
from atlas.validation import body_gate, summary_gate
from atlas.bench import summary as B

DEFAULT_SUMMARY_MODEL = os.environ.get("ATLAS_SUMMARY_MODEL", "qwen/qwen3-235b-a22b-2507|Together")

def biobtree_version():
    try:
        d = json.loads(urllib.request.urlopen("http://127.0.0.1:8000/api/meta", timeout=5).read())
        return d.get("version") or "unknown"
    except Exception:
        return "unknown"

def collect_all(symbol):
    return {s: C.SECTIONS[s](symbol) for s in C.SECTIONS}

def render_all(bundle):
    return "\n\n".join(R.RENDER[s](bundle[s]) for s in R.RENDER)

def _yaml_escape(s):
    return str(s).replace('"', '\\"')

def assemble_page(symbol, summary_text, body_md, meta, bundle=None):
    """Hugo-frontmatter + declarative lead + (optional) AI summary + body.

    bundle: the full {section_id: bundle_dict} from collect_all. When passed,
    the deterministic declarative lead sentence and Updated-date line are
    prepended above the LLM summary. Required for the AI-friendly page shape;
    legacy callers that don't pass it get the prior (no-lead) layout."""
    fm = ["---"]
    for k in ("title", "symbol", "entity_type", "generated_at", "atlas_version", "biobtree_version"):
        fm.append(f'{k}: "{_yaml_escape(meta[k])}"')
    fm.append("---")
    fm.append("")
    head = "\n".join(fm)

    lead = ""
    cancer_overview = ""
    if bundle is not None:
        from atlas.page.declarative import declarative_sentence
        from atlas.page.jsonld import build_jsonld, as_script_tag
        from atlas.gene.render import r_cancer_overview
        sentence = declarative_sentence(bundle)
        # Pull the YYYY-MM-DD from the ISO `generated_at` for a human-visible
        # freshness signal (HTTP Last-Modified is set by Hugo from this same field).
        date = (meta.get("generated_at") or "")[:10]
        updated = f"*Updated: {date}*" if date else ""
        # schema.org Gene JSON-LD — federated-identity signal (sameAs to NCBI/
        # UniProt/Ensembl/HGNC/OMIM). Lives at the top of the body so AI
        # crawlers see it on the rendered page; also written as entity.jsonld
        # sidecar by the publish step for direct machine fetch.
        jsonld_tag = as_script_tag(build_jsonld(bundle))
        lead = jsonld_tag + "\n\n" + sentence + "\n\n" + (updated + "\n\n" if updated else "")
        # Cancer-significance overview block (CIViC paragraph + intOGen
        # driver flag) — placed between the page lead and the LLM Summary
        # so AI agents extract the curated narrative first. Elides for
        # non-cancer genes.
        co = r_cancer_overview(bundle)
        if co:
            cancer_overview = co + "\n\n"

    if summary_text:
        model = meta.get("summary_model", "Qwen3-235B")
        disclosure = (f"*Summary written by {model} from the deterministic data below. "
                      f"Facts in the tables that follow are the authoritative source.*")
        return (head + lead + cancer_overview + "## Summary\n\n" + disclosure + "\n\n"
                + summary_text.strip() + "\n\n" + body_md + "\n")
    return head + lead + cancer_overview + body_md + "\n"

def run_summary(body_md, symbol, model):
    key = B.api_key()
    prompt = B.INSTR.format(g=symbol) + "\n\nBODY:\n" + body_md
    d, dt = B.call(model, prompt, key, max_tokens=600)
    if "choices" not in d:
        raise RuntimeError(f"summary API error: {json.dumps(d)[:200]}")
    txt = str(d["choices"][0]["message"].get("content") or "").strip()
    return txt, dt

def run_gene(symbol, dist_dir, do_summary=True, summary_model=DEFAULT_SUMMARY_MODEL,
             accept_first_run=False, strict_summary=False):
    out_dir = os.path.join(dist_dir, "atlas", "gene", symbol)
    os.makedirs(out_dir, exist_ok=True)
    t0 = time.time()

    print(f"[1/5] collect §1..12 for {symbol}")
    bundle = collect_all(symbol)

    print(f"[2/5] render body")
    body_md = render_all(bundle)

    print(f"[3/5] body_gate")
    bg = body_gate.check(symbol, bundle)
    print(f"      verdict={bg['verdict']} ({bg['summary']})")
    if bg["verdict"] == "regression":
        print("      DIFF (first 15):")
        for d in bg["diff"][:15]:
            if d.get("kind") in {"added", "removed", "changed"}:
                print(f"        §{d['section']:<3} {d['key']:<26} {d['old']} → {d['new']} ({d['kind']})")
        sys.exit(2)
    if bg["verdict"] == "first_run" and not accept_first_run:
        print("      ! first_run — pass --first-run to publish anyway (snapshot still NOT created)")

    summary_text, judge_result = "", None
    if do_summary:
        print(f"[4/5] summary  ({summary_model})")
        summary_text, dt = run_summary(body_md, symbol, summary_model)
        print(f"      {len(summary_text)}c in {dt:.1f}s")
        print(f"[5/5] summary_gate")
        judge_result = summary_gate.check_summary(body_md, summary_text)
        print(f"      verdict={judge_result['verdict']}"
              f" both={len(judge_result.get('both', []))}"
              f" single={len(judge_result.get('single', []))}")
        if strict_summary and judge_result["verdict"] == "flagged":
            print("      STRICT mode: refusing to publish flagged summary")
            for c in judge_result["both"][:5]: print(f"        - {c}")
            sys.exit(3)
    else:
        print(f"[4/5] summary  SKIPPED")
        print(f"[5/5] summary_gate  SKIPPED")

    meta = {
        "title": symbol,
        "symbol": symbol,
        "entity_type": "gene",
        "generated_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "atlas_version": ATLAS_VERSION,
        "biobtree_version": biobtree_version(),
    }
    page_md = assemble_page(symbol, summary_text, body_md, meta)

    with open(os.path.join(out_dir, "page.md"), "w") as f: f.write(page_md)
    with open(os.path.join(out_dir, "bundle.json"), "w") as f:
        json.dump(bundle, f, indent=2, sort_keys=True)
    if summary_text:
        with open(os.path.join(out_dir, "summary.md"), "w") as f: f.write(summary_text + "\n")
    if judge_result is not None:
        with open(os.path.join(out_dir, "judge.json"), "w") as f:
            json.dump(judge_result, f, indent=2)

    print(f"\n✓ {symbol} done in {time.time()-t0:.1f}s -> {out_dir}")
    print(f"   page.md {len(page_md)}c  body_gate={bg['verdict']}"
          + (f"  summary_gate={judge_result['verdict']}" if judge_result else ""))
    return out_dir

def run_disease(name, dist_dir, do_summary=True, summary_model=DEFAULT_SUMMARY_MODEL,
                accept_first_run=False, strict_summary=False):
    """Disease-page pipeline. Mirrors run_gene; key changes:
      - input is a free-text name (or Mondo id); resolved once to anchors + slug
      - bundle key is the slug (filename-safe), not the symbol
      - collect_all/render_all come from atlas.disease, not atlas.gene
      - body_gate snapshot dir is <dist>/snapshots/disease/
    """
    t0 = time.time()

    print(f"[1/5] resolve anchors for {name!r}")
    a = resolve_disease(name)
    slug = slugify(a.canonical_name or name)
    print(f"      → {a.mondo_id} ({a.canonical_name}); cohort={len(a.cohort)}; slug={slug}")

    out_dir = os.path.join(dist_dir, "atlas", "disease", slug)
    os.makedirs(out_dir, exist_ok=True)

    print(f"[2/5] collect §1..14 for {slug}")
    bundle = {sid: DC.REGISTRY[sid].collect_fn(a) for sid in DC.REGISTRY}

    print(f"[3/5] render body")
    body_md = DR.render_all(bundle)

    print(f"[3a/5] body_gate")
    snap_dir = body_gate.snap_dir_for(dist_dir, "disease")
    bg = body_gate.check(slug, bundle, snap_dir=snap_dir)
    print(f"       verdict={bg['verdict']} ({bg['summary']})")
    if bg["verdict"] == "regression":
        print("       DIFF (first 15):")
        for d in bg["diff"][:15]:
            if d.get("kind") in {"added", "removed", "changed"}:
                print(f"         §{d['section']:<3} {d['key']:<28} {d['old']} → {d['new']} ({d['kind']})")
        sys.exit(2)
    if bg["verdict"] == "first_run" and not accept_first_run:
        print("       ! first_run — pass --first-run to publish (snapshot still NOT created)")

    summary_text, judge_result = "", None
    if do_summary:
        print(f"[4/5] summary  ({summary_model})")
        # Use the canonical_name as the entity label for the summary prompt.
        summary_text, dt = run_summary(body_md, a.canonical_name or name, summary_model)
        print(f"       {len(summary_text)}c in {dt:.1f}s")
        print(f"[5/5] summary_gate")
        judge_result = summary_gate.check_summary(body_md, summary_text)
        print(f"       verdict={judge_result['verdict']}"
              f" both={len(judge_result.get('both', []))}"
              f" single={len(judge_result.get('single', []))}")
        if strict_summary and judge_result["verdict"] == "flagged":
            print("       STRICT mode: refusing to publish flagged summary")
            for c in judge_result["both"][:5]: print(f"         - {c}")
            sys.exit(3)
    else:
        print(f"[4/5] summary  SKIPPED")
        print(f"[5/5] summary_gate  SKIPPED")

    meta = {
        "title": a.canonical_name or name,
        "symbol": slug,
        "entity_type": "disease",
        "generated_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "atlas_version": ATLAS_VERSION,
        "biobtree_version": biobtree_version(),
    }
    # No declarative-lead / cancer-overview shim yet — pass bundle=None so the
    # assembler skips the gene-specific lead. We can wire a disease-equivalent
    # lead in a follow-up once we settle on the headline-sentence shape.
    page_md = assemble_page(slug, summary_text, body_md, meta)

    with open(os.path.join(out_dir, "page.md"), "w") as f: f.write(page_md)
    with open(os.path.join(out_dir, "bundle.json"), "w") as f:
        json.dump(bundle, f, indent=2, sort_keys=True, default=str)
    if summary_text:
        with open(os.path.join(out_dir, "summary.md"), "w") as f: f.write(summary_text + "\n")
    if judge_result is not None:
        with open(os.path.join(out_dir, "judge.json"), "w") as f:
            json.dump(judge_result, f, indent=2)

    print(f"\n✓ {slug} done in {time.time()-t0:.1f}s -> {out_dir}")
    print(f"   page.md {len(page_md)}c  body_gate={bg['verdict']}"
          + (f"  summary_gate={judge_result['verdict']}" if judge_result else ""))
    return out_dir

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("entity", choices=["gene", "disease"])
    ap.add_argument("symbol", help="gene symbol OR disease name/Mondo id")
    ap.add_argument("--dist", required=True, help="dist repo root (writes <dist>/atlas/<entity>/<key>/)")
    ap.add_argument("--no-summary", action="store_true", help="skip the LLM summary + summary_gate")
    ap.add_argument("--summary-model", default=DEFAULT_SUMMARY_MODEL)
    ap.add_argument("--first-run", action="store_true", help="proceed despite no body_gate snapshot")
    ap.add_argument("--strict-summary", action="store_true", help="exit non-zero if summary_gate flags claims")
    args = ap.parse_args()
    if args.entity == "gene":
        run_gene(args.symbol, args.dist, do_summary=not args.no_summary,
                 summary_model=args.summary_model, accept_first_run=args.first_run,
                 strict_summary=args.strict_summary)
    else:
        run_disease(args.symbol, args.dist, do_summary=not args.no_summary,
                    summary_model=args.summary_model, accept_first_run=args.first_run,
                    strict_summary=args.strict_summary)

if __name__ == "__main__":
    main()
