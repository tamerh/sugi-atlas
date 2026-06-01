#!/usr/bin/env python3
"""Atlas one-entity end-to-end pipeline.

  collect §1..12  →  render  →  body_gate  →  summary (LLM)  →  summary_gate  →  publish

Writes to <dist>/atlas/<entity>/<symbol>/{page.md, bundle.json, summary.md, judge.json}.

CLI:
  python -m atlas.pipeline gene TP53 --dist /data/sugi-atlas-dist
  python -m atlas.pipeline gene TP53 --dist . --no-summary           # deterministic only
  python -m atlas.pipeline gene TP53 --dist . --first-run            # accept body_gate first_run
"""
import argparse, json, os, re, sys, time, urllib.request
from datetime import datetime, timezone

from atlas import __version__ as ATLAS_VERSION
from atlas.gene import collect as C
from atlas.gene import render as R
from atlas.disease import collect as DC
from atlas.disease import render as DR
from atlas.disease.anchors import resolve as resolve_disease
from atlas.disease.slug import slugify
from atlas.drug import collect as DRC
from atlas.drug import render as DRR
from atlas.drug.anchors import resolve as resolve_drug
from atlas.drug.slug import slugify as drug_slugify
from atlas.validation import body_gate, summary_gate
from atlas.bench import summary as B
from atlas.biobtree import CALLS

DEFAULT_SUMMARY_MODEL = os.environ.get("ATLAS_SUMMARY_MODEL", "qwen/qwen3-235b-a22b-2507|Together")

def biobtree_version():
    """biobtree exposes no version/build field in /api/meta, so we stamp a
    dataset-count fingerprint instead — a real, verifiable signal that changes
    when biobtree's integrated data changes. (`generated_at` is the actual
    reproducibility anchor: replay against biobtree as-of that date.)"""
    try:
        d = json.loads(urllib.request.urlopen("http://127.0.0.1:8000/api/meta", timeout=5).read())
        if d.get("version"):
            return d["version"]
        n = len(d.get("datasets") or {})
        return f"{n} datasets" if n else "unknown"
    except Exception:
        return "unknown"

GENERATED_BY = "Sugi Atlas"  # attribution stamp; details on the /methods page


def datasets_union(registry):
    """Sorted union of every dataset each section DECLARES. Note: undercounts
    for disease/drug, which reach further datasets via cohort fan-out over gene
    collectors — prefer datasets_from_calls() for the true per-page list."""
    seen = set()
    for s in registry.values():
        seen.update(s.datasets or ())
    return sorted(seen)


_DS_TOK = re.compile(r"[a-z0-9_]+")


def datasets_from_calls(calls):
    """The TRUE per-page 'Data sources' list — distinct datasets actually
    queried during collection (the `s` param of search/entry + every chain
    token of map `m`). Captures datasets reached via cohort fan-out that the
    static registry doesn't declare (e.g. disease pages touch depmap / clingen /
    generif through the gene collectors)."""
    seen = set()
    for c in calls:
        p = c.get("params") or {}
        if p.get("s"):
            seen.add(str(p["s"]).lower())
        for tok in (p.get("m") or "").split(">>"):
            tok = tok.strip()
            if tok:
                m = _DS_TOK.match(tok)   # stops before '[' chain filters
                if m:
                    seen.add(m.group(0))
    seen.discard("")
    return sorted(seen)


def collect_all(symbol):
    return {s: C.SECTIONS[s](symbol) for s in C.SECTIONS}

def render_all(bundle):
    # Explicit section order. Expression (§11) is hoisted to right after
    # Transcripts (§2): "where is this expressed" is high-value context that
    # belongs near the top, not buried near the end. Functional genomics +
    # GeneRIFs were carved out of §3 (gene-level, not protein IDs) and render
    # right after it. Each renderer returns "" when it has no data.
    order = ["1", "2", "11", "3", "4", "5", "6", "7", "8", "9", "10", "12"]
    parts = []
    for s in order:
        parts.append(R.RENDER[s](bundle[s]))
        if s == "3":
            b3 = bundle["3"]
            parts.append(R.r_functional_genomics(b3))
            parts.append(R.r_generifs(b3))
    return "\n\n".join(p for p in parts if p)

def _yaml_escape(s):
    return str(s).replace('"', '\\"')

def assemble_page(symbol, summary_text, body_md, meta, bundle=None):
    """Hugo-frontmatter + declarative lead + (optional) AI summary + body.

    bundle: the full {section_id: bundle_dict} from collect_all. When passed,
    the deterministic declarative lead sentence and Updated-date line are
    prepended above the LLM summary. Required for the AI-friendly page shape;
    legacy callers that don't pass it get the prior (no-lead) layout."""
    fm = ["---"]
    for k in ("title", "symbol", "entity_type", "generated_at", "atlas_version",
              "biobtree_version", "generated_by"):
        if meta.get(k) is not None:
            fm.append(f'{k}: "{_yaml_escape(meta[k])}"')
    # datasets: YAML list → theme renders the visible "Data sources" block.
    # The api-call/chain trail stays an internal pipeline artifact (not
    # published) — transparency here is the source list + the generated_by
    # attribution (see the /methods page), not per-page API links.
    datasets = meta.get("datasets") or []
    if datasets:
        fm.append("datasets:")
        fm += [f'  - "{_yaml_escape(d)}"' for d in datasets]
    fm.append("---")
    fm.append("")
    head = "\n".join(fm)

    lead = ""
    cancer_overview = ""
    entity_type = (meta or {}).get("entity_type") or "gene"
    if bundle is not None:
        # No visible "Updated" line — the date lives in frontmatter
        # (generated_at) and the web theme surfaces it in the footer.
        if entity_type == "disease":
            from atlas.page.disease_declarative import declarative_sentence
            from atlas.page.disease_jsonld import build_jsonld, as_script_tag
            sentence = declarative_sentence(bundle)
            # `symbol` here is the disease slug (filename-safe key used by the
            # publish step). disease_jsonld needs it for the page @id URL.
            jsonld_tag = as_script_tag(build_jsonld(bundle, symbol))
            from atlas.page.disease_at_a_glance import at_a_glance
            glance = at_a_glance(bundle)
            if glance:
                sentence += "\n\n" + glance
        elif entity_type == "drug":
            from atlas.page.drug_declarative import declarative_sentence
            from atlas.page.drug_jsonld import build_jsonld, as_script_tag
            sentence = declarative_sentence(bundle)
            jsonld_tag = as_script_tag(build_jsonld(bundle, symbol))
            from atlas.page.drug_at_a_glance import at_a_glance
            glance = at_a_glance(bundle)
            if glance:
                sentence += "\n\n" + glance
        else:
            from atlas.page.declarative import declarative_sentence
            from atlas.page.jsonld import build_jsonld, as_script_tag
            from atlas.gene.render import r_cancer_overview
            sentence = declarative_sentence(bundle)
            jsonld_tag = as_script_tag(build_jsonld(bundle))
            # The opening is one intro region, not a stack of peer ## sections:
            #   lead sentence → RefSeq narrative → At a glance digest.
            # Neither the RefSeq summary nor the digest gets a "## " header (they
            # are part of the beginning, not data sections like Structure etc.).
            #
            # NCBI RefSeq curated summary as an intro paragraph. The trailing
            # "[provided by RefSeq, …]" provenance tag is stripped (our own
            # *Source* line carries attribution).
            b3 = bundle.get("3") or {}
            ncbi = (b3.get("ncbi_summary") or "").strip()
            if ncbi:
                import re as _re
                ncbi = _re.sub(r'\s*\[(?:provided|supplied) by[^\]]*\]\.?\s*$',
                               '', ncbi, flags=_re.I).strip()
                eid = b3.get("entrez_id")
                src = (f"[NCBI Gene {eid}](https://www.ncbi.nlm.nih.gov/gene/{eid})"
                       if eid else "NCBI Gene")
                sentence += f"\n\n{ncbi}\n\n*Source: {src} — RefSeq curated summary.*"
            # "At a glance" itemised digest — intro block (bold label, not a
            # ## section). Elides if nothing qualifies.
            from atlas.page.at_a_glance import at_a_glance
            glance = at_a_glance(bundle)
            if glance:
                sentence += "\n\n" + glance
            # Cancer-overview block is gene-specific (intOGen + CIViC per the
            # canonical gene). Disease pages have their own §4 somatic-driver
            # subblock, so we don't double up.
            co = r_cancer_overview(bundle)
            if co:
                cancer_overview = co + "\n\n"
        # schema.org JSON-LD — federated-identity signal (sameAs to ontology
        # cross-refs). Lives inline at the top of the body so AI crawlers see
        # it on the rendered page; also written as entity.jsonld sidecar by
        # the publish step for direct machine fetch.
        lead = jsonld_tag + "\n\n" + sentence + "\n\n"

    if summary_text:
        model = meta.get("summary_model", "Qwen3-235B")
        disclosure = (f"*Summary written by {model} from the deterministic data below. "
                      f"Facts in the tables that follow are the authoritative source.*")
        return (head + lead + cancer_overview + "## Summary\n\n" + disclosure + "\n\n"
                + summary_text.strip() + "\n\n" + body_md + "\n")
    return head + lead + cancer_overview + body_md + "\n"

def run_summary(body_md, symbol, model, kind="gene"):
    key = B.api_key()
    prompt = B.instruction(symbol, kind) + "\n\nBODY:\n" + body_md
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
    CALLS.clear()  # track datasets actually queried for this page

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
        "generated_by": GENERATED_BY,
        "datasets": datasets_from_calls(CALLS),
    }
    # Pass bundle so assemble_page emits the declarative lead + JSON-LD
    # inline script (parity with the Enju publish task).
    page_md = assemble_page(symbol, summary_text, body_md, meta, bundle=bundle)

    with open(os.path.join(out_dir, "page.md"), "w") as f: f.write(page_md)
    if summary_text:
        with open(os.path.join(out_dir, "summary.md"), "w") as f: f.write(summary_text + "\n")
    if judge_result is not None:
        with open(os.path.join(out_dir, "judge.json"), "w") as f:
            json.dump(judge_result, f, indent=2)

    # Sidecar — schema.org Gene JSON-LD identity card (the inline <script> in
    # the page mirrors it). bundle.json (raw data dump) and provenance.json
    # (api-call trail) are intentionally NOT published — transparency is the
    # frontmatter `datasets:` list + the `generated_by` attribution; the
    # call/chain trail stays an internal pipeline artifact.
    from atlas.page.jsonld import build_jsonld, as_jsonld_string
    with open(os.path.join(out_dir, "entity.jsonld"), "w") as f:
        f.write(as_jsonld_string(build_jsonld(bundle)))

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
    CALLS.clear()  # track datasets actually queried for this page

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
        summary_text, dt = run_summary(body_md, a.canonical_name or name, summary_model, kind="disease")
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
        "generated_by": GENERATED_BY,
        "datasets": datasets_from_calls(CALLS),
    }
    # Disease declarative lead + schema.org/MedicalCondition JSON-LD now
    # flow through assemble_page (entity_type='disease' branch). Same shape
    # as the gene page lead, just disease-shaped sentence + MedicalCondition
    # @type instead of Gene.
    page_md = assemble_page(slug, summary_text, body_md, meta, bundle=bundle)

    with open(os.path.join(out_dir, "page.md"), "w") as f: f.write(page_md)
    if summary_text:
        with open(os.path.join(out_dir, "summary.md"), "w") as f: f.write(summary_text + "\n")
    if judge_result is not None:
        with open(os.path.join(out_dir, "judge.json"), "w") as f:
            json.dump(judge_result, f, indent=2)

    # Sidecar — schema.org MedicalCondition card (bundle.json + provenance.json
    # intentionally not published; see run_gene note).
    from atlas.page.disease_jsonld import build_jsonld, as_jsonld_string
    with open(os.path.join(out_dir, "entity.jsonld"), "w") as f:
        f.write(as_jsonld_string(build_jsonld(bundle, slug)))

    print(f"\n✓ {slug} done in {time.time()-t0:.1f}s -> {out_dir}")
    print(f"   page.md {len(page_md)}c  body_gate={bg['verdict']}"
          + (f"  summary_gate={judge_result['verdict']}" if judge_result else ""))
    return out_dir

def run_drug(name, dist_dir, do_summary=True, summary_model=DEFAULT_SUMMARY_MODEL,
             accept_first_run=False, strict_summary=False):
    """Drug-page pipeline. Mirrors run_disease; key changes:
      - input is a drug name (or ChEMBL id); resolved once to DrugAnchors + slug
      - bundle key is the slug (filename-safe), not the name
      - collect_all/render_all come from atlas.drug
      - body_gate snapshot dir is <dist>/snapshots/drug/
    """
    t0 = time.time()
    CALLS.clear()  # track datasets actually queried for this page

    print(f"[1/5] resolve anchors for {name!r}")
    a = resolve_drug(name)
    slug = drug_slugify(a.canonical_name or name)
    print(f"      → {a.chembl_id} ({a.canonical_name}); type={a.molecule_type}; "
          f"targets={len(a.targets)}; slug={slug}")

    out_dir = os.path.join(dist_dir, "atlas", "drug", slug)
    os.makedirs(out_dir, exist_ok=True)

    print(f"[2/5] collect §1..12 for {slug}")
    bundle = {sid: DRC.REGISTRY[sid].collect_fn(a) for sid in DRC.REGISTRY}

    print(f"[3/5] render body")
    body_md = DRR.render_all(bundle)

    print(f"[3a/5] body_gate")
    snap_dir = body_gate.snap_dir_for(dist_dir, "drug")
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
        summary_text, dt = run_summary(body_md, a.canonical_name or name, summary_model, kind="drug")
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
        "entity_type": "drug",
        "generated_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "atlas_version": ATLAS_VERSION,
        "biobtree_version": biobtree_version(),
        "generated_by": GENERATED_BY,
        "datasets": datasets_from_calls(CALLS),
    }
    page_md = assemble_page(slug, summary_text, body_md, meta, bundle=bundle)

    with open(os.path.join(out_dir, "page.md"), "w") as f: f.write(page_md)
    if summary_text:
        with open(os.path.join(out_dir, "summary.md"), "w") as f: f.write(summary_text + "\n")
    if judge_result is not None:
        with open(os.path.join(out_dir, "judge.json"), "w") as f:
            json.dump(judge_result, f, indent=2)

    # Sidecar — schema.org Drug card (bundle.json + provenance.json
    # intentionally not published; see run_gene note).
    from atlas.page.drug_jsonld import build_jsonld, as_jsonld_string
    with open(os.path.join(out_dir, "entity.jsonld"), "w") as f:
        f.write(as_jsonld_string(build_jsonld(bundle, slug)))

    print(f"\n✓ {slug} done in {time.time()-t0:.1f}s -> {out_dir}")
    print(f"   page.md {len(page_md)}c  body_gate={bg['verdict']}"
          + (f"  summary_gate={judge_result['verdict']}" if judge_result else ""))
    return out_dir

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("entity", choices=["gene", "disease", "drug"])
    ap.add_argument("symbol", help="gene symbol OR disease name/Mondo id OR drug name/ChEMBL id")
    ap.add_argument("--dist", required=True, help="dist repo root (writes <dist>/atlas/<entity>/<key>/)")
    # Dev stage: LLM summary is OFF by default for all entities; opt in with
    # --summary. (--no-summary kept as a no-op alias so old invocations don't break.)
    ap.add_argument("--summary", action="store_true", help="run the LLM summary + summary_gate (off by default)")
    ap.add_argument("--no-summary", action="store_true", help=argparse.SUPPRESS)
    ap.add_argument("--summary-model", default=DEFAULT_SUMMARY_MODEL)
    ap.add_argument("--first-run", action="store_true", help="proceed despite no body_gate snapshot")
    ap.add_argument("--strict-summary", action="store_true", help="exit non-zero if summary_gate flags claims")
    args = ap.parse_args()
    runner = {"gene": run_gene, "disease": run_disease, "drug": run_drug}[args.entity]
    runner(args.symbol, args.dist, do_summary=args.summary,
           summary_model=args.summary_model, accept_first_run=args.first_run,
           strict_summary=args.strict_summary)

if __name__ == "__main__":
    main()
