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
from atlas.atomicio import write_text, write_json
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

_BIOBTREE_META = None  # cached per process — /ws/meta appparams is stable per run


def _biobtree_meta():
    """The `appparams` block from biobtree's /ws/meta (biobtree_version,
    biobtree_commit, biobtree_build_date), fetched once per process. Empty dict
    if the endpoint or block is unavailable (older biobtree)."""
    global _BIOBTREE_META
    if _BIOBTREE_META is None:
        from atlas.biobtree import API
        try:
            d = json.loads(urllib.request.urlopen(f"{API}/ws/meta", timeout=5).read())
            _BIOBTREE_META = d.get("appparams") or {}
        except Exception:
            _BIOBTREE_META = {}
    return _BIOBTREE_META


def biobtree_version():
    """biobtree's own version (appparams.biobtree_version, e.g. 'v2.0.0') — the
    data-snapshot anchor for reproducibility. 'unknown' if unavailable."""
    return _biobtree_meta().get("biobtree_version") or "unknown"


def biobtree_commit():
    """biobtree's build commit (appparams.biobtree_commit), pinning the exact data
    snapshot alongside biobtree_version. 'unknown' if unavailable."""
    return _biobtree_meta().get("biobtree_commit") or "unknown"


_ATLAS_VERSION = None
_ATLAS_COMMIT = None


def _git(*args):
    """Stripped stdout of a git command in the repo, or None on any failure."""
    try:
        import subprocess
        out = subprocess.check_output(["git", *args], cwd=os.path.dirname(__file__),
                                      stderr=subprocess.DEVNULL).decode().strip()
        return out or None
    except Exception:
        return None


def atlas_version():
    """Pipeline version stamped on every page. Prefers the build-time stamp set
    by atlas.sh (ATLAS_BUILD_VERSION = `git describe --tags --always`, captured at
    build START so the corpus is labelled with the code that actually built it);
    falls back to a live git describe, then the packaged __version__. A build off
    a tag reads a clean 'vX.Y.Z'; between tags, 'vX.Y.Z-N-gSHA'. Cached per
    process (every batch worker inherits ATLAS_BUILD_VERSION, so no git per page)."""
    global _ATLAS_VERSION
    if _ATLAS_VERSION is None:
        _ATLAS_VERSION = (os.environ.get("ATLAS_BUILD_VERSION")
                          or _git("describe", "--tags", "--always")
                          or ATLAS_VERSION)
    return _ATLAS_VERSION


def atlas_commit():
    """Short git SHA of the building code (ATLAS_BUILD_COMMIT from atlas.sh,
    captured at build start), else a live lookup, else 'nogit'. Cached."""
    global _ATLAS_COMMIT
    if _ATLAS_COMMIT is None:
        _ATLAS_COMMIT = (os.environ.get("ATLAS_BUILD_COMMIT")
                         or _git("rev-parse", "--short", "HEAD")
                         or "nogit")
    return _ATLAS_COMMIT


GENERATED_BY = "Sugi Atlas"  # attribution stamp; details on the /methods page

# Page filename. "index.md" makes each <entity>/<slug>/ a Hugo leaf page-bundle,
# so biobtree-content can mount the dist as a module directly (no sync wrapper).
# The web team's mount must expect this layout.
PAGE_FILENAME = "index.md"


def build_meta(entity_type, slug, title, datasets, generated_at=None, bundle=None):
    """The single page-frontmatter meta builder — used by run_gene/disease/drug
    AND the batch driver, so the shape can't drift between paths (the m4 fix).
    `slug` is the URL/filename key (gene slug == symbol); `title` is the human
    label (gene=symbol, disease=canonical name, drug=ChEMBL name). When `bundle`
    is passed, the P2/P3 key-facts (identifier, alt_names, tldr, section_defaults)
    are derived and merged in."""
    # Title-case a SHOUTING all-caps title for display (audit #12: drug names
    # like 'IMATINIB'). Genes are exempt — symbols (TP53) are upper by
    # convention; diseases in sentence case get their leading letter capitalized
    # to match the declarative lead. The original ChEMBL name survives in the
    # JSON-LD alternateName.
    if entity_type != "gene" and title:
        from atlas.render_common import display_name
        if title.isupper():
            title = display_name(title)
        elif entity_type == "disease" and title == title.lower():
            title = title[0].upper() + title[1:]
    meta = {
        "title": title,
        "symbol": slug,             # URL slug (legacy); templates key on `identifier`
        "entity_type": entity_type,
        "generated_at": (generated_at
                         or datetime.now(timezone.utc).isoformat(timespec="seconds")),
        "atlas_version": atlas_version(),
        "atlas_commit": atlas_commit(),
        "biobtree_version": biobtree_version(),
        "biobtree_commit": biobtree_commit(),
        "generated_by": GENERATED_BY,
        "datasets": datasets,
    }
    if bundle is not None:
        from atlas.page.meta_facts import entity_facts
        meta.update(entity_facts(entity_type, bundle))
        # Prominence: raw component counts (always) + the 0-100 percentile score
        # when a batch distribution has been loaded (batch render / single-entity
        # rebuild against the last batch). Drives search ranking + featured cards.
        from atlas.page import evidence
        comps = evidence.components(entity_type, bundle)
        meta["evidence_components"] = comps
        score = evidence.lookup(entity_type, slug,
                                  evidence.raw_signal(entity_type, comps))
        if score is not None:
            meta["evidence_score"] = score
    return meta


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
    bundle = {s: C.SECTIONS[s](symbol) for s in C.SECTIONS}
    _scrub_noncoding(bundle)
    return bundle


# Variant/disease/trial fields that biobtree links POSITIONALLY (by genomic
# overlap), so a non-coding gene inherits the overlapping protein-coding
# neighbor's biology (TTN-AS1 → TTN's 6,528 ClinVar variants + disease set).
# Confidently-wrong, worse than empty — so we clear them for non-coding genes
# and the lead / At a glance / mesh / render all reflect "no protein, no
# variant/disease data". (Curated HGNC-keyed links — GenCC/ClinGen/CIViC — are
# cleared defensively too; they're empty for ncRNA anyway.)
_NONCODING_SCRUB = {
    "6":  ("clinvar_total", "clinvar_breakdown", "top_pathogenic", "spliceai_total",
           "top_spliceai", "alphamissense_total", "top_alphamissense", "dbsnp_sample"),
    # disease_trials are positional; civic_evidence_total / molecule_count /
    # is_drug_target also feed the evidence _SPEC + meta_facts, so clear them too
    # (audit #14: non-coding genes carried inherited drug/civic counts in the
    # frontmatter + evidence_score even though the body correctly says "no data").
    "10": ("disease_trials", "disease_trial_count",
           "civic_evidence_total", "molecule_count", "is_drug_target"),
    "12": ("gene_omim", "disease_omim", "gencc", "clingen_validity", "mondo",
           "orphanet", "hpo", "hpo_total", "gwas", "gwas_total", "gwas_studies",
           "efo_traits", "mesh_descriptors", "civic", "intogen"),
}


def _scrub_noncoding(bundle):
    """Clear positionally-inherited variant/disease/trial blocks for non-coding
    genes (biotype != protein_coding). Sets bundle['_noncoding'] = biotype so
    the renderer/At-a-glance can state it affirmatively."""
    b1 = bundle.get("1") or {}
    biotype = ((b1.get("ensembl") or {}).get("biotype") or "").strip()
    if not biotype or biotype == "protein_coding":
        return
    for sid, keys in _NONCODING_SCRUB.items():
        sec = bundle.get(sid)
        if sec:
            for k in keys:
                sec.pop(k, None)
    bundle["_noncoding"] = biotype

from atlas.render_common import demote as _demote, emit_canonical, with_heading_id


def render_all(bundle):
    """Gene page body in the FROZEN canonical H2 order (docs/PAGE_CONTRACT.md):
    Identifiers → Gene structure → Protein → Function → Disease & clinical →
    Drugs & pharmacology. (Summary is wrapped by assemble_page; Related appended
    after.) Section H2s are demoted to H3 and carry stable backend-owned ids."""
    b3 = bundle.get("3") or {}
    noncoding = bundle.get("_noncoding")

    def S(s, anchor):           # registered section, demoted, stable H3 id
        return with_heading_id(_demote(R.RENDER[s](bundle[s])), anchor)

    def D(md, anchor):          # derived renderer, demoted, stable H3 id
        return with_heading_id(_demote(md), anchor)

    def join(*parts):
        return "\n\n".join(p for p in parts if p and p.strip())

    canon = b3.get("canonical_uniprot")
    protein_a = (f'<a id="protein-{canon}"></a>\n\n'
                 if (canon and not noncoding) else "")

    # Non-coding RNA layer (§14) — curated, symbol-keyed (not positional), so it
    # survives the non-coding scrub and is the real content for thin lncRNA/miRNA
    # pages. Gated on data, so coding genes (which rarely carry these) elide them.
    b14 = bundle.get("14") or {}
    nc_function = D(R.r_ncrna_function(b14), "ncrna-function")
    nc_disease = D(R.r_ncrna_disease(b14), "ncrna-disease")
    nc_interactions = D(R.r_ncrna_interactions(b14), "ncrna-interactions")
    nc_drugs = D(R.r_ncrna_drugs(b14), "ncrna-drugs")

    spec = [
        ("Identifiers", "identifiers", S("1", "gene-ids"), None),
        ("Gene structure", "gene-structure",
         join(S("2", "transcripts"), S("11", "expression"),
              D(R.r_hpa_expression(bundle), "hpa-expression"), S("9", "regulation"),
              D(R.r_functional_genomics(b3), "functional-genomics"),
              D(R.r_generifs(b3), "generif"), S("5", "orthologs")), None),
        ("Protein", "protein",
         "" if noncoding else join(S("3", "protein-ids"), S("4", "structure"),
                                   D(R.r_hpa_protein(bundle), "hpa-protein"),
                                   D(R.r_residue_map(b3), "residue-map")),
         "Non-coding RNA — no protein product."),
        ("Function", "function",
         (join(D(R.r_noncoding_genesets(bundle.get("7") or {}), "gene-sets"),
               nc_function, nc_interactions) if noncoding
          else join(S("7", "pathways"), S("8", "interactions"),
                    nc_function, nc_interactions)),
         "No curated pathway, Gene-Ontology, or interaction data."),
        ("Disease & clinical", "disease",
         (nc_disease if noncoding
          else join(D(R.r_cancer_overview(bundle), "cancer"),
                    D(R.r_hpa_cancer(bundle), "hpa-cancer"),
                    S("6", "variants"), S("12", "disease-assoc"), nc_disease)),
         "No curated disease, variant, or cancer-driver associations."),
        ("Drugs & pharmacology", "drugs",
         (nc_drugs if noncoding else join(S("10", "drug-data"), nc_drugs)),
         "No protein-based drug or pharmacology data — Atlas drug coverage is "
         "protein-target-based. Non-coding RNAs may still be targetable by RNA "
         "therapeutics (antisense oligonucleotides, siRNA), which aren't covered here."),
    ]
    return emit_canonical(spec, anchors={"protein": protein_a})

_CTRL = re.compile(r"[\x00-\x08\x0b\x0c\x0e-\x1f\x7f]")  # control chars (keep \t,\n)


def _yaml_escape(s):
    # Make a value safe inside a double-quoted YAML scalar. Strip control chars
    # (e.g. a stray 0x7f from biobtree text), then escape backslash BEFORE quote
    # — a value with a literal/ trailing backslash (e.g. a truncated Mondo
    # synonym `…(formerly cutaneous \"`) would otherwise escape the closing quote
    # and break the whole frontmatter block.
    return _CTRL.sub("", str(s)).replace("\\", "\\\\").replace('"', '\\"')

def assemble_page(symbol, summary_text, body_md, meta, bundle=None):
    """Hugo-frontmatter + declarative lead + (optional) AI summary + body.

    bundle: the full {section_id: bundle_dict} from collect_all. When passed,
    the deterministic declarative lead sentence and Updated-date line are
    prepended above the LLM summary. Required for the AI-friendly page shape;
    legacy callers that don't pass it get the prior (no-lead) layout."""
    fm = ["---"]
    for k in ("title", "identifier", "symbol", "entity_type", "description",
              "generated_at", "atlas_version", "atlas_commit", "biobtree_version",
              "biobtree_commit", "generated_by"):
        if meta.get(k):
            fm.append(f'{k}: "{_yaml_escape(meta[k])}"')
    # Search aliases (P2) — NOT Hugo-reserved `aliases:` (that emits 301s).
    for field in ("alt_names", "tldr"):
        items = [x for x in (meta.get(field) or []) if x]
        if items:
            fm.append(f"{field}:")
            fm += [f'  - "{_yaml_escape(x)}"' for x in items]
    # Prominence (search ranking): a 0-100 percentile within entity type, plus
    # the raw component counts (stable keys per type) for client-side re-weight.
    if meta.get("evidence_score") is not None:
        fm.append(f"evidence_score: {int(meta['evidence_score'])}")
    pc = meta.get("evidence_components") or {}
    if pc:
        fm.append("evidence_components:")
        fm += [f"  {k}: {int(v)}" for k, v in pc.items()]
    # Section open/collapsed hints (P3), keyed by canonical anchor id.
    sd = meta.get("section_defaults") or {}
    if sd:
        fm.append("section_defaults:")
        fm += [f"  {k}: {v}" for k, v in sd.items()]
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
    related_tail = ""   # "## Related Atlas pages" section, appended at page end
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
            # (The medical disclaimer lives in the web theme's existing
            # `entity-footer-about` footer — no body-level block here.)
        elif entity_type == "drug":
            from atlas.page.drug_declarative import declarative_sentence
            from atlas.page.drug_jsonld import build_jsonld, as_script_tag
            sentence = declarative_sentence(bundle)
            jsonld_tag = as_script_tag(build_jsonld(bundle, symbol))
            from atlas.page.drug_at_a_glance import at_a_glance
            glance = at_a_glance(bundle)
            if glance:
                sentence += "\n\n" + glance
        elif entity_type == "pathway":
            # Flat (non-section) bundle; the lead + GO + source ARE the summary.
            from atlas.pathway.render import summary_block
            from atlas.page.pathway_jsonld import build_jsonld, as_script_tag
            sentence = summary_block(bundle)
            jsonld_tag = as_script_tag(build_jsonld(bundle, symbol))
        else:
            from atlas.page.declarative import declarative_sentence
            from atlas.page.jsonld import build_jsonld, as_script_tag
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
            # Cancer-significance (intOGen + CIViC) is no longer a separate
            # pre-body block — it folds into the canonical "Disease & clinical"
            # (#disease) section, emitted by render_all (see PAGE_CONTRACT.md).
        # "Related Atlas pages" — cross-entity navigation, rendered as the final
        # section ("see also" convention; keeps the intro focused on the entity).
        # Machine traversal is handled by the JSON-LD edges in <head>; this is
        # the human-facing block. Elides when nothing is built yet.
        from atlas.page import links
        rel = links.related_block(entity_type, bundle, slug=symbol)
        if rel:
            related_tail = "\n\n" + rel
        # schema.org JSON-LD — federated-identity signal (sameAs to ontology
        # cross-refs). The readable declarative sentence (+ digest) leads, then
        # the inline <script> block (audit #6: it used to precede the lead and
        # bury it under hundreds of lines of JSON; the tag is invisible in
        # rendered HTML, so placement below the lead costs crawlers nothing and
        # restores the one sentence an AI agent should extract first). The
        # complete graph is also written as the entity.jsonld sidecar. The whole
        # intro is one canonical "## Summary {#summary}" section (PAGE_CONTRACT):
        # lead sentence (most-indexable) → RefSeq → At-a-glance → inline JSON-LD.
        lead = "## Summary {#summary}\n\n" + sentence + "\n\n" + jsonld_tag + "\n\n"

    if summary_text:    # legacy LLM-summary path (unused — no LLM summaries)
        model = meta.get("summary_model", "Qwen3-235B")
        disclosure = (f"*Summary written by {model} from the deterministic data below. "
                      f"Facts in the tables that follow are the authoritative source.*")
        return (head + lead + disclosure + "\n\n"
                + summary_text.strip() + "\n\n" + body_md + related_tail + "\n")
    return head + lead + body_md + related_tail + "\n"

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

    # Cross-entity link mesh: load the manifest + register this gene (slug ==
    # symbol), so renderers/JSON-LD can link to already-built Atlas pages.
    from atlas.page import links
    links.load(dist_dir)
    links.upsert(dist_dir, "gene", symbol,
                 id_keys=[symbol, (bundle.get("1") or {}).get("hgnc_id")],
                 canonical=symbol)

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

    meta = build_meta("gene", symbol, symbol, datasets_from_calls(CALLS), bundle=bundle)
    # Pass bundle so assemble_page emits the declarative lead + JSON-LD
    # inline script (parity with the Enju publish task).
    page_md = assemble_page(symbol, summary_text, body_md, meta, bundle=bundle)

    write_text(os.path.join(out_dir, PAGE_FILENAME), page_md)
    if summary_text:
        write_text(os.path.join(out_dir, "summary.md"), summary_text + "\n")
    if judge_result is not None:
        write_json(os.path.join(out_dir, "judge.json"), judge_result, indent=2)

    # Sidecar — schema.org Gene JSON-LD identity card (the inline <script> in
    # the page mirrors it). bundle.json (raw data dump) and provenance.json
    # (api-call trail) are intentionally NOT published — transparency is the
    # frontmatter `datasets:` list + the `generated_by` attribution; the
    # call/chain trail stays an internal pipeline artifact.
    from atlas.page.jsonld import build_jsonld, as_jsonld_string
    write_text(os.path.join(out_dir, "entity.jsonld"), as_jsonld_string(build_jsonld(bundle)))

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

    # Cross-entity link mesh: register this disease under its IDs + name/synonyms.
    from atlas.page import links
    links.load(dist_dir)
    links.upsert(dist_dir, "disease", slug,
                 id_keys=[a.mondo_id, a.efo_id],
                 name_keys=[a.canonical_name, *(a.synonyms or ())],
                 canonical=a.canonical_name)

    print(f"[3/5] render body")
    bundle["_indicated_drugs"] = links.indicated_drugs(dist_dir, slug)
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

    meta = build_meta("disease", slug, a.canonical_name or name, datasets_from_calls(CALLS), bundle=bundle)
    # Disease declarative lead + schema.org/MedicalCondition JSON-LD now
    # flow through assemble_page (entity_type='disease' branch). Same shape
    # as the gene page lead, just disease-shaped sentence + MedicalCondition
    # @type instead of Gene.
    page_md = assemble_page(slug, summary_text, body_md, meta, bundle=bundle)

    write_text(os.path.join(out_dir, PAGE_FILENAME), page_md)
    if summary_text:
        write_text(os.path.join(out_dir, "summary.md"), summary_text + "\n")
    if judge_result is not None:
        write_json(os.path.join(out_dir, "judge.json"), judge_result, indent=2)

    # Sidecar — schema.org MedicalCondition card (bundle.json + provenance.json
    # intentionally not published; see run_gene note).
    from atlas.page.disease_jsonld import build_jsonld, as_jsonld_string
    write_text(os.path.join(out_dir, "entity.jsonld"), as_jsonld_string(build_jsonld(bundle, slug)))

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

    # Cross-entity link mesh: register this drug under its ChEMBL ids + names.
    from atlas.page import links
    links.load(dist_dir)
    from atlas.render_common import display_name
    links.upsert(dist_dir, "drug", slug,
                 id_keys=[a.chembl_id, a.parent_chembl, *(a.child_chembls or ())],
                 name_keys=[a.canonical_name, *(getattr(a, "alt_names", None) or ())],
                 canonical=display_name(a.canonical_name))

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

    meta = build_meta("drug", slug, a.canonical_name or name, datasets_from_calls(CALLS), bundle=bundle)
    page_md = assemble_page(slug, summary_text, body_md, meta, bundle=bundle)

    write_text(os.path.join(out_dir, PAGE_FILENAME), page_md)
    if summary_text:
        write_text(os.path.join(out_dir, "summary.md"), summary_text + "\n")
    if judge_result is not None:
        write_json(os.path.join(out_dir, "judge.json"), judge_result, indent=2)

    # Sidecar — schema.org Drug card (bundle.json + provenance.json
    # intentionally not published; see run_gene note).
    from atlas.page.drug_jsonld import build_jsonld, as_jsonld_string
    write_text(os.path.join(out_dir, "entity.jsonld"), as_jsonld_string(build_jsonld(bundle, slug)))

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
