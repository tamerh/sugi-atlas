#!/usr/bin/env python3
"""Parallel batch builder for the Atlas corpus (gene + disease + drug).

Three phases so the cross-entity mesh resolves completely in ONE collect pass,
with NO shared-manifest writes during the parallel parts (concurrency-safe):

  PHASE A — collect (parallel): per entity, resolve + collect the bundle, cache
            it to .atlas-build/, and return its manifest key-set. Workers touch
            only their own files → no race.
  PHASE B — merge (single): combine every key-set → one complete manifest.json.
  PHASE C — render (parallel): per entity, reload the cached bundle + the
            COMPLETE manifest (read-only) → render + write page.md / entity.jsonld.
            No re-collect; manifest read-only → no race.

Process pool (not threads): each process has its own biobtree.CALLS and
links._MANIFEST, so the module-global state can't race.

  python -m atlas.batch --dist /data/sugi-atlas-dist --workers 12 \
      --genes TP53,KRAS  --drugs imatinib,sotorasib  --diseases cardiomyopathy
  python -m atlas.batch --dist ... --genes @genes.txt --drugs @drugs.txt
"""
import argparse
import json
import os
import sys
import time
from collections import Counter
from datetime import datetime, timezone
from multiprocessing import Pool

from atlas.atomicio import write_text, write_json
from atlas.biobtree import CALLS
from atlas.gene import collect as C
from atlas.disease import collect as DC
from atlas.disease.anchors import resolve as resolve_disease
from atlas.disease.slug import slugify as disease_slug
from atlas.drug import collect as DRC
from atlas.drug.anchors import resolve as resolve_drug
from atlas.drug.slug import slugify as drug_slug
from atlas.disease import render as DR
from atlas.drug import render as DRR
from atlas.validation import body_gate
from atlas.page import links
from atlas import pipeline as P


def _cache_path(cache_dir, etype, slug):
    return os.path.join(cache_dir, etype, f"{slug}.json")


# ---- PHASE A: collect one entity (parallel-safe; own files only) -----------

def collect_one(spec):
    etype, ident, dist_dir, cache_dir = spec
    try:
        CALLS.clear()
        if etype == "gene":
            slug = ident
            bundle = P.collect_all(ident)   # includes the non-coding scrub
            b1 = bundle.get("1") or {}
            id_keys = [ident, b1.get("hgnc_id")]
            name_keys, title = [], ident
        elif etype == "disease":
            a = resolve_disease(ident)
            slug = disease_slug(a.canonical_name or ident)
            bundle = {sid: DC.REGISTRY[sid].collect_fn(a) for sid in DC.REGISTRY}
            id_keys = [a.mondo_id, a.efo_id]
            name_keys = [a.canonical_name, *(a.synonyms or ())]
            title = a.canonical_name or ident
        elif etype == "drug":
            a = resolve_drug(ident)
            slug = drug_slug(a.canonical_name or ident)
            bundle = {sid: DRC.REGISTRY[sid].collect_fn(a) for sid in DRC.REGISTRY}
            id_keys = [a.chembl_id, a.parent_chembl, *(a.child_chembls or ())]
            name_keys = [a.canonical_name, *(getattr(a, "alt_names", None) or ())]
            title = a.canonical_name or ident
        else:
            return {"ok": False, "entity": etype, "ident": ident, "error": "unknown entity"}

        datasets = P.datasets_from_calls(CALLS)
        # body_gate verdict (log only — never blocks the batch). Guarded so a
        # gate hiccup can't fail an otherwise-good collection.
        verdict = "?"
        try:
            snap = body_gate.snap_dir_for(dist_dir, etype)
            verdict = body_gate.check(slug, bundle, snap_dir=snap).get("verdict", "?")
        except Exception:
            pass

        write_json(_cache_path(cache_dir, etype, slug),
                   {"bundle": bundle, "datasets": datasets, "title": title,
                    "entity": etype, "slug": slug})
        # Destination canonical name for the mesh (audit #13). Drug titles are
        # de-SHOUTed for display; gene/disease titles pass through.
        from atlas.render_common import display_name
        canonical = display_name(title) if etype == "drug" else title
        # Raw evidence_score signal (weighted log of evidence counts) — the merge
        # phase ranks these into the 0-100 per-type percentile.
        from atlas.page import evidence
        ev_comps = evidence.components(etype, bundle)
        ev_raw = evidence.raw_signal(etype, ev_comps)
        return {"ok": True, "entity": etype, "slug": slug, "verdict": verdict,
                "canonical": canonical, "evidence_raw": ev_raw,
                "evidence_components": ev_comps,
                "id_keys": [str(k) for k in id_keys if k],
                "name_keys": [k for k in name_keys if k]}
    except (Exception, SystemExit) as e:   # SystemExit too — see harden M1
        return {"ok": False, "entity": etype, "ident": ident, "error": repr(e)}


# ---- PHASE C: render one entity (parallel-safe; reads complete manifest) ----

def render_one(spec):
    etype, slug, dist_dir, cache_dir, generated_at = spec
    try:
        payload = json.load(open(_cache_path(cache_dir, etype, slug)))
        bundle, datasets, title = payload["bundle"], payload["datasets"], payload["title"]
        links.load(dist_dir)                       # the COMPLETE manifest
        links.load_reverse(dist_dir)               # incoming cross-entity edges
        from atlas.page import evidence
        evidence.load(dist_dir)                  # per-type evidence_score scores
        meta = P.build_meta(etype, slug, title, datasets, generated_at, bundle=bundle)
        if etype == "gene":
            body = P.render_all(bundle)
            from atlas.page.jsonld import build_jsonld, as_jsonld_string
            jsonld = as_jsonld_string(build_jsonld(bundle))
        elif etype == "disease":
            bundle["_indicated_drugs"] = links.indicated_drugs(dist_dir, slug)
            # Sparse-subtype fallback: when this term has NO disease-direct
            # indications of its own, surface the PARENT term's indicated drugs
            # (clearly attributed in the render) so the Therapeutics zone isn't an
            # empty placeholder next to a parent that carries real therapeutic data.
            if not bundle["_indicated_drugs"]:
                parent = (bundle.get("1") or {}).get("parent") or {}
                pslug = (links._lookup("disease", parent.get("id"), parent.get("name"))
                         if parent else None)
                if pslug and pslug != slug:
                    pdrugs = links.indicated_drugs(dist_dir, pslug)
                    if pdrugs:
                        bundle["_parent_indicated_drugs"] = pdrugs
                        bundle["_parent_disease"] = {
                            "name": parent.get("name"),
                            "url": links.disease_url(mondo_id=parent.get("id"),
                                                     name=parent.get("name"))}
            body = DR.render_all(bundle)
            from atlas.page.disease_jsonld import build_jsonld, as_jsonld_string
            jsonld = as_jsonld_string(build_jsonld(bundle, slug))
        else:
            body = DRR.render_all(bundle)
            from atlas.page.drug_jsonld import build_jsonld, as_jsonld_string
            jsonld = as_jsonld_string(build_jsonld(bundle, slug))
        page_md = P.assemble_page(slug, "", body, meta, bundle=bundle)
        out_dir = os.path.join(dist_dir, "atlas", etype, slug)
        write_text(os.path.join(out_dir, P.PAGE_FILENAME), page_md)
        write_text(os.path.join(out_dir, "entity.jsonld"), jsonld)
        return {"ok": True, "entity": etype, "slug": slug}
    except (Exception, SystemExit) as e:
        return {"ok": False, "entity": etype, "slug": slug, "error": repr(e)}


def _drain(label, pool, fn, specs):
    """Run fn over specs via imap_unordered with a live per-type progress line
    (rate + ETA + per-entity-type tally) — so a multi-hour collect isn't blind
    and the disease vs gene vs drug throughput is visible as it runs."""
    total = len(specs)
    step = max(1, total // 40)
    out, t, last = [], time.time(), time.time()
    by_type, skip = Counter(), 0
    for i, r in enumerate(pool.imap_unordered(fn, specs, chunksize=1), 1):
        out.append(r)
        if r.get("ok"):
            by_type[r["entity"]] += 1
        else:
            skip += 1
        now = time.time()
        if i % step == 0 or i == total or now - last >= 30:
            rate = i / max(1e-9, now - t)
            eta = (total - i) / rate if rate else 0
            types = " ".join(f"{k}={v}" for k, v in sorted(by_type.items()))
            print(f"[{label}] {i}/{total}  {rate:.1f}/s  eta {eta/60:.0f}m | "
                  f"{types} skip={skip}", flush=True)
            last = now
    return out


def _merge_manifest(collected, dist_dir):
    """PHASE B — one writer builds the whole manifest from every key-set."""
    manifest = {"gene": {}, "disease": {}, "drug": {},
                "canon": {"gene": {}, "disease": {}, "drug": {}}}
    for r in collected:
        bucket = manifest[r["entity"]]
        for k in r["id_keys"]:
            bucket[k] = r["slug"]
        for k in r["name_keys"]:
            nk = links._norm(k)
            if nk:
                bucket[nk] = r["slug"]
        if r.get("canonical"):                       # audit #13 destination name
            manifest["canon"][r["entity"]][r["slug"]] = r["canonical"]
    write_json(os.path.join(dist_dir, "atlas", "manifest.json"),
               manifest, indent=0, sort_keys=True)
    return manifest


def _write_evidence(collected, dist_dir):
    """PHASE B — rank each entity's raw evidence_score signal into a 0-100 percentile
    within its type, and persist the per-type distribution so a later single-
    entity rebuild ranks against the same frozen corpus."""
    from atlas.page import evidence
    raw_by_type = {}
    comp_by_type = {}                    # {etype: {component_key: [counts]}}
    for r in collected:
        raw_by_type.setdefault(r["entity"], {})[r["slug"]] = r.get("evidence_raw", 0.0)
        cd = comp_by_type.setdefault(r["entity"], {})
        for k, v in (r.get("evidence_components") or {}).items():
            try:
                cd.setdefault(k, []).append(int(v or 0))
            except (TypeError, ValueError):
                cd.setdefault(k, []).append(0)
    out = {et: evidence.percentiles(m) for et, m in raw_by_type.items()}
    out["_dist"] = {et: sorted(m.values()) for et, m in raw_by_type.items()}
    # Per-component sorted distributions → corpus-relative "top N%" framing at
    # render (evidence.rank_clause). Same frozen-distribution contract as _dist.
    out["_dist_components"] = {et: {k: sorted(vs) for k, vs in cd.items()}
                               for et, cd in comp_by_type.items()}
    # Per-slug component counts → a render reads a RELATED page's numbers
    # (evidence.components_for): a thin disease quoting its parent's evidence, a
    # lncRNA quoting its sense gene. Drop the all-zero entities to keep it lean.
    out["_components"] = {}
    for r in collected:
        comps = {k: int(v or 0) for k, v in (r.get("evidence_components") or {}).items()}
        if any(comps.values()):
            out["_components"].setdefault(r["entity"], {})[r["slug"]] = comps
    write_json(os.path.join(dist_dir, "atlas", "evidence.json"),
               out, indent=0, sort_keys=True)
    print(f"[B] evidence: scored {sum(len(m) for m in raw_by_type.values())} entities")


def _build_reverse_index(collected, dist_dir, cache_dir):
    """PHASE B (second half) — invert every resolved cross-entity edge so an
    asserted relationship is navigable from both ends (TP53→Venetoclax becomes
    Venetoclax←TP53). Single-pass over the cached bundles, run against the now-
    complete manifest; written as the reverse_edges.json sidecar the render
    reads. Edges are kept directional + labeled at render time (links.py)."""
    links.load(dist_dir)                           # complete manifest, for resolution
    reverse = {}
    for r in collected:
        etype, slug = r["entity"], r["slug"]
        try:
            payload = json.load(open(_cache_path(cache_dir, etype, slug)))
        except (OSError, json.JSONDecodeError):
            continue
        src_url = f"/atlas/{etype}/{slug}/"
        src_label = r.get("canonical") or payload.get("title") or slug
        groups = links.related_targets(etype, payload["bundle"])
        for group, items in groups.items():
            for _label, target_url in items:
                if target_url and target_url != src_url:
                    reverse.setdefault(target_url, []).append(
                        [src_label, src_url, etype, group])
    # Sort each edge list to a stable order (group, label, url) — deterministic
    # regardless of collection order, and groups same-type back-links together.
    for edges in reverse.values():
        edges.sort(key=lambda e: (e[3], (e[0] or "").lower(), e[1]))
    write_json(os.path.join(dist_dir, "atlas", "reverse_edges.json"),
               reverse, indent=0, sort_keys=True)
    return reverse


def _build_indication_index(collected, dist_dir, cache_dir):
    """PHASE B — invert ChEMBL drug indications into a disease→drugs index so a
    disease page surfaces the drugs DIRECTLY indicated for it (not gene-cohort-
    mediated). Phase ≥2 captured (≤1 is noise — 90% fail, Levodopa↔asthma); the
    renderer then TIERS it: phase ≥3 as "indicated", phase 2 as a separate
    investigational sub-block. Keyed by disease URL, written as indicated_drugs.json
    the render reads via links.indicated_drugs(). This is the disease-direct
    therapeutic edge — the only drug content a non-molecular disease can show."""
    links.load(dist_dir)                           # complete manifest, for resolution
    idx = {}                                        # disease_url -> {drug_url: {name, max_phase}}
    for r in collected:
        if r["entity"] != "drug":
            continue
        try:
            payload = json.load(open(_cache_path(cache_dir, r["entity"], r["slug"])))
        except (OSError, json.JSONDecodeError):
            continue
        drug_url = f"/atlas/drug/{r['slug']}/"
        drug_name = r.get("canonical") or payload.get("title") or r["slug"]
        for i in (payload["bundle"].get("4") or {}).get("indications") or []:
            phase = i.get("max_phase") or 0
            if phase < 2:                          # ≤1: exploratory noise, excluded
                continue
            d_url = links.disease_url(mondo_id=i.get("mondo_id"), name=i.get("name"))
            if not d_url or d_url == drug_url:
                continue
            slot = idx.setdefault(d_url, {}).setdefault(
                drug_url, {"name": drug_name, "url": drug_url, "max_phase": 0, "approved": False})
            slot["max_phase"] = max(slot["max_phase"], phase)   # best phase per drug↔disease
            # `approved` (set in s04): phase 4, or an anticancer drug at phase 3 vs
            # a cancer (ChEMBL logs imatinib→CML at phase 3 though approved).
            if i.get("approved"):
                slot["approved"] = True
    # Approved-first, then phase, then name — deterministic ordering.
    out = {d_url: sorted(drugs.values(),
                         key=lambda x: (not x["approved"], -x["max_phase"], x["name"]))
           for d_url, drugs in idx.items()}
    write_json(os.path.join(dist_dir, "atlas", "indicated_drugs.json"),
               out, indent=0, sort_keys=True)
    return out


def _parse_list(val):
    """Comma-separated values, or @file (one entity per line)."""
    if not val:
        return []
    if val.startswith("@"):
        with open(val[1:]) as f:
            return [ln.strip() for ln in f if ln.strip() and not ln.startswith("#")]
    return [x.strip() for x in val.split(",") if x.strip()]


def run(genes, diseases, drugs, dist_dir, cache_dir, workers, limit=None):
    if limit:                       # test slice: first N of each category
        genes, diseases, drugs = genes[:limit], diseases[:limit], drugs[:limit]
    specs_a = ([("gene", g, dist_dir, cache_dir) for g in genes]
               + [("disease", d, dist_dir, cache_dir) for d in diseases]
               + [("drug", x, dist_dir, cache_dir) for x in drugs])
    total = len(specs_a)
    print(f"[batch] {total} entities | {workers} workers | dist={dist_dir}")

    t0 = time.time()
    print(f"[A] collect ({total}) …", flush=True)
    with Pool(workers) as pool:
        collected = _drain("A", pool, collect_one, specs_a)
    # imap_unordered returns in completion order — sort to a stable (entity, slug)
    # order so every phase-B output is byte-identical across runs: manifest
    # name-key collisions resolve deterministically and reverse-edge lists are
    # built in a fixed order (the "deterministic build" guarantee).
    ok = sorted((r for r in collected if r.get("ok")),
                key=lambda r: (r["entity"], r["slug"]))
    failed = [r for r in collected if not r.get("ok")]

    # Retry pass — a TRANSIENT biobtree failure (timeout / 5xx under full-regen
    # contention) usually succeeds on a calm low-concurrency retry once the
    # 30-worker load is gone; this is what dropped TTN/BRCA2 from v1.2.0. TERMINAL
    # failures (no MONDO/HGNC row — the entity genuinely doesn't resolve) are not
    # retryable, so they're excluded. Transient ⇔ a BiobtreeError reached us.
    def _transient(r):
        return "BiobtreeError" in (r.get("error") or "")
    retry = [r for r in failed if _transient(r)]
    if retry:
        print(f"[A] retrying {len(retry)} transient skip(s) at low concurrency …", flush=True)
        with Pool(min(4, workers)) as pool:
            retried = _drain("A-retry", pool, collect_one,
                             [(r["entity"], r["ident"], dist_dir, cache_dir) for r in retry])
        recovered = [r for r in retried if r.get("ok")]
        if recovered:
            ok = sorted(ok + recovered, key=lambda r: (r["entity"], r["slug"]))
        failed = ([r for r in failed if not _transient(r)]      # terminal — keep
                  + [r for r in retried if not r.get("ok")])    # still-failing transient
        print(f"[A] retry recovered {len(recovered)}/{len(retry)}", flush=True)

    # Classify + report the FINAL skips: terminal (expected — no resolution) vs
    # transient (UNEXPECTED — survived a calm retry, so probably data lost to
    # repeated biobtree timeouts; these get a loud flag + the skipped.json report).
    terminal = [r for r in failed if not _transient(r)]
    transient_left = [r for r in failed if _transient(r)]
    write_json(os.path.join(dist_dir, "atlas", "skipped.json"),
               [{"entity": r["entity"], "ident": r.get("ident"),
                 "kind": "transient" if _transient(r) else "terminal",
                 "error": r.get("error")} for r in failed], indent=0, sort_keys=True)
    for r in terminal:
        print(f"    SKIP (terminal) {r['entity']}:{r.get('ident')} — {r.get('error')}")
    for r in transient_left:
        print(f"    ⚠ SKIP (TRANSIENT, unrecovered) {r['entity']}:{r.get('ident')} — {r.get('error')}")
    if transient_left:
        print(f"[A] ⚠ {len(transient_left)} TRANSIENT skip(s) survived retry — likely data loss "
              f"to repeated biobtree timeouts; investigate (atlas/skipped.json).", flush=True)
    print(f"[A] collected {len(ok)}/{total} in {time.time()-t0:.1f}s "
          f"({len(terminal)} terminal, {len(transient_left)} transient skipped)")

    print(f"[B] merge manifest …")
    _merge_manifest(ok, dist_dir)
    print(f"[B] build reverse-edge index …")
    rev = _build_reverse_index(ok, dist_dir, cache_dir)
    print(f"[B] reverse index: {len(rev)} targets with incoming edges")
    ind = _build_indication_index(ok, dist_dir, cache_dir)
    print(f"[B] indication index: {len(ind)} diseases with ≥1 indicated drug")
    _write_evidence(ok, dist_dir)

    gen_at = datetime.now(timezone.utc).isoformat(timespec="seconds")
    specs_c = [(r["entity"], r["slug"], dist_dir, cache_dir, gen_at) for r in ok]
    print(f"[C] render ({len(specs_c)}) …", flush=True)
    t1 = time.time()
    with Pool(workers) as pool:
        rendered = _drain("C", pool, render_one, specs_c)
    rok = [r for r in rendered if r.get("ok")]
    for r in rendered:
        if not r.get("ok"):
            print(f"    RENDER-FAIL {r['entity']}:{r.get('slug')} — {r.get('error')}")
    print(f"[C] rendered {len(rok)}/{len(specs_c)} in {time.time()-t1:.1f}s")
    print(f"[batch] done in {time.time()-t0:.1f}s — {len(rok)} pages published")
    return {"collected": len(ok), "rendered": len(rok), "skipped": len(failed)}


def main(argv=None):
    ap = argparse.ArgumentParser(description="Parallel Atlas corpus builder")
    ap.add_argument("--dist", required=True)
    ap.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 4) - 2))
    ap.add_argument("--cache", default=None,
                    help="bundle cache dir (default: <dist>/cache)")
    ap.add_argument("--genes", default="")
    ap.add_argument("--diseases", default="")
    ap.add_argument("--drugs", default="")
    ap.add_argument("--limit", type=int, default=None,
                    help="build only the first N of each category (test slice)")
    a = ap.parse_args(argv)
    genes, diseases, drugs = _parse_list(a.genes), _parse_list(a.diseases), _parse_list(a.drugs)
    if not (genes or diseases or drugs):
        ap.error("give at least one of --genes / --diseases / --drugs")
    cache = a.cache or os.path.join(a.dist, "cache")   # keep all build artifacts under <dist>/
    run(genes, diseases, drugs, a.dist, cache, a.workers, limit=a.limit)


if __name__ == "__main__":
    main()
