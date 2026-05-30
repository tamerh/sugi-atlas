#!/usr/bin/env python3
"""Biobtree REST client — pure data access layer.

  search(term, source=None)               -> /api/search
  entry(identifier, source)               -> /api/entry
  bbmap(ids, chain, page=None)            -> /api/map (one page; cursor param `p`)
  rows(resp)                              -> search/data rows as dicts keyed by schema
  map_targets(resp)                       -> all map targets as dicts (escape-aware)
  map_all(ids, chain, cap=60)             -> paginated map, deduped, loop-safe

CALLS is a module-level reproducibility log appended on every request. Reset to
[] before a run if you want clean per-run accounting.
"""
import json, urllib.parse, urllib.request

API = "http://127.0.0.1:8000/api"
CALLS = []  # reproducibility log

def _get(path, params):
    qs = urllib.parse.urlencode(params)
    url = f"{API}/{path}?{qs}"
    with urllib.request.urlopen(url, timeout=15) as r:
        body = json.load(r)
    CALLS.append({"path": path, "params": params})
    return body

def search(term, source=None):
    p = {"i": term}
    if source:
        p["s"] = source
    return _get("search", p)

def entry(identifier, source):
    return _get("entry", {"i": identifier, "s": source})

def bbmap(ids, chain, page=None):
    params = {"i": ids, "m": chain}
    if page:
        params["p"] = page  # cursor param is `p` (NOT `page`; FastAPI silently
                            # drops unknown query params so the wrong name fails silently)
    return _get("map", params)

def rows(resp):
    """search/data rows -> list of dicts keyed by schema columns."""
    cols = resp.get("schema", "").split("|")
    return [dict(zip(cols, r.split("|"))) for r in resp.get("data", [])]

def map_targets(resp):
    """map response -> flat list of target dicts keyed by schema columns.

    Schema-aware split: some target fields embed an escaped pipe (\\|),
    e.g. bgee_evidence ids ENSG..\\|UBERON:.. — split on real separators only.
    """
    cols = resp.get("schema", "").split("|")
    out = []
    for m in (resp.get("mappings") or []):
        for t in (m.get("targets") or []):
            parts = [p.replace("\x00", "|") for p in t.replace("\\|", "\x00").split("|")]
            out.append(dict(zip(cols, parts)))
    return out

def map_all(ids, chain, cap=60):
    """All target dicts across pages, deduped. Pagination via the `p=` cursor.
    Dedupe + stop-on-no-new make it loop-safe; `cap` bounds total pages
    (cap*100 rows is the worst case)."""
    out, page, n, seen = [], None, 0, set()
    while True:
        resp = bbmap(ids, chain, page)
        new = 0
        for t in map_targets(resp):
            key = t.get("id") or tuple(t.values())
            if key in seen:
                continue
            seen.add(key)
            out.append(t)
            new += 1
        pg = resp.get("pagination", {}) or {}
        nxt = pg.get("next_token")
        if not pg.get("has_next") or n >= cap or not nxt or new == 0:
            break
        page, n = nxt, n + 1
    return out

def xref_counts(entry_resp):
    """hgnc/ensembl/transcript `entry` carries a dataset|count xref table —
    exact totals for clinvar/spliceai/msigdb/hpo/gwas/collectri/... in 1 call."""
    return {r.split("|")[0]: int(r.split("|")[1])
            for r in entry_resp.get("xrefs", {}).get("data", [])}
