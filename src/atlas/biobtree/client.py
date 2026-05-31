#!/usr/bin/env python3
"""Biobtree client — REST over HTTP/1.1 with a urllib3 connection pool.

Three primitives + four pure parsers + a per-call reproducibility log:
  search(term, source=None)       -> /api/search dict
  entry(identifier, source)       -> /api/entry dict
  bbmap(ids, chain, page=None)    -> /api/map dict
  rows(resp)                      -> dicts keyed by REST schema cols
  map_targets(resp)               -> flat list of target dicts (escape-aware)
  map_all(ids, chain, cap=60)     -> paginated map, deduped, loop-safe
  xref_counts(entry_resp)         -> {dataset: count} from entry's xref table
  CALLS                           -> per-call reproducibility log
"""
import json

import urllib3

API = "http://127.0.0.1:8000/api"
CALLS = []

# One pool, HTTP keep-alive. Retries off — section collectors handle retry
# at their own level when needed. maxsize bounds concurrent connections;
# 16 covers any parallel fan-out we add later.
_POOL = urllib3.PoolManager(
    num_pools=2,
    maxsize=16,
    retries=False,
    timeout=urllib3.Timeout(connect=5.0, read=15.0),
    block=False,
)


def _get(path: str, params: dict) -> dict:
    r = _POOL.request("GET", f"{API}/{path}", fields=params)
    body = json.loads(r.data)
    CALLS.append({"path": path, "params": params})
    return body


def search(term: str, source: str = None) -> dict:
    p = {"i": term}
    if source:
        p["s"] = source
    return _get("search", p)


def entry(identifier: str, source: str) -> dict:
    return _get("entry", {"i": identifier, "s": source})


def bbmap(ids: str, chain: str, page: str = None) -> dict:
    params = {"i": ids, "m": chain}
    if page:
        # cursor param is `p` (NOT `page`; FastAPI silently drops unknown
        # query params so the wrong name fails silently).
        params["p"] = page
    return _get("map", params)


def rows(resp: dict) -> list:
    """search/data rows -> list of dicts keyed by schema columns.

    biobtree returns `"data": null` (literal None, not the missing key) when
    a search has zero results — coerce defensively."""
    cols = (resp.get("schema") or "").split("|")
    return [dict(zip(cols, r.split("|"))) for r in (resp.get("data") or [])]


def map_targets(resp: dict) -> list:
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
