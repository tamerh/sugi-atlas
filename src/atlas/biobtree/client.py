#!/usr/bin/env python3
"""Biobtree client — public API + transport dispatcher.

Public API (stable; never changes regardless of transport):
  search(term, source=None)               -> /api/search-equivalent dict
  entry(identifier, source)               -> /api/entry-equivalent dict
  bbmap(ids, chain, page=None)            -> /api/map-equivalent dict
  rows(resp)                              -> dicts keyed by REST schema cols
  map_targets(resp)                       -> flat list of target dicts (escape-aware)
  map_all(ids, chain, cap=60)             -> paginated map, deduped, loop-safe
  xref_counts(entry_resp)                 -> {dataset: count} from entry's xref table
  CALLS                                   -> per-call reproducibility log (module-level)

Transport selection (env var `ATLAS_BIOBTREE_TRANSPORT`, default `urllib`):
  urllib       — stdlib urllib.request, no connection reuse. Reproducibility baseline.
  urllib_pool  — urllib3 PoolManager with HTTP keep-alive. ~15-20% faster.
  grpc         — (reserved, not implemented yet) gRPC over HTTP/2.

The transport only owns search/entry/bbmap (the wire layer). Everything
else — rows(), map_targets(), map_all(), xref_counts() — is pure parsing
on the returned dict shape and lives here, transport-agnostic.

A/B comparison:
  ATLAS_BIOBTREE_TRANSPORT=urllib       python -m atlas.disease.corpus run ...
  ATLAS_BIOBTREE_TRANSPORT=urllib_pool  python -m atlas.disease.corpus run ...
"""
import os

# Transport selection. Resolved once at import time; flipping the env var
# at runtime won't re-pick a different transport.
_TRANSPORT_NAME = (os.environ.get("ATLAS_BIOBTREE_TRANSPORT") or "urllib").lower()

if _TRANSPORT_NAME == "urllib_pool":
    from atlas.biobtree._transports import urllib_pool as _t
elif _TRANSPORT_NAME == "grpc":
    raise NotImplementedError(
        "gRPC transport not yet implemented. Stubs land at "
        "atlas/biobtree/_pb/ (run `python -m atlas.biobtree._pb.regenerate` "
        "after biobtree's proto changes). Per-dataset field-rename map is "
        "the remaining work. For now use urllib or urllib_pool.")
elif _TRANSPORT_NAME == "urllib":
    from atlas.biobtree._transports import urllib_transport as _t
else:
    raise ValueError(
        f"Unknown ATLAS_BIOBTREE_TRANSPORT={_TRANSPORT_NAME!r}. "
        f"Pick one of: urllib (default), urllib_pool, grpc (reserved).")

# Re-export transport-owned primitives + the shared CALLS log.
API = _t.API
CALLS = _t.CALLS

search = _t.search
entry = _t.entry
bbmap = _t.bbmap


def rows(resp: dict) -> list:
    """search/data rows -> list of dicts keyed by schema columns.

    biobtree returns `"data": null` (literal None, not the missing key) when
    a search has zero results — `.get(key, [])` would still pick that up
    as None. Coerce defensively."""
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
