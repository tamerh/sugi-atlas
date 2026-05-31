"""urllib3 PoolManager transport — HTTP keep-alive variant.

Same wire format and response shape as the urllib transport — same
biobtree REST endpoint, same JSON, same fields. The only difference is
TCP-level: one persistent HTTPConnection is reused across calls instead
of opening a fresh socket per call.

For a disease build (~11k calls/disease) this saves the per-call TCP
setup × 11k (~5-10s/disease). The JSON parsing + biobtree compute costs
stay the same.

Why not requests.Session: requests is a heavier dep with the same
underlying urllib3. Skipping the wrapper.

Limitations vs urllib transport:
- One persistent connection per (scheme, host, port) by default. With
  high concurrency we'd want a larger pool — see PoolManager docstring.
- Slightly different error surface: urllib3.exceptions.HTTPError on
  failures (vs urllib.error.URLError). Callers don't currently catch
  either explicitly so this is fine.

Same exposed primitives: search, entry, bbmap. Returns the same dict shape.
"""
import json

import urllib3

# Shared CALLS log singleton (the dispatcher imports this and exposes it).
CALLS = []
API = "http://127.0.0.1:8000/api"

# One pool, keep-alive, retries off (the section collectors retry at their
# own level when needed). maxsize is the max number of *concurrent*
# connections to one host; 16 covers our parallel-fan use case if we
# add it later.
_POOL = urllib3.PoolManager(
    num_pools=2,
    maxsize=16,
    retries=False,
    timeout=urllib3.Timeout(connect=5.0, read=15.0),
    block=False,
)


def _get(path: str, params: dict) -> dict:
    # urllib3 handles query-string encoding when params is passed via `fields=`.
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
        params["p"] = page
    return _get("map", params)
