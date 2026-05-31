"""urllib-based transport — the reproducibility baseline.

Stdlib only; no connection reuse (each call opens a fresh TCP connection).
Default transport so existing-clone-and-run setups behave exactly as before.
Same wire format and response shape as every other transport.
"""
import json
import urllib.parse
import urllib.request

# Shared with the dispatcher's CALLS log via module-level singleton import.
CALLS = []
API = "http://127.0.0.1:8000/api"


def _get(path: str, params: dict) -> dict:
    qs = urllib.parse.urlencode(params)
    url = f"{API}/{path}?{qs}"
    with urllib.request.urlopen(url, timeout=15) as r:
        body = json.load(r)
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
        params["p"] = page  # cursor param is `p` (NOT `page`; FastAPI silently
                            # drops unknown query params so the wrong name fails silently)
    return _get("map", params)
