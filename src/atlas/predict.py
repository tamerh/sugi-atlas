"""Sugi Predict cross-links — Atlas → the sibling product at sugi.bio/predict.

Predict is chemical k-NN target prediction over ~30M SureChEMBL patent compounds,
grounded to 4,884 human protein targets (shared BioBTree identifier grounding).
Two link directions:

- gene → per-target page (`/predict/target/{SYMBOL}/`): only the 4,884 covered
  targets have a page (a non-target 404s), so gate strictly against the bundled
  target manifest (data/predict/targets.json, refreshed from
  https://sugi.bio/predict/static/targets.json).
- drug → per-compound page (`/predict/compound/{SCHEMBL}`): resolve the drug's
  SMILES to its SureChEMBL id via the internal resolver; a drug not in Predict's
  atlas gets no link, and a down resolver fails the build (no silent fallback).
"""
import functools
import json
import os
from urllib.parse import quote

import urllib3

_BASE = "https://sugi.bio/predict"

# Internal SMILES → SureChEMBL resolver (a drug's structure to its atlas compound
# id). Build-time only; overridable via env. A 303 → SCHEMBL id; a non-303 (200)
# means the structure isn't in the atlas → None. A DOWN resolver RAISES (no silent
# fallback) so the build fails loudly rather than shipping search links.
_RESOLVER = os.environ.get("ATLAS_PREDICT_RESOLVER", "http://127.0.0.1:8012/predict")

# The corpus build forks worker processes (multiprocessing.Pool). A urllib3 pool
# created in the parent and shared across the fork crosses responses between
# processes (drug A getting drug B's SureChEMBL id). Key the pool by PID so every
# process — parent and each forked child — lazily builds its OWN pool.
_pool_by_pid = {}


def _pool():
    pid = os.getpid()
    p = _pool_by_pid.get(pid)
    if p is None:
        p = urllib3.PoolManager(num_pools=2, timeout=urllib3.Timeout(total=8.0))
        _pool_by_pid[pid] = p
    return p


class PredictResolverError(RuntimeError):
    """The Sugi Predict SMILES→SureChEMBL resolver is unreachable."""


def check_resolver():
    """Pre-flight: confirm the resolver answers (any HTTP status = up). Raises
    PredictResolverError if unreachable, so a build with drugs aborts up front
    rather than silently dropping every drug's compound link."""
    try:
        _pool().request("GET", _RESOLVER, fields={"smiles": "C"},
                        redirect=False, retries=False)
    except urllib3.exceptions.HTTPError as e:
        raise PredictResolverError(
            f"Sugi Predict resolver unreachable at {_RESOLVER}: {e}. Start the "
            "resolver (or set ATLAS_PREDICT_RESOLVER) before building drugs.") from e


def resolve_schembl(smiles):
    """SureChEMBL id (e.g. 'SCHEMBL10883') for a drug's SMILES via the internal
    resolver, or None when the structure isn't in the atlas (resolver returns a
    non-303). Raises PredictResolverError if the resolver is unreachable — the
    build must not silently degrade."""
    s = (smiles or "").strip()
    if not s:
        return None
    try:
        r = _pool().request("GET", _RESOLVER, fields={"smiles": s},
                            redirect=False, retries=False)
    except urllib3.exceptions.HTTPError as e:
        raise PredictResolverError(
            f"Sugi Predict resolver unreachable at {_RESOLVER}: {e}. Start the "
            "resolver (or set ATLAS_PREDICT_RESOLVER) before building.") from e
    if r.status == 303:
        loc = r.headers.get("Location") or r.headers.get("location") or ""
        sid = loc.rstrip("/").rsplit("/", 1)[-1]
        return sid if sid.startswith("SCHEMBL") else None
    return None


@functools.lru_cache(maxsize=1)
def _covered():
    """(symbols, accessions) covered by Predict, upper-cased for case-insensitive
    matching. ({}, {}) if the manifest is missing (links simply don't emit)."""
    here = os.path.dirname(os.path.abspath(__file__))
    for path in ("data/predict/targets.json",
                 os.path.join(here, "..", "..", "data", "predict", "targets.json")):
        try:
            with open(path) as f:
                rows = json.load(f)
        except (OSError, json.JSONDecodeError):
            continue
        symbols = {(r.get("gene") or "").upper() for r in rows if r.get("gene")}
        accs = {(r.get("acc") or "").upper() for r in rows if r.get("acc")}
        return symbols, accs
    return set(), set()


def target_url(symbol=None, uniprot=None):
    """Predict per-target page URL for a gene, or None when the gene isn't one of
    Predict's covered targets (avoids linking to a 404). Prefers the gene symbol
    for the URL; falls back to the UniProt accession. Match is by EITHER symbol
    or accession, so a symbol-synonym drift still resolves via the accession."""
    symbols, accs = _covered()
    sym = (symbol or "").strip()
    acc = (uniprot or "").strip()
    covered = (sym and sym.upper() in symbols) or (acc and acc.upper() in accs)
    if not covered:
        return None
    key = sym or acc
    return f"{_BASE}/target/{quote(key)}/"


def compound_url(schembl_id):
    """Predict per-compound page URL for a SureChEMBL id, or None. The only drug
    cross-link — a drug not in Predict's atlas simply gets no link."""
    sid = (schembl_id or "").strip()
    return f"{_BASE}/compound/{quote(sid)}" if sid else None
