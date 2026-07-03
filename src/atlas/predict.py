"""Sugi Predict cross-links — Atlas → the sibling product at sugi.bio/predict.

Predict is chemical k-NN target prediction over ~30M SureChEMBL patent compounds,
grounded to 4,884 human protein targets (shared BioBTree identifier grounding).
Two link directions:

- gene → per-target page (`/predict/target/{SYMBOL}/`): only the 4,884 covered
  targets have a page (a non-target 404s), so gate strictly against the bundled
  target manifest (data/predict/targets.json, refreshed from
  https://sugi.bio/predict/static/targets.json).
- drug → structure prediction (`/predict/?q={SMILES}`): runs k-NN on the drug's
  own structure, so any drug with a SMILES links (no gate; below the 0.3
  confidence floor Predict just reports "novel chemistry").
"""
import functools
import json
import os
from urllib.parse import quote

_BASE = "https://sugi.bio/predict"


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


def smiles_url(smiles):
    """Predict structure-prediction URL for a drug's SMILES, or None when no
    SMILES is available. Uses the dedicated /predict/predict?smiles= endpoint."""
    s = (smiles or "").strip()
    if not s:
        return None
    return f"{_BASE}/predict?smiles={quote(s)}"
