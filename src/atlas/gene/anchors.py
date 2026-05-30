#!/usr/bin/env python3
"""Gene anchors — resolve a symbol to all the IDs the 12 sections need, ONCE.

Every section used to call _resolve(symbol) and pay for the same hgnc + ensembl
+ uniprot lookups (~24 redundant biobtree calls per gene). This module does it
once; sections receive an immutable Anchors record and consume what they need.
"""
import re, sys
from dataclasses import dataclass
from typing import Optional

from atlas.biobtree import search, entry, rows, map_all

@dataclass(frozen=True)
class Anchors:
    symbol: str
    hgnc_id: str
    hgnc_entry: dict
    ensembl_id: Optional[str]
    reviewed_uniprots: tuple   # all reviewed (Swiss-Prot) products; >1 for dual-product genes
    canonical_uniprot: Optional[str]
    canonical_transcript: Optional[str]  # MANE-Select Ensembl ENST, when present

def resolve_hgnc(symbol):
    """Symbol -> (hgnc_id, hgnc_entry) robustly.

    Filter the search to the hgnc dataset (the unfiltered top-50 can omit the
    hgnc row for high-xref genes), then disambiguate ambiguous symbols (e.g.
    'AR' matches amphiregulin AND androgen receptor) by the exact approved
    symbol from each candidate's hgnc entry."""
    cand = [r["id"] for r in rows(search(symbol, source="hgnc"))
            if re.match(r"HGNC:\d+$", r.get("id", ""))]
    if not cand:  # fallback: unfiltered
        cand = [r["id"] for r in rows(search(symbol)) if r.get("dataset") == "hgnc"]
    if not cand:
        sys.exit(f"no HGNC row for {symbol}")
    if len(cand) == 1:
        return cand[0], entry(cand[0], "hgnc")
    for cid in cand[:8]:
        he = entry(cid, "hgnc")
        syms = he.get("Attributes", {}).get("Hgnc", {}).get("symbols", [])
        if symbol.upper() in [s.upper() for s in syms]:
            return cid, he
    return cand[0], entry(cand[0], "hgnc")  # last resort

def _canonical_transcript(ensembl_id):
    """MANE-Select Ensembl ENST for a gene (fallback: first transcript)."""
    if not ensembl_id:
        return None
    mane = map_all(ensembl_id, ">>ensembl>>refseq[is_mane_select==true]")
    mane_nm = next((t["id"] for t in mane if t.get("type") == "mRNA"), None)
    if mane_nm:
        ct = map_all(mane_nm, ">>refseq>>transcript")
        if ct:
            return ct[0]["id"]
    tr = map_all(ensembl_id, ">>ensembl>>transcript")
    return tr[0]["id"] if tr else None

def resolve(symbol):
    """Symbol -> Anchors (all the shared IDs the 12 sections will need)."""
    hgnc_id, he = resolve_hgnc(symbol)
    ens = map_all(hgnc_id, ">>hgnc>>ensembl")
    ensembl_id = ens[0]["id"] if ens else None
    reviewed = tuple(t["id"] for t in map_all(hgnc_id, ">>hgnc>>uniprot"))
    canonical_u = reviewed[0] if reviewed else None
    return Anchors(
        symbol=symbol,
        hgnc_id=hgnc_id,
        hgnc_entry=he,
        ensembl_id=ensembl_id,
        reviewed_uniprots=reviewed,
        canonical_uniprot=canonical_u,
        canonical_transcript=_canonical_transcript(ensembl_id),
    )
