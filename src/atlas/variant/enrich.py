"""Deterministic enrichment for variant pages — the layers that turn a ClinVar
echo into an integrative reference. All computed from primary data:

  - AlphaMissense in-silico pathogenicity (per missense, joined by protein_variant)
  - gnomAD population frequency + a rarity band (joined by rsID)
  - cross-source CONCORDANCE verdict (ClinVar vs AlphaMissense vs gnomAD)
  - submitter-consensus statistics (from the ClinVar submissions[])
  - same-residue hotspot context (from the gene's own ClinVar enumeration)

Nothing here is a clinical call: concordance/consensus are transparent
descriptions of independent evidence, each shown with its source.
"""
import re
from collections import Counter

from atlas.biobtree import map_all

_AA3to1 = {"Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q",
           "Glu": "E", "Gly": "G", "His": "H", "Ile": "I", "Leu": "L", "Lys": "K",
           "Met": "M", "Phe": "F", "Pro": "P", "Ser": "S", "Thr": "T", "Trp": "W",
           "Tyr": "Y", "Val": "V", "Ter": "*"}
_MISSENSE_RE = re.compile(r"p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})$")


def missense_short(hgvs_p):
    """'p.Pro309Ala' → 'P309A' (AlphaMissense key). None for non-missense
    (frameshift/del/dup/nonsense) — AlphaMissense only models missense SNVs."""
    m = _MISSENSE_RE.match((hgvs_p or "").strip())
    if not m:
        return None
    a1, a2 = _AA3to1.get(m.group(1)), _AA3to1.get(m.group(3))
    return f"{a1}{m.group(2)}{a2}" if a1 and a2 else None


def protein_position(hgvs_p):
    """Integer residue position from any p. HGVS ('p.Pro309Ala'→309), or None."""
    m = re.search(r"p\.[A-Za-z]{3}(\d+)", hgvs_p or "")
    return int(m.group(1)) if m else None


# ── per-gene caches (fetch once, reuse across all the gene's variants) ────────
def gene_alphamissense(hgnc_id):
    """{protein_short: (am_class, am_pathogenicity)} for the whole gene."""
    out = {}
    for r in map_all(hgnc_id, ">>hgnc>>uniprot>>alphamissense", cap=200):
        pv = (r.get("protein_variant") or "").strip()
        if pv:
            out[pv] = (r.get("am_class"), r.get("am_pathogenicity"))
    return out


def gnomad_for(rsid):
    """Population-frequency read for an rsID, or None. Absence is a real signal
    (a rare pathogenic allele is expected to be absent from gnomAD)."""
    if not rsid:
        return None
    d = map_all(rsid, ">>dbsnp")
    if not d:
        return None
    freq = (d[0].get("gnomad_frequency") or "").strip()
    absent = freq in ("", "0", "0.0")
    return {"frequency": freq, "absent": absent,
            "is_common": (d[0].get("is_common") == "true"),
            "band": _freq_band(freq, absent)}


def _freq_band(freq, absent):
    if absent:
        return "absent from gnomAD"
    try:
        f = float(freq)
    except ValueError:
        return "present in gnomAD"
    if f < 1e-4:
        return f"ultra-rare (gnomAD MAF {freq})"
    if f < 1e-3:
        return f"rare (gnomAD MAF {freq})"
    if f < 1e-2:
        return f"low-frequency (gnomAD MAF {freq})"
    return f"common (gnomAD MAF {freq})"


# ── derived analyses ─────────────────────────────────────────────────────────
def submitter_consensus(submissions):
    """Agreement among ClinVar submitters — from submissions[]. Returns a small
    dict {n, breakdown, verdict} or None."""
    calls = [(s.get("classification") or "").strip() for s in submissions if s.get("classification")]
    if not calls:
        return None
    bd = Counter(calls)
    top, topn = bd.most_common(1)[0]
    if len(bd) == 1:
        verdict = f"unanimous ({topn}/{len(calls)} {top})"
    else:
        parts = ", ".join(f"{n} {c}" for c, n in bd.most_common())
        verdict = f"conflicting ({parts})"
    return {"n": len(calls), "breakdown": dict(bd), "verdict": verdict}


def concordance(classification, am_entry, gnomad):
    """Cross-source concordance readout: does the independent in-silico +
    population evidence agree with ClinVar? Returns {lines, verdict, flags}.
    Never a clinical determination — a description of evidence agreement."""
    cls = (classification or "").lower()
    clinvar_path = "pathogenic" in cls and "conflict" not in cls
    lines, flags, agree, total = [], [], 0, 0

    if am_entry:
        amclass, amscore = am_entry
        total += 1
        pretty = (amclass or "").replace("_", " ")
        if clinvar_path and amclass == "likely_pathogenic":
            agree += 1
            lines.append(f"AlphaMissense **{pretty}** ({amscore}) — concordant with the pathogenic call")
        elif clinvar_path and amclass == "likely_benign":
            flags.append("AlphaMissense predicts likely-benign")
            lines.append(f"AlphaMissense **{pretty}** ({amscore}) — ⚠ discordant with the pathogenic call")
        else:
            lines.append(f"AlphaMissense **{pretty}** ({amscore})")
    else:
        lines.append("AlphaMissense — not scored (not a missense SNV)")

    if gnomad:
        total += 1
        if gnomad["absent"]:
            if clinvar_path:
                agree += 1
            lines.append("Absent from gnomAD — consistent with a rare pathogenic allele")
        elif gnomad["is_common"]:
            flags.append(f"common in gnomAD ({gnomad['frequency']})")
            lines.append(f"**{gnomad['band']}** — unusual for a pathogenic classification ⚠")
        else:
            lines.append(gnomad["band"].capitalize())

    if total == 0:
        verdict = None
    elif flags:
        verdict = "Evidence sources **disagree** — " + "; ".join(flags)
    elif agree == total and clinvar_path:
        verdict = f"{agree} independent line{'s' if agree != 1 else ''} of evidence **concordant** with the ClinVar classification"
    else:
        verdict = "Mixed / partial evidence (see below)"
    return {"lines": lines, "verdict": verdict, "flags": flags}


def residue_hotspot(hgvs_p, position_index):
    """Other pathogenic ClinVar variants at the same residue, from a prebuilt
    {position: [labels]} index over the gene's P/LP set. Returns a dict or None."""
    pos = protein_position(hgvs_p)
    if pos is None or not position_index:
        return None
    here = [x for x in position_index.get(pos, []) if x.get("hgvs_p") != hgvs_p]
    if not here:
        return None
    return {"position": pos, "others": here}
