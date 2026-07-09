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
    # cap high enough for full per-gene coverage (~2.5k substitutions/small gene)
    # so the percentile stats see the whole distribution.
    for r in map_all(hgnc_id, ">>hgnc>>uniprot>>alphamissense", cap=60):
        pv = (r.get("protein_variant") or "").strip()
        if pv:
            out[pv] = (r.get("am_class"), r.get("am_pathogenicity"))
    return out


_GENOMIC_SNV = re.compile(r"NC_0*\d+\.\d+:g\.(\d+)([ACGT]+)>([ACGT]+)$")


def variant_coordinate(rec):
    """GRCh38 gnomAD-style key 'chr:pos:ref:alt' (no chr prefix). SNVs only.

    The variant carries genomic HGVS in BOTH assemblies (NC_..10 = GRCh37,
    NC_..11 = GRCh38); to pick the GRCh38 one unambiguously we require the HGVS
    position to equal the ClinVar `start` (which ClinVar reports on GRCh38).
    Using the wrong assembly would key gnomAD wrong and read as a false absence.
    """
    start = rec.get("start")
    chrom = str(rec.get("chromosome") or "").strip()
    if not (start and chrom):
        return None
    for e in rec.get("hgvs_expressions") or []:
        m = _GENOMIC_SNV.match(e)
        if m and m.group(1) == str(start):
            return f"{chrom}:{start}:{m.group(2)}:{m.group(3)}"
    return None


def gnomad_frequency(rec):
    """gnomAD v4.1 per-variant frequency — the ACMG BA1/BS1/PM2 layer. Queried by
    coordinate (the rsID reverse-hop isn't live yet). Falls back to the dbSNP
    global frequency when no genomic coordinate can be parsed (e.g. indels)."""
    coord = variant_coordinate(rec)
    if coord:
        g = map_all(coord, ">>gnomad_variant")
        if g:
            af = _f(g[0].get("af"))
            popmax = _f(g[0].get("af_grpmax"))
            anc = g[0].get("grpmax_ancestry")
            return {"af": g[0].get("af"), "popmax": popmax, "ancestry": anc,
                    "absent": False, "is_common": (popmax or 0) >= 0.05,
                    "band": _gnomad_band(af, popmax, anc), "source": "gnomAD v4.1"}
        return {"absent": True, "is_common": False, "popmax": None,
                "band": "absent from gnomAD v4.1", "source": "gnomAD v4.1"}
    # fallback — dbSNP inline global frequency (no coordinate to key gnomAD v4.1)
    return gnomad_for(rec.get("rsid"))


def _f(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def _gnomad_band(af, popmax, ancestry):
    """ACMG-flavoured frequency band from gnomAD popmax (the ACMG-standard
    metric). Descriptive, not a clinical criterion call."""
    p = popmax if popmax is not None else af
    if not p:
        return "absent from gnomAD v4.1"
    anc = f", {ancestry}" if ancestry else ""
    if p >= 0.05:
        return f"common (popmax {p:.1%}{anc}) — too common for a highly-penetrant pathogenic allele"
    if p >= 0.01:
        return f"low-frequency (popmax {p:.2%}{anc})"
    if p >= 1e-3:
        return f"rare (popmax {p:.3%}{anc})"
    if p >= 1e-4:
        return f"ultra-rare (popmax {_sci(p)}{anc})"
    return f"very rare (popmax {_sci(p)}{anc})"


def _sci(x):
    return f"{x:.1e}"


def gnomad_for(rsid):
    """Population-frequency read from the dbSNP inline global frequency (fallback
    when no coordinate is available). Absence is itself a signal."""
    if not rsid:
        return None
    d = map_all(rsid, ">>dbsnp")
    if not d:
        return None
    freq = (d[0].get("gnomad_frequency") or "").strip()
    absent = freq in ("", "0", "0.0")
    return {"frequency": freq, "absent": absent,
            "is_common": (d[0].get("is_common") == "true"),
            "band": _freq_band(freq, absent), "source": "dbSNP/gnomAD"}


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


def concordance(classification, am_entry, gnomad, spliceai=None):
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
        if gnomad.get("absent"):
            if clinvar_path:
                agree += 1
            lines.append("Absent from gnomAD v4.1 — consistent with a rare pathogenic allele")
        elif gnomad.get("is_common"):
            flags.append(gnomad.get("band", "common in gnomAD"))
            lines.append(f"**{gnomad['band']}** ⚠")
        else:
            lines.append(gnomad["band"][0].upper() + gnomad["band"][1:] + " (gnomAD v4.1)")

    if spliceai:
        total += 1
        if clinvar_path:
            agree += 1
        lines.append(f"SpliceAI predicts **{(spliceai.get('effect') or '').replace('_', ' ')}** "
                     f"(Δ {spliceai.get('score')}) — a splice-altering effect consistent with pathogenicity")

    if total == 0:
        verdict = None
    elif flags:
        verdict = "Evidence sources **disagree** — " + "; ".join(flags)
    elif agree == total and clinvar_path:
        verdict = f"{agree} independent line{'s' if agree != 1 else ''} of evidence **concordant** with the ClinVar classification"
    else:
        verdict = "Mixed / partial evidence (see below)"
    return {"lines": lines, "verdict": verdict, "flags": flags}


# ── Batch 3: per-gene context (fetched once per gene, cached in ctx) ─────────
def gene_context(hgnc_id):
    """ACMG-adjacent gene block: gnomAD constraint + ClinGen dosage + gene-disease
    validity + inheritance (GenCC/validity). Fetched once per gene."""
    con = map_all(hgnc_id, ">>hgnc>>gnomad_constraint")
    constraint = None
    if con:
        c = con[0]
        constraint = {"pli": c.get("pli"), "loeuf": c.get("loeuf"), "mis_z": c.get("mis_z")}
    dos = map_all(hgnc_id, ">>hgnc>>clingen_dosage")
    dosage = ({"haplo": dos[0].get("haplo_score"), "triplo": dos[0].get("triplo_score")}
              if dos else None)
    validity = [{"disease": v.get("disease_label"), "moi": v.get("moi"),
                 "classification": v.get("classification")}
                for v in map_all(hgnc_id, ">>hgnc>>clingen_gene_validity") if v.get("disease_label")]
    inh = []
    for g in map_all(hgnc_id, ">>hgnc>>gencc"):
        t = (g.get("moi_title") or "").strip()
        if t and t not in inh:
            inh.append(t)
    return {"constraint": constraint, "dosage": dosage, "validity": validity[:4],
            "inheritance": inh}


# Meaningful UniProt feature types → a human phrase ('{d}' = the description).
_UF_KEEP = {
    "domain": "the {d} domain", "region of interest": "the '{d}' region",
    "binding site": "a binding site", "active site": "an active site",
    "modified residue": "a modified residue ({d})", "helix": "an α-helix",
    "strand": "a β-strand", "turn": "a turn", "motif": "the {d} motif",
    "zinc finger": "a zinc-finger region", "metal binding": "a metal-binding site",
    "dna-binding region": "a DNA-binding region", "site": "a functional site ({d})",
    "nucleotide binding region": "a nucleotide-binding region",
    "cross-link": "a cross-link site", "disulfide bond": "a disulfide bond",
    "sequence variant": "a UniProt-annotated variant site ({d})",
}


def gene_structure(hgnc_id):
    """PDB structures + AlphaFold confidence + parsed UniProt feature intervals
    (for per-residue structural context). Fetched once per gene."""
    pdb = [{"id": p.get("id"), "title": p.get("title"), "method": p.get("method"),
            "resolution": p.get("resolution")}
           for p in map_all(hgnc_id, ">>hgnc>>uniprot>>pdb")]
    af = map_all(hgnc_id, ">>hgnc>>uniprot>>alphafold")
    alphafold = ({"plddt": af[0].get("global_metric"),
                  "frac_high": af[0].get("fraction_plddt_very_high")} if af else None)
    intervals = []
    for f in map_all(hgnc_id, ">>hgnc>>uniprot>>ufeature", cap=60):
        t = (f.get("type") or "").strip().lower()
        if t not in _UF_KEEP:
            continue
        try:
            b, e = int(f.get("location_begin")), int(f.get("location_end"))
        except (TypeError, ValueError):
            continue
        intervals.append({"type": t, "desc": (f.get("description") or "").strip(),
                          "begin": b, "end": e})
    return {"pdb": pdb, "alphafold": alphafold, "intervals": intervals}


def gene_spliceai(hgnc_id):
    """{coordinate: {effect, score}} of the gene's SpliceAI splice-impact
    predictions (chr:pos:ref:alt keys), fetched once per gene."""
    out = {}
    for r in map_all(hgnc_id, ">>hgnc>>spliceai", cap=60):
        cid = r.get("id")
        if cid and r.get("score"):
            out[cid] = {"effect": r.get("effect"), "score": r.get("score")}
    return out


def spliceai_for(coord, cache):
    """SpliceAI prediction for a variant's coordinate, if it has a meaningful
    (>=0.2) delta score — SpliceAI only annotates splice-relevant positions."""
    if not coord or not cache:
        return None
    hit = cache.get(coord)
    if not hit:
        return None
    try:
        if float(hit["score"]) < 0.2:
            return None
    except (TypeError, ValueError):
        return None
    return hit


def gene_pharmgkb(hgnc_id):
    """{rsID: {annotations, clinical}} pharmacogenomics for the gene, fetched once.
    Empty for the vast majority of genes (only pharmacogenes carry PGx)."""
    out = {}
    for a in map_all(hgnc_id, ">>hgnc>>pharmgkb_var_annotation", cap=20):
        rs = a.get("variant")
        if rs and a.get("drugs"):
            out.setdefault(rs, {"annotations": [], "clinical": []})["annotations"].append(
                {"drugs": a.get("drugs"), "category": a.get("phenotype_category"),
                 "significance": a.get("significance"), "sentence": a.get("sentence"),
                 "pmid": a.get("pmid")})
    for c in map_all(hgnc_id, ">>hgnc>>pharmgkb_clinical", cap=20):
        rs = c.get("variant")
        if rs and c.get("chemicals"):
            out.setdefault(rs, {"annotations": [], "clinical": []})["clinical"].append(
                {"chemicals": c.get("chemicals"), "level": c.get("level_of_evidence"),
                 "type": c.get("type"), "phenotypes": c.get("phenotypes")})
    return out


def pharmgkb_for(rsid, cache):
    """PGx for a variant's rsID from the per-gene cache, deduped."""
    if not rsid or rsid not in (cache or {}):
        return None
    e = cache[rsid]
    drugs = sorted({a["drugs"] for a in e["annotations"] if a.get("drugs")})
    return {"drugs": drugs[:8], "n": len(e["annotations"]),
            "clinical": sorted(e["clinical"], key=lambda c: (c.get("level") or "9"))[:5]}


def gene_has_civic(hgnc_id):
    """Whether the gene has any CIViC curation — gates the per-variant CIViC
    lookup so non-oncology genes cost one probe, not one call per variant."""
    return bool(map_all(hgnc_id, ">>hgnc>>civic", cap=1))


def civic_for(variation_id, has_civic):
    """CIViC predictive/prognostic evidence for a variant (cancer genes). Joined
    via the clinvar↔civic_variant xref. None unless the gene has CIViC data."""
    if not has_civic:
        return None
    cv = map_all(variation_id, ">>clinvar>>civic_variant")
    if not cv:
        return None
    ev = map_all(cv[0]["id"], ">>civic_variant>>civic_evidence")
    if not ev:
        return None
    return {"name": cv[0].get("name"),
            "evidence": [{"disease": e.get("disease"), "therapies": e.get("therapies"),
                          "type": e.get("evidence_type"), "level": e.get("evidence_level"),
                          "significance": e.get("significance")} for e in ev[:8]]}


def gene_variant_landscape(recs):
    """Aggregate profile of the gene's built variants — for the per-gene index.
    Zero new calls: everything from the in-memory record list."""
    from collections import Counter
    by_cls = Counter(r.get("classification") for r in recs)
    by_type = Counter(r.get("variant_type") for r in recs if r.get("variant_type"))
    pos = Counter()
    for r in recs:
        p = protein_position(r.get("hgvs_p"))
        if p is not None and "pathogenic" in (r.get("classification") or "").lower():
            pos[p] += 1
    recurrent = [(p, n) for p, n in pos.most_common(8) if n > 1]
    span = None
    ps = [p for p in pos]
    if ps:
        span = (min(ps), max(ps))
    return {"by_class": dict(by_cls), "by_type": dict(by_type.most_common(6)),
            "recurrent": recurrent, "n": len(recs), "residue_span": span}


def gene_condition_digest(hgnc_id):
    """Patient-facing digest of the gene's PRIMARY condition — routed via the
    gene's best-annotated Orphanet disease (NOT the narrow ClinVar MONDO, which
    is usually unannotated). Inheritance / onset / prevalence / HPO-with-frequency."""
    diseases = [d for d in map_all(hgnc_id, ">>hgnc>>orphanet")
                if (d.get("disorder_type") or "") == "Disease"]
    if not diseases:
        return None
    best = max(diseases, key=lambda d: int(d.get("phenotype_count") or 0))
    o = _orphanet_entry(best.get("id"))
    if not o:
        return None
    phen = sorted((o.get("phenotypes") or []),
                  key=lambda p: -(p.get("frequency_value") or 0))
    prev = (o.get("prevalences") or [{}])[0]
    return {
        "name": o.get("name") or best.get("name"),
        "definition": (o.get("definition") or "").strip(),
        "inheritance": o.get("inheritance") or [],
        "onset": o.get("onset") or [],
        "prevalence": (f"{prev.get('prevalence_class')} ({prev.get('geographic')})"
                       if prev.get("prevalence_class")
                       and prev.get("prevalence_class").lower() != "unknown" else None),
        "phenotypes": [{"term": p.get("hpo_term"), "freq": p.get("frequency")}
                       for p in phen[:8] if p.get("hpo_term")],
        "n_phenotypes": len(o.get("phenotypes") or []),
    }


def _orphanet_entry(oid):
    from atlas.biobtree import entry
    a = (entry(oid, "orphanet") or {}).get("Attributes") or {}
    return a.get("Orphanet") or (next(iter(a.values()), {}) if a else {})


def gene_panelapp(hgnc_id):
    """Green (diagnostic-grade) Genomics England panels the gene is on."""
    return [p.get("panel_name") for p in map_all(hgnc_id, ">>hgnc>>panelapp_gene")
            if (p.get("confidence") or "").lower() == "green" and p.get("panel_name")]


def condition_links(mondo_id, cache):
    """GARD registry + clinical trials for a condition, climbing one MONDO parent
    level when the exact term has none. Cached per MONDO id."""
    if not mondo_id:
        return {}
    if mondo_id in cache:
        return cache[mondo_id]
    gard = [g.get("id") for g in map_all(mondo_id, ">>mondo>>gard") if g.get("id")]
    trials = map_all(mondo_id, ">>mondo>>clinical_trials")
    if not trials:
        for par in map_all(mondo_id, ">>mondo>>mondoparent")[:2]:
            trials = map_all(par.get("id"), ">>mondo>>clinical_trials")
            if trials:
                break
    out = {"gard": gard[0] if gard else None, "trial_count": len(trials)}
    cache[mondo_id] = out
    return out


# ── Batch 3: per-variant derived (mostly in-memory, no new calls) ────────────
def structural_context(hgvs_p, intervals):
    """UniProt features overlapping the variant's residue → human phrases."""
    pos = protein_position(hgvs_p)
    if pos is None or not intervals:
        return None
    out = []
    for f in intervals:
        if f["begin"] <= pos <= f["end"]:
            phrase = _UF_KEEP[f["type"]].replace("{d}", f["desc"] or f["type"])
            out.append(phrase)
    return {"position": pos, "features": out} if out else None


def am_percentile(score, am_map):
    """The variant's AlphaMissense percentile within the gene's full modeled set
    ('top 3% most-pathogenic-predicted of N'), or None."""
    try:
        s = float(score)
    except (TypeError, ValueError):
        return None
    vals = []
    for _, sc in am_map.values():
        try:
            vals.append(float(sc))
        except (TypeError, ValueError):
            pass
    if len(vals) < 20:
        return None
    n_below = sum(1 for v in vals if v < s)
    pct = round(100 * (len(vals) - n_below) / len(vals))
    return {"top_pct": max(pct, 1), "n": len(vals)}


def similar_variants(rec, recs):
    """Up to 6 sibling variant pages most related to this one — same residue,
    then same condition, then same variant type. Internal-mesh aid for both
    personas + crawlers."""
    pos = protein_position(rec.get("hgvs_p"))
    conds = {c.get("mondo_id") for c in (rec.get("conditions") or [])}
    vtype = rec.get("variant_type")
    me = rec["canonical_slug"]
    scored = []
    for r in recs:
        if r["canonical_slug"] == me:
            continue
        s = 0
        if pos is not None and protein_position(r.get("hgvs_p")) == pos:
            s += 3
        if conds & {c.get("mondo_id") for c in (r.get("conditions") or [])}:
            s += 2
        if vtype and r.get("variant_type") == vtype:
            s += 1
        if s:
            scored.append((s, r))
    scored.sort(key=lambda x: (-x[0], x[1]["canonical_slug"]))
    return [{"label": f"{r['gene_symbol']} {r.get('hgvs_p') or r.get('hgvs_c')}",
             "slug": r["canonical_slug"]} for _, r in scored[:6]]


def submission_timeline(submissions):
    """First/last submission year + whether classifications diverged over time."""
    dated = [s for s in submissions if s.get("date")]
    years = sorted(int(s["date"][:4]) for s in dated if (s.get("date") or "")[:4].isdigit())
    if not years:
        return None
    calls = {(s.get("classification") or "").split("/")[0].strip().lower() for s in submissions}
    return {"first": years[0], "last": years[-1], "n": len(submissions),
            "stable": len(calls) <= 1}


def plain_summary(rec, digest):
    """Deterministic plain-language one-liner for patients (NOT an LLM)."""
    cls = (rec.get("classification") or "").lower()
    if "conflict" in cls:
        meaning = "of uncertain / conflicting significance (experts disagree)"
    elif "pathogenic" in cls:
        meaning = "considered disease-causing"
    else:
        meaning = f"classified as {rec.get('classification')}"
    gene = rec.get("gene_symbol")
    n = rec.get("submitter_count") or 0
    who = (f" reported by {n} clinical lab" + ("s" if n != 1 else "")) if n else ""
    cond = (digest or {}).get("name") or (rec.get("conditions") or [{}])[0].get("name")
    link = f", and is linked to {cond}" if cond else ""
    return f"This is a change in the {gene} gene that is {meaning}{who}{link}."


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
