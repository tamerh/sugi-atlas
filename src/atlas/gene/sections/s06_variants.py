"""§6 — variants: ClinVar (per-class breakdown + best-reviewed pathogenic), ClinGen
expert panels, SpliceAI, AlphaMissense (counts + hotspot residues).

Per-variant depth (every rsID / every AlphaMissense score) is Sugi Variant's job;
this section summarises. The dbSNP sample table was dropped in v1.11.8: an
arbitrary ~300-rsID slice (8% of gene-page bytes) that missed the famous SNPs."""
import re

from atlas.biobtree import entry, map_all, xref_counts
from atlas.gene.sections.base import Section
from atlas.render_common import clinvar_stars, clinvar_id_num

TOP_PATHOGENIC_CAP = 15
SPLICEAI_CAP = 10
AM_HOTSPOT_CAP = 10


def _dedup_disease_names(names):
    """Collapse condition labels that differ only by a trailing numeric subtype
    suffix — 'Li-Fraumeni syndrome' subsumes 'Li-Fraumeni syndrome 1' (two
    ontology granularities for the same entity). Keep the shortest representative
    per base; order preserved by base name."""
    by_base = {}
    for n in names:
        n = (n or "").strip()
        if not n:
            continue
        base = re.sub(r"\s+\d+$", "", n).lower()      # drop a trailing " 1"/" 2"…
        if base not in by_base or len(n) < len(by_base[base]):
            by_base[base] = n
    return sorted(by_base.values())

_PROT_VAR = re.compile(r"^(?:p\.)?([A-Z])(\d+)([A-Z])$")


def am_hotspots(rows):
    """Group AlphaMissense likely-pathogenic substitutions by residue → hotspot
    list [{residue: 'R175', n: 19, max: 1.0}], most-intolerant first (substitution
    count, then max score, then position). Rows whose protein_variant isn't a
    simple missense ('R175H') are skipped."""
    by = {}
    for t in rows:
        m = _PROT_VAR.match((t.get("protein_variant") or "").strip())
        if not m:
            continue
        ref, pos = m.group(1), int(m.group(2))
        try:
            sc = float(t.get("am_pathogenicity"))
        except (TypeError, ValueError):
            sc = 0.0
        h = by.setdefault(pos, {"residue": f"{ref}{pos}", "n": 0, "max": 0.0, "_pos": pos})
        h["n"] += 1
        h["max"] = max(h["max"], sc)
    out = sorted(by.values(), key=lambda h: (-h["n"], -h["max"], h["_pos"]))
    for h in out:
        h.pop("_pos")
    return out


CHAINS = (
    '>>hgnc>>clinvar[germline_classification=="<class>"]',  # 5 classes
    ">>hgnc>>spliceai",
    '>>transcript>>alphamissense[am_class=="likely_pathogenic"]',
)
DATASETS = ("clinvar", "spliceai", "alphamissense", "transcript", "hgnc", "clingen_variant")

def collect(a):
    bundle = {"section": "06_variants", "symbol": a.symbol, "hgnc_id": a.hgnc_id}
    xc = xref_counts(a.hgnc_entry)
    bundle["clinvar_total"] = xc.get("clinvar", 0)
    bundle["spliceai_total"] = xc.get("spliceai", 0)

    # ClinVar/SpliceAI records carry the gene_symbol they map to. For a non-coding
    # gene these positional records belong to the OVERLAPPING protein-coding
    # gene(s) (e.g. lncRNA CAHM's variants are QKI's). Capture those distinct
    # symbols (minus self) so the non-coding render can orient the reader to the
    # overlapping gene. Survives the non-coding scrub (not a scrubbed key).
    self_sym = (a.symbol or "").upper()
    overlap = set()
    def _add_overlap(field):
        for g in re.split(r"[;,]", field or ""):
            g = g.strip()
            if g and g.upper() != self_sym:
                overlap.add(g)

    classes = ["Pathogenic", "Likely pathogenic", "Uncertain significance",
               "Likely benign", "Benign"]
    breakdown, plp = {}, []
    for cls in classes:
        rs = map_all(a.hgnc_id, f'>>hgnc>>clinvar[germline_classification=="{cls}"]')
        breakdown[cls] = len(rs)
        for t in rs:
            _add_overlap(t.get("gene_symbol"))
        if cls in ("Pathogenic", "Likely pathogenic"):
            plp += [{"id": t["id"], "hgvs": t.get("name"),
                     "classification": t.get("germline_classification"),
                     "review_status": t.get("review_status")} for t in rs]
    bundle["clinvar_breakdown"] = breakdown
    # Best-reviewed first (ClinVar stars), Pathogenic before Likely pathogenic,
    # then earliest ClinVar id (established hallmark variants before recent ones). Capped — the full set is Sugi Variant's.
    plp.sort(key=lambda v: (-clinvar_stars(v["review_status"]),
                            v["classification"] != "Pathogenic", clinvar_id_num(v["id"])))
    bundle["top_pathogenic_total"] = len(plp)
    bundle["top_pathogenic"] = plp[:TOP_PATHOGENIC_CAP]

    sp = sorted(map_all(a.hgnc_id, ">>hgnc>>spliceai"),
                key=lambda t: float(t.get("score") or 0), reverse=True)
    for t in sp:
        _add_overlap(t.get("gene_symbol"))
    bundle["top_spliceai"] = [{"id": t["id"], "effect": t.get("effect"),
                               "score": t.get("score")} for t in sp[:SPLICEAI_CAP]]
    bundle["overlap_genes"] = sorted(overlap)

    ct = a.canonical_transcript
    bundle["canonical_transcript"] = ct
    if ct:
        bundle["alphamissense_total"] = xref_counts(entry(ct, "transcript")).get("alphamissense", 0)
        am = sorted(map_all(ct, '>>transcript>>alphamissense[am_class=="likely_pathogenic"]'),
                    key=lambda t: float(t.get("am_pathogenicity") or 0), reverse=True)
        bundle["am_lp_total"] = len(am)
        hot = am_hotspots(am)
        bundle["am_lp_residues"] = len(hot)
        bundle["am_hotspots"] = hot[:AM_HOTSPOT_CAP]

    # ClinGen VCEP expert-panel interpretations — ACMG calls reviewed by a
    # Variant Curation Expert Panel, a higher authority tier than raw ClinVar
    # submissions (provenance already advertised clingen_variant; this is the
    # collector that was missing). Schema: id|gene_symbol|disease|assertion|vcep.
    from collections import Counter
    cg = map_all(a.hgnc_id, ">>hgnc>>clingen_variant")
    if cg:
        order = ["Pathogenic", "Likely Pathogenic", "Uncertain Significance",
                 "Likely Benign", "Benign"]
        cnt = Counter((r.get("assertion") or "").strip() for r in cg if r.get("assertion"))
        bundle["clingen_variant_total"] = len(cg)
        bundle["clingen_variant_breakdown"] = (
            [(k, cnt[k]) for k in order if cnt.get(k)]
            + [(k, n) for k, n in cnt.items() if k not in order])
        bundle["clingen_variant_vceps"] = sorted(
            {(r.get("vcep") or "").strip() for r in cg if r.get("vcep")})
        bundle["clingen_variant_diseases"] = _dedup_disease_names(
            r.get("disease") for r in cg if r.get("disease"))

    return bundle

SECTION = Section(
    id="6", name="variants",
    description="ClinVar variants (per-class breakdown), SpliceAI splice impact, AlphaMissense hotspot residues",
    needs=("hgnc_id", "hgnc_entry", "canonical_transcript"),
    produces=("clinvar_total", "clinvar_breakdown", "top_pathogenic", "top_pathogenic_total",
              "top_spliceai", "alphamissense_total", "am_lp_total", "am_lp_residues", "am_hotspots",
              "clingen_variant_total", "clingen_variant_breakdown",
              "clingen_variant_vceps", "clingen_variant_diseases", "overlap_genes"),
    datasets=DATASETS, chains=CHAINS, collect_fn=collect,
)
