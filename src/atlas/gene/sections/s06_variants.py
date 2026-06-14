"""§6 — variants: ClinVar (per-class breakdown + top pathogenic), SpliceAI, AlphaMissense, dbSNP sample."""
import re

from atlas.biobtree import entry, map_all, xref_counts
from atlas.gene.sections.base import Section


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

CHAINS = (
    '>>hgnc>>clinvar[germline_classification=="<class>"]',  # 5 classes
    ">>hgnc>>spliceai",
    '>>transcript>>alphamissense[am_class=="likely_pathogenic"]',
    ">>hgnc>>entrez>>dbsnp",  # via entrez (hgnc>>dbsnp empty — BIOBTREE_ISSUES.md)
)
DATASETS = ("clinvar", "spliceai", "alphamissense", "dbsnp", "entrez", "transcript",
            "hgnc", "clingen_variant")

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
    breakdown, top_path = {}, []
    for cls in classes:
        rs = map_all(a.hgnc_id, f'>>hgnc>>clinvar[germline_classification=="{cls}"]')
        breakdown[cls] = len(rs)
        for t in rs:
            _add_overlap(t.get("gene_symbol"))
        if cls in ("Pathogenic", "Likely pathogenic"):
            for t in rs:
                if len(top_path) >= 30:
                    break
                top_path.append({"id": t["id"], "hgvs": t.get("name"),
                                 "classification": t.get("germline_classification")})
    bundle["clinvar_breakdown"] = breakdown
    bundle["top_pathogenic"] = top_path

    sp = sorted(map_all(a.hgnc_id, ">>hgnc>>spliceai"),
                key=lambda t: float(t.get("score") or 0), reverse=True)
    for t in sp:
        _add_overlap(t.get("gene_symbol"))
    bundle["top_spliceai"] = [{"id": t["id"], "effect": t.get("effect"),
                               "score": t.get("score")} for t in sp[:30]]
    bundle["overlap_genes"] = sorted(overlap)

    ct = a.canonical_transcript
    bundle["canonical_transcript"] = ct
    if ct:
        bundle["alphamissense_total"] = xref_counts(entry(ct, "transcript")).get("alphamissense", 0)
        am = sorted(map_all(ct, '>>transcript>>alphamissense[am_class=="likely_pathogenic"]'),
                    key=lambda t: float(t.get("am_pathogenicity") or 0), reverse=True)
        bundle["top_alphamissense"] = [{"id": t["id"], "variant": t.get("protein_variant"),
                                        "am_pathogenicity": t.get("am_pathogenicity")} for t in am[:30]]

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

    # dbSNP rsIDs via ENTREZ (direct hgnc>>dbsnp unbacked; see BIOBTREE_ISSUES.md).
    dbs = map_all(a.hgnc_id, ">>hgnc>>entrez>>dbsnp", cap=2)
    bundle["dbsnp_sample"] = [{"id": t["id"], "pos": f"{t.get('chromosome')}:{t.get('position')}",
                               "change": f"{t.get('ref_allele')}>{t.get('alt_allele')}"} for t in dbs[:30]]
    bundle["dbsnp_sampled"] = len(dbs)
    return bundle

SECTION = Section(
    id="6", name="variants",
    description="ClinVar variants (per-class breakdown), SpliceAI splice impact, AlphaMissense pathogenicity, dbSNP sample",
    needs=("hgnc_id", "hgnc_entry", "canonical_transcript"),
    produces=("clinvar_total", "clinvar_breakdown", "top_pathogenic", "top_spliceai",
              "alphamissense_total", "top_alphamissense", "dbsnp_sample",
              "clingen_variant_total", "clingen_variant_breakdown",
              "clingen_variant_vceps", "clingen_variant_diseases", "overlap_genes"),
    datasets=DATASETS, chains=CHAINS, collect_fn=collect,
)
