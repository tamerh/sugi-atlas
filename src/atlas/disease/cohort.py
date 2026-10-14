"""Cohort fan-out helper — runs a gene section's `collect_fn(GeneAnchors)`
over every gene in the disease cohort and returns the bundles per gene.

Disease sections §5–§12 and §14 don't define new biobtree chains; they
*aggregate* over the per-gene bundles produced by the existing gene
collectors. This module is the single place where the fan-out happens so the
disease section bodies stay tiny ("call cohort.fan(gene_section, anchors),
then aggregate") and we avoid re-implementing per-gene logic on the disease
side.

Cost: serial today (50 genes × N chains). Threaded fan-out would shrink wall
clock at the cost of contention on biobtree's process-local cache —
introduce only if a section becomes the bottleneck.
"""
from typing import Callable, List, Optional, Dict, Any, Tuple
from atlas.gene.anchors import Anchors as GeneAnchors
from atlas.biobtree import map_all

import re

# Generic descriptors that would over-match when comparing two disease names
# (every cancer shares "cancer", "carcinoma"…). Dropping them keeps the
# substantive tokens that actually identify the disease.
_GENCC_STOP = {"disease", "diseases", "syndrome", "cancer", "carcinoma", "tumor",
               "tumour", "neoplasm", "disorder", "disorders", "type", "familial",
               "hereditary", "susceptibility", "predisposition", "complementation",
               "group", "deficiency", "autosomal", "dominant", "recessive",
               "with", "without"}


def disease_tokens(s):
    """Substantive tokens of a disease name for on-disease matching (drops the
    generic descriptors that would over-match). Shared by the GenCC dedup and
    the dual-evidence on-disease filter so both use identical matching."""
    return {t for t in re.findall(r"[a-z0-9]+", (s or "").lower())
            if len(t) >= 4 and t not in _GENCC_STOP}


# Definitive(6) + Strong(5) — the GenCC classifications strong enough to call a
# gene "causal". Moderate/Limited stay out of the causal claim (still in §4).
_CAUSAL_GENCC_FLOOR = 5


# Chromosomal (contiguous-gene) syndromes: a ClinVar record on a 2q13 deletion
# lists every gene in the region, so "caused by RGPD8" is an artefact of how the
# CNV was submitted, not a causal gene.
_CHROMOSOMAL = re.compile(
    r"micro(?:deletion|duplication)|\b(?:deletion|duplication|trisomy|monosomy|"
    r"tetrasomy|triploidy|uniparental disomy|ring chromosome)\b|"
    r"\b(?:chromosome )?(?:[1-9]|1[0-9]|2[0-2]|x|y)[pq]\d", re.I)


def is_chromosomal_anomaly(name):
    return bool(_CHROMOSOMAL.search(name or ""))


# Antisense / divergent / intronic transcripts overlapping a coding gene pick up
# its ClinVar records ("CFTR-related disorder … CFTR, CFTR-AS1, and CFTR-AS2").
_NONCODING_NEIGHBOR = re.compile(r"-(AS|DT|IT|OT)\d*$|^(MIR\d|LINC\d|SNOR[AD]\d)|^SNHG\d")


def _clinvar_tier_applies(bundle):
    """The ClinVar-only causal tier is for gene-defined rare diseases. Not for a
    chromosomal syndrome (region CNVs), a common/complex trait (≥20 GWAS
    associations: asthma → CARD11, MDD → TDO2 were ClinVar artefacts), or an
    umbrella term (its genes belong to the subtypes)."""
    b1 = bundle.get("1") or {}
    if is_chromosomal_anomaly(b1.get("canonical_name") or b1.get("name")):
        return False
    # Cancers: ClinVar germline records on a tumour type are predisposition or
    # somatic noise, not "the" gene (AML with maturation → MIR181A1HG).
    if b1.get("is_cancer"):
        return False
    # A disease with no OMIM / Orphanet record has no Mendelian anchor — the
    # common/acquired tail (first-degree AV block → KCNJ16).
    if not (b1.get("omim_ids") or b1.get("orphanet_ids")):
        return False
    if ((bundle.get("2") or {}).get("assoc_total") or 0) >= 20:
        return False
    if (b1.get("child_count") or 0) >= 5:
        return False
    return True


def causal_genes(bundle):
    """High-confidence causal gene(s) for a disease page, strongest-first:
    genes GenCC-curated Definitive/Strong for THIS disease's exact Mondo id, or
    an OMIM Mendelian overlap (gene MIM ∩ disease MIM). Returns [(symbol, label)];
    empty for polygenic / GWAS-only diseases (most cancers, common disease).

    Uses the exact-id §4 `disease_gencc` (>>mondo>>gencc), NOT disease-name
    token matching — so a neonatal-diabetes / MODY gene whose GenCC record merely
    contains "diabetes" never gets claimed as causing type-2 diabetes. A separate,
    purely additive signal for the lead + at-a-glance; no cohort table changes."""
    from atlas.render_common import gencc_rank
    b4 = bundle.get("4") or {}

    best = {}                                # symbol -> (rank, classification)
    for g in (b4.get("disease_gencc") or []):
        r = gencc_rank(g.get("classification"))
        if r < _CAUSAL_GENCC_FLOOR:          # Definitive/Strong only
            continue
        sym = g.get("symbol")
        if sym and r > best.get(sym, (0, ""))[0]:
            best[sym] = (r, g.get("classification"))

    out = [(sym, f"GenCC {cls}") for sym, (r, cls) in best.items()]
    seen = set(best)
    # OMIM Mendelian overlap — strong causal signal even without a GenCC record.
    for g in (b4.get("omim_genes") or []):
        sym = g.get("symbol")
        if sym and sym not in seen:
            out.append((sym, "OMIM Mendelian"))
            seen.add(sym)

    out.sort(key=lambda t: (-best.get(t[0], (0,))[0], t[0]))  # rank desc, then symbol

    # ClinVar fallback tier — for a gene-defined Mendelian disease with no GenCC/
    # OMIM record (most rare ClinVar-only syndromes), the genes that reached the
    # cohort via the ClinVar germline route FOR THIS DISEASE are the causal
    # candidates. Bounded to a focused set (≤5) so a polygenic / cancer cohort with
    # many ClinVar genes never claims them all causal; labelled weaker than GenCC.
    # Orphanet tier — "Disease-causing germline mutation(s) in" associations for a
    # disease whose cohort came only from Orphanet (anchors fallback).
    if not out:
        b5 = bundle.get("5") or {}
        orph = sorted(g.get("symbol") for g in (b5.get("genes") or [])
                      if g.get("symbol") and (g.get("evidence") or {}).get("orphanet_causing"))
        if 0 < len(orph) <= 5:
            out = [(s, "Orphanet disease-causing") for s in orph]
    if not out and _clinvar_tier_applies(bundle):
        b5 = bundle.get("5") or {}
        # ≥2 (likely) pathogenic non-CNV records for THIS disease (v1.12.2 audit:
        # 16/25 correct with any single record), coding genes only.
        plp = (bundle.get("3") or {}).get("clinvar_plp_by_gene")
        clinvar_genes = sorted(
            g.get("symbol") for g in (b5.get("genes") or [])
            if g.get("symbol") and (g.get("evidence") or {}).get("clinvar")
            and not _NONCODING_NEIGHBOR.search(g["symbol"])
            and (plp is None or plp.get(g["symbol"], 0) >= 2))
        if 0 < len(clinvar_genes) <= 5:
            out = [(s, "ClinVar-linked") for s in clinvar_genes]
    return out


def enrichment_fan(enrichment_cohort: "tuple[tuple[str, str], ...]",
                   chain: str, cap: int = 10) -> List[Tuple[str, str, list]]:
    """Run ONE cheap chain over the WIDE enrichment cohort (hgnc_id, symbol
    pairs), returning [(hgnc_id, symbol, rows), ...]. This is the breadth track:
    aggregate-only sections (pathway enrichment, druggability) fan a single
    chain over ~250 genes for sharper statistics, instead of the full per-gene
    gene-plan depth the display cohort gets. One bad gene → empty rows, never
    fatal."""
    out: List[Tuple[str, str, list]] = []
    for hgnc, sym in enrichment_cohort:
        try:
            rows = map_all(hgnc, chain, cap=cap)
        except Exception:
            rows = []
        out.append((hgnc, sym, rows))
    return out


def fan(collect_fn: Callable[[GeneAnchors], dict],
        cohort: "tuple[GeneAnchors, ...]",
        on_error: str = "skip") -> List[dict]:
    """Run `collect_fn(gene_anchors)` for every gene in the cohort.

    Returns one bundle per gene, in the same order as `cohort`. The bundle
    includes the gene's symbol + hgnc_id at minimum (assumed by every gene
    collector), so disease sections can group/sort/filter without a second
    lookup.

    on_error:
      'skip'  — log and drop the gene from the result (default; one bad gene
                shouldn't fail a 50-gene fan-out).
      'raise' — propagate; for debugging.
    """
    out: List[dict] = []
    for ga in cohort:
        try:
            b = collect_fn(ga)
            # Defensive — make sure symbol/hgnc_id always present in result.
            b.setdefault("symbol", ga.symbol)
            b.setdefault("hgnc_id", ga.hgnc_id)
            out.append(b)
        except Exception as e:
            if on_error == "raise":
                raise
            # Keep going — partial coverage beats a hard failure on rare HGNC.
            out.append({"symbol": ga.symbol, "hgnc_id": ga.hgnc_id,
                        "_error": str(e)})
    return out


def by_symbol(bundles: List[dict]) -> Dict[str, dict]:
    """Index a fan-out result by gene symbol for fast lookup in render code."""
    return {b["symbol"]: b for b in bundles if "symbol" in b}


def filter_evidence(cohort: "tuple[GeneAnchors, ...]",
                    cohort_evidence: Dict[str, Dict[str, bool]],
                    route: str) -> "tuple[GeneAnchors, ...]":
    """Subset the cohort to genes hit by a specific evidence route
    ('gwas' / 'gencc' / 'clinvar' / 'civic_evidence'). Used by §4 Mendelian
    overlap and §16 druggability slicing."""
    return tuple(g for g in cohort
                 if (cohort_evidence.get(g.hgnc_id) or {}).get(route))
