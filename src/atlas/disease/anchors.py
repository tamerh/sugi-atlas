#!/usr/bin/env python3
"""Disease anchors — resolve a disease name (or Mondo id) to the ID set + gene
cohort the 14 disease sections will need, ONCE.

Disease anchor is more involved than gene's because:
- the disease has multiple ontology IDs (Mondo + EFO + MeSH + OMIM + Orphanet),
- the canonical *gene cohort* is the union of four routes
  (GWAS, GenCC Mendelian, ClinVar germline-rare, CIViC-evidence somatic),
- each cohort gene needs its own GeneAnchors resolved so disease §5–§12/§14
  can reuse the existing gene collectors via the cohort fan-out helper.

Cancer somatic data flows in through the CIViC-evidence route — the
2026-05-30 biobtree refresh shipped `civic_evidence` on the mondo edge, which
partially lifts the somatic-coverage blocker noted in
/data/biobtree-content/biobtree/disease/NOTES.md.
"""
import re, sys
from dataclasses import dataclass, field
from typing import Optional, Tuple, Dict, Set

from atlas.biobtree import search, entry, rows, map_all, xref_counts
from atlas.gene.anchors import Anchors as GeneAnchors, resolve as resolve_gene_anchors

# Heuristic — used when neither mondo definition nor parent ancestry helps.
_CANCER_RX = re.compile(
    r"\b(cancer|carcinoma|leukemia|leukaemia|lymphoma|sarcoma|"
    r"glioma|melanoma|myeloma|neoplasm|tumou?r|adenoma)\b", re.I)

# Cap for the per-disease gene cohort. The migrated pages all rendered
# "TOP 50" sections; we mirror that. Cohort genes feed gene-collector
# fan-outs (§5–§12/§14) so a tight cap bounds wall-clock cost
# (50 genes × ~12 collectors × ~5 chains each ≈ 3k biobtree calls per disease).
COHORT_CAP = 50

# Evidence-route definitions for the gene cohort union. Each entry:
# (flag_name, chain). Ranking later prefers genes hit by more routes.
_COHORT_ROUTES = (
    ("gwas",           ">>mondo>>gwas>>hgnc"),
    ("gencc",          ">>mondo>>gencc>>hgnc"),
    ("clinvar",        ">>mondo>>clinvar>>hgnc"),
    ("civic_evidence", ">>mondo>>civic_evidence>>hgnc"),
)

@dataclass(frozen=True)
class DiseaseAnchors:
    name: str                       # input name as given by the user
    mondo_id: str                   # MONDO:NNNNNNN
    mondo_entry: dict               # full /api/entry for traceability
    canonical_name: Optional[str]   # mondo's preferred name (often more formal)
    xref_counts: Dict[str, int]     # one-shot counts from mondo entry xrefs
    efo_id: Optional[str]
    mesh_ids: Tuple[str, ...]
    omim_ids: Tuple[str, ...]
    orphanet_ids: Tuple[str, ...]
    is_cancer: bool
    # Gene cohort — top COHORT_CAP genes by evidence-route count.
    cohort: Tuple[GeneAnchors, ...]
    # Parallel mapping hgnc_id -> evidence-flags dict {route: bool}
    cohort_evidence: Dict[str, Dict[str, bool]]
    # Full cohort hgnc ids before capping (useful for §4 Mendelian-overlap stats
    # and "X of Y genes in cohort" counts).
    cohort_full: Tuple[str, ...]
    # GWAS study count (top-level, from mondo xrefs — cheap signal for §2).
    gwas_study_count: int

def _is_cancer(mondo_id: str, mondo_entry: dict, canonical_name: Optional[str]) -> bool:
    """Cheap is-cancer heuristic. Mondo's parent ancestry would be more rigorous
    but isn't surfaced via single-call traversal; the name regex catches the
    common cancer terms (carcinoma, leukemia, sarcoma, ...). Refine when biobtree
    exposes mondo ancestry."""
    if canonical_name and _CANCER_RX.search(canonical_name):
        return True
    # Definition string sometimes has 'malignant neoplasm of ...' — check it too.
    attrs = (mondo_entry.get("Attributes") or {}).get("Mondo") or {}
    definition = attrs.get("definition") or ""
    if _CANCER_RX.search(definition):
        return True
    return False

def resolve_mondo(name_or_id: str) -> Tuple[str, dict, Optional[str]]:
    """name_or_id -> (mondo_id, full_entry, canonical_name). If the input is
    already a MONDO:NNNN id we fetch the entry directly; otherwise we search
    the mondo dataset and pick the highest-xref-count match (the canonical
    Mondo node beats deprecated/obsolete duplicates on xref count)."""
    if re.match(r"^MONDO:\d+$", name_or_id):
        en = entry(name_or_id, "mondo")
        canonical = (((en.get("Attributes") or {}).get("Mondo") or {}).get("name"))
        return name_or_id, en, canonical

    resp = search(name_or_id, source="mondo")
    cand = [r for r in rows(resp) if r.get("id", "").startswith("MONDO:")]
    if not cand:
        # Fallback — unfiltered search, take any MONDO row.
        cand = [r for r in rows(search(name_or_id))
                if r.get("id", "").startswith("MONDO:")]
    if not cand:
        sys.exit(f"no MONDO row for {name_or_id!r}")
    # Prefer non-obsolete (the canonical row has the highest xref_count by far,
    # the search results already rank by relevance + xref_count).
    cand.sort(key=lambda r: int(r.get("xref_count") or 0), reverse=True)
    chosen = cand[0]
    en = entry(chosen["id"], "mondo")
    canonical = (((en.get("Attributes") or {}).get("Mondo") or {}).get("name")
                 or chosen.get("name"))
    return chosen["id"], en, canonical

def _build_cohort(mondo_id: str):
    """Union the four evidence routes into a per-gene evidence map.
    Returns (cohort_hgnc_ids_ranked, evidence_map).

    Ranking: descending by number of evidence routes (gwas+gencc+clinvar+civic),
    ties broken by the order returned by biobtree (already relevance-ranked
    for GWAS routes at least). This puts dual/triple-evidence genes first so
    the COHORT_CAP retains the strongest-evidence subset.
    """
    evidence: Dict[str, Dict[str, bool]] = {}
    for flag, chain in _COHORT_ROUTES:
        for t in map_all(mondo_id, chain):
            h = t.get("id")
            if not h or not h.startswith("HGNC:"):
                continue
            evidence.setdefault(h, {f: False for f, _ in _COHORT_ROUTES})[flag] = True

    def _score(hgnc):
        return sum(1 for v in evidence[hgnc].values() if v)

    ranked = sorted(evidence.keys(), key=lambda h: (-_score(h), h))
    return ranked, evidence

def resolve(name_or_id: str) -> DiseaseAnchors:
    """Disease name (or Mondo id) → DiseaseAnchors.

    Resolution steps:
      1. name → mondo_id + entry (search if name; direct fetch if Mondo id).
      2. xref_counts from the mondo entry (the same shortcut gene anchors use
         on hgnc — exact totals for gwas / gwas_study / clinical_trials /
         clinvar / antibody / cellosaurus / civic_evidence in a single call).
      3. One-hop sibling ontology ids (efo / mesh / omim / orphanet) from the
         mondo entry's xrefs section. Surfaced as explicit fields because §1
         renders them as a federated-identifier table.
      4. Union the 4-route gene cohort (gwas / gencc / clinvar / civic_evidence),
         rank by evidence breadth, cap to COHORT_CAP, pre-resolve each gene to
         a GeneAnchors so disease §5–§12/§14 reuse existing gene collectors
         without re-paying anchor cost.
    """
    mondo_id, mondo_entry, canonical_name = resolve_mondo(name_or_id)
    xc = xref_counts(mondo_entry)

    # Sibling ontology ids — the mondo entry's xrefs section lists per-dataset
    # *counts* but not the actual IDs. So we call the cheap per-hop maps to
    # capture the IDs themselves; both are bounded (each disease has at most
    # a handful of cross-ontology mappings).
    efo_rows  = map_all(mondo_id, ">>mondo>>efo")
    mesh_rows = map_all(mondo_id, ">>mondo>>mesh")
    mim_rows  = map_all(mondo_id, ">>mondo>>mim")
    orph_rows = map_all(mondo_id, ">>mondo>>orphanet")

    cohort_full, evidence = _build_cohort(mondo_id)

    # Pre-resolve the top-N gene anchors so downstream cohort sections don't
    # re-pay anchor cost. Each gene resolve is ~4 biobtree calls so this
    # adds ~200 calls upfront for a 50-gene cap — once, not per section.
    cohort: list = []
    for hgnc in cohort_full[:COHORT_CAP]:
        # resolve_gene_anchors takes a symbol; pull it from the hgnc entry.
        # We already have the hgnc id; fetch its entry once to get the symbol.
        try:
            he = entry(hgnc, "hgnc")
            symbol = ((he.get("Attributes") or {}).get("Hgnc") or {}).get("symbols", [None])[0]
            if not symbol:
                continue
            cohort.append(resolve_gene_anchors(symbol))
        except Exception:
            # A single gene anchor failure (rare HGNC, retired symbol) shouldn't
            # block the whole disease build — skip and continue.
            continue

    return DiseaseAnchors(
        name=name_or_id,
        mondo_id=mondo_id,
        mondo_entry=mondo_entry,
        canonical_name=canonical_name,
        xref_counts=xc,
        efo_id=(efo_rows[0]["id"] if efo_rows else None),
        mesh_ids=tuple(r["id"] for r in mesh_rows),
        omim_ids=tuple(r["id"] for r in mim_rows),
        orphanet_ids=tuple(r["id"] for r in orph_rows),
        is_cancer=_is_cancer(mondo_id, mondo_entry, canonical_name),
        cohort=tuple(cohort),
        cohort_evidence=evidence,
        cohort_full=tuple(cohort_full),
        gwas_study_count=int(xc.get("gwas_study") or 0),
    )
