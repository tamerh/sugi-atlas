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

def _name_overlap(a, b):
    """Token-overlap of two disease names (significant words, length ≥4), as a
    fraction of the shorter name's tokens. Order-insensitive, so "Chondrodysplasia
    punctata, brachytelephalangic, autosomal" and Orphanet's "Brachytelephalangic
    chondrodysplasia punctata" score 1.0 — but an unrelated sibling form scores
    low. Used to gate the Orphanet-via-OMIM fallback."""
    def toks(s):
        return {w for w in re.findall(r"[a-z0-9]+", (s or "").lower()) if len(w) >= 4}
    ta, tb = toks(a), toks(b)
    if not ta or not tb:
        return 0.0
    return len(ta & tb) / min(len(ta), len(tb))


# Heuristic — used when neither mondo definition nor parent ancestry helps.
_CANCER_RX = re.compile(
    r"\b(cancer|carcinoma|leukemia|leukaemia|lymphoma|sarcoma|"
    r"glioma|melanoma|myeloma|neoplasm|tumou?r|adenoma)\b", re.I)

# Cap for the per-disease gene cohort. The migrated pages all rendered
# "TOP 50" sections; we mirror that. Cohort genes feed gene-collector
# fan-outs (§5–§12/§14) so a tight cap bounds wall-clock cost
# (50 genes × ~12 collectors × ~5 chains each ≈ 3k biobtree calls per disease).
# Cohort selection. The fan-out reuses the full gene plan per gene, so cohort
# size is the dominant cost of a disease build — but a flat top-N both drops
# strong genes on well-studied diseases and admits GWAS-tail noise on others.
# So: an EVIDENCE FLOOR keeps every gene with curated evidence (GenCC/CIViC) or
# >=2 routes regardless of the cap (they're never noise), then COHORT_CAP fills
# the rest from the single-route tail; COHORT_MAX is a hard ceiling so a
# pathological disease can't fan unbounded.
COHORT_CAP = 75
COHORT_MAX = 150

# Enrichment cohort — a WIDER, evidence-ranked gene set used only to feed the
# aggregate sections (pathway enrichment, druggability breadth). It is fanned
# with one cheap chain per gene (not the full per-gene plan), so "wide" is
# cheap; the bound is about signal (the GWAS tail dilutes enrichment), not cost.
ENRICH_CAP = 250

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
    synonyms: Tuple[str, ...]       # Mondo Ontology synonyms (e.g. "AT", "Lou Gehrig's…")
    xref_counts: Dict[str, int]     # one-shot counts from mondo entry xrefs
    efo_id: Optional[str]
    mesh_ids: Tuple[str, ...]
    omim_ids: Tuple[str, ...]
    orphanet_ids: Tuple[str, ...]
    # Additional cross-ontology IDs from Mondo's OBO — exposed by biobtree
    # 2026-06-01. Keys: doid, sctid, umls, ncit, medgen, icd10cm, icd11,
    # gard, meddra, nord. Per-disease coverage varies (rare diseases get
    # the full set; cancers get most; common diseases get less).
    obo_xrefs: Dict[str, Tuple[str, ...]]
    # UBERON anatomy from Mondo's `disease_has_location` axiom — surfaces
    # the affected tissue (e.g. UBERON:0000310 = "breast" for breast cancer).
    # Drives schema.org `associatedAnatomy`. Empty for systemic diseases.
    anatomy_uberon_ids: Tuple[str, ...]
    # Primary Orphanet entry's full attrs — used for clinical phenotypes
    # (HPO list with frequencies) and prevalences (epidemiology).
    # Picked as the first orphanet_id whose disorder_type is "Disease"
    # (vs "Clinical subtype" / "Group"). {} for non-rare-disease entries.
    orphanet_attrs: dict
    is_cancer: bool
    # Gene cohort — all strong-evidence genes (evidence floor) + single-route
    # tail filled to COHORT_CAP, bounded by COHORT_MAX (see _select_cohort).
    cohort: Tuple[GeneAnchors, ...]
    # Parallel mapping hgnc_id -> evidence-flags dict {route: bool}
    cohort_evidence: Dict[str, Dict[str, bool]]
    # Full cohort hgnc ids before capping (useful for §4 Mendelian-overlap stats
    # and "X of Y genes in cohort" counts).
    cohort_full: Tuple[str, ...]
    # Enrichment cohort — (hgnc_id, symbol) for the top ENRICH_CAP evidence-ranked
    # genes. NOT resolved to GeneAnchors; aggregate sections fan one cheap chain
    # over it for breadth (pathway enrichment, druggability) without per-gene depth.
    enrichment_cohort: Tuple[Tuple[str, str], ...]
    # GWAS study count (top-level, from mondo xrefs — cheap signal for §2).
    gwas_study_count: int

def _mondo_attrs(mondo_entry: dict) -> dict:
    """biobtree's mondo entry now exposes name/synonyms/definition under
    Attributes.Ontology; older payloads used Attributes.Mondo. Accept either
    so the resolve survives the schema rename."""
    attrs = (mondo_entry.get("Attributes") or {})
    return attrs.get("Ontology") or attrs.get("Mondo") or {}


def _is_cancer(mondo_id: str, mondo_entry: dict, canonical_name: Optional[str]) -> bool:
    """Cheap is-cancer heuristic. Mondo's parent ancestry would be more rigorous
    but isn't surfaced via single-call traversal; the name regex catches the
    common cancer terms (carcinoma, leukemia, sarcoma, ...). Refine when biobtree
    exposes mondo ancestry."""
    if canonical_name and _CANCER_RX.search(canonical_name):
        return True
    # Definition string sometimes has 'malignant neoplasm of ...' — check it too.
    definition = _mondo_attrs(mondo_entry).get("definition") or ""
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
        canonical = _mondo_attrs(en).get("name")
        if not canonical:
            # Defensive fallback — biobtree's name-by-id might fail (deprecated
            # term, schema in flux). Search by the id to recover the canonical
            # name from the index (search returns name even when entry doesn't).
            for r in rows(search(name_or_id, source="mondo")):
                if r.get("id") == name_or_id:
                    canonical = r.get("name")
                    break
        return name_or_id, en, canonical

    resp = search(name_or_id, source="mondo")
    cand = [r for r in rows(resp) if r.get("id", "").startswith("MONDO:")]
    if not cand:
        # Fallback — unfiltered search, take any MONDO row.
        cand = [r for r in rows(search(name_or_id))
                if r.get("id", "").startswith("MONDO:")]
    if not cand:
        raise LookupError(f"no MONDO row for {name_or_id!r}")
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


def _is_strong(ev_row: Dict[str, bool]) -> bool:
    """Evidence floor: a gene is 'strong' if it carries curated evidence
    (GenCC or CIViC) or is hit by >=2 routes. Strong genes are never dropped by
    the cap; the single-route GWAS/ClinVar tail is what the cap limits."""
    return bool(ev_row.get("gencc") or ev_row.get("civic_evidence")
                or sum(1 for v in ev_row.values() if v) >= 2)


def _select_cohort(ranked, evidence, cap=COHORT_CAP, hard_max=COHORT_MAX):
    """Pick the cohort to fan out from the ranked union. Keep ALL strong genes
    (evidence floor, even past `cap`), then fill remaining slots up to `cap` from
    the single-route tail in rank order; bounded by `hard_max`. `ranked` is
    already route-count-descending, so order is preserved within each group."""
    strong = [h for h in ranked if _is_strong(evidence[h])]
    tail   = [h for h in ranked if not _is_strong(evidence[h])]
    return (strong + tail[:max(0, cap - len(strong))])[:hard_max]

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

    # Orphanet-via-OMIM fallback (BIOBTREE_ISSUES #41): a Mondo node often lacks
    # the Orphanet xref while the SAME disease's Orphanet record — the rich HPO
    # phenotype + prevalence bundle — is reachable through its OMIM id. Pull those
    # candidates when the direct route is empty. The OMIM→Orphanet hop can land on
    # a sibling form, so these candidates are name-gated below before their attrs
    # are trusted (the direct >>mondo>>orphanet rows are curator-asserted, so they
    # are not gated).
    orph_from_omim = not orph_rows and bool(mim_rows)
    if orph_from_omim:
        seen = set()
        for mr in mim_rows:
            for r in map_all(mr["id"], ">>mim>>orphanet"):
                if r.get("id") and r["id"] not in seen:
                    seen.add(r["id"])
                    orph_rows.append(r)

    # Cross-ontology xrefs from biobtree's Mondo OBO ingest (2026-06-01).
    # Each is typically 0-1 row per Mondo term; total <30 calls per disease.
    # Only fetched when xref_counts says the dataset has at least one row,
    # avoiding ~half the round-trips for sparser terms.
    OBO_XREF_DATASETS = ("doid", "sctid", "umls", "ncit", "medgen",
                          "icd10cm", "icd11", "gard", "meddra", "nord")
    obo_xrefs = {}
    for ds in OBO_XREF_DATASETS:
        if xc.get(ds, 0) > 0:
            ids = tuple(r["id"] for r in map_all(mondo_id, f">>mondo>>{ds}"))
            if ids:
                obo_xrefs[ds] = ids

    # UBERON anatomy — drives schema.org `associatedAnatomy` in the JSON-LD.
    anatomy_uberon_ids = tuple(r["id"] for r in
                                map_all(mondo_id, ">>mondo>>uberon")) if xc.get("uberon", 0) > 0 else ()

    # Pick the primary Orphanet entry (disorder_type=="Disease"). Pull its
    # full attrs once — carries the per-disease HPO phenotype list with
    # frequencies + multi-region prevalence data. Both drive new content
    # in §1 + JSON-LD (signOrSymptom, epidemiology).
    def _name_ok(attrs):
        # Only gate the OMIM-fallback candidates: their name must overlap the
        # disease's so we never attach a sibling form's phenotypes. Direct Mondo
        # xrefs are curator-asserted → always accepted.
        if not orph_from_omim:
            return True
        return _name_overlap(canonical_name, attrs.get("name")) >= 0.6

    orphanet_attrs = {}
    for r in orph_rows:
        try:
            e = entry(r["id"], "orphanet")
            attrs = (e.get("Attributes") or {}).get("Orphanet") or {}
            if attrs.get("disorder_type") == "Disease" and _name_ok(attrs):
                orphanet_attrs = attrs
                break
        except Exception:
            continue
    # Fallback: if no "Disease" found but there are entries, take the first
    # name-accepted subtype's attrs — better than empty for terms Mondo only
    # links to a subtype.
    if not orphanet_attrs and orph_rows:
        for r in orph_rows:
            try:
                attrs = (entry(r["id"], "orphanet").get("Attributes") or {}).get("Orphanet") or {}
            except Exception:
                continue
            if _name_ok(attrs):
                orphanet_attrs = attrs
                break

    cohort_full, evidence = _build_cohort(mondo_id)

    # Resolve the symbol for the top ENRICH_CAP evidence-ranked genes ONCE
    # (one hgnc entry each). This serves both the wide enrichment cohort and the
    # deep display cohort (a subset), so symbols aren't resolved twice.
    symbol_of: Dict[str, str] = {}
    for hgnc in cohort_full[:ENRICH_CAP]:
        try:
            he = entry(hgnc, "hgnc")
            s = ((he.get("Attributes") or {}).get("Hgnc") or {}).get("symbols", [None])[0]
            if s:
                symbol_of[hgnc] = s
        except Exception:
            continue
    enrichment_cohort = tuple((h, symbol_of[h]) for h in cohort_full[:ENRICH_CAP]
                              if h in symbol_of)

    # Pre-resolve the DISPLAY gene anchors (full per-gene plan) so downstream
    # per-gene sections don't re-pay anchor cost. Selection applies the evidence
    # floor + cap (see _select_cohort); symbols come from the map built above.
    cohort: list = []
    for hgnc in _select_cohort(cohort_full, evidence):
        symbol = symbol_of.get(hgnc)
        if not symbol:
            continue
        try:
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
        synonyms=tuple(_mondo_attrs(mondo_entry).get("synonyms") or ()),
        xref_counts=xc,
        efo_id=(efo_rows[0]["id"] if efo_rows else None),
        mesh_ids=tuple(r["id"] for r in mesh_rows),
        omim_ids=tuple(r["id"] for r in mim_rows),
        orphanet_ids=tuple(r["id"] for r in orph_rows),
        obo_xrefs=obo_xrefs,
        anatomy_uberon_ids=anatomy_uberon_ids,
        orphanet_attrs=orphanet_attrs,
        is_cancer=_is_cancer(mondo_id, mondo_entry, canonical_name),
        cohort=tuple(cohort),
        cohort_evidence=evidence,
        cohort_full=tuple(cohort_full),
        enrichment_cohort=enrichment_cohort,
        gwas_study_count=int(xc.get("gwas_study") or 0),
    )
