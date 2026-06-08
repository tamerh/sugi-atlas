"""§10 — drugs: ChEMBL targets + phased molecules, PharmGKB, GtoPdb curated
pharmacology, potency-ranked BindingDB, clinical trials via the DISEASE route
(gene → MONDO → trials)."""
import html
import re
from collections import Counter
from atlas.biobtree import map_all, entry, xref_counts
from atlas.civic import aggregate_predictive
from atlas.gene.sections.base import Section

# Binding-affinity unit → nanomolar, for ranking heterogeneous BindingDB values
# (strings like "1.0 nM", "0.5 µM") onto one comparable scale.
_NM_PER_UNIT = {"m": 1e9, "mm": 1e6, "um": 1e3, "nm": 1.0, "pm": 1e-3, "fm": 1e-6}
_AFFINITY_RE = re.compile(r'\s*[<>~=]*\s*([\d.]+)\s*([a-z]+)', re.I)
_HUMAN_ORG = {"homo sapiens", "human"}


def _clean_ligand(name):
    """BindingDB ligand_name is a '::'-joined synonym dump (IUPAC, CHEMBL ids,
    patent refs, analog labels). Pick the first non-CHEMBL segment and return it
    in FULL — no baked-in truncation; the frontend clamps long IUPAC names for
    display. Patent-extracted compounds carry only a patent reference here (e.g.
    'US8524722, 5'); recovering a chemical name for those needs a structure join
    to ChEMBL/PubChem upstream in BioBTree."""
    if not name:
        return name
    parts = [p.strip() for p in str(name).split("::") if p.strip()]
    named = [p for p in parts if not p.upper().startswith("CHEMBL")]
    return (named or parts)[0]


_PATENT_ID_RX = re.compile(r"^((?:US|WO|EP|CN|JP)\d+)", re.I)


def _patent_id(name):
    """Leading patent number of a patent-reference ligand name
    ('US8524722, 5' -> 'US8524722'), else None."""
    m = _PATENT_ID_RX.match((name or "").strip())
    return m.group(1).upper() if m else None


def _is_patent_ref(name):
    """True if a ligand name is only a patent reference (a patent-extracted
    compound with no chemical name in the BindingDB record)."""
    return _patent_id(name) is not None


_PATENT_KINDS = ("B2", "B1", "A1", "A2", "B", "A")


def _patent_title(pn):
    """Invention title of a patent from biobtree's Patent dataset. Our
    patent_number is plain (US8524722) but the dataset is keyed 'US-8524722-B2'
    (hyphens + a kind code not in our data), so reformat and try the common kind
    codes (B2 dominates). None if none resolve."""
    m = re.match(r"^([A-Z]{2})(\d+)$", (pn or "").strip())
    if not m:
        return None
    base = f"{m.group(1)}-{m.group(2)}"
    for kc in _PATENT_KINDS:
        try:
            a = entry(f"{base}-{kc}", "patent").get("Attributes", {}).get("Patent", {})
            if a.get("title"):
                return a["title"]
        except Exception:
            pass
    return None


# Numeric character references in IUPAC ligand names — both the valid form
# (&#8243; = ″) and BindingDB's mangled form ($#8243;, where the & is corrupted
# upstream). Decode both to the actual character; html.unescape handles named
# entities (&amp; etc.) on top.
_NUM_ENTITY_RX = re.compile(r"[$&]#(\d+);")


def _decode_entities(s):
    if not s:
        return s
    s = _NUM_ENTITY_RX.sub(lambda m: chr(int(m.group(1))), s)
    return html.unescape(s)


_MEASURE_NUM = re.compile(r"^(\s*[<>~=]*\s*)([\d.]+)(.*)$")


def _round_measure(s):
    """Round the numeric part of a measure string to 3 sig figs for display —
    BindingDB ships spurious precision (0.031623 nM = 0.0316 nM) on heterogeneous
    assays. Keeps qualifier (>/~) + unit."""
    if not s:
        return s
    m = _MEASURE_NUM.match(str(s))
    if not m:
        return s
    try:
        v = float(m.group(2))
    except ValueError:
        return s
    return f"{m.group(1)}{float(f'{v:.3g}'):g}{m.group(3)}"


def _affinity_nm(s):
    """Parse a BindingDB affinity string to nanomolar; None if unparseable."""
    if not s:
        return None
    m = _AFFINITY_RE.match(str(s).strip().replace("µ", "u").replace("μ", "u"))
    if not m:
        return None
    try:
        val = float(m.group(1))
    except ValueError:
        return None
    factor = _NM_PER_UNIT.get(m.group(2).lower())
    return val * factor if (factor and val > 0) else None

CHAINS = (
    ">>uniprot>>chembl_target",
    '>>chembl_target>>chembl_molecule[highestDevelopmentPhase>=1]',
    ">>uniprot>>chembl_activity",
    ">>uniprot>>chembl_target>>chembl_assay",
    ">>chembl_assay>>chembl_document",
    ">>chembl_molecule>>patent_compound",
    ">>hgnc>>cellosaurus",
    ">>hgnc>>pharmgkb_gene",
    ">>hgnc>>pharmgkb_clinical",
    ">>hgnc>>pharmgkb_variant",
    ">>hgnc>>pharmgkb_guideline",
    ">>uniprot>>bindingdb",
    ">>uniprot>>gtopdb",
    ">>uniprot>>gtopdb>>gtopdb_interaction",
    ">>uniprot>>pubchem_activity",
    ">>hgnc>>entrez>>ctd_gene_interaction",
    ">>hgnc>>gencc>>mondo>>clinical_trials",
    ">>hgnc>>clinvar>>mondo>>clinical_trials",
    ">>hgnc>>civic_evidence",
)
DATASETS = ("chembl_target", "chembl_molecule", "chembl_activity", "chembl_assay",
            "chembl_document", "patent_compound", "cellosaurus",
            "pharmgkb_gene", "pharmgkb_clinical", "pharmgkb_variant",
            "pharmgkb_guideline",
            "bindingdb", "gtopdb", "pubchem_activity",
            "ctd_gene_interaction", "entrez",
            "clinical_trials", "mondo", "gencc", "clinvar", "uniprot", "hgnc",
            "civic_evidence")

# ChEMBL assay type codes — single-letter classification we surface as a
# breakdown ("how heavily this target is profiled, and by what experimental
# style"). Source: ChEMBL data dictionary / BAO.
_ASSAY_TYPE_NAMES = {
    "B": "Binding",
    "F": "Functional",
    "A": "ADMET",
    "P": "Physicochemical",
    "T": "Toxicity",
    "U": "Unclassified",
}

# Activity types we sort by potency (low value = potent). Other outcomes carry
# qualifier flags rather than affinity values.
_AFFINITY_TYPES = {"ki", "ic50", "kd", "ec50"}

def _phase(d):
    return int(d) if (d or "").isdigit() else 0

def _f(x):
    try: return float(x)
    except (TypeError, ValueError): return float("inf")

def collect(a):
    bundle = {"section": "10_drugs", "symbol": a.symbol}
    uni = a.canonical_uniprot

    targets = map_all(uni, ">>uniprot>>chembl_target") if uni else []
    bundle["chembl_targets"] = [{"id": t["id"], "title": t.get("title"), "type": t.get("type")}
                                for t in targets]

    # Phased drugs only (highestDevelopmentPhase>=1) — the raw chembl_molecule edge
    # returns thousands of ID-ordered screening compounds; filter + rank for the real drugs.
    drugs = {}
    for t in targets:
        for m in map_all(t["id"], '>>chembl_target>>chembl_molecule[highestDevelopmentPhase>=1]', cap=10):
            drugs[m["id"]] = {"id": m["id"], "name": m.get("name"), "type": m.get("type"),
                              "phase": m.get("highestDevelopmentPhase")}
    bundle["molecules"] = sorted(drugs.values(), key=lambda d: _phase(d["phase"]), reverse=True)
    bundle["molecule_count"] = len(drugs)

    # Patent literature coverage per phased molecule — chemistry IP intensity.
    # Each chembl_molecule maps to 0-N patent_compound records (PubChem CIDs);
    # each carries an xref count of patents that mention the compound.
    # Approved/late-phase drugs typically score in the thousands (Erlotinib
    # ≈2358, Gefitinib ≈2154); early-phase have few or zero. Capped at top 20
    # by phase to bound entry-fetch cost (~3 entries/molecule worst case).
    for mol in bundle["molecules"][:20]:
        pcids = [t.get("id") for t in map_all(mol["id"], ">>chembl_molecule>>patent_compound", cap=2)
                 if t.get("id")]
        total = 0
        for pcid in pcids:
            try:
                total += xref_counts(entry(pcid, "patent_compound")).get("patent", 0)
            except Exception:
                pass
        mol["patent_count"] = total
        mol["patent_compound_ids"] = pcids
    bundle["patent_total"] = sum(m.get("patent_count", 0) for m in bundle["molecules"][:20])

    # is_vip is dropped on purpose: it's broken upstream (always true, even for
    # ACTB/GAPDH/TTN), so it's noise. has_cpic_guideline is real — keep that.
    bundle["pharmgkb"] = [{"id": t["id"],
                           "cpic_guideline": t.get("has_cpic_guideline")}
                          for t in map_all(a.hgnc_id, ">>hgnc>>pharmgkb_gene")]

    # PharmGKB clinical annotations + variant pages — biobtree #13 fix
    # landed for these two (guideline still 0-rows). Clinical annotations
    # carry the variant + drug + phenotype + level_of_evidence tuple that
    # makes per-gene PGx context concrete (vs the bare existence flag).
    bundle["pharmgkb_clinical"] = [{
        "id": t.get("id"),
        "variant": t.get("variant"),
        "type": t.get("type"),
        "level_of_evidence": t.get("level_of_evidence"),
        "chemicals": t.get("chemicals"),
        "phenotypes": t.get("phenotypes"),
    } for t in map_all(a.hgnc_id, ">>hgnc>>pharmgkb_clinical")]
    bundle["pharmgkb_variant"] = [{
        "id": t.get("id"),
        "name": t.get("variant_name"),
        "gene_symbols": t.get("gene_symbols"),
        "level_of_evidence": t.get("level_of_evidence"),
        "score": t.get("score"),
        "clinical_annotation_count": t.get("clinical_annotation_count"),
        "associated_drugs": t.get("associated_drugs"),
    } for t in map_all(a.hgnc_id, ">>hgnc>>pharmgkb_variant")]
    # PharmGKB guidelines — CPIC / DPWG / CPNDS dosing guidance per
    # gene-drug pair. Many for canonical pharmacogenes (CYP2C19 has 37,
    # CYP2D6 has 69); zero for non-pharmacogenes (most cancer genes).
    bundle["pharmgkb_guideline"] = [{
        "id": t.get("id"),
        "name": t.get("name"),
        "source": t.get("source"),  # CPIC / DPWG / CPNDS / etc.
        "gene_symbols": t.get("gene_symbols"),
        "chemical_names": t.get("chemical_names"),
        "has_dosing_info": t.get("has_dosing_info") == "true",
        "has_recommendation": t.get("has_recommendation") == "true",
    } for t in map_all(a.hgnc_id, ">>hgnc>>pharmgkb_guideline")]

    # GtoPdb / IUPHAR — hand-curated pharmacology (TIER 1). gene-keyed
    # (uniprot→gtopdb), structurally clean. Target classification + the most
    # potent curated ligand interactions (affinity = pIC50/pKi, higher = tighter).
    gtopdb = map_all(uni, ">>uniprot>>gtopdb") if uni else []
    if gtopdb:
        g0 = gtopdb[0]
        bundle["gtopdb_target"] = {"id": g0.get("id"), "type": g0.get("type"),
                                   "family": g0.get("family_name")}
        inter = map_all(uni, ">>uniprot>>gtopdb>>gtopdb_interaction", cap=200)

        def _aff(r):
            try:
                return float(r.get("affinity") or 0)
            except (TypeError, ValueError):
                return 0.0
        inter = [r for r in inter if _aff(r) > 0]
        inter.sort(key=_aff, reverse=True)
        bundle["gtopdb_interactions"] = [{
            "ligand": r.get("ligand_name"), "type": r.get("type"),
            "action": r.get("action"), "affinity": round(_aff(r), 2),
            "parameter": r.get("affinity_parameter"),
        } for r in inter[:25]]
        bundle["gtopdb_interaction_count"] = len(inter)

    # BindingDB — measured binding affinities (TIER 2: heterogeneous assays, not
    # directly comparable). Promote from a 2-row sample to a potency-ranked
    # slice: filter to the human target, rank by the tightest available measure
    # (Ki/IC50/Kd/EC50 normalized to nM).
    bd = map_all(uni, ">>uniprot>>bindingdb") if uni else []
    # Recover a chemical name for patent-extracted compounds, whose ligand_name is
    # only a patent reference. biobtree #35 now projects pubchem_cids on the
    # bindingdb map, and one >>pubchem call gives cid -> title (IUPAC), so we join
    # on the cid with no per-row entry() fan-out. The source patent itself is kept
    # and surfaced in its own column (parsed from the original reference).
    pc_title = ({p["id"]: p.get("title") for p in map_all(uni, ">>uniprot>>bindingdb>>pubchem")
                 if p.get("id")} if uni else {})

    def _best_measure(r):
        best = None
        for m in ("ki", "ic50", "kd", "ec50"):
            nm = _affinity_nm(r.get(m))
            if nm is not None and (best is None or nm < best[0]):
                best = (nm, m.upper(), r.get(m))
        return best
    human = [r for r in bd
             if (r.get("target_source_organism") or "").lower() in _HUMAN_ORG]
    ranked = []
    for r in human:
        bm = _best_measure(r)
        if not bm:
            continue
        name = _clean_ligand(r.get("ligand_name"))
        # Authoritative patent source: the projected patent_number (populated for
        # named AND patent-extracted compounds), else the patent parsed from the
        # ligand reference. So a blank Patent cell means genuinely no patent.
        patent = (r.get("patent_number") or "").strip() or _patent_id(name)
        # Recover a chemical name when ligand_name is only a patent reference.
        if _is_patent_ref(name):
            cid = (r.get("pubchem_cids") or "").split(",")[0].strip()
            recovered = pc_title.get(cid)
            if recovered and not _is_patent_ref(recovered):
                name = recovered                  # IUPAC from PubChem
        ranked.append({"ligand": _decode_entities(name), "patent": patent, "measure": bm[1],
                       "value": _round_measure(bm[2]), "nm": round(bm[0], 4)})
    ranked.sort(key=lambda x: x["nm"])
    top = ranked[:50]
    # Resolve each distinct patent's title (invention name — more readable than the
    # IUPAC), cached so each patent is fetched once across the shown rows.
    pat_titles = {}
    for x in top:
        pn = x.get("patent")
        if pn and pn not in pat_titles:
            pat_titles[pn] = _patent_title(pn)
        x["patent_title"] = pat_titles.get(pn)
    bundle["bindingdb_ranked"] = top
    bundle["bindingdb_total"] = len(bd)
    bundle["bindingdb_human"] = len(human)
    bundle["bindingdb_measured"] = len(ranked)

    # PubChem BioAssay activities — Active outcomes with a real affinity value,
    # sorted by potency (low IC50/Ki/Kd/EC50 = most potent). The activity_id
    # itself encodes CID_AID_VERSION so we surface those without per-row entry
    # fetches (PMID enrichment would require entries — deferred to scale-out).
    pa = map_all(uni, ">>uniprot>>pubchem_activity") if uni else []
    # Compound names for the CIDs (PubChem title / IUPAC), one bulk hop. PubChem
    # BioAssay is largely COMPLEMENTARY to BindingDB (tiny overlap — it contributes
    # ~25-30% of the compound union that BindingDB lacks), so it adds real chemical
    # matter and is worth showing identifiably rather than as bare CIDs.
    pc_names = ({p["id"]: p.get("title") for p in map_all(uni, ">>uniprot>>pubchem_activity>>pubchem")
                 if p.get("id") and p.get("title")} if uni else {})
    # Assay names (what each AID measured, e.g. "Inhibition of EGFR") — one bulk
    # hop, keyed by AID, for readable experimental context alongside the compound.
    assay_names = ({a["id"]: a.get("name") for a in map_all(uni, ">>uniprot>>pubchem_activity>>pubchem_assay")
                    if a.get("id") and a.get("name")} if uni else {})
    # Keep value==0 rows: PubChem rounds to 4-decimal µM, so the MOST potent
    # (sub-0.1 nM) round to 0.0000 — excluding them dropped the actual top.
    actives = [r for r in pa if r.get("activity_outcome") == "Active"
               and (r.get("activity_type") or "").lower() in _AFFINITY_TYPES
               and r.get("value") not in (None, "")
               and _f(r.get("value")) != float("inf")]
    actives.sort(key=lambda r: _f(r.get("value")))     # most potent first
    bioassays, seen = [], set()
    for r in actives:
        if len(bioassays) >= 50:
            break
        cid = (r.get("id") or "").split("_")[0] or None
        if cid in seen:                                # distinct compounds
            continue
        seen.add(cid)
        parts = (r.get("id") or "").split("_")
        aid = parts[1] if len(parts) > 1 else None
        v = r.get("value")
        bioassays.append({
            "id": r["id"], "cid": cid, "aid": aid,
            "name": pc_names.get(cid),
            "assay_name": assay_names.get(aid),
            "activity_type": r.get("activity_type"),
            "value": "<0.0001" if _f(v) == 0.0 else v,   # below the 4-decimal µM floor
            "unit": r.get("unit"),
        })
    bundle["pubchem_bioassay"] = bioassays
    bundle["pubchem_bioassay_total"] = len(pa)
    bundle["pubchem_bioassay_active_count"] = len(actives)

    # ChEMBL bioactivities for the target — pchembl is the gold metric
    # (-log10(M) standardized potency: 10 = 0.1 nM, 9 = 1 nM, 6 = 1 µM).
    # Sort by pchembl desc, keep rows with pchembl ≥ 5 (≤10 µM, "real binding").
    # Closes the gap for targets where PubChem activity is empty (e.g. KRAS,
    # see BIOBTREE_ISSUES.md #12).
    def _pchembl(r):
        try: return float(r.get("pchembl") or 0)
        except (TypeError, ValueError): return 0.0
    ca = map_all(uni, ">>uniprot>>chembl_activity") if uni else []
    potent = [r for r in ca if _pchembl(r) >= 5.0]
    potent.sort(key=_pchembl, reverse=True)
    # Resolve each displayed activity to its compound. The molecule is not stored
    # on the activity record (only an xref edge — see BIOBTREE_ISSUES #36), so we
    # traverse it per row for the top 50. The activity id (CHEMBL_ACT_…) is opaque;
    # the molecule (named when known, e.g. Mobocertinib, else its ChEMBL id) is the
    # identity that makes the row interpretable.
    acts, seen = [], set()
    for r in potent[:100]:            # cap the per-activity molecule lookups
        if len(acts) >= 50:
            break
        mol = map_all(r["id"], ">>chembl_activity>>chembl_molecule")
        m = mol[0] if mol else {}
        key = (r.get("pchembl"), r.get("standard_type"), r.get("standard_value"),
               r.get("standard_units"), m.get("id"))
        if key in seen:               # drop duplicate assay records (same compound + measurement)
            continue
        seen.add(key)
        acts.append({
            "id": r["id"],
            "type": r.get("standard_type"),
            "value": r.get("standard_value"),
            "unit": r.get("standard_units"),
            "pchembl": r.get("pchembl"),
            "molecule_id": m.get("id"),
            "molecule_name": (m.get("name") or "").strip() or None,
        })
    bundle["chembl_activities"] = acts
    bundle["chembl_activity_total"] = len(ca)
    bundle["chembl_activity_potent_count"] = len(potent)

    # ChEMBL assays — aggregated across all ChEMBL targets for this protein.
    # Per-target chembl_assay is capped (≈100/target on heavily-screened
    # targets), so totals are screening-floor signals; the type breakdown
    # captures the experimental mix. Sample descriptions surface the biology
    # behind the numbers (compound × target × mechanism).
    seen_assays = {}
    type_counts = Counter()
    for t in targets:
        for r in map_all(t["id"], ">>chembl_target>>chembl_assay", cap=100):
            aid = r.get("id")
            if not aid or aid in seen_assays:
                continue
            ty = (r.get("type") or "U").upper()
            type_counts[ty] += 1
            seen_assays[aid] = {"id": aid, "type": ty,
                                "desc": (r.get("desc") or "").strip()}
    bundle["chembl_assay_total"] = len(seen_assays)
    bundle["chembl_assay_type_counts"] = {_ASSAY_TYPE_NAMES.get(k, k): v
                                          for k, v in type_counts.most_common()}
    # Three representative assays — prefer rows with descriptions, biased
    # toward type variety so the sample illustrates the breadth.
    samples = []
    by_type = {}
    for r in seen_assays.values():
        if not r["desc"]:
            continue
        by_type.setdefault(r["type"], []).append(r)
    for ty in type_counts:
        if by_type.get(ty):
            samples.append(by_type[ty][0])
        if len(samples) >= 3:
            break
    # Source publication for each sample assay — one chembl_document per
    # assay edge. Adds title + journal context to the assay sample table so
    # readers can trace the screen to its paper. Tiny cost: 3 entry calls.
    enriched_samples = []
    for s in samples:
        rec = {"id": s["id"], "type": _ASSAY_TYPE_NAMES.get(s["type"], s["type"]),
               "desc": s["desc"][:240]}
        try:
            for d in map_all(s["id"], ">>chembl_assay>>chembl_document", cap=1):
                doc_id = d.get("id")
                if doc_id:
                    de = entry(doc_id, "chembl_document")
                    doc = (((de.get("Attributes") or {}).get("Chembl") or {})
                           .get("doc") or {})
                    rec["doc_id"] = doc_id
                    rec["doc_title"] = (doc.get("title") or "").strip()
                    rec["doc_journal"] = doc.get("journal")
                    break
        except Exception:
            pass
        enriched_samples.append(rec)
    bundle["chembl_assay_samples"] = enriched_samples

    # Cellosaurus — cell lines associated with this gene (mutated, deficient,
    # expressed-in, model-of, etc.). Totals run into the thousands for
    # heavily-studied tumor suppressors / oncogenes (TP53 ≈10k, KRAS ≈5k);
    # the category breakdown is the load-bearing signal. Surfaces "what
    # experimental resources exist for this gene" for downstream wet-lab
    # choices. Uncapped — full pagination.
    cell_rows = map_all(a.hgnc_id, ">>hgnc>>cellosaurus")
    cell_cats = Counter()
    cells = []
    for r in cell_rows:
        cat = (r.get("category") or "Unclassified").strip()
        cell_cats[cat] += 1
        cells.append({"id": r.get("id"), "name": r.get("name"),
                      "category": cat, "sex": r.get("sex")})
    bundle["cellosaurus_total"] = len(cells)
    bundle["cellosaurus_category_counts"] = dict(cell_cats.most_common())
    bundle["cellosaurus_samples"] = cells[:10]

    # CTD literature-mined chemical-gene interactions. Each row is a
    # PubMed-evidenced interaction with CV-coded action verbs (e.g.
    # 'increases^expression', 'decreases^activity'). Sort by pubmed_count
    # (top hits are the most-cited / best-supported); filter to Homo sapiens
    # since the CTD index has rodent/cell-line rows too. CTD's chemical_id
    # (e.g. C000228) is the head of the activity_id; we surface it for
    # clickable links to ctdbase.org.
    def _pmc(r):
        try: return int(r.get("pubmed_count") or 0)
        except (TypeError, ValueError): return 0
    ctd = map_all(a.hgnc_id, ">>hgnc>>entrez>>ctd_gene_interaction")
    ctd_human = [r for r in ctd if (r.get("organism") or "") == "Homo sapiens"]
    ctd_human.sort(key=_pmc, reverse=True)
    bundle["ctd_interactions"] = [{
        "id": r["id"],
        "chemical_id": (r.get("id") or "").split("_", 1)[0],
        "chemical": r.get("chemical_name"),
        "actions": [act.strip() for act in (r.get("interaction_actions") or "").split(";") if act.strip()],
        "pmids": r.get("pubmed_count"),
    } for r in ctd_human[:30]]
    bundle["ctd_interaction_total"] = len(ctd_human)

    # Clinical trials via the DISEASE route — chembl_molecule>>clinical_trials would
    # pollute with off-target drugs (ChEMBL target↔molecule is bioactivity-based).
    # Disease-level is biobtree's intended pattern.
    trials = {}
    for chain in (">>hgnc>>gencc>>mondo>>clinical_trials",
                  ">>hgnc>>clinvar>>mondo>>clinical_trials"):
        for t in map_all(a.hgnc_id, chain, cap=2):
            trials[t["id"]] = {"id": t["id"], "title": t.get("brief_title"),
                               "phase": t.get("phase"), "status": t.get("overall_status")}
    def _ph(t):
        p = (t.get("phase") or "").upper().replace("PHASE", "").strip()
        return int(p) if p.isdigit() else 0
    bundle["disease_trials"] = sorted(trials.values(), key=_ph, reverse=True)
    bundle["disease_trial_count"] = len(trials)

    # CIViC clinical evidence — the drug × variant × indication triple, the
    # precision-medicine narrative ("drug X for variant Y in disease Z"). The
    # civic_evidence dataset carries curated predictive / prognostic /
    # diagnostic associations per molecular profile, already normalized
    # (xrefs to hgnc/mondo/chembl_molecule/pubmed). We surface the PREDICTIVE
    # subset — the drug-actionable view — deduped to unique (variant, therapy,
    # indication, effect) and ranked by CIViC evidence level (A validated →
    # E inferential), keeping a representative evidence id + supporting-item
    # count per association. The Effect column (significance) distinguishes
    # Sensitivity/Response from Resistance — opposite clinical meaning, both
    # actionable. Reached via >>hgnc>>civic_evidence; empty for non-cancer
    # genes (TTN etc.) so the block cleanly elides on the long tail.
    civic = map_all(a.hgnc_id, ">>hgnc>>civic_evidence", cap=15)
    ranked, cstats = aggregate_predictive(civic, top=30)
    bundle["civic_evidence"] = ranked
    bundle["civic_evidence_total"] = cstats["evidence_total"]
    bundle["civic_predictive_total"] = cstats["predictive_total"]
    bundle["civic_association_total"] = cstats["association_total"]
    bundle["civic_evidence_type_counts"] = cstats["evidence_type_counts"]

    bundle["is_drug_target"] = bool(targets)
    return bundle

SECTION = Section(
    id="10", name="drugs",
    description=("ChEMBL targets + phased targeting molecules, PharmGKB, BindingDB "
                 "affinities, PubChem BioAssay actives (sorted by potency, with "
                 "clickable CID/AID), clinical trials via disease route, CIViC "
                 "clinical evidence (drug × variant × indication precision triple)"),
    needs=("hgnc_id", "canonical_uniprot"),
    produces=("chembl_targets", "molecules", "chembl_activities",
              "chembl_assay_total", "chembl_assay_type_counts",
              "chembl_assay_samples", "patent_total",
              "cellosaurus_total", "cellosaurus_category_counts",
              "cellosaurus_samples", "pharmgkb", "pharmgkb_clinical", "pharmgkb_variant",
              "pharmgkb_guideline",
              "gtopdb_target", "gtopdb_interactions", "gtopdb_interaction_count",
              "bindingdb_ranked", "bindingdb_total", "bindingdb_human",
              "bindingdb_measured", "pubchem_bioassay", "ctd_interactions",
              "disease_trials", "civic_evidence", "civic_predictive_total",
              "civic_evidence_total", "is_drug_target"),
    datasets=DATASETS, chains=CHAINS, collect_fn=collect,
)
