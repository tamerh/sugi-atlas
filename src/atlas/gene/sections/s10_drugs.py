"""§10 — drugs: ChEMBL targets + phased molecules, PharmGKB, BindingDB sample,
clinical trials via the DISEASE route (gene → MONDO → trials)."""
from collections import Counter
from atlas.biobtree import map_all, entry, xref_counts
from atlas.civic import aggregate_predictive
from atlas.gene.sections.base import Section

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
    ">>uniprot>>bindingdb",
    ">>uniprot>>pubchem_activity",
    ">>hgnc>>entrez>>ctd_gene_interaction",
    ">>hgnc>>gencc>>mondo>>clinical_trials",
    ">>hgnc>>clinvar>>mondo>>clinical_trials",
    ">>hgnc>>civic_evidence",
)
DATASETS = ("chembl_target", "chembl_molecule", "chembl_activity", "chembl_assay",
            "chembl_document", "patent_compound", "cellosaurus",
            "pharmgkb_gene", "pharmgkb_clinical", "pharmgkb_variant",
            "bindingdb", "pubchem_activity",
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

    bundle["pharmgkb"] = [{"id": t["id"], "vip": t.get("is_vip"),
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

    bd = map_all(uni, ">>uniprot>>bindingdb", cap=2) if uni else []
    bundle["bindingdb_sample"] = [{"ligand": t.get("ligand_name"), "ki": t.get("ki"),
                                   "ic50": t.get("ic50")} for t in bd[:30]]
    bundle["bindingdb_sampled"] = len(bd)

    # PubChem BioAssay activities — Active outcomes with a real affinity value,
    # sorted by potency (low IC50/Ki/Kd/EC50 = most potent). The activity_id
    # itself encodes CID_AID_VERSION so we surface those without per-row entry
    # fetches (PMID enrichment would require entries — deferred to scale-out).
    pa = map_all(uni, ">>uniprot>>pubchem_activity") if uni else []
    actives = [r for r in pa if r.get("activity_outcome") == "Active"
               and (r.get("activity_type") or "").lower() in _AFFINITY_TYPES
               and _f(r.get("value")) not in (0.0, float("inf"))]
    actives.sort(key=lambda r: _f(r.get("value")))
    bioassays = []
    for r in actives[:30]:
        parts = (r.get("id") or "").split("_")
        cid = parts[0] if parts else None
        aid = parts[1] if len(parts) > 1 else None
        bioassays.append({
            "id": r["id"], "cid": cid, "aid": aid,
            "activity_type": r.get("activity_type"),
            "value": r.get("value"), "unit": r.get("unit"),
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
    bundle["chembl_activities"] = [{
        "id": r["id"],
        "type": r.get("standard_type"),
        "value": r.get("standard_value"),
        "unit": r.get("standard_units"),
        "pchembl": r.get("pchembl"),
    } for r in potent[:30]]
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
              "bindingdb_sample", "pubchem_bioassay", "ctd_interactions",
              "disease_trials", "civic_evidence", "civic_predictive_total",
              "civic_evidence_total", "is_drug_target"),
    datasets=DATASETS, chains=CHAINS, collect_fn=collect,
)
