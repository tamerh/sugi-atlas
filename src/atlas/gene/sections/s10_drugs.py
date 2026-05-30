"""§10 — drugs: ChEMBL targets + phased molecules, PharmGKB, BindingDB sample,
clinical trials via the DISEASE route (gene → MONDO → trials)."""
from collections import Counter
from atlas.biobtree import map_all
from atlas.gene.sections.base import Section

CHAINS = (
    ">>uniprot>>chembl_target",
    '>>chembl_target>>chembl_molecule[highestDevelopmentPhase>=1]',
    ">>uniprot>>chembl_activity",
    ">>uniprot>>chembl_target>>chembl_assay",
    ">>hgnc>>pharmgkb_gene",
    ">>uniprot>>bindingdb",
    ">>uniprot>>pubchem_activity",
    ">>hgnc>>entrez>>ctd_gene_interaction",
    ">>hgnc>>gencc>>mondo>>clinical_trials",
    ">>hgnc>>clinvar>>mondo>>clinical_trials",
)
DATASETS = ("chembl_target", "chembl_molecule", "chembl_activity", "chembl_assay",
            "pharmgkb_gene", "bindingdb", "pubchem_activity",
            "ctd_gene_interaction", "entrez",
            "clinical_trials", "mondo", "gencc", "clinvar", "uniprot", "hgnc")

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

    bundle["pharmgkb"] = [{"id": t["id"], "vip": t.get("is_vip"),
                           "cpic_guideline": t.get("has_cpic_guideline")}
                          for t in map_all(a.hgnc_id, ">>hgnc>>pharmgkb_gene")]

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
    bundle["chembl_assay_samples"] = [
        {"id": s["id"], "type": _ASSAY_TYPE_NAMES.get(s["type"], s["type"]),
         "desc": s["desc"][:240]} for s in samples
    ]

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
    bundle["is_drug_target"] = bool(targets)
    return bundle

SECTION = Section(
    id="10", name="drugs",
    description=("ChEMBL targets + phased targeting molecules, PharmGKB, BindingDB "
                 "affinities, PubChem BioAssay actives (sorted by potency, with "
                 "clickable CID/AID), clinical trials via disease route"),
    needs=("hgnc_id", "canonical_uniprot"),
    produces=("chembl_targets", "molecules", "chembl_activities",
              "chembl_assay_total", "chembl_assay_type_counts",
              "chembl_assay_samples", "pharmgkb",
              "bindingdb_sample", "pubchem_bioassay", "ctd_interactions",
              "disease_trials", "is_drug_target"),
    datasets=DATASETS, chains=CHAINS, collect_fn=collect,
)
