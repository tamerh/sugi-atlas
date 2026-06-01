"""Schema.org JSON-LD emitter for Atlas disease pages.

Produces a `MedicalCondition` JSON-LD blob whose primary value is the
`sameAs` cross-ref to Mondo / EFO / MeSH / OMIM / Orphanet — the federated
disease-identity signal AI agents and KG pipelines rely on to disambiguate
disease references.

Two output forms (mirroring atlas.page.jsonld for genes):
  - as_script_tag(jsonld)    → <script type="application/ld+json">…</script>,
                                inlined at the top of the page body
  - as_jsonld_string(jsonld) → pretty-printed JSON for the entity.jsonld sidecar
"""
import json
from atlas.page import links

BASE_URL = "https://sugi.bio/atlas"
_HOST = BASE_URL.rsplit("/atlas", 1)[0]  # "https://sugi.bio" — prefix for internal links

# Per-ontology URL templates. Mondo's OLS4 entry is the most stable.
def _mondo_url(mid: str) -> str:
    return (f"https://www.ebi.ac.uk/ols4/ontologies/mondo/classes/"
            f"http%253A%252F%252Fpurl.obolibrary.org%252Fobo%252F"
            f"{mid.replace(':', '_')}")

def _efo_url(eid: str) -> str:
    return (f"https://www.ebi.ac.uk/ols4/ontologies/efo/classes/"
            f"http%253A%252F%252Fwww.ebi.ac.uk%252Fefo%252F{eid.replace(':', '_')}")

def _mesh_url(mid: str) -> str:
    return f"https://meshb.nlm.nih.gov/record/ui?ui={mid}"

def _omim_url(oid: str) -> str:
    return f"https://www.omim.org/entry/{oid}"

def _orphanet_url(oid: str) -> str:
    # Orphanet ids in biobtree look like "279" — bare number.
    return f"https://www.orpha.net/en/disease/detail/{oid}"

# Mondo OBO cross-ontology URL templates. Selected to land on the canonical
# entry page for each authority (no redirect chains).
_OBO_URL = {
    "doid":     lambda i: f"https://disease-ontology.org/term/DOID:{i}" if not i.startswith("DOID:") else f"https://disease-ontology.org/term/{i}",
    "ncit":     lambda i: f"https://ncithesaurus.nci.nih.gov/ncitbrowser/ConceptReport.jsp?dictionary=NCI_Thesaurus&code={i}",
    "umls":     lambda i: f"https://uts.nlm.nih.gov/uts/umls/concept/{i}",
    "medgen":   lambda i: f"https://www.ncbi.nlm.nih.gov/medgen/{i}",
    "sctid":    lambda i: f"https://browser.ihtsdotools.org/?perspective=full&conceptId1={i}",
    "icd10cm":  lambda i: f"https://www.icd10data.com/ICD10CM/Codes/{i}",
    "icd11":    lambda i: f"https://icd.who.int/browse11/l-m/en#/{i}",
    "gard":     lambda i: f"https://rarediseases.info.nih.gov/diseases/{i}",
}

def _uberon_url(uid: str) -> str:
    # uid looks like 'UBERON:0000310'
    return f"https://www.ebi.ac.uk/ols4/ontologies/uberon/classes/http%253A%252F%252Fpurl.obolibrary.org%252Fobo%252F{uid.replace(':', '_')}"


def same_as_urls(bundle: dict) -> list:
    """Federated-identity URLs to emit under `sameAs`. Order intentional:
    Mondo (the disease ontology authority) first, then the cross-ontology
    references."""
    b1 = bundle.get("1") or {}
    out = []
    mondo_id = b1.get("mondo_id")
    if mondo_id:
        out.append(_mondo_url(mondo_id))
    efo_id = b1.get("efo_id")
    if efo_id:
        out.append(_efo_url(efo_id))
    for mesh in (b1.get("mesh_ids") or []):
        out.append(_mesh_url(mesh))
    for mim in (b1.get("omim_ids") or []):
        out.append(_omim_url(mim))
    for orph in (b1.get("orphanet_ids") or []):
        out.append(_orphanet_url(orph))
    # Mondo OBO cross-ontology xrefs (DOID, NCIT, UMLS, etc.).
    obo = b1.get("obo_xrefs") or {}
    for ds in ("doid", "ncit", "umls", "medgen", "sctid", "icd10cm", "icd11", "gard"):
        builder = _OBO_URL.get(ds)
        if not builder:
            continue
        for ident in (obo.get(ds) or []):
            out.append(builder(ident))
    return out


def _associated_anatomy(b1: dict) -> list:
    """schema.org `associatedAnatomy` — AnatomicalStructure nodes from
    Mondo's `disease_has_location` UBERON axioms. Empty for systemic
    diseases (no body-site annotation in Mondo)."""
    uids = b1.get("anatomy_uberon_ids") or []
    out = []
    for uid in uids:
        out.append({
            "@type": "AnatomicalStructure",
            "identifier": uid,
            "url": _uberon_url(uid),
        })
    return out


def _drugs(bundle: dict) -> list:
    """schema.org `drug` array — top 5 trial drugs by disease-scoped trial
    count. Source: §13 trial_drugs (already deduplicated via the chembl
    parent/child fold). Each entry is a shallow Drug node with name +
    ChEMBL link. AI agents asking 'drugs used for X' key off this."""
    b13 = bundle.get("13") or {}
    drugs = b13.get("trial_drugs") or []
    out = []
    for d in drugs[:5]:
        name = d.get("name")
        mid = d.get("molecule_id")
        if not name and not mid:
            continue
        rec = {"@type": "Drug", "name": name or mid}
        sameas = []
        if mid:
            rec["identifier"] = mid
            sameas.append(f"https://www.ebi.ac.uk/chembl/compound_report_card/{mid}/")
        # Canonical url = the Atlas drug page when built, else the ChEMBL card;
        # the authority link is always preserved under sameAs.
        internal = links.drug_url(chembl_id=mid, name=name)
        if internal:
            rec["url"] = _HOST + internal
        elif sameas:
            rec["url"] = sameas[0]
        if sameas:
            rec["sameAs"] = sameas
        out.append(rec)
    return out


def _associated_genes(bundle: dict) -> list:
    """Compact MedicalCondition.associatedGene list for the cohort. We surface
    the cohort genes' symbols + HGNC id + their canonical UniProt — enough for
    an AI agent to follow back to NCBI/UniProt without paginating the page."""
    b5 = bundle.get("5") or {}
    genes = b5.get("genes") or []
    out = []
    for g in genes[:50]:  # full cohort
        sym = g.get("symbol")
        hgnc = g.get("hgnc_id")
        uni = g.get("canonical_uniprot")
        if not sym:
            continue
        rec = {"@type": "Gene", "name": sym}
        sameas = []
        if hgnc:
            rec["identifier"] = hgnc
            sameas.append(f"https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/{hgnc}")
        if uni:
            sameas.append(f"https://www.uniprot.org/uniprotkb/{uni}")
        internal = links.gene_url(symbol=sym, hgnc_id=hgnc)
        if internal:
            rec["url"] = _HOST + internal
        elif sameas:
            rec["url"] = sameas[0]
        if sameas:
            rec["sameAs"] = sameas
        out.append(rec)
    return out


def _summary_description(bundle: dict) -> str:
    """One-sentence canonical description from §1 + §2 counts (no LLM)."""
    b1 = bundle.get("1") or {}
    b2 = bundle.get("2") or {}
    b5 = bundle.get("5") or {}
    name = b1.get("canonical_name") or b1.get("name") or "?"
    bits = []
    if b2.get("assoc_total"):
        bits.append(f"{b2['assoc_total']:,} GWAS associations")
    if b2.get("study_total"):
        bits.append(f"{b2['study_total']:,} studies")
    if b5.get("gene_count"):
        bits.append(f"{b5['gene_count']} cohort genes")
    suffix = (" (" + ", ".join(bits) + ")") if bits else ""
    cancer = ", a cancer" if b1.get("is_cancer") else ""
    return f"{name}{cancer}{suffix}."


def build_jsonld(bundle: dict, slug: str, base_url: str = BASE_URL) -> dict:
    """Compose the schema.org MedicalCondition JSON-LD dict from a full
    collector bundle. Skips fields whose value is None/empty."""
    b1 = bundle.get("1") or {}
    canonical_page = f"{base_url}/disease/{slug}/"
    name = b1.get("canonical_name") or b1.get("name") or slug

    out = {
        "@context": "https://schema.org",
        "@type": "MedicalCondition",
        "@id": canonical_page,
        "name": name,
        "identifier": b1.get("mondo_id"),
        "url": canonical_page,
        "alternateName": (b1.get("synonyms") or [])[:20] or None,
        "description": _summary_description(bundle),
        "sameAs": same_as_urls(bundle) or None,
        "associatedGene": _associated_genes(bundle) or None,
        "drug": _drugs(bundle) or None,
        "epidemiology": _epidemiology(b1) or None,
        "signOrSymptom": _sign_or_symptom(b1) or None,
        "associatedAnatomy": _associated_anatomy(b1) or None,
    }
    # MedicalCode entries for each ontology (more granular than sameAs).
    codes = []
    if b1.get("mondo_id"):
        codes.append({"@type": "MedicalCode", "codingSystem": "Mondo",
                      "codeValue": b1["mondo_id"]})
    if b1.get("efo_id"):
        codes.append({"@type": "MedicalCode", "codingSystem": "EFO",
                      "codeValue": b1["efo_id"]})
    for mesh in (b1.get("mesh_ids") or []):
        codes.append({"@type": "MedicalCode", "codingSystem": "MeSH",
                      "codeValue": mesh})
    for mim in (b1.get("omim_ids") or []):
        codes.append({"@type": "MedicalCode", "codingSystem": "OMIM",
                      "codeValue": mim})
    for orph in (b1.get("orphanet_ids") or []):
        codes.append({"@type": "MedicalCode", "codingSystem": "Orphanet",
                      "codeValue": orph})
    # Mondo OBO cross-ontology xrefs as additional MedicalCode entries.
    _CODING_SYSTEM = {"doid": "DOID", "ncit": "NCI Thesaurus",
                      "umls": "UMLS", "medgen": "MedGen",
                      "sctid": "SNOMED CT", "icd10cm": "ICD-10-CM",
                      "icd11": "ICD-11", "gard": "GARD",
                      "meddra": "MedDRA", "nord": "NORD"}
    for ds, ids in (b1.get("obo_xrefs") or {}).items():
        sys_name = _CODING_SYSTEM.get(ds, ds)
        for ident in ids:
            codes.append({"@type": "MedicalCode", "codingSystem": sys_name,
                          "codeValue": ident})
    if codes:
        out["code"] = codes if len(codes) > 1 else codes[0]
    return {k: v for k, v in out.items() if v not in (None, [], "")}


def _epidemiology(b1: dict) -> str:
    """schema.org `epidemiology` is a free-text field. Prefer the Validated
    point-prevalence row with the widest geography (Worldwide > continental >
    country). Falls back to the first prevalence row, or empty."""
    prevs = b1.get("prevalences") or []
    if not prevs:
        return ""
    def _rank(p):
        geo = p.get("geographic") or ""
        return (0 if p.get("validation_status") == "Validated" else 1,
                0 if p.get("prevalence_type") == "Point prevalence" else 1,
                0 if geo == "Worldwide" else (1 if geo in ("Europe", "Americas") else 2))
    best = sorted(prevs, key=_rank)[0]
    parts = [best.get("prevalence_type") or "Prevalence",
             best.get("prevalence_class") or ""]
    if best.get("geographic"):
        parts.append(f"({best['geographic']})")
    if best.get("validation_status") == "Validated":
        parts.append("[Orphanet-validated]")
    return " ".join(p for p in parts if p)


def _sign_or_symptom(b1: dict) -> list:
    """schema.org `signOrSymptom` — MedicalSignOrSymptom nodes per HPO
    phenotype (top 20 by frequency). Each carries the HPO id as identifier
    and `frequency` as natural-language label so consumers can interpret."""
    phs = b1.get("phenotypes") or []
    if not phs:
        return []
    out = []
    for p in phs[:20]:
        hpo_id = p.get("hpo_id")
        node = {
            "@type": "MedicalSignOrSymptom",
            "name": p.get("hpo_term"),
            "identifier": hpo_id,
        }
        if hpo_id:
            node["url"] = f"https://hpo.jax.org/app/browse/term/{hpo_id}"
        if p.get("frequency"):
            node["frequency"] = p["frequency"]
        out.append(node)
    return out


def as_script_tag(jsonld: dict) -> str:
    """JSON-LD as an inline <script> block for the page body."""
    body = json.dumps(jsonld, indent=2)
    return f'<script type="application/ld+json">\n{body}\n</script>'


def as_jsonld_string(jsonld: dict) -> str:
    """Pretty-printed JSON-LD for the entity.jsonld sidecar."""
    return json.dumps(jsonld, indent=2) + "\n"
