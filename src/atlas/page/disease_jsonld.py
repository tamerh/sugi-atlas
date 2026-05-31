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

BASE_URL = "https://sugi.bio/atlas"

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
        if mid:
            rec["identifier"] = mid
            rec["url"] = f"https://www.ebi.ac.uk/chembl/compound_report_card/{mid}/"
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
        if hgnc:
            rec["identifier"] = hgnc
            rec["url"] = f"https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/{hgnc}"
        if uni:
            rec["sameAs"] = [f"https://www.uniprot.org/uniprotkb/{uni}"]
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
        "description": _summary_description(bundle),
        "sameAs": same_as_urls(bundle) or None,
        "associatedGene": _associated_genes(bundle) or None,
        "drug": _drugs(bundle) or None,
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
    if codes:
        out["code"] = codes if len(codes) > 1 else codes[0]
    return {k: v for k, v in out.items() if v not in (None, [], "")}


def as_script_tag(jsonld: dict) -> str:
    """JSON-LD as an inline <script> block for the page body."""
    body = json.dumps(jsonld, indent=2)
    return f'<script type="application/ld+json">\n{body}\n</script>'


def as_jsonld_string(jsonld: dict) -> str:
    """Pretty-printed JSON-LD for the entity.jsonld sidecar."""
    return json.dumps(jsonld, indent=2) + "\n"
