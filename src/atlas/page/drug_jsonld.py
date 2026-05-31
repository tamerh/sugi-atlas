"""Schema.org JSON-LD emitter for Atlas drug pages.

Produces a `Drug` JSON-LD blob whose primary value is the `sameAs` cross-ref to
ChEMBL / PubChem / ChEBI / ATC — the federated chemical-identity signal — plus
`target` (Gene nodes) and `treats` (MedicalCondition nodes) so AI agents can
answer "what does drug X target / treat" directly from the structured head.

Two output forms (mirror atlas.page.jsonld / disease_jsonld):
  - as_script_tag(jsonld)    → inline <script> at the top of the page body
  - as_jsonld_string(jsonld) → pretty-printed JSON for the entity.jsonld sidecar
"""
import json

BASE_URL = "https://sugi.bio/atlas"


def same_as_urls(b1: dict) -> list:
    out = []
    if b1.get("chembl_id"):
        out.append(f"https://www.ebi.ac.uk/chembl/compound_report_card/{b1['chembl_id']}/")
    if b1.get("pubchem_cid"):
        out.append(f"https://pubchem.ncbi.nlm.nih.gov/compound/{b1['pubchem_cid']}")
    if b1.get("chebi_id"):
        out.append(f"https://www.ebi.ac.uk/chebi/searchId.do?chebiId={b1['chebi_id']}")
    for atc in (b1.get("atc_codes") or []):
        out.append(f"https://www.whocc.no/atc_ddd_index/?code={atc}")
    return out


def _targets(b2: dict) -> list:
    out = []
    for t in (b2.get("primary_targets") or []):
        sym = t.get("gene_symbol")
        if not sym:
            continue
        rec = {"@type": "Gene", "name": sym}
        if t.get("hgnc_id"):
            rec["identifier"] = t["hgnc_id"]
            rec["url"] = (f"https://www.genenames.org/data/gene-symbol-report/"
                          f"#!/hgnc_id/{t['hgnc_id']}")
        if t.get("uniprot"):
            rec["sameAs"] = [f"https://www.uniprot.org/uniprotkb/{t['uniprot']}"]
        out.append(rec)
    return out


def _treats(b4: dict) -> list:
    out = []
    for i in (b4.get("indications") or [])[:10]:
        name = i.get("name")
        mondo = i.get("mondo_id")
        if not name and not mondo:
            continue
        rec = {"@type": "MedicalCondition", "name": name or mondo}
        if mondo:
            rec["identifier"] = mondo
        out.append(rec)
    return out


def _description(b1: dict, b2: dict, b6: dict) -> str:
    name = b1.get("canonical_name") or b1.get("chembl_id") or "?"
    roles = (b6 or {}).get("chebi_roles") or []
    klass = roles[0] if roles else (b1.get("molecule_type") or "drug")
    genes = [t.get("gene_symbol") for t in (b2.get("primary_targets") or [])
             if t.get("gene_symbol")][:3]
    tgt = (" targeting " + ", ".join(genes)) if genes else ""
    return f"{name} — {klass}{tgt}."


def build_jsonld(bundle: dict, slug: str, base_url: str = BASE_URL) -> dict:
    b1 = bundle.get("1") or {}
    b2 = bundle.get("2") or {}
    b4 = bundle.get("4") or {}
    b6 = bundle.get("6") or {}
    page = f"{base_url}/drug/{slug}/"
    name = b1.get("canonical_name") or slug

    roles = (b6 or {}).get("chebi_roles") or []
    out = {
        "@context": "https://schema.org",
        "@type": "Drug",
        "@id": page,
        "name": name,
        "identifier": b1.get("chembl_id"),
        "url": page,
        "description": _description(b1, b2, b6),
        "alternateName": list(b1.get("alt_names") or [])[:12] or None,
        "drugClass": roles[0] if roles else None,
        "sameAs": same_as_urls(b1) or None,
        "target": _targets(b2) or None,
        "treats": _treats(b4) or None,
    }
    # Chemistry descriptors (small molecules only)
    if b1.get("inchi_key"):
        out["inChIKey"] = b1["inchi_key"]
    if b1.get("smiles"):
        out["smiles"] = b1["smiles"]
    if b1.get("molecular_formula"):
        out["molecularFormula"] = b1["molecular_formula"]
    if b1.get("molecular_weight"):
        out["molecularWeight"] = b1["molecular_weight"]
    return {k: v for k, v in out.items() if v not in (None, [], "")}


def as_script_tag(jsonld: dict) -> str:
    return f'<script type="application/ld+json">\n{json.dumps(jsonld, indent=2)}\n</script>'


def as_jsonld_string(jsonld: dict) -> str:
    return json.dumps(jsonld, indent=2) + "\n"
