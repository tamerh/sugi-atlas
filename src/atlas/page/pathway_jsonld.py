#!/usr/bin/env python3
"""schema.org JSON-LD for a pathway page. A Reactome pathway is modelled as a
`DefinedTerm` (a term in the Reactome vocabulary), with sameAs to its Reactome
detail page and the GO term it maps to — the federated-identity signal AI agents
and search use. Mirrors the gene/disease/drug jsonld helper shape.
"""
import json

BASE_URL = "https://sugi.bio/atlas"


def build_jsonld(bundle: dict, slug: str, base_url: str = BASE_URL) -> dict:
    rid = bundle.get("reactome_id")
    j = {
        "@context": "https://schema.org",
        "@type": "DefinedTerm",
        "@id": f"{base_url}/pathway/{slug}/",
        "identifier": rid,
        "name": bundle.get("name"),
        "inDefinedTermSet": "https://reactome.org/",
    }
    same = []
    if rid:
        same.append(f"https://reactome.org/content/detail/{rid}")
    if bundle.get("go_id"):
        same.append("http://purl.obolibrary.org/obo/" + bundle["go_id"].replace(":", "_"))
    if same:
        j["sameAs"] = same
    # member genes as hasPart (capped in the inline view by compact_for_inline)
    members = [m.get("symbol") for m in (bundle.get("members") or []) if m.get("symbol")]
    if members:
        j["hasPart"] = [{"@type": "Gene", "name": s} for s in members]
    return {k: v for k, v in j.items() if v is not None}


def as_script_tag(jsonld: dict) -> str:
    from atlas.page.jsonld_inline import compact_for_inline
    return (f'<script type="application/ld+json">\n'
            f'{json.dumps(compact_for_inline(jsonld), indent=2)}\n</script>')


def as_jsonld_string(jsonld: dict) -> str:
    return json.dumps(jsonld, indent=2) + "\n"
