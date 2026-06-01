"""Schema.org JSON-LD emitter for Atlas gene pages.

Produces a `Gene` JSON-LD blob whose primary value is the `sameAs` cross-ref
to NCBI Gene / UniProt / Ensembl / HGNC / OMIM — the federated-identity
signal AI agents and Google Knowledge-Graph use to ground "is this the same
gene I just heard about?" decisions. Per the landscape report none of the
13 surveyed gene-info competitors emit this; shipping it is the single
highest-leverage AI-friendliness move.

Two output forms:
  - as_script_tag(jsonld)   → <script type="application/ld+json">…</script>,
                              inlined at the top of the page body markdown
  - as_jsonld_string(jsonld) → pretty-printed JSON for the entity.jsonld
                              sidecar (machine-fetchable; can be advertised
                              via <link rel="alternate" type="application/ld+json">)

Both consume the same dict from build_jsonld(bundle).
"""
import json
from atlas.page.declarative import declarative_sentence

BASE_URL = "https://sugi.bio/atlas"

def _strip_md(s):
    # The declarative lead uses **symbol** for visual emphasis; in JSON-LD's
    # description field we want plain text.
    return s.replace("**", "")

def same_as_urls(bundle):
    """Federated-identity URLs to emit under `sameAs`. Order is intentional:
    HGNC (the symbol authority) first, then the data-rich repositories."""
    b1 = bundle.get("1") or {}
    b3 = bundle.get("3") or {}
    out = []
    hgnc_id = b1.get("hgnc_id")
    if hgnc_id:
        out.append(f"https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/{hgnc_id}")
    entrez = (b1.get("entrez") or [None])[0]
    if entrez:
        out.append(f"https://www.ncbi.nlm.nih.gov/gene/{entrez}")
    canon = b3.get("canonical_uniprot")
    if canon:
        out.append(f"https://www.uniprot.org/uniprotkb/{canon}")
    ens = b1.get("ensembl_id")
    if ens:
        out.append(f"https://www.ensembl.org/Homo_sapiens/Gene/Summary?g={ens}")
    mims = b1.get("mim") or []
    if mims:
        out.append(f"https://www.omim.org/entry/{mims[0]}")
    return out

def _encodes(bundle):
    """encodesBioChemEntity for each reviewed UniProt product. Dual-product
    genes (CDKN2A) get a list; single-product gets a dict; ncRNA gets None."""
    b3 = bundle.get("3") or {}
    rev = b3.get("reviewed_uniprot") or []
    if not rev:
        return None
    proteins = [{
        "@type": "Protein",
        "name": u,
        "identifier": f"UniProtKB:{u}",
        "url": f"https://www.uniprot.org/uniprotkb/{u}",
    } for u in rev]
    return proteins[0] if len(proteins) == 1 else proteins

def build_jsonld(bundle, base_url=BASE_URL):
    """Compose the schema.org Gene JSON-LD dict from a full collector bundle.

    Skips fields whose value would be None/empty so the emitted JSON stays
    clean. The shape mirrors schema.org/Gene, with `sameAs` and
    `encodesBioChemEntity` doing the most semantic work."""
    b1 = bundle.get("1") or {}
    hgnc = b1.get("hgnc") or {}
    sym = b1.get("symbol") or "?"
    name = hgnc.get("name") or ""
    aliases = list(hgnc.get("aliases") or [])
    canonical_page = f"{base_url}/gene/{sym}/"

    alt_names = []
    if name and name != sym:
        alt_names.append(name)
    alt_names.extend(aliases)

    out = {
        "@context": "https://schema.org",
        "@type": "Gene",
        "@id": canonical_page,
        "name": sym,
        "identifier": b1.get("hgnc_id"),
        "url": canonical_page,
        "description": _strip_md(declarative_sentence(bundle)),
        "alternateName": alt_names or None,
        "sameAs": same_as_urls(bundle) or None,
        "encodesBioChemEntity": _encodes(bundle),
        "taxonomicRange": "https://www.ncbi.nlm.nih.gov/taxonomy/9606",
    }
    # Surface UniProt CC FUNCTION as a structured `disambiguatingDescription`
    # alongside the short declarative `description` — AI agents that want the
    # curated function paragraph (not just identifier facts) can pick it
    # without scraping the body. Source: biobtree uniprot.comments.function.
    function_cc = ((bundle.get("3") or {}).get("cc") or {}).get("function")
    if function_cc:
        out["disambiguatingDescription"] = function_cc
    loc = hgnc.get("location")
    if loc:
        out["isPartOfBioChemEntity"] = {"@type": "Chromosome", "name": loc}
    # Gene → disease / drug edges via JSON-LD @reverse: schema.org has no
    # forward Gene→MedicalCondition/Drug predicate, but a disease's
    # `associatedGene` and a drug's `target` both point *to* this gene — so we
    # emit them under @reverse with those real predicates. Targets are the
    # built Atlas pages (internal traversal); the block elides if none built.
    rev = _reverse_edges(bundle, base_url)
    if rev:
        out["@reverse"] = rev
    # Drop None/empty values for clean machine-readable output.
    return {k: v for k, v in out.items() if v not in (None, [], "")}


def _reverse_edges(bundle, base_url):
    """{@reverse: {associatedGene: [MedicalCondition…], target: [Drug…]}} —
    diseases/drugs whose edge points at this gene, limited to built Atlas
    pages with their internal URL."""
    from atlas.page import links
    host = base_url.rsplit("/atlas", 1)[0]
    groups = links.related_targets("gene", bundle)
    rev = {}
    dz = [{"@type": "MedicalCondition", "name": n, "url": host + p}
          for n, p in groups.get("Diseases", [])[:20]]
    dr = [{"@type": "Drug", "name": n, "url": host + p}
          for n, p in groups.get("Drugs", [])[:20]]
    if dz:
        rev["associatedGene"] = dz if len(dz) > 1 else dz[0]
    if dr:
        rev["target"] = dr if len(dr) > 1 else dr[0]
    return rev or None

def as_script_tag(jsonld):
    """JSON-LD as an inline <script> block for the page body."""
    body = json.dumps(jsonld, indent=2)
    return f'<script type="application/ld+json">\n{body}\n</script>'

def as_jsonld_string(jsonld):
    """Pretty-printed JSON-LD for the entity.jsonld sidecar file."""
    return json.dumps(jsonld, indent=2) + "\n"
