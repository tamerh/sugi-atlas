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
    if canon and not bundle.get("_noncoding"):   # ncRNA: a stray inherited UniProt isn't this gene's product (audit #12)
        out.append(f"https://www.uniprot.org/uniprotkb/{canon}")
    ens = b1.get("ensembl_id")
    if ens:
        out.append(f"https://www.ensembl.org/Homo_sapiens/Gene/Summary?g={ens}")
    mims = b1.get("mim") or []
    if mims:
        out.append(f"https://www.omim.org/entry/{mims[0]}")
    return out

# UniProt feature types that map to a druggable/functional residue — the
# Bioschemas `hasSequenceAnnotation` seed for the residue map (layer A). Order
# is display priority; the full set stays in the §3/§4 body (the residue map).
_RESIDUE_TYPES = ("active site", "binding site", "site", "modified residue",
                  "disulfide bond", "glycosylation site",
                  "lipid moiety-binding region", "mutagenesis site")
_ANN_CAP = 12      # representative residues in the inline graph, not the full list
_PDB_CAP = 6


def _to_int(v):
    try:
        return int(str(v).strip())
    except (TypeError, ValueError):
        return None


def _seq_annotations(ufeatures, acc, cap=_ANN_CAP):
    """Bioschemas SequenceAnnotation nodes for an accession's functional
    residues. `ufeatures` are already accession-stamped by s03. Sorted by
    druggability priority (_RESIDUE_TYPES order) so a gene with dozens of
    mutagenesis sites still surfaces its active/binding/PTM residues first."""
    prio = {t: i for i, t in enumerate(_RESIDUE_TYPES)}
    matched = [f for f in ufeatures
               if f.get("uniprot") == acc and f.get("type") in prio]
    matched.sort(key=lambda f: (prio[f["type"]], _to_int(f.get("begin")) or 0))
    out = []
    for f in matched:
        node = {"@type": "SequenceAnnotation", "name": f.get("type")}
        if f.get("description"):
            node["description"] = f["description"]
        start = _to_int(f.get("begin"))
        if start is not None:
            end = _to_int(f.get("end"))
            node["sequenceLocation"] = {"@type": "SequenceRange",
                                        "rangeStart": start,
                                        "rangeEnd": end if end is not None else start}
        out.append(node)
        if len(out) >= cap:
            break
    return out


def _representations(b4, acc, is_canonical, cap=_PDB_CAP):
    """Bioschemas hasRepresentation — AlphaFold (accession-tagged) + a sample of
    PDB structures. PDB is a flat list not yet accession-tagged, so it attaches
    to the canonical product (single-product genes: correct; dual-product: a
    documented approximation pending per-accession PDB tagging in s04)."""
    out = []
    for af in (b4.get("alphafold") or []):
        if af.get("uniprot") == acc and af.get("present"):
            out.append({"@type": "3DModel", "name": f"AlphaFold {af.get('id')}",
                        "url": f"https://alphafold.ebi.ac.uk/entry/{acc}"})
    if is_canonical:
        for p in (b4.get("pdb") or [])[:cap]:
            pid = p.get("id")
            if pid:
                out.append({"@type": "3DModel", "name": f"PDB {pid}",
                            "url": f"https://www.rcsb.org/structure/{pid}"})
    return out


def _protein_nodes(bundle, gene_page):
    """encodesBioChemEntity — one fully-typed, addressable Protein node per
    reviewed UniProt product (layer A). Each carries its own `@id`
    (…/gene/SYM/#protein-<acc>), the reciprocal `isEncodedByBioChemEntity` edge
    back to the Gene, and a Bioschemas Protein-profile seed (sequence
    annotations = the residue map, structural representations). Dual-product
    genes (CDKN2A) get a list; single-product a dict; ncRNA None."""
    b3 = bundle.get("3") or {}
    b4 = bundle.get("4") or {}
    # ncRNA has no protein product; a non-coding gene can carry a stray,
    # positionally-inherited reviewed UniProt (SCP2D1-AS1 → SCP2D1's Q9BR46).
    # Emitting a Protein node then contradicts the body's "non-coding, no protein"
    # (audit #12) — so gate on the same _noncoding flag the body renders from.
    if bundle.get("_noncoding"):
        return None
    rev = b3.get("reviewed_uniprot") or []
    if not rev:
        return None
    canon = b3.get("canonical_uniprot")
    ufeatures = b3.get("ufeatures") or []
    nodes = []
    for acc in rev:
        is_canon = (acc == canon) or (len(rev) == 1)
        # Per-accession protein names aren't in the bundle yet (only the
        # canonical's protein_name); non-canonical products fall back to the
        # accession until s03 carries per-product names.
        pname = (b3.get("protein_name") if is_canon else None) or acc
        node = {
            "@type": "Protein",
            "@id": f"{gene_page}#protein-{acc}",
            "name": pname,
            "identifier": f"UniProtKB:{acc}",
            "url": f"https://www.uniprot.org/uniprotkb/{acc}",
            "sameAs": f"https://www.uniprot.org/uniprotkb/{acc}",
            "taxonomicRange": "https://www.ncbi.nlm.nih.gov/taxonomy/9606",
            "isEncodedByBioChemEntity": {"@id": gene_page},
            "hasSequenceAnnotation": _seq_annotations(ufeatures, acc) or None,
            "hasRepresentation": _representations(b4, acc, is_canon) or None,
        }
        nodes.append({k: v for k, v in node.items() if v not in (None, [], "")})
    return nodes[0] if len(nodes) == 1 else nodes

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
        # Each reviewed product is a fully-typed Protein node with its own @id —
        # the molecular layer (residue map, structure, binding) hangs off the
        # Protein, not the Gene (layer A; see docs/MOLECULAR_ENRICHMENT.md).
        "encodesBioChemEntity": _protein_nodes(bundle, canonical_page),
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
    dr = [{"@type": "Substance", "name": n, "url": host + p}  # not Drug: Drug ⊂ Product trips Google's product validator
          for n, p in groups.get("Drugs", [])[:20]]
    if dz:
        rev["associatedGene"] = dz if len(dz) > 1 else dz[0]
    if dr:
        rev["target"] = dr if len(dr) > 1 else dr[0]
    return rev or None

def as_script_tag(jsonld):
    """JSON-LD as an inline <script> block for the page body.

    Compacted (audit #6): over-long arrays (sameAs, @reverse disease/drug
    edges) are capped for the inline copy; the full graph stays in the
    entity.jsonld sidecar."""
    from atlas.page.jsonld_inline import compact_for_inline
    body = json.dumps(compact_for_inline(jsonld), indent=2)
    return f'<script type="application/ld+json">\n{body}\n</script>'

def as_jsonld_string(jsonld):
    """Pretty-printed JSON-LD for the entity.jsonld sidecar file."""
    return json.dumps(jsonld, indent=2) + "\n"
