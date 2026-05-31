"""Per-page provenance sidecar for drug pages — mirror of
atlas.page.provenance (gene) / disease_provenance. Uses the shared UPSTREAM
dict and the drug REGISTRY for section metadata."""
from atlas.drug.sections import REGISTRY
from atlas.page.provenance import UPSTREAM, BASE_URL, _section_provenance


def build_provenance(bundle: dict, slug: str, meta: dict = None,
                     base_url: str = BASE_URL) -> dict:
    """Compose the provenance Dataset blob for a drug page.

    bundle: full collector bundle ({section_id: bundle_dict}).
    slug:   URL-safe drug slug used in the page URL.
    meta:   optional dict with generated_at / atlas_version / biobtree_version.
    """
    meta = meta or {}
    b1 = bundle.get("1") or {}
    b2 = bundle.get("2") or {}
    name = b1.get("canonical_name") or slug
    page_url = f"{base_url}/drug/{slug}/"

    sections = [_section_provenance(REGISTRY[sid]) for sid in REGISTRY]
    primary_targets = [t.get("gene_symbol") for t in (b2.get("primary_targets") or [])
                       if t.get("gene_symbol")]

    return {
        "@context": "https://schema.org",
        "@type": "Dataset",
        "name": f"{name} — provenance trail",
        "description": (
            "Per-section provenance for Atlas's drug page: datasets touched, "
            "biobtree chains used, upstream-source URLs. Every fact in the "
            "page can be traced back through this trail to its primary source."
        ),
        "isPartOf": page_url,
        "url": f"{page_url}provenance.json",
        "generated_at": meta.get("generated_at"),
        "atlas_version": meta.get("atlas_version"),
        "biobtree_version": meta.get("biobtree_version"),
        "anchors": {
            "name": name,
            "slug": slug,
            "chembl_id": b1.get("chembl_id"),
            "pubchem_cid": b1.get("pubchem_cid"),
            "chebi_id": b1.get("chebi_id"),
            "atc_codes": b1.get("atc_codes") or [],
            "molecule_type": b1.get("molecule_type"),
            "max_phase": b1.get("max_phase"),
            "primary_targets": primary_targets,
            "parent_chembl": b1.get("parent_chembl"),
        },
        "data_access": {
            "biobtree_project": "https://biobtree.org",
            "chain_syntax": "https://biobtree.org",
            "note": ("Chains in each section's `chains` field can be replayed "
                     "against any biobtree instance to verify the source data."),
        },
        "sections": sections,
    }


def as_provenance_string(prov: dict) -> str:
    import json
    return json.dumps(prov, indent=2) + "\n"
