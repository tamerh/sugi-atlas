"""Per-page provenance sidecar for disease pages — mirror of
atlas.page.provenance for gene. Uses the shared UPSTREAM dict and the
disease REGISTRY for section metadata."""
from atlas.disease.sections import REGISTRY
from atlas.page.provenance import UPSTREAM, BASE_URL, _section_provenance


def build_provenance(bundle: dict, slug: str, meta: dict = None,
                     base_url: str = BASE_URL) -> dict:
    """Compose the provenance Dataset blob for a disease page.

    bundle: full collector bundle ({section_id: bundle_dict}).
    slug:   URL-safe disease slug used in the page URL.
    meta:   optional dict with generated_at / atlas_version / biobtree_version.
    """
    meta = meta or {}
    b1 = bundle.get("1") or {}
    name = b1.get("canonical_name") or b1.get("name") or slug
    page_url = f"{base_url}/disease/{slug}/"

    sections = [_section_provenance(REGISTRY[sid]) for sid in REGISTRY]

    return {
        "@context": "https://schema.org",
        "@type": "Dataset",
        "name": f"{name} — provenance trail",
        "description": (
            "Per-section provenance for Atlas's disease page: datasets touched, "
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
            "mondo_id": b1.get("mondo_id"),
            "efo_id": b1.get("efo_id"),
            "mesh_ids": b1.get("mesh_ids") or [],
            "omim_ids": b1.get("omim_ids") or [],
            "orphanet_ids": b1.get("orphanet_ids") or [],
            "is_cancer": b1.get("is_cancer"),
            "cohort_size": len((bundle.get("5") or {}).get("genes") or []),
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
    """Pretty-printed JSON for the provenance.json sidecar."""
    import json
    return json.dumps(prov, indent=2) + "\n"
