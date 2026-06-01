"""§1 — disease_ids: federated identifier set (Mondo + EFO + MeSH + OMIM +
Orphanet) + canonical name + the per-dataset xref count table.

NEW collector (not a gene fanout) — operates directly off the DiseaseAnchors
record, which already pre-resolved the ID set during resolve(). The work
here is mostly shaping; no new biobtree calls beyond the anchor."""
from atlas.section import Section

CHAINS   = (">>mondo>>efo", ">>mondo>>mesh", ">>mondo>>mim", ">>mondo>>orphanet",
            ">>mondo>>doid", ">>mondo>>sctid", ">>mondo>>umls", ">>mondo>>ncit",
            ">>mondo>>medgen", ">>mondo>>icd10cm", ">>mondo>>icd11",
            ">>mondo>>gard", ">>mondo>>meddra", ">>mondo>>nord",
            ">>mondo>>uberon")
DATASETS = ("mondo", "efo", "mesh", "mim", "orphanet",
            "doid", "sctid", "umls", "ncit", "medgen",
            "icd10cm", "icd11", "gard", "meddra", "nord", "uberon")

def collect(a):
    # HPO phenotypes from primary Orphanet entry — frequency-sorted desc so
    # render can slice the most clinically-relevant features first.
    oa = a.orphanet_attrs or {}
    phenotypes = list(oa.get("phenotypes") or [])
    phenotypes.sort(key=lambda p: float(p.get("frequency_value") or 0), reverse=True)

    bundle = {
        "section": "01_disease_ids",
        "name": a.name,
        "canonical_name": a.canonical_name,
        "mondo_id": a.mondo_id,
        "efo_id": a.efo_id,
        "mesh_ids": list(a.mesh_ids),
        "omim_ids": list(a.omim_ids),
        "orphanet_ids": list(a.orphanet_ids),
        # Cross-ontology xrefs from Mondo OBO ingest. {prefix: [ids,...]} —
        # only keys with data present.
        "obo_xrefs": {k: list(v) for k, v in (a.obo_xrefs or {}).items()},
        # UBERON anatomy ids (drives schema.org `associatedAnatomy`).
        "anatomy_uberon_ids": list(a.anatomy_uberon_ids or ()),
        # Orphanet primary-entry attrs — resolved once at anchor time.
        # prevalences = multi-geography epidemiology rows (drives JSON-LD
        # `epidemiology`). phenotypes = HPO list with both label and
        # numeric frequency_value (drives JSON-LD `signOrSymptom`).
        # Empty for non-rare-disease conditions (most cancers, common dz).
        "orphanet_name": oa.get("name") or "",
        "orphanet_disorder_type": oa.get("disorder_type") or "",
        "prevalences": list(oa.get("prevalences") or []),
        "phenotypes": phenotypes,
        "phenotype_count": oa.get("phenotype_count") or len(phenotypes),
        "is_cancer": a.is_cancer,
        "xref_counts": dict(a.xref_counts),
    }
    return bundle

SECTION = Section(
    id="1", name="disease_ids",
    description=("Federated disease identifiers (Mondo, EFO, MeSH, OMIM, "
                 "Orphanet) + canonical Mondo name + per-dataset xref counts "
                 "+ Orphanet epidemiology (prevalences) and clinical features "
                 "(HPO phenotype list with frequencies)."),
    needs=("mondo_id", "canonical_name", "efo_id", "mesh_ids", "omim_ids",
           "orphanet_ids", "orphanet_attrs", "obo_xrefs", "anatomy_uberon_ids",
           "xref_counts", "is_cancer"),
    produces=("mondo_id", "canonical_name", "efo_id", "mesh_ids", "omim_ids",
              "orphanet_ids", "obo_xrefs", "anatomy_uberon_ids",
              "orphanet_name", "orphanet_disorder_type",
              "prevalences", "phenotypes", "phenotype_count",
              "xref_counts", "is_cancer"),
    datasets=DATASETS, chains=CHAINS, collect_fn=collect,
)
