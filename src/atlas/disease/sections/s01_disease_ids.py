"""§1 — disease_ids: federated identifier set (Mondo + EFO + MeSH + OMIM +
Orphanet) + canonical name + the per-dataset xref count table.

NEW collector (not a gene fanout) — operates directly off the DiseaseAnchors
record, which already pre-resolved the ID set during resolve(). The work
here is mostly shaping; no new biobtree calls beyond the anchor."""
from atlas.section import Section

CHAINS   = (">>mondo>>efo", ">>mondo>>mesh", ">>mondo>>mim", ">>mondo>>orphanet")
DATASETS = ("mondo", "efo", "mesh", "mim", "orphanet")

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
           "orphanet_ids", "orphanet_attrs", "xref_counts", "is_cancer"),
    produces=("mondo_id", "canonical_name", "efo_id", "mesh_ids", "omim_ids",
              "orphanet_ids", "orphanet_name", "orphanet_disorder_type",
              "prevalences", "phenotypes", "phenotype_count",
              "xref_counts", "is_cancer"),
    datasets=DATASETS, chains=CHAINS, collect_fn=collect,
)
