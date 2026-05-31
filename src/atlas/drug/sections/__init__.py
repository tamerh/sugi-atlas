"""Per-section registry for drug pages — mirrors atlas.gene.sections /
atlas.disease.sections.

Each s<NN>_*.py module exposes a SECTION (atlas.section.Section). The
orchestrator in atlas.drug.collect imports REGISTRY and runs each section
against a single shared atlas.drug.anchors.DrugAnchors record.

Sections land incrementally per the spec's sequencing. Present so far:
  §1  drug_ids            (NEW, anchor read + chemistry)
  §2  targets             (NEW, GtoPdb curated + chembl bioactivity)
  §3  bioactivity         (NEW, chembl_activity)
  §4  indications         (NEW, anchor read, efo/mesh→mondo)
  §5  clinical_trials     (NEW, chembl_molecule→clinical_trials)
  §6  pharmacology        (NEW, ChEBI roles + ATC)
  §10 clinical_evidence   (NEW, CIViC via chembl_molecule→civic_evidence)
  §11 patent_literature   (NEW, patent_compound)
  §12 salt_forms          (NEW, anchor read, parent/child nav)
Still to add: §7 related_molecules, §8 target_pathways (cohort fan-out),
§9 pharmacogenomics (partly blocked on biobtree #13).
"""
from atlas.section import Section
from atlas.drug.sections import (
    s01_drug_ids, s02_targets, s03_bioactivity, s04_indications,
    s05_clinical_trials, s06_pharmacology, s10_clinical_evidence,
    s11_patent_literature, s12_salt_forms,
)

_MODULES = (
    s01_drug_ids, s02_targets, s03_bioactivity, s04_indications,
    s05_clinical_trials, s06_pharmacology, s10_clinical_evidence,
    s11_patent_literature, s12_salt_forms,
)

REGISTRY: "dict[str, Section]" = {m.SECTION.id: m.SECTION for m in _MODULES}

__all__ = ["Section", "REGISTRY"]
