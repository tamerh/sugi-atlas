"""Per-section registry for drug pages — mirrors atlas.gene.sections /
atlas.disease.sections.

Each s<NN>_*.py module exposes a SECTION (atlas.section.Section). The
orchestrator in atlas.drug.collect imports REGISTRY and runs each section
against a single shared atlas.drug.anchors.DrugAnchors record.

Sections land incrementally per the spec's sequencing. Present so far:
  §1  drug_ids            (NEW, anchor read + chemistry)
  §3  bioactivity         (NEW, chembl_activity)
  §5  clinical_trials     (NEW, chembl_molecule→clinical_trials)
  §6  pharmacology        (NEW, ChEBI roles + ATC)
  §11 patent_literature   (NEW, patent_compound)
Still to add: §2 targets, §4 indications, §7 related_molecules,
§8 target_pathways, §9 pharmacogenomics, §10 clinical_evidence (CIViC),
§12 salt_forms_and_parent.
"""
from atlas.section import Section
from atlas.drug.sections import (
    s01_drug_ids, s03_bioactivity, s05_clinical_trials,
    s06_pharmacology, s11_patent_literature,
)

_MODULES = (
    s01_drug_ids, s03_bioactivity, s05_clinical_trials,
    s06_pharmacology, s11_patent_literature,
)

REGISTRY: "dict[str, Section]" = {m.SECTION.id: m.SECTION for m in _MODULES}

__all__ = ["Section", "REGISTRY"]
