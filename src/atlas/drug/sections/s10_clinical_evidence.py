"""§10 — clinical evidence (CIViC), drug-anchored. The third face of the same
drug × variant × indication triple gene §10 and disease §13 render: here the
drug is fixed and we surface which variant, in which indication, with what
effect. Reached via the ID-join `>>chembl_molecule>>civic_evidence` (NOT
therapy-name matching), fed through the shared aggregator. Heavy for targeted
therapies (Imatinib ~263 evidence rows); elides for non-precision drugs."""
from atlas.biobtree import map_all
from atlas.civic import aggregate_predictive
from atlas.section import Section


def collect(a):
    civic = map_all(a.chembl_id, ">>chembl_molecule>>civic_evidence", cap=15)
    ranked, stats = aggregate_predictive(civic, top=30)
    return {
        "section": "10_clinical_evidence",
        "civic_evidence": ranked,
        "civic_evidence_total": stats["evidence_total"],
        "civic_predictive_total": stats["predictive_total"],
        "civic_association_total": stats["association_total"],
        "civic_evidence_type_counts": stats["evidence_type_counts"],
    }


SECTION = Section(
    id="10", name="clinical_evidence",
    description=("CIViC clinical evidence, drug-anchored (variant × indication × "
                 "effect) via chembl_molecule→civic_evidence; predictive subset, "
                 "ranked by evidence level"),
    needs=("chembl_id",),
    produces=("civic_evidence", "civic_predictive_total", "civic_association_total"),
    datasets=("chembl_molecule", "civic_evidence"),
    chains=(">>chembl_molecule>>civic_evidence",),
    collect_fn=collect,
)
