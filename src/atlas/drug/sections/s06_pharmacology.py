"""§6 — pharmacology. Drug-class semantics from ChEBI roles (open-licensed,
decoded to names like "tyrosine kinase inhibitor", "antineoplastic agent") +
the raw ATC code(s) (link-out only — biobtree exposes the code, not the
WHO-licensed level names). Both come from the anchor; no extra calls."""
from atlas.section import Section


def collect(a):
    return {
        "section": "06_pharmacology",
        "molecule_type": a.molecule_type,
        "chebi_roles": list(a.chebi_roles),
        "atc_codes": list(a.atc_codes),
    }


SECTION = Section(
    id="6", name="pharmacology",
    description=("Pharmacological class: ChEBI roles (decoded drug-class semantics) "
                 "+ raw ATC code(s) with whocc.no link-out (WHO ATC names are "
                 "licensing-restricted, not reproduced)"),
    needs=("chebi_roles", "atc_codes", "molecule_type"),
    produces=("molecule_type", "chebi_roles", "atc_codes"),
    datasets=("chebi", "chembl_molecule"),
    chains=(">>chembl_molecule>>chebi",),
    collect_fn=collect,
)
