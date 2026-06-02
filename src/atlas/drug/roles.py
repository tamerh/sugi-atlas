"""Pharmacological-class selection from ChEBI roles (audit #9).

A ChEBI role list is arbitrarily ordered and mixes pharmacology with
toxicology / environmental / metabolic annotations — atorvastatin's roles are
``['environmental contaminant', 'xenobiotic']``; caffeine leads with
``'mutagen'`` before ``'central nervous system stimulant'``. Leading the
drug-class / MoA with ``roles[0]`` therefore prints nonsense.

``pharma_class`` skips the non-pharmacological roles and prefers a clearly
pharmacological one (a mechanism/class keyword), falling back to the molecule
type when ChEBI offers nothing usable. The full role list is still shown in the
pharmacology section; this only governs the single headline class/MoA label.
"""

# Roles that are never a drug class: environmental / metabolic / lab-reagent
# (same spirit as the corpus gate's _NONPHARMA_ROLE denylist) plus toxicology
# hazard descriptors that are true of many drugs but aren't a therapeutic class.
_NON_CLASS = (
    "solvent", "metabolite", "greenhouse gas", "fertilis", "fertiliz", "fuel",
    "reference compound", "nmr", "chemical shift", "food", "nutrient",
    "contaminant", "pollutant", "reagent", "dye", "stain", "buffer",
    "cosmetic", "flavour", "flavor", "fragrance",
    "xenobiotic", "mutagen", "teratogen", "carcinogen", "toxin", "poison",
    "allergen", "irritant",
)

# Signals that a role IS a pharmacological class / mechanism of action.
_PHARMA_KW = (
    "inhibitor", "antagonist", "agonist", "blocker", "modulator", "agent",
    "stimulant", "depressant", "drug", "antibiotic", "antineoplastic",
    "anti-inflammatory", "analgesic", "anaesthetic", "anesthetic",
    "vasodilator", "vasoconstrictor", "diuretic", "antihypertensive",
    "anticoagulant", "antidepressant", "antipsychotic", "anticonvulsant",
    "chelator", "ligand", "potentiator", "activator", "antibody",
    "vaccine", "hormone",
)


def _displayable(role: str) -> bool:
    r = (role or "").strip().lower()
    return bool(r) and not any(k in r for k in _NON_CLASS)


def pharma_class(roles, fallback=None):
    """Best pharmacological-class label from a ChEBI role list.

    Drops non-pharmacological roles; among the rest prefers one carrying a
    pharmacological keyword (mechanism/class), else the first displayable role.
    Returns ``fallback`` (molecule type / 'drug' / None) when nothing qualifies
    — better to omit the class than to print 'environmental contaminant'."""
    disp = [r for r in (roles or []) if _displayable(r)]
    if not disp:
        return fallback
    pharma = [r for r in disp if any(k in r.lower() for k in _PHARMA_KW)]
    return (pharma or disp)[0]
