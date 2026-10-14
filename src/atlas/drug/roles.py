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


# ── ATC class (DrugCentral, biobtree v2.17) ─────────────────────────────────
#
# ChEBI roles often name a side activity or an enzyme class rather than what
# the drug is for: aripiprazole "H1-receptor antagonist", fluconazole "P450
# inhibitor", ibrutinib "EC 2", amodiaquine "NSAID". WHO ATC is the therapeutic
# classification, so it decides unless ChEBI's class agrees with it.

import re as _re

# A route/site-restricted ATC group ranks after a systemic one (fluconazole:
# J02 "antimycotics for systemic use" before D01 "antifungals for topical use").
_LOCAL_ATC = _re.compile(r"TOPICAL|DERMATOLOG|OPHTHALM|OTOLOG|STOMATOLOG|GYNECOLOG|LOCAL ORAL|"
                         r"HEMORRHOIDS|NASAL|THROAT|DIAGNOSTIC", _re.I)
_KEEP_CAPS = {"I", "II", "III", "IV", "A", "B", "C", "D", "E", "K", "ACE", "HIV", "ADHD", "HMG",
              "GABA", "MAO", "DNA", "RNA", "COX", "CNS"}
_MECH = _re.compile(r"inhibitor|antagonist|agonist|blocker|modulator", _re.I)
_CLASS_STOP = {"agent", "agents", "drug", "drugs", "inhibitor", "inhibitors", "antagonist",
               "antagonists", "agonist", "agonists", "other", "and", "for", "of", "the", "use",
               "systemic", "plain", "excl", "incl", "with", "derivatives", "preparations",
               "products", "receptor", "selective", "non"}


def _caps(c):
    """An upper-case (level 1-3) ATC name — parentheticals ignored ("… (ARBs)")."""
    c = _re.sub(r"\(.*?\)", "", c)
    return c == c.upper()


def _atc_groups(classes):
    """Split DrugCentral's flattened ATC names into per-code groups: upper-case
    level 2/3 names, each group closed by its mixed-case level-4 name."""
    groups, cur = [], []
    for c in classes or ():
        cur.append(c)
        if not _caps(c):
            groups.append(cur)
            cur = []
    if cur:
        groups.append(cur)
    return groups


def _singular(w):
    lo = w.lower()
    if lo.endswith("ies"):
        return w[:-3] + ("Y" if w.isupper() else "y")
    if lo.endswith("s") and not lo.endswith(("ss", "us", "is")):
        return w[:-1]
    return w


def _phrase(name):
    """'OTHER ANTINEOPLASTIC AGENTS' → 'antineoplastic agent'; 'ANGIOTENSIN II
    RECEPTOR BLOCKERS (ARBs), PLAIN' → 'angiotensin II receptor blocker'."""
    n = _re.sub(r"\s*\(.*?\)", "", name)
    n = n.split(",")[0]
    n = _re.sub(r"\s+(FOR|USED|IN|EXCL\.?)\b.*$", "", n, flags=_re.I)
    n = _re.sub(r"^(ALL\s+)?OTHER\s+", "", n, flags=_re.I).strip()
    n = _re.sub(r"\s+AND\s+RELATED\s+\w+$", "", n, flags=_re.I)       # "… and related agents"
    parts = []
    for part in _re.split(r"\s+AND\s+", n, flags=_re.I):
        ws = part.split()
        if ws:
            ws[-1] = _singular(ws[-1])
        parts.append(" ".join(w if w.upper() in _KEEP_CAPS or (w != w.upper()) else w.lower()
                              for w in ws))
    return " and ".join(p for p in parts if p)


def atc_class(atc_classes):
    """(class phrase, mechanism phrase|None) — see _atc_pick."""
    return _atc_pick(atc_classes)[:2]


def _atc_pick(atc_classes):
    """(class phrase, mechanism phrase|None, route-restricted?) from the drug's
    ATC names; (None, None, False) when there are none. Class = the most specific upper-case (level-3)
    name of the first systemic code; mechanism = its level-4 name when that is
    a mechanism ("Protein kinase inhibitors"), not a chemical family."""
    groups = _atc_groups(atc_classes)
    if not groups:
        return None, None, False
    # Systemic before route-restricted; a code with a specific level-4 before an
    # "Other …" catch-all (celecoxib: M01AH "Coxibs" before L01XX "Other
    # antineoplastic agents"). Stable sort keeps DrugCentral's order otherwise.
    groups.sort(key=lambda g: (bool(_LOCAL_ATC.search(" ".join(g))),
                               any(c.lower().startswith("other") for c in g if not _caps(c))))
    for g in groups:
        caps = [c for c in g if _caps(c) and not _re.match(r"ALL OTHER", c, _re.I)]
        klass = None
        for c in reversed(caps):            # most specific informative upper-case name
            k = _ATC_NAME_OVERRIDE.get(c.upper(), _phrase(c))
            if k and not _re.search(r"\b[A-Z] and [A-Z]\b", k) and not _GENERIC_CLASS.search(k):
                klass = k
                break
        if not klass:
            continue
        l4 = next((c for c in g if not _caps(c)), None)
        mech = None
        if l4 and _MECH.search(l4) and not l4.lower().startswith("other"):
            mech = _phrase(l4)
            mech = mech[0].lower() + mech[1:] if mech[:2] != mech[:2].upper() else mech
        return klass, mech, bool(_LOCAL_ATC.search(" ".join(g)))
    return None, None, False


# ATC names whose mechanical singular reads badly.
_ATC_NAME_OVERRIDE = {
    "ANTIINFLAMMATORY AND ANTIRHEUMATIC PRODUCTS, NON-STEROIDS": "non-steroidal anti-inflammatory drug",
    "ANTIINFLAMMATORY AND ANTIRHEUMATIC PRODUCTS": "anti-inflammatory and antirheumatic drug",
}

# Phrases that name no class: "drug" (DRUGS FOR OBSTRUCTIVE AIRWAY DISEASES),
# "systemic drug", "alimentary tract and metabolism product".
_GENERIC_CLASS = _re.compile(r"^(systemic\s+)?(drug|agent|product|preparation|substance)$|"
                             r"\bproduct$|\btract\b", _re.I)


def _stems(s):
    return {w[:6] for w in _re.findall(r"[a-z0-9]+", (s or "").lower())
            if w not in _CLASS_STOP and len(w) > 2}


_JUNK_CLASS = _re.compile(r"^EC \d|geroprotector|prodrug|P450|metaboli|^drug$", _re.I)


def drug_class(roles, atc_classes, fallback=None):
    """Headline drug class: ChEBI's pharmacological role when it agrees with the
    ATC classification (shares a distinctive word: imatinib "tyrosine kinase
    inhibitor" ~ "Protein kinase inhibitors"), else the ATC class (amodiaquine
    "antimalarial", not ChEBI's "non-steroidal anti-inflammatory drug"), else
    ChEBI as before."""
    chebi = pharma_class(roles)
    klass, mech, local = _atc_pick(atc_classes)
    if not klass:
        return chebi or fallback
    if local and chebi and not _JUNK_CLASS.search(chebi):
        return chebi                      # only a topical/local ATC use (amlexanox)
    if chebi and not _JUNK_CLASS.search(chebi) and \
            _stems(chebi) & _stems(" ".join(atc_classes or ())):
        return chebi
    return f"{klass} ({mech})" if mech and _stems(mech) - _stems(klass) else klass
