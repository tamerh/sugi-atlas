"""Deterministic one-sentence lead for an Atlas drug page.

Mirror of atlas.page.declarative (gene) / disease_declarative. Pure function
of the drug bundle — every fact from biobtree, no LLM. Sits above the LLM exec
summary so it dominates the first-extracted text AI agents pick up.

Shape (sub-clauses drop silently when data is absent):

  **Imatinib** (CHEMBL941) is an approved small-molecule tyrosine kinase
  inhibitor (ATC L01EA01) targeting ABL1, DDR1, and DDR2; indicated across
  52 conditions including chronic myelogenous leukemia and gastrointestinal
  stromal tumor; with CIViC clinical evidence for 203 variant-indication
  associations (e.g. BCR::ABL1 fusion in chronic myeloid leukemia).
"""


def _join(names):
    names = [n for n in names if n]
    if not names:
        return ""
    if len(names) == 1:
        return names[0]
    if len(names) == 2:
        return " and ".join(names)
    return ", ".join(names[:-1]) + ", and " + names[-1]


def _class_clause(b1, b6):
    """'approved small-molecule tyrosine kinase inhibitor (ATC L01EA01)'."""
    bits = []
    if b1.get("is_fda_approved") or (b1.get("max_phase") == 4):
        bits.append("approved")
    mtype = (b1.get("molecule_type") or "").lower()
    roles = (b6 or {}).get("chebi_roles") or []
    # Lead with the first ChEBI role as the functional class; molecule type as
    # adjective (small-molecule / antibody). Avoid duplicating "agent".
    klass = roles[0] if roles else (mtype or "drug")
    if mtype and mtype not in klass:
        bits.append(mtype.replace("small molecule", "small-molecule"))
    bits.append(klass)
    s = " ".join(bits)
    atc = (b1.get("atc_codes") or [])
    if atc:
        s += f" (ATC {atc[0]})"
    return s


def _targets_clause(b2):
    prim = (b2 or {}).get("primary_targets") or []
    genes = [t.get("gene_symbol") for t in prim if t.get("gene_symbol")]
    if not genes:
        return ""
    return " targeting " + _join(genes[:3])


def _indications_clause(b4):
    n = (b4 or {}).get("indication_count") or 0
    inds = (b4 or {}).get("indications") or []
    if not n:
        return ""
    names = [(i.get("name") or "").strip().lower() for i in inds[:2]]
    names = [x for x in names if x and not x.startswith("mondo:")]
    s = f"; indicated across {n} condition" + ("s" if n != 1 else "")
    if names:
        s += " including " + _join(names)
    return s


def _civic_clause(b10):
    n = (b10 or {}).get("civic_association_total") or 0
    rows = (b10 or {}).get("civic_evidence") or []
    if not n:
        return ""
    s = (f"; with CIViC clinical evidence for {n} variant-indication "
         f"association" + ("s" if n != 1 else ""))
    if rows:
        top = rows[0]
        prof, dis = top.get("profile"), top.get("disease")
        if prof and dis:
            s += f" (e.g. {prof} in {dis.lower()})"
    return s


def declarative_sentence(bundle):
    """Compose the drug lead from a full bundle dict {section_id: bundle_dict}."""
    b1 = bundle.get("1") or {}
    b2 = bundle.get("2") or {}
    b4 = bundle.get("4") or {}
    b6 = bundle.get("6") or {}
    b10 = bundle.get("10") or {}

    name = b1.get("canonical_name") or "?"
    display = name.title() if name.isupper() else name
    chembl = b1.get("chembl_id")
    head = f"**{display}**" + (f" ({chembl})" if chembl else "")

    sentence = f"{head} is a{'n' if _class_clause(b1, b6)[:1].lower() in 'aeiou' else ''} "
    sentence += _class_clause(b1, b6)
    sentence += _targets_clause(b2)
    sentence += _indications_clause(b4)
    sentence += _civic_clause(b10)
    sentence += "."
    return sentence
