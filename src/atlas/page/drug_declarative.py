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


def _class_clause(b1, b6, b5=None):
    """'approved small-molecule tyrosine kinase inhibitor (ATC L01EA01)'."""
    from atlas.indication import molecule_approved
    bits = []
    # Don't call a drug "approved" in the lead when the FDA flag explicitly says no
    # (audit: 926 phase-4 pages had FDA=no — non-US/historical/withdrawn, e.g. the
    # withdrawn MAOI phenoxypropazine). Omit the approval adjective in that case; the
    # nuanced status (reached phase 4, not FDA-approved) lives in At-a-glance + IDs.
    if molecule_approved(b1, b5) and b1.get("is_fda_approved") is not False:
        bits.append("approved")
    elif b1.get("max_phase") == 3:
        bits.append("phase-3 clinical-stage")
    mtype = (b1.get("molecule_type") or "").lower()
    roles = (b6 or {}).get("chebi_roles") or []
    # Lead with the best PHARMACOLOGICAL ChEBI role as the functional class
    # (audit #9: roles[0] can be 'environmental contaminant'/'mutagen'); molecule
    # type as adjective (small-molecule / antibody). Avoid duplicating "agent".
    from atlas.drug.roles import pharma_class
    klass = pharma_class(roles, fallback=(mtype or "drug"))
    if mtype and mtype not in klass:
        bits.append(mtype.replace("small molecule", "small-molecule"))
    bits.append(klass)
    s = " ".join(bits)
    atc = (b1.get("atc_codes") or [])
    if atc:
        s += f" (ATC {atc[0]})"
    return s


def _curated_target_genes(b2):
    """Target genes for the lead / at-a-glance: ChEMBL curated mechanism-of-action
    genes FIRST (the drug's actual mechanism targets — imatinib → ABL1, PDGFRB,
    KIT; adalimumab → TNF), then any GtoPdb primary targets not already listed.
    Previously GtoPdb led and MOA genes were only a fallback, so imatinib read
    "targeting DDR1, DDR2, and ABL1" and adalimumab named no target at all."""
    b2 = b2 or {}
    out = []
    for g in (b2.get("mechanism_genes") or []):
        sym = g.get("gene_symbol")
        if sym and sym not in out:
            out.append(sym)
    def _aff(t):
        try:
            return -float(t.get("affinity"))
        except (TypeError, ValueError):
            return 0.0
    # GtoPdb rows arrive unordered; strongest affinity first (diazepam led with
    # TRHR at pKi 5.2 ahead of its GABA-A α1-3 targets at ~7.8).
    for t in sorted(b2.get("primary_targets") or [], key=_aff):
        sym = t.get("gene_symbol")
        if sym and sym not in out:
            out.append(sym)
    return out


def _targets_clause(b2):
    genes = _curated_target_genes(b2)
    if not genes:
        return ""
    return " targeting " + _join(genes[:3])


def _near(a, b):
    """True if one of a/b is the other followed by a word boundary — i.e. a
    parent and its qualified child ("osteoarthritis" vs "osteoarthritis, knee",
    "anxiety" vs "anxiety disorder"). Used to keep the lead's two named diseases
    visibly distinct."""
    lo, hi = sorted((a, b), key=len)
    return hi.startswith(lo) and (len(hi) == len(lo) or not hi[len(lo)].isalnum())


def _indications_clause(b4):
    n = (b4 or {}).get("indication_count") or 0
    inds = (b4 or {}).get("indications") or []
    if not n:
        return ""
    # Dedup by name (audit #11: indications cross-walked from both EFO and MeSH
    # can resolve to the same disease → "including neoplasm and neoplasm"); skip
    # blank and raw ontology-id names. Also skip a name that is a prefix-relative
    # of one already chosen (e.g. "anxiety" vs "anxiety disorder", "bacterial
    # infectious disease" vs "…with sepsis") — the lead should name two visibly
    # DISTINCT diseases, not a parent and its child.
    from atlas.render_common import is_ontology_id
    seen, names = set(), []
    for i in inds:
        nm = (i.get("name") or "").strip().lower()
        if not nm or is_ontology_id(nm) or nm in seen:
            continue
        if any(_near(nm, c) for c in names):     # "osteoarthritis" vs
            continue                             # "osteoarthritis, knee" etc.
        seen.add(nm)
        names.append(nm)
        if len(names) == 2:
            break
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
    b5 = bundle.get("5") or {}
    b6 = bundle.get("6") or {}
    b10 = bundle.get("10") or {}

    name = b1.get("canonical_name") or "?"
    display = name.title() if name.isupper() else name
    chembl = b1.get("chembl_id")
    head = f"**{display}**" + (f" ({chembl})" if chembl else "")

    sentence = f"{head} is a{'n' if _class_clause(b1, b6, b5)[:1].lower() in 'aeiou' else ''} "
    sentence += _class_clause(b1, b6, b5)
    sentence += _targets_clause(b2)
    sentence += _indications_clause(b4)
    sentence += _civic_clause(b10)
    sentence += "."
    return sentence
