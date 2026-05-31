"""Build the deterministic one-sentence lead for an Atlas disease page.

Mirror of atlas.page.declarative for the disease entity type. Pure function
of the disease bundle — every fact sourced from biobtree; no LLM, no network
beyond the collect step.

Sits above the LLM exec summary so it dominates the first-extracted text
that AI agents pick up.

Shape (filled with whichever sub-clauses have data — missing fields drop
silently rather than emitting "None"):

  **Endometrial carcinoma** (MONDO:0011962) is a cancer with 91 cohort
  genes (113 GWAS associations across 26 studies, plus 33 CIViC-evidence
  somatic drivers and 4 ClinVar predisposition genes) and 1,034 clinical
  trials. Top therapeutic interventions include megestrol acetate,
  dostarlimab, and medroxyprogesterone acetate.
"""

_CANCER_PHRASE = "cancer"


def _disease_class(b1):
    """Plain-English class label. Cancers get a more specific phrase; everything
    else stays as 'disease' (Mondo doesn't expose a coarser-grained category in
    its entry attrs, so refining further would require ontology walks)."""
    return _CANCER_PHRASE if b1.get("is_cancer") else "disease"


def _format_int(n):
    try:
        return f"{int(n):,}"
    except (TypeError, ValueError):
        return None


def _evidence_clause(b1, b2, b4, b5):
    """Compose the parenthetical 'with X cohort genes (Y GWAS / Z somatic / ...)'
    clause from §1/§2/§4/§5. Returns '' when no countable evidence is present."""
    cohort_n = (b5 or {}).get("gene_count") or 0
    sub_clauses = []

    gwas_assoc = (b2 or {}).get("assoc_total")
    gwas_studies = (b2 or {}).get("study_total")
    if gwas_assoc:
        s = f"{_format_int(gwas_assoc)} GWAS associations"
        if gwas_studies:
            s += f" across {_format_int(gwas_studies)} studies"
        sub_clauses.append(s)

    somatic_n = len((b4 or {}).get("somatic_driver_genes") or [])
    if somatic_n:
        sub_clauses.append(f"{somatic_n} CIViC-evidence somatic driver"
                           + ("s" if somatic_n != 1 else ""))

    clinvar_n = (b1 or {}).get("xref_counts", {}).get("clinvar") or 0
    # Predisposition vs somatic — for cancers, clinvar is the germline-rare
    # predisposition set distinct from civic_evidence somatic.
    if clinvar_n and (b1 or {}).get("is_cancer") and somatic_n:
        sub_clauses.append(f"{_format_int(clinvar_n)} ClinVar predisposition record"
                           + ("s" if clinvar_n != 1 else ""))

    if not cohort_n and not sub_clauses:
        return ""
    parts = []
    if cohort_n:
        parts.append(f"{cohort_n} cohort gene" + ("s" if cohort_n != 1 else ""))
    if sub_clauses:
        parens = "; ".join(sub_clauses)
        return " with " + (parts[0] + f" ({parens})" if parts else parens)
    return " with " + parts[0]


def _trials_clause(b1, b13):
    n = (b1 or {}).get("xref_counts", {}).get("clinical_trials") or \
        (b13 or {}).get("trial_count") or 0
    if not n:
        return ""
    return f" and {_format_int(n)} clinical trial" + ("s" if n != 1 else "")


def _drugs_clause(b13):
    drugs = (b13 or {}).get("trial_drugs") or []
    if not drugs:
        return ""
    # Pick top 3 by trial_count (already sorted in the bundle).
    names = [(d.get("name") or "").strip() for d in drugs[:3]]
    names = [n.lower() for n in names if n]
    if not names:
        return ""
    if len(names) == 1:
        joined = names[0]
    elif len(names) == 2:
        joined = " and ".join(names)
    else:
        joined = ", ".join(names[:-1]) + ", and " + names[-1]
    return f" Top therapeutic interventions include {joined}."


def _pathway_clause(b14):
    """When there's a clearly-dominant pathway (≥3 cohort genes), name it."""
    tp = (b14 or {}).get("top_pathways") or []
    if not tp:
        return ""
    top = tp[0]
    if (top.get("gene_count") or 0) < 3:
        return ""
    name = (top.get("name") or "").strip()
    if not name:
        return ""
    # Avoid Reactome's broad-umbrella names that aren't informative.
    NOISE = {"disease", "signal transduction", "metabolism", "metabolism of proteins"}
    if name.lower() in NOISE:
        return ""
    return f" The dominant Reactome pathway is *{name}* ({top['gene_count']} cohort genes)."


def declarative_sentence(bundle):
    """Compose the disease lead from a full bundle dict {section_id: bundle_dict}.

    Sources used:
      bundle["1"]  — canonical_name, mondo_id, is_cancer, xref_counts (clinvar, clinical_trials)
      bundle["2"]  — assoc_total, study_total
      bundle["4"]  — somatic_driver_genes (cancer only)
      bundle["5"]  — gene_count, cohort_function_summary
      bundle["13"] — trial_drugs (top therapeutic interventions)
      bundle["14"] — top_pathways (dominant pathway, when clearly so)

    Returns a single Markdown sentence (no trailing newline)."""
    b1  = bundle.get("1")  or {}
    b2  = bundle.get("2")  or {}
    b4  = bundle.get("4")  or {}
    b5  = bundle.get("5")  or {}
    b13 = bundle.get("13") or {}
    b14 = bundle.get("14") or {}

    name = b1.get("canonical_name") or b1.get("name") or "?"
    mondo_id = b1.get("mondo_id")
    klass = _disease_class(b1)

    parens = [mondo_id] if mondo_id else []
    head = f"**{name.capitalize() if name == name.lower() else name}**"
    if parens:
        head += f" ({', '.join(parens)})"

    sentence = f"{head} is a {klass}"
    sentence += _evidence_clause(b1, b2, b4, b5)
    sentence += _trials_clause(b1, b13)
    sentence += "."
    sentence += _pathway_clause(b14)
    sentence += _drugs_clause(b13)
    return sentence
