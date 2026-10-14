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

import re

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


def _trials_clause(b1, b13, joined=False):
    # §13 trial_count is the authoritative total off the now exact-match
    # mondo→clinical_trials edge. It may legitimately be 0 for a rare disease,
    # which is the correct answer — so don't treat 0 as "missing".
    # `joined`: a preceding evidence clause ("… with N cohort genes") is present,
    # so continue its list with "and". Otherwise open a new "with …" — writing
    # " and N trials" straight after "is a disease" reads as a category error
    # ("a disease and 193 trials").
    if b13 is None or "trial_count" not in b13:
        return ""
    n = b13.get("trial_count") or 0
    if not n:
        return ""
    lead = " and " if joined else " with "
    return f"{lead}{_format_int(n)} registered clinical trial" + ("s" if n != 1 else "")


# Perioperative / supportive / diagnostic agents that top trial counts for
# procedures done ON patients with the disease (tonsillectomy → dexmedetomidine),
# not treatments OF it. Name-matched: trial_drugs carry no ATC code.
_PROCEDURAL = re.compile(
    r"^(dexmedetomidine|propofol|midazolam|ketamine|esketamine|fentanyl|remifentanil|"
    r"sufentanil|alfentanil|morphine|sevoflurane|desflurane|isoflurane|lidocaine|bupivacaine|"
    r"levobupivacaine|ropivacaine|mepivacaine|ondansetron|granisetron|palonosetron|"
    r"droperidol|sugammadex|rocuronium|cisatracurium|vecuronium|neostigmine|"
    r"sodium chloride|dextrose|glucose|water|placebo|fluorodeoxyglucose.*|"
    r"florbetapir.*|flutemetamol.*|florbetaben.*|gadolinium.*|gado\w+|iohexol|iodixanol|"
    r"technetium.*|indocyanine green|fluorescein.*|methylene blue)\b", re.I)


def _name_list(names):
    if len(names) == 1:
        return names[0]
    if len(names) == 2:
        return " and ".join(names)
    return ", ".join(names[:-1]) + ", and " + names[-1]


def _drugs_clause(b13, indicated=None):
    """Treatment sentence. Approved drugs first (disease-direct phase-4 ChEMBL
    indications), ranked by how many of THIS disease's trials use them; else
    approved drugs with phase-3 evidence here; else the drugs most tested in its
    trials — labelled as such, never as "therapeutic interventions" (v1.12 put an
    anaesthetic first on tonsillitis and listed trial drugs as treatments on
    3,896 pages)."""
    from atlas.render_common import display_name
    trials = {(d.get("name") or "").strip().lower(): d.get("trial_count") or 0
              for d in (b13 or {}).get("trial_drugs") or []}
    rows = [r for r in (indicated or []) if r.get("name")]

    def _top(grp, k=3):
        grp = sorted(grp, key=lambda r: (-trials.get(r["name"].strip().lower(), 0), r["name"]))
        return [display_name(r["name"]).lower() for r in grp[:k]]

    appr = [r for r in rows if r.get("approved")]
    late = [r for r in rows if not r.get("approved") and r.get("late")]
    if appr:
        return f" Approved drugs (ChEMBL) include {_name_list(_top(appr))}."
    if late:
        return (f" Approved drugs with phase-3 evidence for it include "
                f"{_name_list(_top(late))}.")
    drugs = [d for d in (b13 or {}).get("trial_drugs") or []
             if (d.get("name") or "").strip() and not _PROCEDURAL.match(d["name"].strip())]
    names = [display_name(d["name"].strip()).lower() for d in drugs[:3]]
    if not names:
        return ""
    return f" Drugs most tested in its clinical trials include {_name_list(names)}."


# Polygenic test, calibrated on real bundles (biobtree v2.12.0): Mendelian
# diseases carry up to ~60 GWAS associations (sickle cell 60, DMD 47, CF 37) but
# are anchored by a GenCC-Definitive gene; complex traits start at ~60 (epilepsy
# 61 → SHROOM4, ALS 169, SLE 1,562). So: ≥100 GWAS hits is polygenic outright;
# ≥20 is polygenic unless a causal gene is GenCC Definitive.
_POLYGENIC_GWAS = 100
_POLYGENIC_GWAS_SOFT = 20
_UMBRELLA_SUBTYPES = 5    # ≥ this many Mondo subtypes → an umbrella/grouping term


def _gene_list(syms, cg):
    if len(syms) == 1:
        return syms[0]
    if len(syms) == 2:
        return f"{syms[0]} and {syms[1]}"
    if len(syms) == 3:
        return f"{syms[0]}, {syms[1]}, and {syms[2]}"
    n = len(syms) - 3
    return f"{syms[0]}, {syms[1]}, {syms[2]}, and {n} other gene{'s' if n != 1 else ''}"


def causal_kind(bundle):
    """(kind, causal genes) — how the curated genes relate to THIS disease:
      "caused"    a specific, non-polygenic disease: "caused by CFTR"
      "forms"     a complex trait whose curated genes are its rare monogenic
                  (or, for cancers, hereditary) forms
      "subtypes"  an umbrella term (≥5 Mondo subtypes)
      "clinvar"   pathogenic ClinVar submissions only (no curated validity call)
      None        no causal genes
    One decision shared by the lead, the inheritance fallback and the FAQ, so a
    polygenic disease is never "caused by" / "inherited" anywhere on the page."""
    from atlas.disease.cohort import causal_genes
    cg = causal_genes(bundle)
    if not cg:
        return None, cg
    b1 = bundle.get("1") or {}
    n_gwas = (bundle.get("2") or {}).get("assoc_total") or 0
    definitive = any("definitive" in (ev or "").lower() for _, ev in cg)
    if n_gwas >= _POLYGENIC_GWAS or (n_gwas >= _POLYGENIC_GWAS_SOFT and not definitive):
        return "forms", cg
    if (b1.get("child_count") or 0) >= _UMBRELLA_SUBTYPES:
        return "subtypes", cg
    if all(lbl == "ClinVar-linked" for _, lbl in cg):
        return "clinvar", cg
    return "caused", cg


def _causal_clause(bundle):
    """Causal-gene clause for the lead, driven by causal_kind(). '' when none.

    v1.11.7 read "SLE … caused by variants in SAT1 and TREX1"; a complex trait's
    curated genes are its rare monogenic forms, an umbrella's belong to its
    subtypes, and a ClinVar-only gene is reported as evidence, not cause."""
    kind, cg = causal_kind(bundle)
    if not kind:
        return _subtype_gene_clause(bundle)
    syms = [s for s, _ in cg]
    genes = _gene_list(syms, cg)
    if kind == "forms":
        # Cancers: the curated genes are germline predisposition genes →
        # "hereditary forms" (breast cancer → BRCA1/BARD1, prostate → HOXB13).
        k = "hereditary" if (bundle.get("1") or {}).get("is_cancer") else "monogenic"
        return f" whose {k} forms are linked to {genes}"
    if kind == "subtypes":
        # Prefer the genes curated across the subtypes (ranked by how many they
        # cover) over the few curated on the umbrella term itself (inherited
        # retinal dystrophy: PRPH2/BEST1/RHO, not CLN3/UBAP1L/IMPG1).
        return _subtype_gene_clause(bundle) or f" whose subtypes are caused by variants in {genes}"
    if kind == "clinvar":
        return f" with pathogenic variants reported in {genes} (ClinVar)"
    if len(syms) == 1:
        return f" caused by {syms[0]} ({cg[0][1]})"
    return f" caused by variants in {genes}"


def _subtype_gene_clause(bundle):
    """Umbrella term with no curated gene of its own: name the genes GenCC
    curates across its subtypes ("whose subtypes are caused by variants in
    PRPH2, BEST1, RHO, and 147 other genes")."""
    b1 = bundle.get("1") or {}
    b4 = bundle.get("4") or {}
    st = [g.get("symbol") for g in (b4.get("subtype_genes") or []) if g.get("symbol")]
    if (b1.get("child_count") or 0) < _UMBRELLA_SUBTYPES or not st:
        return ""
    total = max(b4.get("subtype_gene_total") or 0, len(st))
    syms = st[:3] + [""] * (total - min(3, len(st)))
    if b1.get("is_cancer"):          # germline predisposition, not "the cause"
        return f" whose hereditary forms are linked to {_gene_list(syms, None)}"
    return f" whose subtypes are caused by variants in {_gene_list(syms, None)}"


def _same_mode(a, b):
    """'Autosomal recessive' ~ 'Autosomal recessive'; GenCC's 'X-linked' matches
    Orphanet's 'X-linked recessive/dominant'."""
    a, b = (a or "").lower().replace(" inheritance", ""), (b or "").lower().replace(" inheritance", "")
    return a == b or a.startswith(b) or b.startswith(a)


def inheritance_modes(bundle):
    """(modes, source): Orphanet's disease-level inheritance, else — only for a
    disease its causal genes actually cause — their GenCC mode (type 2 diabetes
    must not inherit GCK-MODY's "autosomal dominant")."""
    from atlas.disease.summary import inheritance, orphanet_is_parent
    b1 = bundle.get("1") or {}
    kind, cg = causal_kind(bundle)
    gm = inheritance({}, bundle.get("4"), cg) if kind == "caused" else []
    orph = [] if orphanet_is_parent(b1) else inheritance(b1)
    if orph:
        # Orphanet lists every mode seen across a disease's forms; when the causal
        # gene's own curated modes are known, keep only the modes they share
        # (Fanconi anemia group C: Orphanet "AR or X-linked", FANCC is AR).
        if gm and len(orph) > 1:
            shared = [m for m in orph if any(_same_mode(m, g) for g in gm)]
            if shared:
                return shared, "Orphanet, consistent with GenCC"
        return orph, "Orphanet"
    return gm, ("GenCC" if gm else "")


def _with_source(text, src):
    """'… short stature.' + Orphanet → '… short stature (Orphanet).'"""
    if not text:
        return ""
    if text.endswith("."):
        return f"{text[:-1]} ({src})."
    return f"{text} ({src})."


def _features_clause(b1):
    from atlas.disease.summary import key_features
    terms, src = key_features(b1)
    if not terms:
        return ""
    lst = ", ".join(terms[:-1]) + (", and " if len(terms) > 2 else " and ") + terms[-1] \
        if len(terms) > 1 else terms[0]
    if src == "HPO":
        return f" Its most frequent clinical features (HPO) include {lst}."
    return f" Commonly reported symptoms (Wikidata) include {lst}."


def _evidence_sentence(b1, b2, b4, b5, b13):
    """'This page aggregates 43 associated genes (37 GWAS associations across 4
    studies) and 1,441 registered clinical trials.' — the Atlas evidence counts,
    kept out of the defining first sentence."""
    ev = _evidence_clause(b1, b2, b4, b5)            # " with N cohort genes (…)"
    ev = ev[len(" with "):] if ev.startswith(" with ") else ev
    ev = ev.replace("cohort gene", "associated gene")
    tr = _trials_clause(b1, b13, joined=bool(ev))     # " and N trials" / " with N trials"
    if not ev and not tr:
        return ""
    body = ev + (tr if ev else tr[len(" with "):])
    return f" This page aggregates {body}."


def declarative_sentence(bundle):
    """Compose the disease lead from a full bundle dict {section_id: bundle_dict}.

    Order: what the disease is (class, causal gene, inheritance), its curated
    definition, its key clinical features, then the Atlas evidence counts. The
    first sentence/definition is what search snippets and AI answers quote, so
    it carries the facts readers ask for ("is X hereditary", "what is X")
    rather than cohort counts.

    Sources used:
      bundle["1"]  — name, mondo_id, is_cancer, definitions, inheritance, HPO, Wikidata
      bundle["2"]  — assoc_total, study_total
      bundle["4"]  — GenCC causal genes + mode of inheritance, somatic drivers
      bundle["5"]  — gene_count
      bundle["13"] — trial_count, trial_drugs, CIViC
      bundle["14"] — top_pathways (dominant pathway, when clearly so)

    Returns Markdown (no trailing newline)."""
    from atlas.disease import summary as S
    b1  = bundle.get("1")  or {}
    b2  = bundle.get("2")  or {}
    b4  = bundle.get("4")  or {}
    b5  = bundle.get("5")  or {}
    b13 = bundle.get("13") or {}
    b14 = bundle.get("14") or {}

    name = b1.get("canonical_name") or b1.get("name") or "?"
    mondo_id = b1.get("mondo_id")
    klass = _disease_class(b1)

    head = f"**{name.capitalize() if name == name.lower() else name}**"
    if mondo_id:
        head += f" ({mondo_id})"

    n_sub = b1.get("child_count") or 0       # umbrella / grouping term?
    umbrella = f" (an umbrella term covering {n_sub} Mondo subtypes)" if n_sub >= 5 else ""
    causal = _causal_clause(bundle)          # " caused by CFTR (GenCC Definitive)"
    inh = S.inheritance_phrase(inheritance_modes(bundle)[0])
    defn, src = S.definition(b1)
    # A parent-group Orphanet definition ("A rare group of …") never becomes the
    # predicate of "<this subtype> is …"; it's quoted with its attribution.
    pred = "" if src.startswith("Orphanet, for") else S.definition_predicate(defn, name)
    form_of = S.templated_class(defn)        # "Any <class> in which … the X gene"
    if form_of:
        defn = pred = ""                     # restates the causal gene: use the class

    if causal or inh:
        what = f"a form of {form_of}" if form_of else f"a {klass}"
        sentence = f"{head} is {what}{umbrella}{causal}"
        if inh:
            sentence += ("," if causal else "") + f" {inh}"
        sentence += "."
        if defn:
            sentence += " " + _with_source(f"It is {pred}" if pred else defn, src)
    elif pred:
        sentence = f"{head} is " + _with_source(pred, src)
        if umbrella:
            sentence += f" It is an umbrella term covering {n_sub} Mondo subtypes."
    else:
        sentence = (f"{head} is a form of {form_of}{umbrella}." if form_of
                    else f"{head} is a {klass}{umbrella}.")
        if defn:
            sentence += " " + _with_source(defn, src)
    sentence += _features_clause(b1)
    sentence += _evidence_sentence(b1, b2, b4, b5, b13)
    # (No "dominant Reactome pathway" sentence: on the summary it read as a
    # finding, yet it was cohort-enrichment noise on many pages — Huntington and
    # sickle cell → "Olfactory Signaling" from clustered OR genes in the GWAS tail.
    # The enrichment table stays in the body with its statistics.)
    sentence += _civic_subtype_clause(b13)
    sentence += _drugs_clause(b13, bundle.get("_indicated_drugs"))
    sentence += _parent_clause(b1, b5)
    return sentence


def _parent_clause(b1, b5):
    """When this page's gene cohort is empty but it sits under a broader Mondo
    term, point the reader to the parent — granular subtype terms carry little
    direct evidence while the parent holds the cohort (IDH-wildtype glioblastoma
    → glioblastoma). Plain text; the linked navigation is in the Disease family
    section. '' for content-rich diseases."""
    parent = (b1 or {}).get("parent") or None
    if not parent or ((b5 or {}).get("gene_count") or 0) > 0:
        return ""
    nm = parent.get("name")
    if not nm:
        return ""
    return (f" A subtype of {nm} — broader associated-gene and molecular evidence "
            "is on the parent page (see Disease family below).")


def _civic_subtype_clause(b13):
    """Answer-first precision-medicine verdict from the top CIViC predictive
    association on the disease's subtype map (§13). '' for non-cancer diseases
    or those with no curated evidence."""
    from atlas.civic import predictive_verdict
    rows = (b13 or {}).get("civic_evidence") or []
    verdict = predictive_verdict(rows)
    if not verdict:
        return ""
    total = (b13 or {}).get("civic_association_total") or 0
    more = (f"; {total - 1} further subtype–drug associations are mapped below"
            if total > 1 else "")
    return f" Molecularly, {verdict}{more}."
