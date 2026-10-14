"""Intent-specific FAQ questions — the questions people (and AI search) actually
ask about a gene / disease / drug, answered from the page's own data.

Bing queries that reached Atlas pages: "how many exons in the mfsd8 gene",
"is bronchiectasis … autosomal recessive or dominant?", "sotagliflozin molecular
targets", "what is the pde1c gene? are there any drugs … that target this gene?".
Each question is emitted only when the bundle answers it; every answer names
its source. Pure functions of the bundle; [] when nothing applies.
"""
import re


def _join(xs):
    xs = [x for x in xs if x]
    if len(xs) <= 1:
        return "".join(xs)
    if len(xs) == 2:
        return " and ".join(xs)
    return ", ".join(xs[:-1]) + ", and " + xs[-1]


def _more(n_total, n_shown):
    return f" (+{n_total - n_shown} more)" if n_total > n_shown else ""


# ── disease ─────────────────────────────────────────────────────────────────

def _disease(bundle, title):
    from atlas.page.disease_declarative import (causal_kind, inheritance_modes,
                                                inheritance_scope, inheritance_predicate,
                                                hereditary_named)
    from atlas.disease import summary as S
    from atlas.render_common import display_name
    b1 = bundle.get("1") or {}
    qas = []

    kind, cg = causal_kind(bundle)
    modes, src = inheritance_modes(bundle)
    src = {"Orphanet, consistent with GenCC": "Orphanet and GenCC"}.get(src, src)
    phrase = S.inheritance_phrase(modes)
    genes = _join([f"{s} ({l})" for s, l in cg[:5]]) + _more(len(cg), 5)
    cause = {"caused": f"{title} is caused by variants in {genes}.",
             "forms": (f"{title} is a cancer; its hereditary forms are linked to germline "
                       f"variants in {genes}." if b1.get("is_cancer") and not hereditary_named(b1) else
                       f"{title} is a complex (multifactorial) disease; its rare monogenic "
                       f"forms are linked to {genes}."),
             "subtypes": f"{title} is an umbrella term; its subtypes are caused by "
                         f"variants in {genes}."}.get(kind)
    st = [g for g in ((bundle.get("4") or {}).get("subtype_genes") or []) if g.get("symbol")]
    if st and kind in ("subtypes", None) and (b1.get("child_count") or 0) >= 5:
        tot = max((bundle.get("4") or {}).get("subtype_gene_total") or 0, len(st))
        top = _join([f"{g['symbol']} ({g['subtypes']} subtype{'s' if g['subtypes'] != 1 else ''})"
                     for g in st[:5]])
        cause = (f"{title} is an umbrella term covering {b1.get('child_count')} Mondo subtypes. "
                 f"GenCC curates {tot} gene{'s' if tot != 1 else ''} as Definitive or Strong causes of its subtypes; "
                 f"the most frequent are {top}.")
    scope = inheritance_scope(bundle)
    if cause:
        if phrase:
            if scope != "disease":
                cause += (f" Its {'subtypes' if scope == 'subtypes' else 'familial forms'} "
                          f"{inheritance_predicate(phrase)} ({src}).")
            else:
                cause += (f" It has {phrase[5:]} ({src})." if phrase.startswith("with ")
                          else f" It is {phrase} ({src}).")
        qas.append((f"What causes {title}?", cause))
    if phrase:
        if scope == "forms":
            # A complex disease or sporadic cancer: Orphanet's modes describe its
            # familial forms (ALS, moyamoya, cataract answered a flat "Yes" in v1.12).
            a = (f"Mostly not: {title} is usually sporadic or multifactorial. Its familial "
                 f"forms {inheritance_predicate(phrase)} ({src}).")
        elif scope == "subtypes":
            a = f"{title} is an umbrella term; its subtypes {inheritance_predicate(phrase)} ({src})."
        elif phrase.startswith("with "):
            a = f"According to {src}, {title} has {phrase[5:]}."
        else:
            a = f"Yes. According to {src}, {title} is {phrase}."
        qas.append((f"Is {title} inherited?", a))

    terms, fsrc = S.key_features(b1, k=8)
    if terms:
        if fsrc == "HPO":
            a = (f"The most frequent clinical features of {title} (HPO, ranked by "
                 f"Orphanet frequency) are {_join(terms)}.")
        else:
            a = f"Commonly reported symptoms of {title} (Wikidata) include {_join(terms)}."
        qas.append((f"What are the symptoms of {title}?", a))

    rows = [r for r in (bundle.get("_indicated_drugs") or []) if r.get("name")]
    appr = [r for r in rows if r.get("approved")]
    late = [r for r in rows if not r.get("approved") and r.get("late")]
    # Only phase-4 rows answer "approved for" (v1.13 QA: IPF listed ambrisentan,
    # whose IPF trial stopped for harm); phase-3 uses of approved drugs get their
    # own, hedged question.
    if appr:
        names = [display_name(r["name"]) for r in appr[:8]]
        qas.append((f"What drugs are approved for {title}?",
                    f"Drugs with an approved indication for {title} in ChEMBL (phase 4): "
                    f"{_join(names)}{_more(len(appr), len(names))}."))
    if late:
        names = [display_name(r["name"]) for r in late[:8]]
        qas.append((f"Which approved drugs have been tested in phase 3 for {title}?",
                    f"Drugs approved for other uses that reached phase 3 for {title} in ChEMBL "
                    f"(phase-3 evidence, not an approval for this disease): "
                    f"{_join(names)}{_more(len(late), len(names))}."))
    b13 = bundle.get("13") or {}
    n, ids = recruiting(b13.get("status_counts"), b13.get("top_trials"))
    if n:
        qas.append((f"Are there recruiting clinical trials for {title}?",
                    f"Yes — at least {n:,} registered trial{'s are' if n != 1 else ' is'} recruiting or about "
                    f"to recruit (ClinicalTrials.gov)" + (f", e.g. {', '.join(ids[:3])}" if ids else "") + "."))
    sub = [r for r in (bundle.get("_subtype_approved_drugs") or []) if r.get("name")]
    if sub and not appr:
        names = [f"{display_name(r['name'])} ({r['subtypes'][0]})" for r in sub[:8]]
        extra = (f"No drug is recorded as approved for {title} as a whole; drugs approved "
                 f"for its subtypes (ChEMBL) include {_join(names)}{_more(len(sub), len(names))}.")
        q = f"What drugs are approved for {title}?"
        prev = [i for i, (qq, _) in enumerate(qas) if qq == q]
        if prev:
            qas[prev[0]] = (q, qas[prev[0]][1] + " " + extra)
        else:
            qas.append((q, extra))
    return qas


# ── drug ────────────────────────────────────────────────────────────────────

def _drug(bundle, title):
    from atlas.page.drug_declarative import _curated_target_genes, _indication_names
    b1 = bundle.get("1") or {}
    b2 = bundle.get("2") or {}
    b4 = bundle.get("4") or {}
    qas = []

    moa = []
    for m in b2.get("mechanisms") or []:
        t = (m.get("mechanism_of_action") or "").strip().rstrip(".")
        if t.lower() in ("unknown", "not applicable", "n/a", "none"):
            continue
        if t and t.lower() not in {x.lower() for x in moa}:
            moa.append(t)
    genes = (list(dict.fromkeys(g.get("gene_symbol") for g in b2.get("mechanism_genes") or []
                                if g.get("gene_symbol")))
             if moa else _curated_target_genes(b2))       # same rule as the lead
    orgs = b2.get("mechanism_organisms") or []
    if moa or genes:
        parts = []
        if moa:
            parts.append(f"{title} acts as: {'; '.join(moa[:4])}{_more(len(moa), 4)} "
                         "(ChEMBL curated mechanism of action).")
        if genes:
            parts.append(f"Its target genes are {_join(genes[:6])}{_more(len(genes), 6)}.")
        if orgs:
            parts.append(f"Target organism: {_join(orgs[:3])}.")
        qas.append((f"What is the mechanism of action of {title}?", " ".join(parts)))

    # Approval status is the MOLECULE's (FDA flag / ChEMBL phase 4 / phase-4
    # trial) — the per-indication phases can't say "not approved" for an approved
    # drug (v1.12: 77 FDA-approved drugs answered "no approved indication").
    from atlas.indication import indication_tiers
    from atlas.page.drug_at_a_glance import _status
    appr, late, trial = indication_tiers(b4)
    status = _status(b1, bundle.get("5"), bundle.get("14"))
    from atlas.indication import molecule_approved
    # ChEMBL phase 4 is an approval somewhere, even when DrugCentral says "not
    # FDA" (cefamandole: "No" beside "approved for 1 indication").
    approved_mol = status.startswith(("FDA-approved", "Approved")) or \
        molecule_approved(b1, bundle.get("5"))
    if appr:
        names = _indication_names(appr, k=8)
        a = (f"Approved indications of {title} (ChEMBL phase 4): {_join(names)}"
             f"{_more(len(appr), len(names))}.")
        if late:
            a += (f" ChEMBL also records phase-3 evidence for {len(late)} other condition"
                  + ("s" if len(late) != 1 else "") + ".")
        qas.append((f"What is {title} approved for?", a))
    elif late and approved_mol:
        names = _indication_names(late, k=6)
        qas.append((f"What is {title} used for?",
                    f"{title} is an approved drug ({re.sub(r' [(](.*)[)]$', r'; \1', status)}). ChEMBL records no phase-4 "
                    f"indication for it; its most advanced recorded uses (phase 3) include "
                    f"{_join(names)}{_more(len(late), len(names))}. Check the drug label for "
                    "the approved indications."))
    if status:
        if approved_mol and status.startswith("Not FDA-approved"):
            a = (f"Yes — ChEMBL records {title} as an approved drug (max phase 4), "
                 "though DrugCentral lists no FDA approval.")
        elif approved_mol:
            a = f"Yes. {title} status: {status}."
        else:
            a = f"No. {title} status: {status}."
            if trial or late:
                rows = late + trial
                names = _indication_names(rows, k=3)
                a += (f" It has been investigated in clinical trials for {len(rows)} condition"
                      + ("s" if len(rows) != 1 else "")
                      + (f", including {_join(names)}" if names else "") + ".")
        qas.append((f"Is {title} an approved drug?", a))

    labels, guide, pgenes = pgx_summary(bundle.get("9"))
    if labels or guide:
        parts = []
        if labels:
            parts.append("Drug labels with pharmacogenomic content (PharmGKB): " + "; ".join(
                f"{src} — {lvl}" + (f" ({', '.join(gs[:3])})" if gs else "") for src, lvl, gs in labels[:4]) + ".")
        if guide:
            parts.append(f"Dosing guidelines: {', '.join(guide)}"
                         + (f" ({', '.join(pgenes[:4])})" if pgenes else "") + ".")
        top = labels[0][1].lower() if labels else ""
        lead = ("Yes — genetic testing is required by at least one label. " if top == "testing required"
                else "Genetic testing is recommended by at least one label. " if top == "testing recommended"
                else "Pharmacogenomic guidance exists. ")
        qas.append((f"Is genetic testing recommended before taking {title}?", lead + " ".join(parts)))

    ae = adverse_events(bundle.get("13"))
    tot = (bundle.get("13") or {}).get("total_reports") or 0
    if ae:
        qas.append((f"What adverse events are reported for {title}?",
                    f"Most frequently reported reactions in FDA FAERS spontaneous reports "
                    f"({tot:,} reports; reports do not establish causation): "
                    f"{_join([r['reaction'] for r in ae])}."))

    n, ids = recruiting((bundle.get("5") or {}).get("status_counts"), (bundle.get("5") or {}).get("top_trials"))
    if n:
        qas.append((f"Are there recruiting clinical trials for {title}?",
                    f"Yes — at least {n:,} registered trial{'s are' if n != 1 else ' is'} recruiting or about "
                    f"to recruit (ClinicalTrials.gov)" + (f", e.g. {', '.join(ids[:3])}" if ids else "") + "."))
    return qas


# ── gene ────────────────────────────────────────────────────────────────────

def _phase(p):
    try:
        return float(p)
    except (TypeError, ValueError):
        return 0.0


def _gene(bundle, title):
    from atlas.render_common import display_name
    b2 = bundle.get("2") or {}
    b10 = bundle.get("10") or {}
    qas = []

    n_ex = b2.get("canonical_exon_count")
    tx = b2.get("canonical_transcript")
    if n_ex:
        a = (f"The canonical Ensembl transcript of {title}"
             + (f" ({tx})" if tx else "") + f" has {n_ex} exon" + ("s" if n_ex != 1 else "") + ".")
        nt = b2.get("ensembl_transcript_count")
        if nt:
            a += f" {title} has {nt} annotated Ensembl transcript" + ("s" if nt != 1 else "") + "."
        qas.append((f"How many exons does {title} have?", a))

    # Where is it expressed? (HPA RNA tissue specificity + top tissues + location)
    hpa = (bundle.get("13") or {}).get("hpa") or {}
    tops = [t.split("|")[0] for t in hpa.get("top_tissues") or [] if t][:4]
    if hpa.get("rna_tissue_specificity") or tops:
        a = (f"Human Protein Atlas RNA tissue specificity: {hpa.get('rna_tissue_specificity') or 'not classified'}"
             + (f"; highest in {_join([t.lower() for t in tops])}" if tops else "") + ".")
        loc = hpa.get("subcellular_main") or []
        if loc:
            a += f" Main subcellular location: {_join(loc[:3])}."
        qas.append((f"Where is {title} expressed?", a))

    # Which diseases? curated validity first, then pathogenic-ClinVar-supported.
    b12 = bundle.get("12") or {}
    from atlas.render_common import gencc_rank
    cur = []
    for g in (b12.get("clingen_validity") or []) + (b12.get("gencc") or []):
        d = (g.get("disease") or "").strip()
        if d and gencc_rank(g.get("classification")) >= gencc_rank("Strong") and d.lower() not in {x.lower() for x in cur}:
            cur.append(d)
    # Other = ≥2 pathogenic ClinVar records or GenCC Moderate+ (a Limited /
    # Refuted GenCC call alone is not "linked").
    other = [m.get("name") for m in b12.get("mondo") or [] if m.get("name")
             and (m.get("clinvar_records") or gencc_rank(m.get("gencc")) >= gencc_rank("Moderate"))
             and m["name"].lower() not in {x.lower() for x in cur}]
    if cur or other:
        a = ""
        if cur:
            a = (f"Diseases with Definitive or Strong curated gene–disease validity (GenCC / ClinGen): "
                 f"{_join(cur[:8])}{_more(len(cur), 8)}.")
        if other:
            a += ((" " if a else "") + f"{'Other diseases' if cur else 'Diseases'} linked through curated or "
                  f"pathogenic ClinVar evidence: {_join(other[:6])}{_more(len(other), 6)}.")
        qas.append((f"What diseases is {title} associated with?", a))

    mp = (bundle.get("5") or {}).get("mouse_phenotypes") or []
    if mp:
        n = (bundle.get("5") or {}).get("mouse_phenotype_total") or len(mp)
        qas.append((f"Is there a mouse model for {title}?",
                    f"Yes — MGI records {n:,} phenotype{'s' if n != 1 else ''} for mutant alleles of the "
                    f"mouse orthologue, including {_join([p['statement'] for p in mp[:4] if p.get('statement')])}."))

    moa = [d for d in (b10.get("moa_drugs") or []) if d.get("name")]
    if moa and any(_phase(d.get("phase")) > 0 for d in moa):
        appr = [d for d in moa if _phase(d.get("phase")) >= 4]
        clin = [d for d in moa if 0 < _phase(d.get("phase")) < 4]
        parts = []
        if appr:
            names = [display_name(d["name"]) for d in appr[:8]]
            parts.append(f"Approved drugs whose curated ChEMBL mechanism of action targets "
                         f"{title}: {_join(names)}{_more(len(appr), len(names))}.")
        if clin:
            names = [display_name(d["name"]) for d in clin[:5]]
            parts.append(("Clinical-stage" if appr else f"No approved drug targets {title} directly; clinical-stage")
                         + f" drugs with a curated mechanism on {title} include "
                           f"{_join(names)}{_more(len(clin), len(names))}.")
        if parts:
            qas.append((f"What drugs target {title}?", " ".join(parts)))
    cx = [d for d in (b10.get("complex_moa_drugs") or []) if d.get("name") and _phase(d.get("phase")) >= 4]
    if cx:
        names = [display_name(d["name"]) for d in cx[:8]]
        a = (f"Approved drugs acting on a protein complex or family that includes {title} "
             f"({'; '.join(list(dict.fromkeys(d.get('target') for d in cx if d.get('target')))[:2])}): "
             f"{_join(names)}{_more(len(cx), len(names))}.")
        if qas and qas[-1][0] == f"What drugs target {title}?":
            qas[-1] = (qas[-1][0], qas[-1][1] + " " + a)
        else:
            qas.append((f"What drugs target {title}?", a))
    return qas


# ── shared: PGx, adverse events, recruiting trials ──────────────────────────

_PGX_RANK = {"testing required": 4, "testing recommended": 3, "actionable pgx": 2,
             "informative pgx": 1}


def pgx_summary(b9):
    """(labels ≥ Actionable as [(agency, level, genes)], guideline sources, genes)
    from PharmGKB drug-label + guideline annotations; ([], [], []) if none."""
    b9 = b9 or {}
    labels = []
    for l in b9.get("drug_labels") or []:
        lvl = (l.get("level") or "").strip()
        if _PGX_RANK.get(lvl.lower(), 0) >= 2:
            labels.append((l.get("source") or "", lvl, [g for g in l.get("genes") or [] if g]))
    labels.sort(key=lambda x: (-_PGX_RANK.get(x[1].lower(), 0), x[0] != "FDA", x[0]))
    guide = sorted({(g.get("source") or "").strip() for g in b9.get("guidelines") or []} - {""})
    genes = []
    for _, _, gs in labels:
        genes += [g for g in gs if g not in genes]
    for g in b9.get("guidelines") or []:
        genes += [x.strip() for x in str(g.get("genes") or "").split(",") if x.strip() and x.strip() not in genes]
    return labels, guide, genes


# FAERS terms that are not adverse reactions (36% of drugs led with one).
_NON_REACTION = re.compile(
    r"drug ineffective|off label|product |no adverse event|death|inappropriate|incorrect|"
    r"wrong |dose omission|intentional|unapproved|medication error|exposure|"
    r"condition aggravated|disease progression|therapeutic (response|product)|"
    r"treatment (failure|noncompliance)|drug interaction|drug intolerance|"
    r"not effective|lack of efficacy|underdose|overdose|poor quality|"
    r"device|malaise$|general physical health deterioration|"
    r"^(pain|fall|illness)$", re.I)


def adverse_events(b13, k=8):
    rows = [r for r in (b13 or {}).get("most_reported") or []
            if r.get("reaction") and not _NON_REACTION.search(r["reaction"])]
    return rows[:k]


def recruiting(status_counts, trials, strict=False):
    """(n, example ids): trials recruiting or about to (strict: recruiting only,
    for a bare "N recruiting" label)."""
    want = ("RECRUITING",) if strict else ("RECRUITING", "NOT_YET_RECRUITING")
    n = sum(v for k, v in (status_counts or {}).items() if str(k).upper() in want)
    ids = [t.get("id") for t in (trials or []) if str(t.get("status") or "").upper() == "RECRUITING"
           and t.get("id")]
    return n, ids


_BUILDERS = {"disease": _disease, "drug": _drug, "gene": _gene}


def intent_qas(entity_type, bundle, title):
    """[(question, answer)] for this page; [] when nothing applies or on error
    (the FAQ is additive — it must never break a page build)."""
    fn = _BUILDERS.get(entity_type)
    title = (title or "").strip()
    if not fn or not title:
        return []
    try:
        return [(q, re.sub(r"\s+", " ", a).strip()) for q, a in fn(bundle or {}, title) if a]
    except Exception:
        return []
