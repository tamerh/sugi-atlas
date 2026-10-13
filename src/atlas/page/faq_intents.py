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
    from atlas.page.disease_declarative import causal_kind, inheritance_modes
    from atlas.disease import summary as S
    from atlas.render_common import display_name
    b1 = bundle.get("1") or {}
    qas = []

    kind, cg = causal_kind(bundle)
    modes, src = inheritance_modes(bundle)
    phrase = S.inheritance_phrase(modes)
    genes = _join([f"{s} ({l})" for s, l in cg[:5]]) + _more(len(cg), 5)
    cause = {"caused": f"{title} is caused by variants in {genes}.",
             "forms": (f"{title} is a cancer; its hereditary forms are linked to germline "
                       f"variants in {genes}." if b1.get("is_cancer") else
                       f"{title} is a complex (multifactorial) disease; its rare monogenic "
                       f"forms are linked to {genes}."),
             "subtypes": f"{title} is an umbrella term; its subtypes are caused by "
                         f"variants in {genes}."}.get(kind)
    if cause:
        if phrase:
            cause += f" It is {phrase} ({src})."
        qas.append((f"What causes {title}?", cause))
    if phrase:
        if phrase.startswith("with "):
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

    appr = [r for r in (bundle.get("_indicated_drugs") or []) if r.get("approved") and r.get("name")]
    if appr:
        names = [display_name(r["name"]) for r in appr[:8]]
        qas.append((f"What drugs are approved for {title}?",
                    f"Drugs with an approved indication for {title} in ChEMBL: "
                    f"{_join(names)}{_more(len(appr), len(names))}."))
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
        if t and t.lower() not in {x.lower() for x in moa}:
            moa.append(t)
    genes = _curated_target_genes(b2)
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

    inds = b4.get("indications") or []
    appr = [i for i in inds if i.get("approved")]
    trial = [i for i in inds if not i.get("approved")]
    if appr:
        names = _indication_names(appr, k=8)
        qas.append((f"What is {title} approved for?",
                    f"Approved indications of {title} (ChEMBL): {_join(names)}"
                    f"{_more(len(appr), len(names))}."
                    + (f" It is also in clinical trials for {len(trial)} other condition"
                       + ("s" if len(trial) != 1 else "") + "." if trial else "")))
    else:
        phase = b1.get("max_phase")
        try:
            phase = int(float(phase))
        except (TypeError, ValueError):
            phase = None
        if phase is not None and phase < 4:
            a = f"No approved indication is recorded for {title} in ChEMBL; its highest development phase is {phase}."
            if trial:
                names = _indication_names(trial, k=3)
                a += (f" It has been investigated in clinical trials for {len(trial)} condition"
                      + ("s" if len(trial) != 1 else "")
                      + (f", including {_join(names)}" if names else "") + ".")
            qas.append((f"Is {title} an approved drug?", a))
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

    moa = [d for d in (b10.get("moa_drugs") or []) if d.get("name")]
    if moa:
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
    return qas


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
