"""Plain-language summary facts for a disease page: what the disease IS
(curated definition), how it is inherited, and its key clinical features.

Shared by the lead sentence, the clinical-description block and the JSON-LD so
the three can't disagree. Pure functions of the §1/§4 bundles; '' / [] when a
fact is absent.
"""
import re

# Definition sources, most disease-specific first. Orphanet is written for the
# rare-disease reader; Mondo covers every disease; MeSH and DOID fill the
# common-disease gap.
_DEF_SOURCES = (("orphanet_definition", "Orphanet"), ("mondo_definition", "Mondo"),
                ("mesh_scope_note", "MeSH"), ("doid_definition", "Disease Ontology"))
_DEF_MAX = 260

# MeSH scope notes upper-case cross-referenced descriptors ("a degenerative disease
# of the BRAIN characterized by … DEMENTIA"). Lower-case those runs; short
# stand-alone acronyms (DNA, HIV, AIDS) are left alone.
_CAPS_RUN = re.compile(r"\b[A-Z][A-Z-]{1,}(?:[ ,][A-Z][A-Z-]{1,})+\b|\b[A-Z][A-Z-]{4,}\b")
_ACRONYMS = {"DNA", "RNA", "HIV", "AIDS", "MHC", "HLA", "CNS", "ATP", "NADH", "NADPH",
             "IGE", "IGG", "IGA", "IGM", "MRNA", "TRNA", "RRNA", "COVID-19", "SARS"}


def _uncap_mesh(text):
    def fix(m):
        words = m.group(0)
        if words.upper() in _ACRONYMS or all(w in _ACRONYMS for w in re.split(r"[ ,]+", words)):
            return words
        return words.lower()
    return _CAPS_RUN.sub(fix, text)


# Disease Ontology logical-definition relations leak into its text
# ("… that has_material_basis_in mutation of …", 740 pages).
_DO_RELATIONS = {"has_material_basis_in": "is caused by", "located_in": "is located in",
                 "transmitted_by": "is transmitted by", "realized_by": "is realized by",
                 "has_material_basis_in_mutation_in": "is caused by mutation in"}
_REL = re.compile(r"\b[a-z]+(?:_[a-z]+)+\b")


def clean_text(text):
    """Source-definition cleanup: HTML entities / tags (Orphanet '&lt;i&gt;' —
    143 pages), DO relation tokens, whitespace and space-before-punctuation."""
    import html
    t = text or ""
    for _ in range(2):                         # double-escaped '&amp;lt;i&amp;gt;'
        t = html.unescape(t)
    t = re.sub(r"</?[a-zA-Z][^>]{0,40}>", "", t)
    t = _REL.sub(lambda m: _DO_RELATIONS.get(m.group(0), m.group(0).replace("_", " ")), t)
    t = " ".join(t.split())
    return re.sub(r"\s+([,;:.])", r"\1", t)


def _first_sentences(text, limit=_DEF_MAX):
    """Leading sentence(s) of a definition, within `limit` characters. A sentence
    ends at '. ' followed by a capital (abbreviations like "e.g. the" don't). An
    over-long first sentence is cut at a clause boundary (', which', '; ', ', ')
    rather than mid-phrase."""
    text = clean_text(text)
    text = re.sub(r"\s+([,;:.])", r"\1", text)   # source text: "fatal , epileptic"
    parts = re.split(r"(?<=[.!?])\s+(?=[A-Z(])", text)
    out = ""
    for p in parts:
        if out and len(out) + 1 + len(p) > limit:
            break
        out = f"{out} {p}".strip()
        if len(out) >= 120:
            break
    if len(out) > limit:                       # one over-long sentence
        head = out[:limit]
        cut = max(head.rfind(", which"), head.rfind("; "), head.rfind(", "), head.rfind(" (")) \
            if len(head) > 120 else -1
        out = (head[:cut] if cut >= 120 else head.rsplit(" ", 1)[0]).rstrip(",;: (") + "…"
    if out and out[-1] not in ".!?…":
        out += "."
    return out


_MARKER = re.compile(r"\b(?:\d+[a-z]?|(?:type|group|form|class)\s+[a-z0-9]{1,4})\b", re.I)


def _markers(name):
    """Subtype markers, normalised so "type 11" == "11"."""
    return {re.sub(r"^(?:type|form|class)\s+(?=\d)", "", re.sub(r"\s+", " ", m.group(0).lower()))
            for m in _MARKER.finditer(name or "")}


def orphanet_is_parent(b1):
    """True when the page's Orphanet record is the broader GROUP, not this
    subtype: the Mondo name carries a subtype marker (a number, "type 4",
    "group C") the Orphanet name lacks — long QT syndrome 2 → Orphanet "Congenital
    long QT syndrome"; Stickler syndrome type 4 → "Stickler syndrome". Its
    definition ("A rare group of …") and union of inheritance modes ("AD or AR")
    then describe the group, not this disease (~650 pages)."""
    b1 = b1 or {}
    on, cn = b1.get("orphanet_name") or "", b1.get("canonical_name") or b1.get("name") or ""
    if not on:
        return False
    # Authoritative when known: does the Orphanet record map back to THIS Mondo
    # id? (The name diff misfired on renamed nomenclature — LGMD 2A vs
    # "Calpain-3-related LGMD R1", HSP 11 vs "… type 11": ~150 pages.)
    back = b1.get("orphanet_mondo")
    if back and b1.get("mondo_id") and b1["mondo_id"] in back:
        return False
    # Otherwise a subtype marker this name has and Orphanet's lacks (a mapping to
    # another Mondo id alone isn't enough: Birt-Hogg-Dubé maps to a sibling id).
    return bool(cn and (_markers(cn) - _markers(on)))


def definition(b1):
    """(text, source) — the best curated definition, trimmed to its opening
    sentence(s); ('', '') when none. A parent-group Orphanet definition is used
    only as a last resort, attributed to the group."""
    b1 = b1 or {}
    parent = orphanet_is_parent(b1)
    cands = []
    for key, src in _DEF_SOURCES:
        if parent and key == "orphanet_definition":
            continue
        t = (b1.get(key) or "").strip()
        if len(t) < 20:
            continue
        if src == "MeSH":
            t = _uncap_mesh(t)
        cands.append((_first_sentences(t), src))
    # The first source whose opening sentence fits whole; else the first one
    # (clause-cut). v1.12: 3,177 leads ended mid-sentence with "…".
    for txt, src in cands:
        if not txt.endswith("…"):
            return txt, src
    if cands:
        return cands[0]
    t = (b1.get("orphanet_definition") or "").strip()
    if parent and len(t) >= 20:
        return _first_sentences(t), f"Orphanet, for the broader {b1.get('orphanet_name')}"
    return "", ""


def definition_predicate(text, name):
    """The definition as the predicate of "<name> is …", or '' when it doesn't
    read that way. "A rare, genetic pulmonary disorder …" → "a rare, genetic
    pulmonary disorder …"; "Asthma is a bronchial disease …" → "a bronchial
    disease …"."""
    t = (text or "").strip()
    nm = (name or "").strip()
    if nm and t.lower().startswith(nm.lower() + " is "):
        return _article(t[len(nm) + 4:])
    if re.match(r"(A|An|The)\s", t):
        return _article(t[0].lower() + t[1:])
    return ""


def _article(p):
    """Fix a/an agreement the source got clearly wrong ("An disease or disorder
    …"). Only unambiguous cases change: letter-only rules broke correct text
    ("a urea cycle disorder" → "an urea", "an X-linked" → "a X-linked")."""
    m = re.match(r"(a|an) (\w[\w-]*)", p)
    if not m:
        return p
    art, w = m.group(1), m.group(2)
    lo = w.lower()
    if art == "a" and lo[0] in "aeio" and not lo.startswith(("one", "once")):
        want = "an"
    elif art == "an" and w == lo and lo[0] not in "aeiouh":
        want = "a"           # lowercase consonant word (acronyms keep "an": an RNA, an X-linked)
    else:
        return p
    return want + p[len(art):]


_INHERITANCE_SKIP = {"not applicable", "unknown", "not yet documented", "no data available", ""}


def inheritance(b1, b4=None, causal=None):
    """Modes of inheritance: Orphanet's disease-level annotation, else the GenCC
    mode recorded for the given causal genes' Definitive/Strong records. Callers
    gate the GenCC fallback on the disease actually being caused by those genes
    (disease_declarative.inheritance_modes)."""
    out = []
    if orphanet_is_parent(b1):
        return [] if not causal else inheritance({}, b4, causal)
    for m in (b1 or {}).get("orphanet_inheritance") or []:
        if m and m.strip().lower() not in _INHERITANCE_SKIP and m not in out:
            out.append(m.strip())
    if out:
        return out
    keep = {s for s, _ in (causal or [])}
    for g in (b4 or {}).get("disease_gencc") or []:
        if g.get("symbol") not in keep:
            continue
        for m in g.get("moi") or []:
            if m.strip().lower() not in _INHERITANCE_SKIP and m not in out:
                out.append(m.strip())
    return out


def inheritance_phrase(modes):
    """'inherited in an autosomal recessive pattern' / '… autosomal dominant or
    autosomal recessive pattern' / 'with mitochondrial inheritance'."""
    if not modes:
        return ""
    pattern, other = [], []
    for m in modes:
        lo = m.lower()
        if re.search(r"autosomal|x-linked|y-linked|semidominant|semi-dominant", lo):
            pattern.append(re.sub(r"\b([xy])-linked", lambda m: m.group(1).upper() + "-linked",
                                  lo.replace(" inheritance", "")))
        else:
            other.append(lo.replace(" inheritance", ""))
    if pattern:
        p3 = pattern[:3]
        joined = " or ".join(p3) if len(p3) < 3 else f"{p3[0]}, {p3[1]}, or {p3[2]}"
        art = "an" if joined[:1].lower() in "aeiox" else "a"     # "an X-linked"
        return f"inherited in {art} {joined} pattern"
    return "with " + " or ".join(other[:2]) + " inheritance"


def key_features(b1, k=3):
    """Top clinical features: Orphanet-frequency HPO terms (only when frequencies
    exist — an unranked OMIM list would be arbitrary), else Wikidata symptoms for
    a disease HPO barely covers. Returns (terms, source)."""
    phs = (b1 or {}).get("phenotypes") or []
    ranked = [p for p in phs if p.get("frequency_value")]
    if ranked:
        ranked.sort(key=lambda p: -(p.get("frequency_value") or 0))
        return [p.get("hpo_term").lower() if _plain_lower(p.get("hpo_term")) else p.get("hpo_term")
                for p in ranked[:k] if p.get("hpo_term")], "HPO"
    if len(phs) < 10:
        # Short labels only: some P780 values are whole phrases ("difficulty
        # completing familiar tasks at home, at work or at leisure").
        ws = [s.get("name") for s in ((b1 or {}).get("wikidata_symptoms") or [])
              if s.get("name") and len(s["name"]) <= 40]
        if ws:
            return ws[:k], "Wikidata"
    return [], ""


def _plain_lower(term):
    """Lower-case an HPO label for running text unless it holds an acronym
    ("Abnormal EEG")."""
    return bool(term) and not re.search(r"\b[A-Z]{2,}\b", term[1:])


# Mondo's templated gene-defined subtype definition: "Any autosomal dominant
# nonsyndromic deafness in which the cause of the disease is a mutation in the
# DIAPH1 gene." It only restates the causal gene; the useful part is the class.
_MONDO_TEMPLATE = re.compile(
    r"^Any (.+?) in which the cause of the disease is (?:a mutation|mutations|a variation) in the \S+ gene",
    re.I)


def templated_class(text):
    """'autosomal dominant nonsyndromic deafness' from Mondo's templated
    definition, else ''."""
    m = _MONDO_TEMPLATE.match((text or "").strip())
    return m.group(1).strip() if m else ""
