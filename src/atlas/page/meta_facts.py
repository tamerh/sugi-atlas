"""Frontmatter key-facts (web-team P2/P3, docs/internal/research/NEXT.md) — derived
deterministically from the bundle, structured for the frontend "Key facts" card
+ Pagefind, separate from the prose body:

  identifier        — a TYPED, trustworthy id per entity (gene HGNC symbol,
                      disease MONDO id, drug ChEMBL id); `symbol` stays the URL
                      slug. Templates should key on `identifier`, not `symbol`.
  alt_names         — search aliases (NOT Hugo-reserved `aliases:`): gene
                      prev/alias symbols, disease synonyms, drug brand/INN names.
  tldr              — 3–5 key-fact bullets (reuses the At-a-glance digest).
  section_defaults  — open/collapsed hints keyed by canonical anchor id (P3).
"""
import importlib
import re

_AT_A_GLANCE = {
    "gene": "atlas.page.at_a_glance",
    "disease": "atlas.page.disease_at_a_glance",
    "drug": "atlas.page.drug_at_a_glance",
}


def _identifier(entity_type, b1):
    if entity_type == "gene":
        return b1.get("symbol") or b1.get("hgnc_id")
    if entity_type == "disease":
        return b1.get("mondo_id")
    return b1.get("chembl_id")          # drug


def _clean_aliases(names, cap=20):
    """Strip + case-insensitively dedup synonyms (keep first-seen form). Does NOT
    split on commas: comma-splitting shatters chemical names (gemcitabine's IUPAC
    name → "2'", "5R)-3") and Mondo synonyms with internal commas ("AT,
    complementation group A"). Drug brand/code pairs ('GLEEVEC,STI-571') are
    split upstream in the drug anchor, where a separator comma can be told apart
    from a chemistry comma; here we just normalize."""
    out, seen = [], set()
    for raw in names or []:
        p = str(raw).strip()
        k = p.lower()
        # Drop malformed aliases with unbalanced parentheses — a Mondo source
        # quirk ("…type V)") or an encoding artifact ("ChC)diak-Higashi", a
        # mangled "Chédiak"). The well-formed name survives as a separate synonym.
        if p and k not in seen and p.count("(") == p.count(")"):
            seen.add(k)
            out.append(p)
        if len(out) >= cap:
            break
    return out


def _alt_names(entity_type, b1):
    if entity_type == "gene":
        src = (b1.get("hgnc") or {}).get("aliases")
        title = (b1.get("hgnc") or {}).get("symbol")
    elif entity_type == "disease":
        src = b1.get("synonyms")
        title = b1.get("canonical_name")
    else:
        src = b1.get("alt_names")        # drug brand / INN synonyms
        title = b1.get("canonical_name")
    # Drop the entity's own title (case-insensitive) — a self-reference chip
    # ("GEMCITABINE" alt-name on the Gemcitabine page).
    t = (title or "").strip().lower()
    return [a for a in _clean_aliases(src) if a.strip().lower() != t]


def _tldr(entity_type, bundle):
    """Reuse the At-a-glance bullets (already curated key facts) as the
    structured TL;DR; prepend the encoded-protein identity for genes."""
    mod = importlib.import_module(_AT_A_GLANCE[entity_type])
    md = mod.at_a_glance(bundle) or ""
    out = [b.replace("**", "").strip() for b in re.findall(r"^- (.+)$", md, re.M)]
    if entity_type == "gene":
        b3 = bundle.get("3") or {}
        pn, uni = b3.get("protein_name"), b3.get("canonical_uniprot")
        if pn and uni:
            out.insert(0, f"Encodes {pn} (UniProt {uni})")
    return out[:5]


def _section_defaults(entity_type, bundle):
    """P3 — open the lede + the section most likely to match intent. Keyed by
    canonical anchor id (PAGE_CONTRACT). The frontend may override."""
    d = {"summary": "open"}
    if entity_type == "gene" and (bundle.get("10") or {}).get("is_drug_target"):
        d["drugs"] = "open"
    elif entity_type == "disease":
        # Clinical presentation is the headline — open it when present.
        if (bundle.get("1") or {}).get("phenotypes"):
            d["clinical"] = "open"
        d["trials" if (bundle.get("13") or {}).get("trial_count") else "genes"] = "open"
    elif entity_type == "drug":
        d["indications"] = "open"
    return d


_DECLARATIVE = {"gene": "atlas.page.declarative",
                "disease": "atlas.page.disease_declarative",
                "drug": "atlas.page.drug_declarative"}


def _meta_description(entity_type, bundle, limit=155):
    """Per-page <meta description>: the declarative lead, stripped of markdown and
    truncated at a word boundary. Each page gets a distinct description instead of
    the generic site-wide fallback (audit: all 52k shared one)."""
    try:
        mod = importlib.import_module(_DECLARATIVE[entity_type])
        s = mod.declarative_sentence(bundle) or ""
    except Exception:
        return ""
    s = re.sub(r"\[([^\]]+)\]\([^)]*\)", r"\1", s)   # link → label
    s = re.sub(r"[*`_]", "", s)                      # drop emphasis
    s = re.sub(r"\s+", " ", s).strip()
    if len(s) > limit:
        s = s[:limit].rsplit(" ", 1)[0].rstrip(" ,;.") + "…"
    return s


def entity_facts(entity_type, bundle):
    """{identifier, alt_names, tldr, description, section_defaults} for the
    frontmatter."""
    if entity_type == "pathway":
        # Flat (non-section) bundle; description = the declarative lead, markdown-
        # stripped, for the meta/SEO description.
        from atlas.pathway.render import declarative_sentence
        desc = re.sub(r"[*`\[\]]|\((?:/|http)[^)]*\)", "", declarative_sentence(bundle)).strip()
        return {"identifier": bundle.get("reactome_id"), "alt_names": [], "tldr": [],
                "description": desc[:300], "section_defaults": {}}
    b1 = bundle.get("1") or {}
    return {
        "identifier": _identifier(entity_type, b1),
        "alt_names": _alt_names(entity_type, b1),
        "tldr": _tldr(entity_type, bundle),
        "description": _meta_description(entity_type, bundle),
        "section_defaults": _section_defaults(entity_type, bundle),
    }
