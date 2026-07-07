"""Corpus FAQPage JSON-LD — built from title + description + TL;DR, using the
real tldr shapes from gene/disease pages. Pure; no biobtree."""
from atlas.page.faq_jsonld import build_faq, as_script_tag

_TLDR = [  # verbatim from a built TP53 page
    "Encodes Cellular tumor antigen p53 (UniProt P04637)",
    "GWAS associations: 50 (top 2% of genes corpus-wide)",
    "Clinical variants (ClinVar): 3,923 total — 771 pathogenic",
]


def test_build_faq_two_questions():
    faq = build_faq("TP53", "TP53 encodes a tumor suppressor.", _TLDR,
                    "https://sugi.bio/atlas/gene/TP53/")
    assert faq["@type"] == "FAQPage" and faq["@id"].endswith("#faq")
    qs = {q["name"]: q["acceptedAnswer"]["text"] for q in faq["mainEntity"]}
    assert "What is TP53?" in qs
    kf = "What are the key facts about TP53?"
    assert kf in qs and "GWAS associations: 50" in qs[kf] and "3,923 total" in qs[kf]


def test_build_faq_description_only_and_facts_only():
    d_only = build_faq("X", "X is a thing.", [], "u")
    assert len(d_only["mainEntity"]) == 1 and "What is X?" in d_only["mainEntity"][0]["name"]
    f_only = build_faq("Y", "", ["Classification: Cancer"], "u")
    assert len(f_only["mainEntity"]) == 1 and "key facts" in f_only["mainEntity"][0]["name"]


def test_build_faq_empty_returns_none():
    assert build_faq("Z", "", [], "u") is None
    assert build_faq("", "desc", ["fact"], "u") is None
    assert as_script_tag("Z", "", [], "u") == ""


def test_script_tag_is_valid_json():
    import json, re
    tag = as_script_tag("TP53", "desc", _TLDR, "https://sugi.bio/atlas/gene/TP53/")
    body = re.search(r">\n(.*)\n</script>", tag, re.S).group(1)
    assert json.loads(body)["@type"] == "FAQPage"
