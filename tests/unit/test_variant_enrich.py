"""Deterministic variant-enrichment logic: missense parsing, frequency banding,
submitter consensus, and the cross-source concordance verdict. Pure functions —
no biobtree. These gate the derived claims before any regen."""
from atlas.variant.enrich import (missense_short, protein_position, _freq_band,
                                  submitter_consensus, concordance, residue_hotspot)


def test_missense_short_and_position():
    assert missense_short("p.Pro309Ala") == "P309A"
    assert missense_short("p.Arg14Ter") == "R14*"
    assert missense_short("p.Lys328del") is None      # not a missense → AM n/a
    assert missense_short("p.Asn117fs") is None
    assert protein_position("p.Pro309Ala") == 309
    assert protein_position("p.Lys328del") == 328


def test_freq_band():
    assert _freq_band("", True) == "absent from gnomAD"
    assert _freq_band("0.00005", False).startswith("ultra-rare")
    assert _freq_band("0.0005", False).startswith("rare")
    assert _freq_band("0.27", False).startswith("common")


def test_submitter_consensus():
    uni = submitter_consensus([{"classification": "Pathogenic"}] * 3)
    assert uni["n"] == 3 and "unanimous" in uni["verdict"]
    conf = submitter_consensus([{"classification": "Pathogenic"}] * 2
                               + [{"classification": "Uncertain significance"}])
    assert "conflicting" in conf["verdict"] and "2 Pathogenic" in conf["verdict"]
    assert submitter_consensus([]) is None


def test_concordance_agree_and_flag():
    # ClinVar Pathogenic + AlphaMissense likely-pathogenic + absent from gnomAD → concordant
    ok = concordance("Pathogenic", ("likely_pathogenic", "0.99"),
                     {"absent": True, "is_common": False, "frequency": "", "band": "absent from gnomAD"})
    assert "concordant" in ok["verdict"].lower() and not ok["flags"]
    # ClinVar Pathogenic but common in gnomAD → disagreement flagged
    bad = concordance("Pathogenic", None,
                      {"absent": False, "is_common": True, "frequency": "0.2", "band": "common (gnomAD MAF 0.2)"})
    assert bad["flags"] and "disagree" in bad["verdict"].lower()
    # AlphaMissense likely-benign vs ClinVar Pathogenic → discordant flag
    disc = concordance("Pathogenic", ("likely_benign", "0.1"), None)
    assert disc["flags"]


def test_residue_hotspot():
    idx = {309: [{"hgvs_p": "p.Pro309Ala", "label": "ACTA1 p.Pro309Ala", "slug": "a", "classification": "Pathogenic"},
                 {"hgvs_p": "p.Pro309Leu", "label": "ACTA1 p.Pro309Leu", "slug": "b", "classification": "Pathogenic"}]}
    hs = residue_hotspot("p.Pro309Ala", idx)
    assert hs["position"] == 309 and len(hs["others"]) == 1        # excludes itself
    assert residue_hotspot("p.Pro309Ala", {309: [{"hgvs_p": "p.Pro309Ala", "slug": "a"}]}) is None


# ── Batch 3 derived logic ────────────────────────────────────────────────────
def test_am_percentile():
    from atlas.variant.enrich import am_percentile
    am = {f"V{i}": ("x", str(i / 100)) for i in range(100)}   # scores 0.00..0.99
    p = am_percentile("0.98", am)
    assert p["n"] == 100 and p["top_pct"] <= 3          # 0.98 is near the top
    assert am_percentile("0.5", {"a": ("x", "0.1")}) is None   # too few to rank


def test_structural_context():
    from atlas.variant.enrich import structural_context
    intervals = [{"type": "helix", "desc": "", "begin": 300, "end": 320},
                 {"type": "region of interest", "desc": "actin-binding", "begin": 100, "end": 150}]
    sc = structural_context("p.Pro309Ala", intervals)
    assert sc["position"] == 309 and any("helix" in f for f in sc["features"])
    assert structural_context("p.Pro200Ala", intervals) is None   # 200 in no interval


def test_similar_variants_ranks_same_residue_first():
    from atlas.variant.enrich import similar_variants
    me = {"canonical_slug": "g-p-pro9ala", "hgvs_p": "p.Pro9Ala", "gene_symbol": "G",
          "variant_type": "single nucleotide variant", "conditions": [{"mondo_id": "MONDO:1"}]}
    recs = [me,
            {"canonical_slug": "g-p-pro9leu", "hgvs_p": "p.Pro9Leu", "gene_symbol": "G",
             "variant_type": "x", "conditions": []},                         # same residue (score 3)
            {"canonical_slug": "g-p-arg5his", "hgvs_p": "p.Arg5His", "gene_symbol": "G",
             "variant_type": "y", "conditions": [{"mondo_id": "MONDO:1"}]}]   # same condition (score 2)
    sim = similar_variants(me, recs)
    assert sim[0]["slug"] == "g-p-pro9leu"          # same-residue ranks first
    assert me["canonical_slug"] not in [s["slug"] for s in sim]


def test_submission_timeline():
    from atlas.variant.enrich import submission_timeline
    tl = submission_timeline([{"classification": "Pathogenic", "date": "2002-01-01"},
                              {"classification": "Pathogenic", "date": "2018-06-01"}])
    assert tl["first"] == 2002 and tl["last"] == 2018 and tl["stable"]
    tl2 = submission_timeline([{"classification": "Pathogenic", "date": "2010-01-01"},
                               {"classification": "Uncertain significance", "date": "2020-01-01"}])
    assert not tl2["stable"]


def test_plain_summary():
    from atlas.variant.enrich import plain_summary
    s = plain_summary({"gene_symbol": "ACTA1", "classification": "Pathogenic",
                       "submitter_count": 3, "conditions": [{"name": "nemaline myopathy"}]}, None)
    assert "ACTA1" in s and "disease-causing" in s and "3 clinical labs" in s
    conf = plain_summary({"gene_symbol": "X", "classification": "Conflicting classifications of pathogenicity",
                          "submitter_count": 0, "conditions": []}, None)
    assert "experts disagree" in conf
