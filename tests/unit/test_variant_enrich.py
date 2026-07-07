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
