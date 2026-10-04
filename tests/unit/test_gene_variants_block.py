"""Gene #variants block (v1.11.8): summary lead, Sugi Variant link right under the
ClinVar sample, AlphaMissense as hotspot residues, no dbSNP table, PharmGKB shell
rows dropped."""
from atlas.gene import render as R
from atlas.gene.sections.s06_variants import am_hotspots


def test_am_hotspots_grouped_and_ranked():
    rows = [{"protein_variant": v, "am_pathogenicity": s} for v, s in [
        ("R175H", "0.99"), ("R175G", "0.98"), ("R175P", "1.0"),
        ("G245S", "0.97"), ("G245D", "0.99"),
        ("K13E", "1.0"), ("bogus", "1.0"), ("p.Y220C", "0.95")]]
    h = am_hotspots(rows)
    assert [x["residue"] for x in h] == ["R175", "G245", "K13", "Y220"]
    assert h[0] == {"residue": "R175", "n": 3, "max": 1.0}


B6 = {"symbol": "TP53", "clinvar_total": 3000,
      "clinvar_breakdown": {"Pathogenic": 900, "Likely pathogenic": 300, "Benign": 10},
      "top_pathogenic": [{"id": "1", "hgvs": "x", "classification": "Pathogenic",
                          "review_status": "reviewed by expert panel"}],
      "top_pathogenic_total": 1200, "clingen_variant_total": 50,
      "alphamissense_total": 4000, "am_lp_total": 1500, "am_lp_residues": 300,
      "am_hotspots": [{"residue": "R175", "n": 6, "max": 1.0}]}


def test_summary_lead_links_to_evidence_blocks():
    s = R.r_variants_summary({"6": B6,
                              "10": {"civic_variant_total": 9, "pharmgkb_guideline": [{}, {}]},
                              "12": {"gwas_total": 39}})
    assert s.startswith("**Variant evidence for TP53:**")
    for frag in ("[3,000 ClinVar records](#clinvar) (1,200 pathogenic / likely pathogenic)",
                 "(#clingen-variants)", "[9 CIViC cancer variants](#civic-variants)",
                 "[2 PharmGKB dosing guidelines](#pharmgkb-guidelines)",
                 "[39 GWAS associations](#gwas-assoc)", "Sugi Variant"):
        assert frag in s
    assert R.r_variants_summary({"6": {"symbol": "X"}}) == ""


def test_block_order_and_trims():
    md = R.r_variants(B6, summary="**Variant evidence for TP53:** …")
    assert md.startswith("## Genetic variants")
    assert md.index("Variant evidence") < md.index("{#clinvar}")
    assert md.index("{#top-pathogenic}") < md.index("{#sugi-variant}") < md.index("{#clingen-variants}")
    assert "dbsnp" not in md.lower()
    assert "| R175 | 6 | 1 |" in md and "1,500 missense substitutions" in md
    assert "/19" not in md


def test_no_sugi_link_without_plp():
    b = dict(B6, clinvar_breakdown={"Benign": 3})
    assert "{#sugi-variant}" not in R.r_variants(b)


def test_pharmgkb_shell_rows_dropped():
    b = {"symbol": "EGFR", "molecules": [], "molecule_count": 0, "pharmgkb_variant": [
        {"name": "rs712829"},                                            # empty shell
        {"name": "rs2227983", "level_of_evidence": "3", "associated_drugs": "gefitinib"}]}
    md = R.r_drugs(b)
    assert "rs2227983" in md and "rs712829" not in md
    md2 = R.r_drugs(dict(b, pharmgkb_variant=[{"name": "rs712829"}]))
    assert "{#pharmgkb-variants}" not in md2
