"""GWAS p-value display tidying — "8.000000e-11" → "8e-11"."""
from atlas.render_common import pval


def test_pval_tidies_exponent_form():
    assert pval("8.000000e-11") == "8e-11"
    assert pval("5.000000e-08") == "5e-8"
    assert pval("1.500000e-08") == "1.5e-8"
    assert pval("1.000000e-08") == "1e-8"


def test_pval_passes_through_non_numeric_and_empty():
    assert pval("") == ""
    assert pval(None) == ""
    assert pval("NS") == "NS"
    assert pval("0") == "0"          # zero/garbage isn't a real p-value
