"""FAERS adverse-event aggregation (§13) — pure grouping of faers_reaction rows
across a drug's FAERS name records. The live collect() hits biobtree and is
exercised by the corpus build; here we pin the deterministic aggregation."""
from atlas.drug.sections.s13_faers import _aggregate_reactions


# Three name-record rows for the same drug: "diarrhoea" appears under two
# records (791@1.90 and 1766@2.09), "oral pigmentation" once with a high PRR
# off few reports, "rash" once below the PRR report floor.
_ROWS = [
    {"reaction": "diarrhoea", "report_count": "791", "prr": "1.90", "serious_count": "692"},
    {"reaction": "diarrhoea", "report_count": "1766", "prr": "2.09", "serious_count": "1500"},
    {"reaction": "oral pigmentation", "report_count": "18", "prr": "536.0", "serious_count": "3"},
    {"reaction": "rash", "report_count": "5", "prr": "300.0", "serious_count": "5"},
]


def test_sums_report_count_across_name_records():
    most, _disp, distinct = _aggregate_reactions(_ROWS)
    assert distinct == 3                                  # diarrhoea merged
    diarrhoea = next(r for r in most if r["reaction"] == "diarrhoea")
    assert diarrhoea["report_count"] == 791 + 1766
    assert diarrhoea["serious_count"] == 692 + 1500       # serious summed too
    assert diarrhoea["records"] == 2


def test_prr_is_report_count_weighted():
    most, _disp, _ = _aggregate_reactions(_ROWS)
    diarrhoea = next(r for r in most if r["reaction"] == "diarrhoea")
    # (791*1.90 + 1766*2.09) / (791+1766) = 2.031...
    assert diarrhoea["prr"] == round((791 * 1.90 + 1766 * 2.09) / (791 + 1766), 2)


def test_most_reported_sorted_by_volume():
    most, _disp, _ = _aggregate_reactions(_ROWS)
    assert [r["reaction"] for r in most][0] == "diarrhoea"   # 2557 reports


def test_disproportionate_respects_report_floor():
    # rash (PRR 300 but only 5 reports) is below the default floor of 10 and must
    # be excluded; oral pigmentation (18 reports) survives and leads by PRR.
    _most, disp, _ = _aggregate_reactions(_ROWS, prr_min_reports=10)
    names = [r["reaction"] for r in disp]
    assert "rash" not in names
    assert names[0] == "oral pigmentation"


def test_empty_and_zero_rows():
    most, disp, distinct = _aggregate_reactions(
        [{"reaction": "", "report_count": "5", "prr": "2"},
         {"reaction": "x", "report_count": "0", "prr": "2"}])
    assert most == [] and disp == [] and distinct == 0
