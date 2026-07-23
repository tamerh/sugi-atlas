"""Disease §18 single-cell (CZ CELLxGENE) — pure dedup/sort/cap in collect() and
the renderer. The live map_all is exercised by the corpus build; here we pin the
deterministic shaping and the honest-framing render (count + largest datasets,
NO cell-type list — that edge is multi-disease-atlas contaminated)."""
from types import SimpleNamespace

from atlas.disease.sections import s18_single_cell as SC
from atlas.disease.render import r_single_cell


_ROWS = [
    {"id": "d1", "title": "Global Atlas",     "organism": "Homo sapiens", "cell_count": "621200"},
    {"id": "d2", "title": "Immune",           "organism": "Homo sapiens", "cell_count": "274555"},
    {"id": "d1", "title": "Global Atlas dup",  "organism": "Homo sapiens", "cell_count": "621200"},  # dup id
    {"id": "d3", "title": "No count",         "organism": "Homo sapiens", "cell_count": None},        # → 0
    {"id": "d4", "title": "Garbage count",    "organism": "Homo sapiens", "cell_count": "n/a"},       # → 0
]


def _collect(monkeypatch, rows, n_pad=0):
    padded = rows + [{"id": f"p{i}", "title": f"pad{i}", "organism": "Homo sapiens",
                      "cell_count": str(i + 1)} for i in range(n_pad)]
    monkeypatch.setattr(SC, "map_all", lambda *a, **k: padded)
    return SC.collect(SimpleNamespace(mondo_id="MONDO:0000001"))


def test_collect_dedupes_sorts_and_sums(monkeypatch):
    b = _collect(monkeypatch, _ROWS)
    assert b["dataset_count"] == 4                     # d1 deduped (5 rows → 4)
    assert b["total_cells"] == 621200 + 274555         # garbage/None count as 0
    titles = [d["title"] for d in b["top_datasets"]]
    assert titles[0] == "Global Atlas" and titles[1] == "Immune"   # sorted by cells desc
    assert b["top_datasets"][0]["cells"] == 621200


def test_collect_caps_display_but_counts_all(monkeypatch):
    b = _collect(monkeypatch, _ROWS, n_pad=20)         # 24 distinct datasets
    assert b["dataset_count"] == 24                    # full count preserved
    assert len(b["top_datasets"]) == SC._TOP_DATASETS  # display capped at 8


def test_render_empty_elides():
    assert r_single_cell({"dataset_count": 0}) == ""   # bonus resource, no placeholder


def test_render_count_framing_and_source(monkeypatch):
    b = _collect(monkeypatch, _ROWS)
    md = r_single_cell(b)
    assert "4 single-cell RNA-seq datasets" in md      # plural + count
    assert "annotated with this disease" in md         # honest framing, not "studies of"
    assert "multi-disease atlases" in md               # the contamination caveat
    assert "CELLxGENE Census (CC BY 4.0)" in md        # attribution
    assert "cell type" not in md.lower()               # cell-type edge NOT surfaced


def test_render_singular_grammar(monkeypatch):
    b = _collect(monkeypatch, _ROWS[:1])
    md = r_single_cell(b)
    assert "1 single-cell RNA-seq dataset**" in md and " is annotated" in md
