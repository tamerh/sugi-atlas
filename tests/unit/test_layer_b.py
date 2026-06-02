"""Layer B — Gene/Protein/Clinical body zoning + the functional-residue map
(docs/MOLECULAR_ENRICHMENT.md)."""
import atlas.pipeline as P
from atlas.gene import render as R


# ── zoning ──────────────────────────────────────────────────────────────────
def _stub_renderers(monkeypatch):
    monkeypatch.setattr(R, "RENDER",
                        {s: (lambda b, _s=s: f"## Section {_s}\n\nbody{_s}")
                         for s in (str(i) for i in range(1, 13))})
    monkeypatch.setattr(R, "r_functional_genomics", lambda b: "## Functional genomics\n\nfg")
    monkeypatch.setattr(R, "r_generifs", lambda b: "## GeneRIF\n\ngr")
    monkeypatch.setattr(R, "r_residue_map", lambda b: "## Functional residue map\n\nrm")


def _h2(md):
    return [l for l in md.splitlines() if l.startswith("## ")]


def test_three_zones_for_coding_gene(monkeypatch):
    _stub_renderers(monkeypatch)
    bundle = {s: {} for s in (str(i) for i in range(1, 13))}
    bundle["3"] = {"reviewed_uniprot": ["P1"], "canonical_uniprot": "P1"}
    md = P.render_all(bundle)
    assert _h2(md) == ["## Gene — the locus", "## Protein product(s)",
                       "## Clinical & disease"]
    assert '<a id="protein-P1"></a>' in md
    # sections demoted under their zone
    assert "### Section 1" in md
    assert "\n## Section 1" not in md
    # protein-layer content lands in the protein zone
    prot = md.split("## Protein product(s)")[1].split("## Clinical")[0]
    assert "### Functional residue map" in prot
    assert "### Section 10" in prot          # drugs
    # locus content in the gene zone
    gene = md.split("## Gene — the locus")[1].split("## Protein")[0]
    assert "### Functional genomics" in gene


def test_ncrna_collapses_to_gene_zone(monkeypatch):
    _stub_renderers(monkeypatch)
    bundle = {s: {} for s in (str(i) for i in range(1, 13))}
    bundle["3"] = {"reviewed_uniprot": [], "canonical_uniprot": None}
    bundle["_noncoding"] = "lncRNA"
    md = P.render_all(bundle)
    assert _h2(md) == ["## Gene — the locus"]
    assert "protein-" not in md


def test_demote_bumps_headings():
    assert P._demote("## A\ntext\n### B") == "### A\ntext\n#### B"
    assert P._demote("no heading") == "no heading"


# ── residue map ───────────────────────────────────────────────────────────────
def _feat(u, t, b, e=None, d=None):
    return {"uniprot": u, "type": t, "begin": str(b),
            "end": str(e if e is not None else b), "description": d}


def test_residue_map_categorizes_and_orders():
    b = {"reviewed_uniprot": ["P0"], "canonical_uniprot": "P0", "ufeatures": [
        _feat("P0", "active site", 837, d="proton acceptor"),
        _feat("P0", "binding site", 718, 726),
        _feat("P0", "modified residue", 693),
        _feat("P0", "disulfide bond", 31, 58),
        _feat("P0", "glycosylation site", 56),
        _feat("P0", "mutagenesis site", 275, d="reduced activity"),
        _feat("P0", "strand", 1, 5),          # structural — excluded
    ]}
    md = R.r_residue_map(b)
    assert "## Functional residue map" in md
    assert "**Catalytic / active sites (1):** **837** (proton acceptor)" in md
    assert "**Ligand- & substrate-binding residues (1):** **718–726**" in md
    assert "Post-translational modifications (1)" in md
    assert "Disulfide bonds (1)" in md
    assert "Mutagenesis-validated functional residues (1)" in md
    assert "strand" not in md                  # secondary structure not surfaced


def test_residue_map_per_product_skips_empty():
    b = {"reviewed_uniprot": ["P0", "P1", "P2"], "canonical_uniprot": "P0",
         "ufeatures": [
            _feat("P0", "active site", 10),
            _feat("P1", "binding site", 20, 25),
            _feat("P2", "chain", 1, 99),        # only a non-residue feature
         ]}
    md = R.r_residue_map(b)
    assert "**P0 (canonical)**" in md           # two products with residues → headers
    assert "**P1**" in md
    assert "**P2**" not in md                    # empty product gets no orphan header


def test_residue_map_single_product_no_header():
    """One product with residues → no redundant per-product header."""
    b = {"reviewed_uniprot": ["P0", "P9"], "canonical_uniprot": "P0", "ufeatures": [
        _feat("P0", "active site", 10),
        _feat("P9", "chain", 1, 50),
    ]}
    md = R.r_residue_map(b)
    assert "**P0 (canonical)**" not in md
    assert "Catalytic / active sites (1)" in md


def test_residue_map_empty_when_no_features():
    assert R.r_residue_map({"reviewed_uniprot": ["P0"], "ufeatures": []}) == ""
