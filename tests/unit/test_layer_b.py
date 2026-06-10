"""Gene page canonical H2 taxonomy + the functional-residue map
(docs/PAGE_CONTRACT.md, docs/MOLECULAR_ENRICHMENT.md)."""
import atlas.pipeline as P
from atlas.gene import render as R

# The FROZEN gene H2 set (order matters) — Summary/#related are added by
# assemble_page/related_block, so render_all emits these six.
GENE_H2 = [
    "## Identifiers {#identifiers}",
    "## Gene structure {#gene-structure}",
    "## Protein {#protein}",
    "## Function {#function}",
    "## Disease & clinical {#disease}",
    "## Drugs & pharmacology {#drugs}",
]


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


def test_canonical_h2_set_and_order_coding(monkeypatch):
    _stub_renderers(monkeypatch)
    bundle = {s: {} for s in (str(i) for i in range(1, 13))}
    bundle["3"] = {"reviewed_uniprot": ["P1"], "canonical_uniprot": "P1"}
    md = P.render_all(bundle)
    assert _h2(md) == GENE_H2                      # exact set + frozen order
    assert '<a id="protein-P1"></a>' in md         # JSON-LD @id linkage kept
    # sub-sections demoted to H3 under their canonical H2
    assert "### Section 1" in md and "\n## Section 1" not in md
    prot = md.split("{#protein}")[1].split("{#function}")[0]
    assert "### Functional residue map" in prot
    gs = md.split("{#gene-structure}")[1].split("{#protein}")[0]
    assert "### Functional genomics" in gs
    assert "### Section 10" in md.split("{#drugs}")[1]   # drugs


def test_ncrna_emits_all_sections_with_placeholders(monkeypatch):
    """Emit-even-if-empty: the TOC is identical across genes; ncRNA shows the
    same H2 set, with informative placeholders for the protein-layer sections."""
    _stub_renderers(monkeypatch)
    bundle = {s: {} for s in (str(i) for i in range(1, 13))}
    bundle["3"] = {"reviewed_uniprot": [], "canonical_uniprot": None}
    bundle["_noncoding"] = "lncRNA"
    md = P.render_all(bundle)
    assert _h2(md) == GENE_H2                      # same set, not collapsed
    assert "*Non-coding RNA — no protein product" in md.split("{#protein}")[1]
    assert "protein-" not in md                    # no JSON-LD @id anchor for ncRNA


DISEASE_H2 = [
    "## Clinical features {#clinical}",
    "## Identifiers {#identifiers}",
    "## Disease family {#family}",
    "## Genetics & variants {#genetics}",
    "## Genes & proteins {#genes}",
    "## Function {#function}",
    "## Therapeutics {#drugs}",
    "## Clinical trials & evidence {#trials}",
]
DRUG_H2 = [
    "## Identifiers {#identifiers}",
    "## Targets {#targets}",
    "## Indications & clinical {#indications}",
    "## Pharmacology {#pharmacology}",
    "## Related molecules {#related-molecules}",
]


def test_disease_canonical_h2_set_and_order():
    """Every canonical section always emits, in the frozen order
    (docs/PAGE_CONTRACT.md) — sub-renderers carry their own empty-state text, so
    the H2 set is identical across every disease page even with no data."""
    from atlas.disease import render as DR
    md = DR.render_all({str(i): {} for i in range(1, 15)})
    assert _h2(md) == DISEASE_H2


def test_drug_canonical_h2_set_and_order():
    from atlas.drug import render as DRR
    md = DRR.render_all({str(i): {} for i in range(1, 13)})
    assert _h2(md) == DRUG_H2


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
    # each residue category is now its own H4 sub-section (count kept in heading),
    # with the residue body underneath (web renderer wraps each into a collapsible)
    assert "### Catalytic / active sites (1)" in md
    assert "**837** (proton acceptor)" in md
    assert "### Ligand- & substrate-binding residues (1)" in md
    assert "**718–726**" in md
    assert "### Post-translational modifications (1)" in md
    assert "### Disulfide bonds (1)" in md
    assert "### Mutagenesis-validated functional residues (1)" in md
    assert "strand" not in md                  # secondary structure not surfaced


def test_residue_map_per_product_skips_empty():
    b = {"reviewed_uniprot": ["P0", "P1", "P2"], "canonical_uniprot": "P0",
         "ufeatures": [
            _feat("P0", "active site", 10),
            _feat("P1", "binding site", 20, 25),
            _feat("P2", "chain", 1, 99),        # only a non-residue feature
         ]}
    md = R.r_residue_map(b)
    # two products with residues → each category H4 carries a product suffix
    # (no separate per-product bold header — keeps the heading hierarchy flat)
    assert "### Catalytic / active sites (1) — P0 (canonical)" in md
    assert "### Ligand- & substrate-binding residues (1) — P1" in md
    assert "P2" not in md                        # residue-less product gets no heading


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
