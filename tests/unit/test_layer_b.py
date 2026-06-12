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


def test_clinical_description_from_mesh_scope_note():
    """The MeSH scope note fills the Clinical-features zone for common diseases
    where HPO is empty; '' when there is no scope note."""
    from atlas.disease import render as DR
    md = DR.r_clinical_description({"mesh_scope_note": "Persistently high systemic "
                                    "arterial BLOOD PRESSURE."})
    assert "## Clinical description" in md
    assert "Persistently high systemic arterial BLOOD PRESSURE." in md
    assert "MeSH descriptor scope note" in md
    assert DR.r_clinical_description({}) == ""
    assert DR.r_clinical_description({"mesh_scope_note": "  "}) == ""


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
    # with a residue table underneath (web renderer wraps each into a collapsible)
    assert "### Catalytic / active sites (1)" in md
    assert "| 837 | proton acceptor |" in md    # perres: Position | Role table
    assert "### Ligand- & substrate-binding residues (1)" in md
    assert "| 718–726 |" in md
    assert "### Post-translational modifications (1)" in md
    assert "### Disulfide bonds (1)" in md
    assert "### Mutagenesis-validated functional residues (1)" in md
    assert "strand" not in md                  # secondary structure not surfaced


def test_residue_map_labels_modifications_by_type():
    """PTM / glycosylation positions are grouped by modification identity (from
    the UniProt feature description) into a `Type | Positions` table instead of a
    bare, meaningless number list. Unlabeled features (most intrachain disulfide
    bonds) collapse onto one row tagged with the group's singular noun."""
    b = {"reviewed_uniprot": ["P0"], "canonical_uniprot": "P0", "ufeatures": [
        _feat("P0", "modified residue", 229, d="phosphoserine"),
        _feat("P0", "modified residue", 695, d="phosphoserine"),
        _feat("P0", "modified residue", 869, d="phosphotyrosine; by src"),
        _feat("P0", "glycosylation site", 56, d="n-linked (glcnac...) asparagine"),
        _feat("P0", "disulfide bond", 31, 58),     # no description
    ]}
    md = R.r_residue_map(b)
    assert "| Phosphoserine | 229, 695 |" in md    # grouped, ordered, no kinase tail
    assert "| Phosphotyrosine | 869 |" in md        # '; by src' tail dropped
    assert "| N-linked (glcnac...) asparagine | 56 |" in md
    assert "### Disulfide bonds (1)" in md
    # the disulfide bond carries no description → falls onto a 'Disulfide bond' row
    assert "| Disulfide bond | 31–58 |" in md


def test_residue_map_binding_site_shows_ligand():
    """Binding-site rows show the ligand (ATP, Mg(2+)…) from the clean `ligand`
    field, not the lowercased description; active sites keep their description."""
    b = {"reviewed_uniprot": ["P0"], "canonical_uniprot": "P0", "ufeatures": [
        {"uniprot": "P0", "type": "binding site", "begin": "718", "end": "726",
         "description": "atp", "ligand": "ATP"},
        {"uniprot": "P0", "type": "active site", "begin": "837", "end": "837",
         "description": "proton acceptor", "ligand": None},
    ]}
    md = R.r_residue_map(b)
    assert "| 718–726 | ATP |" in md            # clean ligand casing, not 'atp'
    assert "| 837 | proton acceptor |" in md     # active site falls back to description


def test_hpa_protein_shows_antibody_reliability():
    """HPA protein block surfaces the antibody-staining reliability tier (IH/IF),
    omitting whichever field is absent."""
    md = R.r_hpa_protein({"13": {"hpa": {
        "subcellular_main": ["Nucleoplasm"],
        "reliability_ih": "enhanced", "reliability_if": "supported"}}})
    assert "**Antibody reliability (HPA):** IH enhanced · IF supported" in md
    # IF absent → only IH shown, no dangling separator
    md2 = R.r_hpa_protein({"13": {"hpa": {
        "subcellular_main": ["Cytosol"], "reliability_ih": "enhanced"}}})
    assert "**Antibody reliability (HPA):** IH enhanced" in md2
    assert "·" not in md2.split("Antibody reliability")[1]


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


def test_sense_gene_orientation_bullet(monkeypatch):
    """An antisense lncRNA orients to its sense gene, quoting the SENSE gene's
    biology (attributed), manifest-gated."""
    from atlas.page import at_a_glance as AAG
    from atlas.page import evidence, links
    monkeypatch.setattr(links, "gene_url",
                        lambda symbol=None, **k: f"/atlas/gene/{symbol}/" if symbol == "FER" else None)
    evidence._COMP = {"gene": {"FER": {"gwas_count": 33, "drug_count": 56}}}
    b = AAG._sense_gene_bullet("FER-AS1")
    assert "antisense to" in b and "[FER](/atlas/gene/FER/)" in b
    assert "33 GWAS associations" in b and "56 ChEMBL molecules" in b
    assert "not this transcript" in b                       # honest attribution
    assert AAG._sense_gene_bullet("TP53") == ""             # not antisense
    assert AAG._sense_gene_bullet("XXX-AS1") == ""          # sense gene not built
    evidence.reset()


def test_parent_evidence_clause(monkeypatch):
    """A thin subtype quantifies its parent's evidence from the parent's frozen
    components — attributed, never folded into this page's score."""
    from atlas.disease import render as DR
    from atlas.page import evidence, links
    monkeypatch.setattr(links, "_lookup", lambda et, *a: "heavy-metal-poisoning")
    evidence._COMP = {"disease": {"heavy-metal-poisoning":
                                  {"trial_count": 4, "gwas_count": 1, "drug_count": 1}}}
    c = DR._parent_evidence_clause({"id": "MONDO:1", "name": "heavy metal poisoning"})
    assert "4 clinical trials" in c and "1 GWAS association" in c and "1 approved drug" in c
    evidence.reset()


def test_mechanism_synthesis(monkeypatch):
    """Cross-entity join: indicated drugs that also target a cohort gene, matched
    by manifest URL (never by fuzzy name)."""
    from atlas.disease import render as DR
    from atlas.page import links
    monkeypatch.setattr(links, "drug_url",
                        lambda chembl_id=None, name=None: f"/atlas/drug/{(chembl_id or name or '').lower()}/")
    bundles = {
        "5": {"gene_count": 8},
        "10": {"drugs": [{"id": "CHEMBL1", "name": "Vemurafenib", "gene_targets": ["BRAF"]},
                         {"id": "CHEMBL2", "name": "OtherDrug", "gene_targets": ["XYZ"]}]},
        "_indicated_drugs": [{"name": "Vemurafenib", "url": "/atlas/drug/chembl1/"},
                             {"name": "Aspirin", "url": "/atlas/drug/aspirin/"}],
    }
    md = DR.r_mechanism_synthesis(bundles)
    assert "Mechanistic alignment" in md
    assert "1 of the 2 drugs indicated" in md
    assert "Vemurafenib" in md and "BRAF" in md
    assert "Aspirin" not in md                     # not a cohort drug → not aligned
    assert DR.r_mechanism_synthesis({"5": {}, "10": {}, "_indicated_drugs": []}) == ""


def test_clingen_disease_name_dedup():
    """ClinGen condition labels differing only by a trailing numeric subtype
    collapse to the shortest representative (Li-Fraumeni syndrome subsumes
    'Li-Fraumeni syndrome 1'); genuinely distinct names are kept."""
    from atlas.gene.sections.s06_variants import _dedup_disease_names as dd
    assert dd(["Li-Fraumeni syndrome", "Li-Fraumeni syndrome 1"]) == ["Li-Fraumeni syndrome"]
    assert dd(["Li-Fraumeni syndrome 1"]) == ["Li-Fraumeni syndrome 1"]   # nothing to subsume
    assert dd(["BRCA1-related cancer", "Lynch syndrome"]) == ["BRCA1-related cancer", "Lynch syndrome"]
    assert dd(["", None, "  "]) == []
