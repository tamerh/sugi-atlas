"""Drug-page MeSH pharmacology (§15) + PharmGKB drug labels / pathways (§9).
biobtree is monkeypatched with small fixtures shaped like the live responses
(search rows "id|dataset|name|xref_count"; mesh / pharmgkb / chembl_molecule
entry Attributes; map_all target dicts)."""
from types import SimpleNamespace

import pytest

from atlas.drug import render as R
from atlas.drug.sections import s09_pharmacogenomics as P
from atlas.drug.sections import s15_mesh_pharmacology as M


# ───── MeSH note cleaning ─────────────────────────────────────────────────────

@pytest.mark.parametrize("note", [
    "structure in first source",
    "structure given in first source",
    "RN & structure given in first source",
    "RN given refers to parent cpd; structure",
    "was minor descriptor; RN given refers to parent cpd",
    "N Engl J Med.2022 Mar 17;386(11):1013-1025",
    "cerivastatin is the ((E)-(+))-isomer; structure given in first source",
    "Minor descriptor (66-83); on-line & Index Medicus search PROPIOPHENONES (66-83)",
    "InChIKey: LIRYPHYGHXZJBZ-UHFFFAOYSA-N",
    "Synonyms 3024 CERM & respilene refers to di-HCl",
    "Chloride was minor descriptor (75-80)",
    "File maintained to MORPHOLINES (66-86)",
    "Note: tradenames that start with Histex refer to more than one drug.",
    "Was EP to DIPHENOXYLATE (77-81)",
    "N1 in Chemline is same as synonym 8",
    "",
    None,
])
def test_clean_note_indexing_chatter_is_empty(note):
    assert M.clean_note(note) == ""


def test_clean_note_strips_junk_clauses_and_capitalises_scr_note():
    assert (M.clean_note("a HER2 inhibitor with antineoplastic activity; structure in first source")
            == "A HER2 inhibitor with antineoplastic activity.")
    assert (M.clean_note("potent analgesic; RN given refers to ((5R)))-isomer (dezocin); structure")
            == "Potent analgesic.")
    assert (M.clean_note("Phosphate ester of betamethasone; structure in Negwer,5th ed, 4975")
            == "Phosphate ester of betamethasone.")
    assert (M.clean_note("An aldose reductase inhibitor; AS-3201 and SX-3202 are the "
                         "(R-(-))- and (S-(+))-isomers, respectively")
            == "An aldose reductase inhibitor.")
    assert (M.clean_note("Structure & RN given in first source; a non-ergot dopamine D(2)-agonist")
            == "A non-ergot dopamine D(2)-agonist.")
    assert (M.clean_note("A recombinant form of tissue factor pathway inhibitor, RN & amino acid "
                         "sequence in first source; MW about 32 kDa")
            == "A recombinant form of tissue factor pathway inhibitor; MW about 32 kDa.")
    assert (M.clean_note("Surface anesthetic used mainly orally in chronic gastritis; "
                         "INDEX MEDICUS search ETHANOLAMINES (72)")
            == "Surface anesthetic used mainly orally in chronic gastritis.")
    assert (M.clean_note("antiscabies, antipruritic drug; active ingredient in Eurax Cream "
                         "& Lotion (Geigy); request from searcher 3/76")
            == "Antiscabies, antipruritic drug; active ingredient in Eurax Cream & Lotion (Geigy).")
    assert (M.clean_note("was MH 1975-92 (see under SULFONYLUREA COMPOUNDS 1975-90); use "
                         "SULFONYLUREA COMPOUNDS to search GLIBORNURIDE 1975-92; an oral, "
                         "sulfonylurea hypoglycemic agent")
            == "An oral, sulfonylurea hypoglycemic agent.")


def test_clean_note_keeps_descriptor_scope_note_verbatim():
    note = ("An element with atomic symbol O, atomic number 8, and atomic weight "
            "[15.99903; 15.99977]. It is the most abundant element on earth.")
    assert M.clean_note(note, supplementary=False) == note
    # descriptor notes are never clause-filtered (SCR junk patterns would eat this)
    iso = "Dextromethorphan is the d-isomer of the codeine analog of levorphanol."
    assert M.clean_note(iso, supplementary=False) == iso
    # mixed-case first word (mAb, c-Met) is not re-cased
    assert M.clean_note("mAb against CD20") == "mAb against CD20."


# ───── MeSH fixtures ──────────────────────────────────────────────────────────

MESH = {
    "D014859": {"descriptor_ui": "D014859", "descriptor_name": "Warfarin",
                "tree_numbers": ["D03.383.663.283.446.520.914"],
                "scope_note": "An anticoagulant that acts by inhibiting the synthesis of "
                              "vitamin K-dependent coagulation factors.",
                "pharmacological_actions": ["Anticoagulants", "Rodenticides"]},
    "C536683": {"descriptor_ui": "C536683", "descriptor_name": "Warfarin syndrome",
                "is_supplementary": True, "heading_mapped_to": "D000014"},
    "D000014": {"descriptor_name": "Abnormalities, Drug-Induced", "tree_numbers": ["C16.131.077"]},
    "C106538": {"descriptor_ui": "C106538", "descriptor_name": "abacavir", "is_supplementary": True,
                "scope_note": "a carbocyclic nucleoside with potent selective anti-HIV activity",
                "pharmacological_actions": ["D018894", "D019380"], "heading_mapped_to": "D003521"},
    "D003521": {"descriptor_name": "Cyclopropanes", "tree_numbers": ["D02.455.326.271.210"]},
    "D018894": {"descriptor_name": "Reverse Transcriptase Inhibitors", "tree_numbers": ["D27.505"]},
    "D019380": {"descriptor_name": "Anti-HIV Agents", "tree_numbers": ["D27.505.954"]},
    "D000068877": {"descriptor_name": "Imatinib Mesylate", "tree_numbers": ["D03.383.606.405"],
                   "scope_note": "A tyrosine kinase inhibitor and ANTINEOPLASTIC AGENT.",
                   "pharmacological_actions": ["Antineoplastic Agents"]},
    # a same-name disease SCR — must fail the chemical-category guard
    "C999001": {"descriptor_name": "Fooxin", "is_supplementary": True, "heading_mapped_to": "D000014",
                "scope_note": "a rare syndrome named after Dr Fooxin"},
    # same-name SCR + descriptor: descriptor wins
    "C999002": {"descriptor_name": "Barazol", "is_supplementary": True, "heading_mapped_to": "D003521",
                "scope_note": "an SCR note"},
    "D999002": {"descriptor_name": "Barazol", "tree_numbers": ["D02.1"],
                "scope_note": "The Barazol descriptor scope note."},
    # two same-name descriptors: ambiguous → nothing
    "D999003": {"descriptor_name": "Quuxin", "tree_numbers": ["D02.2"], "scope_note": "First quuxin."},
    "D999004": {"descriptor_name": "Quuxin", "tree_numbers": ["D02.3"], "scope_note": "Second quuxin."},
    # SCR with neither a usable note nor a pharmacological action
    "C999005": {"descriptor_name": "zorbamide", "is_supplementary": True,
                "heading_mapped_to": "D003521", "scope_note": "structure in first source"},
}
# search keyword → MeSH ids whose name or entry term carries it
SEARCH = {
    "WARFARIN": ["C536683", "D014859"],
    "ABACAVIR": ["C106538"],
    "IMATINIB": ["D000068877"],            # entry term only — not an exact name
    "IMATINIB MESYLATE": ["D000068877"],
    "FOOXIN": ["C999001"],
    "BARAZOL": ["C999002", "D999002"],
    "QUUXIN": ["D999003", "D999004"],
    "ZORBAMIDE": ["C999005"],
}
PARENTS = {"D014859": [{"id": "D015110", "descriptor_name": "4-Hydroxycoumarins"}],
           "D000068877": [{"id": "D001549", "descriptor_name": "Benzamides"},
                          {"id": "D010879", "descriptor_name": "Piperazines"}]}
CHEMBL_NAMES = {"CHEMBL1642": "IMATINIB MESYLATE", "CHEMBL2386595": "IMATINIB HYDROCHLORIDE",
                "CHEMBL1464": "WARFARIN"}


@pytest.fixture
def mesh(monkeypatch):
    def search(term, source=None):
        assert source == "mesh"
        ids = SEARCH.get(term.upper(), [])
        return {"schema": "id|dataset|name|xref_count",
                "data": [f"{i}|mesh|{MESH[i]['descriptor_name']}|1" for i in ids]}

    def entry(i, source):
        if source == "mesh":
            if i not in MESH:
                raise RuntimeError("entry → HTTP 400")
            return {"Attributes": {"Mesh": MESH[i]}}
        if source == "chembl_molecule":
            return {"Attributes": {"Chembl": {"molecule": {"name": CHEMBL_NAMES.get(i, "")}}}}
        raise AssertionError(source)

    def map_all(i, chain, cap=60):
        assert chain == ">>mesh>>meshparent"
        return PARENTS.get(i, [])

    monkeypatch.setattr(M, "search", search)
    monkeypatch.setattr(M, "entry", entry)
    monkeypatch.setattr(M, "map_all", map_all)


def _anchor(name, chembl="CHEMBL1", parent=None, childs=()):
    return SimpleNamespace(canonical_name=name, name=chembl, chembl_id=chembl,
                           parent_chembl=parent, child_chembls=tuple(childs))


# ───── MeSH collector ─────────────────────────────────────────────────────────

def test_descriptor_exact_name_match_skips_same_keyword_records(mesh):
    b = M.collect(_anchor("WARFARIN", "CHEMBL1464"))
    assert b["mesh_id"] == "D014859"                 # not "Warfarin syndrome"
    assert b["mesh_record_type"] == "descriptor"
    assert b["pharmacological_actions"] == ["Anticoagulants", "Rodenticides"]
    assert b["broader_headings"] == ["4-Hydroxycoumarins"]
    assert b["definition"].startswith("An anticoagulant that acts")
    assert b["matched_via"] == "name"


def test_supplementary_concept_resolves_action_ids_and_mapped_heading(mesh):
    b = M.collect(_anchor("ABACAVIR"))
    assert b["mesh_id"] == "C106538" and b["mesh_record_type"] == "supplementary"
    assert b["pharmacological_actions"] == ["Reverse Transcriptase Inhibitors", "Anti-HIV Agents"]
    assert b["broader_headings"] == ["Cyclopropanes"]
    assert b["definition"] == "A carbocyclic nucleoside with potent selective anti-HIV activity."


def test_entry_term_hit_is_not_an_exact_match_but_salt_form_is(mesh):
    # "IMATINIB" only hits the descriptor via an entry term — the own-name pass
    # must miss; the salt form's exact name ("Imatinib Mesylate") then matches.
    b = M.collect(_anchor("IMATINIB", "CHEMBL941", childs=("CHEMBL2386595", "CHEMBL1642")))
    assert b["mesh_id"] == "D000068877"
    assert b["matched_via"] == "salt" and b["matched_chembl"] == "CHEMBL1642"


def test_salt_page_falls_back_to_parent_name(mesh):
    b = M.collect(_anchor("WARFARIN SODIUM", "CHEMBL1200879", parent="CHEMBL1464"))
    assert b["mesh_id"] == "D014859"
    assert b["matched_via"] == "parent" and b["matched_chembl"] == "CHEMBL1464"


def test_same_name_disease_record_fails_category_guard(mesh):
    assert M.collect(_anchor("FOOXIN"))["mesh_id"] is None


def test_descriptor_beats_same_name_supplementary_concept(mesh):
    assert M.collect(_anchor("BARAZOL"))["mesh_id"] == "D999002"


def test_two_same_name_descriptors_are_ambiguous(mesh):
    assert M.collect(_anchor("QUUXIN"))["mesh_id"] is None


def test_no_match_renders_nothing(mesh):
    b = M.collect(_anchor("NOSUCHDRUG-123"))
    assert b == {"section": "15_mesh_pharmacology", "mesh_id": None}
    assert R.r_mesh_pharmacology(b) == ""
    assert R.r_mesh_pharmacology({}) == ""


def test_record_without_note_or_action_renders_nothing(mesh):
    b = M.collect(_anchor("ZORBAMIDE"))
    assert b["mesh_id"] == "C999005" and b["definition"] == ""
    assert R.r_mesh_pharmacology(b) == ""


# ───── MeSH renderer ──────────────────────────────────────────────────────────

def test_render_mesh_block(mesh):
    md = R.r_mesh_pharmacology(M.collect(_anchor("WARFARIN", "CHEMBL1464")))
    assert md.startswith("## MeSH pharmacological classification\n\n**Definition (MeSH):** An anticoagulant")
    assert "**Pharmacological action:** Anticoagulants; Rodenticides." in md
    assert "**Broader MeSH heading:** 4-Hydroxycoumarins." in md
    assert ("*MeSH descriptor [D014859](https://meshb.nlm.nih.gov/record/ui?ui=D014859) "
            "“Warfarin” — exact-name match.*") in md


def test_render_mesh_block_salt_provenance(mesh):
    md = R.r_mesh_pharmacology(
        M.collect(_anchor("IMATINIB", "CHEMBL941", childs=("CHEMBL1642",))))
    assert "exact-name match on this drug's salt/hydrate form `CHEMBL1642`" in md
    assert "**Broader MeSH headings:** Benzamides; Piperazines." in md


def test_render_mesh_scr_note_label(mesh):
    md = R.r_mesh_pharmacology(M.collect(_anchor("ABACAVIR")))
    assert "**MeSH note:** A carbocyclic nucleoside" in md and "Definition" not in md
    # comma-bearing MeSH names stay unambiguous in the list
    md = R.r_mesh_pharmacology({"mesh_id": "D003061", "mesh_name": "Codeine",
                                "definition": "An opioid analgesic.",
                                "pharmacological_actions": ["Analgesics, Opioid", "Narcotics"]})
    assert "**Pharmacological action:** Analgesics, Opioid; Narcotics." in md


def test_render_mesh_actions_only():
    md = R.r_mesh_pharmacology({"mesh_id": "C582435", "mesh_name": "pembrolizumab",
                                "mesh_record_type": "supplementary", "definition": "",
                                "pharmacological_actions": ["Immune Checkpoint Inhibitors"],
                                "broader_headings": [], "matched_via": "name"})
    assert "Definition" not in md and "Broader" not in md
    assert "**Pharmacological action:** Immune Checkpoint Inhibitors." in md
    assert "MeSH supplementary concept [C582435]" in md


def test_render_all_places_mesh_block_first_in_pharmacology():
    b15 = {"mesh_id": "D014859", "mesh_name": "Warfarin", "mesh_record_type": "descriptor",
           "definition": "An anticoagulant.", "pharmacological_actions": ["Anticoagulants"],
           "broader_headings": [], "matched_via": "name"}
    md = R.render_all({"15": b15, "9": {}})
    i_ph = md.index("## Pharmacology {#pharmacology}")
    i_mesh = md.index("### MeSH pharmacological classification {#mesh-pharmacology}")
    i_pgx = md.index("### Pharmacogenomics {#pharmacogenomics}")
    assert i_ph < i_mesh < i_pgx
    # no match → no H3 at all
    assert "mesh-pharmacology" not in R.render_all({"15": {"mesh_id": None}, "9": {}})


# ───── PharmGKB drug labels + pathways (§9) ───────────────────────────────────

LABELS = [
    {"label_id": "PA1", "name": "Annotation of HCSC Label for x and CYP2D6", "source": "HCSC",
     "testing_level": "Actionable PGx", "genes": ["CYP2D6"], "has_prescribing_info": True},
    {"label_id": "PA2", "name": "Annotation of FDA Label for x and HLA-B", "source": "FDA",
     "testing_level": "Testing Required", "genes": ["HLA-B"], "has_prescribing_info": True,
     "has_alternate_drug": True, "variants": ["HLA-B*57:01"]},
    {"label_id": "PA3", "name": "Annotation of EMA Label for x", "source": "EMA"},   # empty → dropped
    {"label_id": "PA4", "name": "Annotation of EMA Label for x and CYP3A4", "source": "EMA",
     "genes": ["CYP3A4"]},                                                           # no level
    {"label_id": "PA5", "name": "Annotation of FDA Label for x and CYP2D6", "source": "FDA",
     "testing_level": "Actionable PGx", "genes": ["CYP2D6"], "has_dosing_info": True},
]
PATHWAYS = [
    {"id": "PA145011114", "name": "Warfarin Pathway, Pharmacodynamics",
     "is_pharmacokinetic": "false", "is_pharmacodynamic": "true"},
    {"id": "PA145011113", "name": "Warfarin Pathway, Pharmacokinetics",
     "is_pharmacokinetic": "true", "is_pharmacodynamic": "false"},
    {"id": "PA0", "name": ""},                                                       # nameless → dropped
]


@pytest.fixture
def pgkb(monkeypatch):
    def entry(i, source):
        assert (i, source) == ("PA451906", "pharmgkb")
        return {"Attributes": {"Pharmgkb": {"drug_labels": LABELS}}}

    def map_all(i, chain, cap=60):
        assert (i, chain) == ("PA451906", ">>pharmgkb>>pharmgkb_pathway")
        return PATHWAYS

    monkeypatch.setattr(P, "entry", entry)
    monkeypatch.setattr(P, "map_all", map_all)


def test_drug_labels_collector(pgkb):
    labels = P._drug_labels("PA451906")
    assert [l["label_id"] for l in labels] == ["PA1", "PA2", "PA4", "PA5"]   # PA3 dropped
    pa2 = labels[1]
    assert pa2 == {"label_id": "PA2", "source": "FDA", "level": "Testing Required",
                   "genes": ["HLA-B"], "has_prescribing_info": True,
                   "has_dosing_info": False, "has_alternate_drug": True}
    assert labels[2]["level"] == ""
    assert P._drug_labels(None) == []


def test_pathways_collector(pgkb):
    pws = P._pathways("PA451906")
    assert [p["id"] for p in pws] == ["PA145011114", "PA145011113"]          # name-sorted
    assert pws[0]["pd"] and not pws[0]["pk"]
    assert pws[1]["pk"] and not pws[1]["pd"]
    assert P._pathways(None) == []


def test_render_labels_and_pathways(pgkb):
    b = {"pharmgkb_chemical_id": "PA451906",
         "drug_labels": P._drug_labels("PA451906"), "pathways": P._pathways("PA451906")}
    md = R.r_pharmacogenomics(b)
    rows = [ln for ln in md.splitlines() if ln.startswith("| ") and "---" not in ln][1:]
    # Testing Required first; within a level FDA before HCSC; no-level last
    assert [r.split(" | ")[0].lstrip("| ") for r in rows] == ["FDA", "FDA", "HCSC", "EMA"]
    assert rows[0].split(" | ")[1] == "Testing Required"
    assert "prescribing, alternate drug" in rows[0]
    assert "[PA2](https://www.pharmgkb.org/labelAnnotation/PA2)" in rows[0]
    assert rows[3].split(" | ")[1] == "—"
    assert "4 label annotations:" in md
    assert ("- [Warfarin Pathway, Pharmacokinetics](https://www.pharmgkb.org/pathway/"
            "PA145011113) — pharmacokinetic") in md
    assert "**PharmGKB pathways**" in md
    # labels/pathways are PGx data → the "nothing in PharmGKB" line must not show
    assert "No CPIC/DPWG dosing guideline" not in md


def test_render_without_labels_or_pathways_is_unchanged():
    md = R.r_pharmacogenomics({"pharmgkb_chemical_id": "PA1", "drug_labels": [], "pathways": []})
    assert "PGx drug labels" not in md and "PharmGKB pathway" not in md
    assert "No CPIC/DPWG dosing guideline" in md
    assert R.r_pharmacogenomics({}) == ("## Pharmacogenomics\n\n"
                                        "*No PharmGKB pharmacogenomic data curated for this drug.*")


def test_in_progress_pathways_dropped(monkeypatch):
    from atlas.drug.sections import s09_pharmacogenomics as S9
    monkeypatch.setattr(S9, "map_all", lambda r, c, **k: [
        {"id": "PA1", "name": "In Progress: Drug X Pathway"},
        {"id": "PA2", "name": "Warfarin Pathway, Pharmacodynamics", "is_pharmacodynamic": "true"}])
    assert [p["id"] for p in S9._pathways("PA449082")] == ["PA2"]
