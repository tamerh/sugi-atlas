"""§8 — interactions: PPIs (STRING/IntAct/BioGRID via interaction records,
score-bearing), SIGNOR signaling, ESM2/Diamond similarity."""
from atlas.biobtree import map_all
from atlas.gene.sections.base import Section

CHAINS = (
    ">>uniprot>>string_interaction",
    ">>uniprot>>intact",
    ">>uniprot>>biogrid_interaction",
    ">>uniprot>>signor",
    ">>uniprot>>esm2_similarity>>uniprot",
    ">>uniprot>>diamond_similarity>>uniprot",
)
DATASETS = ("uniprot", "string_interaction", "intact", "biogrid_interaction",
            "signor", "esm2_similarity", "diamond_similarity")

def _f(x):
    try: return float(x)
    except (TypeError, ValueError): return 0.0

def collect(a):
    bundle = {"section": "08_interactions", "symbol": a.symbol}
    uni = a.canonical_uniprot

    # PPIs via the interaction RECORDS (NOT >>...>>uniprot, which collapses to
    # bare partner ids) — these carry the per-edge scores/evidence.
    st = map_all(uni, ">>uniprot>>string_interaction") if uni else []
    st.sort(key=lambda t: _f(t.get("score")), reverse=True)
    bundle["string"] = [{"partner": t.get("uniprot_b"), "score": t.get("score")} for t in st[:30]]
    bundle["string_count"] = len(st)

    ia = map_all(uni, ">>uniprot>>intact") if uni else []
    ia.sort(key=lambda t: _f(t.get("confidence_score")), reverse=True)
    bundle["intact"] = [{"a": t.get("protein_a_gene"), "b": t.get("protein_b_gene"),
                         "type": t.get("interaction_type"),
                         "score": t.get("confidence_score")} for t in ia[:30]]
    bundle["intact_count"] = len(ia)

    bg = map_all(uni, ">>uniprot>>biogrid_interaction") if uni else []
    bundle["biogrid"] = [{"partner": t.get("interactor_b_symbol"),
                          "method": t.get("experimental_system")} for t in bg[:30]]
    bundle["biogrid_count"] = len(bg)

    def partners(ds):
        return [t["id"] for t in (map_all(uni, f">>uniprot>>{ds}>>uniprot") if uni else [])]
    bundle["esm2_similar"] = partners("esm2_similarity")
    bundle["diamond_similar"] = partners("diamond_similarity")

    sig = map_all(uni, ">>uniprot>>signor") if uni else []
    bundle["signor"] = [{"a": t.get("entity_a"), "b": t.get("entity_b"),
                         "effect": t.get("effect"), "mechanism": t.get("mechanism")} for t in sig]
    bundle["signor_count"] = len(sig)
    return bundle

SECTION = Section(
    id="8", name="interactions",
    description="Protein-protein interactions (STRING/IntAct/BioGRID with scores), SIGNOR signaling, ESM2/Diamond structural similarity",
    needs=("canonical_uniprot",),
    produces=("string", "intact", "biogrid", "signor", "esm2_similar", "diamond_similar"),
    datasets=DATASETS, chains=CHAINS, collect_fn=collect,
)
