"""§8 — protein interactions among cohort genes (STRING / IntAct / BioGRID).
REUSE wrapper over gene §8 + inter-cohort edge graph.

Key value-add: surface UNDRUGGED cohort genes that interact with DRUGGED ones
(indirect druggability path — §17 feeds off this)."""
from atlas.section import Section

CHAINS   = (">>uniprot>>string_interaction", ">>uniprot>>intact",
            ">>uniprot>>biogrid_interaction")
DATASETS = ("uniprot", "string_interaction", "intact", "biogrid_interaction")

def collect(a):
    bundle = {"section": "08_protein_interactions", "mondo_id": a.mondo_id,
              "_todo": "Implement: fan gene §8 over cohort + cross-cohort edges"}
    return bundle

SECTION = Section(
    id="8", name="protein_interactions",
    description=("Cohort-internal interaction graph (STRING / IntAct / BioGRID): "
                 "hub genes, undrugged-↔-drugged bridges for indirect "
                 "druggability."),
    needs=("cohort",),
    produces=("per_gene_interactions", "cohort_edges", "hub_genes",
              "undrugged_to_drugged_bridges"),
    datasets=DATASETS, chains=CHAINS, collect_fn=collect,
)
