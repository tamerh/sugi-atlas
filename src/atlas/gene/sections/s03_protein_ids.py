"""§3 — protein_ids: UniProt + RefSeq proteins, InterPro/Pfam domains,
antibodies, UniProt sequence features (with per-accession species filter)."""
from collections import Counter
from atlas.biobtree import map_all
from atlas.gene.sections.base import Section

CHAINS = (
    ">>ensembl>>uniprot",
    '>>ensembl>>refseq[type=="protein"]',
    ">>uniprot>>interpro",
    ">>uniprot>>pfam",
    ">>uniprot>>antibody",
    ">>uniprot>>ufeature",
)
DATASETS = ("uniprot", "ensembl", "refseq", "interpro", "pfam", "antibody", "ufeature")

def collect(a):
    bundle = {
        "section": "03_protein_ids", "symbol": a.symbol,
        "hgnc_id": a.hgnc_id,
        "reviewed_uniprot": list(a.reviewed_uniprots),
        "canonical_uniprot": a.canonical_uniprot,
        "ensembl_id": a.ensembl_id,
    }

    allu = map_all(a.ensembl_id, ">>ensembl>>uniprot") if a.ensembl_id else []
    bundle["uniprot_all"] = [t["id"] for t in allu]
    bundle["uniprot_count"] = len(allu)

    if a.ensembl_id:
        nps = map_all(a.ensembl_id, '>>ensembl>>refseq[type=="protein"]')
        bundle["refseq_protein"] = [{"id": t["id"], "mane": t.get("is_mane_select") == "true"}
                                    for t in nps]
        bundle["refseq_protein_count"] = len(nps)

    # domains/families + antibody + UniProt sequence features, unioned across
    # all reviewed products. (ufeature was previously leaking ortholog features
    # — fixed upstream, BIOBTREE_ISSUES.md #11 RESOLVED — so no post-filter
    # needed anymore.)
    interpro, pfam, antibody = {}, set(), 0
    ufeatures = []
    for u in a.reviewed_uniprots:
        for t in map_all(u, ">>uniprot>>interpro"):
            interpro[t["id"]] = {"id": t["id"], "name": t.get("short_name"),
                                 "type": t.get("type")}
        pfam.update(t["id"] for t in map_all(u, ">>uniprot>>pfam"))
        antibody += len(map_all(u, ">>uniprot>>antibody"))
        for t in map_all(u, ">>uniprot>>ufeature", cap=100):
            ufeatures.append({"uniprot": u, "id": t["id"], "type": t.get("type"),
                              "description": t.get("description"),
                              "begin": t.get("location_begin"),
                              "end": t.get("location_end")})
    bundle["interpro"] = list(interpro.values())
    bundle["pfam"] = sorted(pfam)
    bundle["antibody_count"] = antibody
    bundle["ufeature_counts"] = dict(Counter(f["type"] for f in ufeatures))
    bundle["ufeatures"] = ufeatures
    return bundle

SECTION = Section(
    id="3", name="protein_ids",
    description="UniProt accessions (canonical+all), RefSeq proteins, InterPro/Pfam domains, antibodies, UniProt features",
    needs=("hgnc_id", "ensembl_id", "reviewed_uniprots", "canonical_uniprot"),
    produces=("reviewed_uniprot", "uniprot_all", "refseq_protein", "interpro",
              "pfam", "antibody_count", "ufeatures", "ufeature_counts"),
    datasets=DATASETS, chains=CHAINS, collect_fn=collect,
)
