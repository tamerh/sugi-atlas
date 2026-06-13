"""§5 — orthologs: cross-species orthologs (Ensembl Compara) + paralogs, plus
UniProt-wide cross-species homologs (ESM2/Diamond similarity) that reach species
beyond Compara's model-organism set."""
from atlas.biobtree import map_all
from atlas.gene.sections.base import Section

CHAINS = (">>ensembl>>ortholog", ">>ensembl>>paralog",
          ">>uniprot>>esm2_similarity", ">>uniprot>>diamond_similarity",
          ">>uniprot>>taxonomy")
DATASETS = ("ensembl", "ortholog", "paralog",
            "uniprot", "esm2_similarity", "diamond_similarity", "taxonomy")

_HOMOLOG_SPECIES_CAP = 15   # distinct non-model species to surface
_HOMOLOG_PROBE_CAP = 40     # bound the per-accession taxonomy resolutions

# biobtree's >>ensembl>>ortholog edge leaves name/genome EMPTY for some model-
# organism-database namespaces (WormBase C. elegans confirmed; FlyBase/zebrafish/
# mouse/rat come back populated). Recover the organism from the id prefix so the
# row isn't half-blank; the symbol isn't resolvable via biobtree, so it stays
# empty (the id remains in its own column). See docs/internal/BIOBTREE_ISSUES.md.
_NS_ORGANISM = (("WBGene", "caenorhabditis_elegans"),
                ("FBgn", "drosophila_melanogaster"))


def _organism_from_id(gid):
    g = (gid or "").upper()
    return next((org for pre, org in _NS_ORGANISM if g.startswith(pre.upper())), "")


def collect(a):
    bundle = {"section": "05_orthologs", "symbol": a.symbol, "ensembl_id": a.ensembl_id}
    orths = map_all(a.ensembl_id, ">>ensembl>>ortholog") if a.ensembl_id else []
    bundle["orthologs"] = [{"id": t["id"], "symbol": t.get("name"),
                            "organism": t.get("genome") or _organism_from_id(t["id"])}
                           for t in orths]
    bundle["ortholog_count"] = len(orths)
    paras = map_all(a.ensembl_id, ">>ensembl>>paralog") if a.ensembl_id else []
    bundle["paralogs"] = [{"id": t["id"], "symbol": t.get("name")} for t in paras]
    bundle["paralog_count"] = len(paras)

    # Cross-species homologs (ESM2/Diamond similarity) — UniProt-wide, so they
    # reach species the Ensembl Compara orthologs above (model organisms only)
    # miss. The similar protein's gene symbol / Ensembl id aren't in biobtree for
    # non-model species, so we present organism (via >>uniprot>>taxonomy) +
    # UniProt accession + similarity score, deduped to one row per species. Human
    # hits are skipped (paralogs already cover them).
    # Species already covered by the Compara orthologs table above — so the
    # homolog block extends BEYOND the model set, not duplicating it. Compara
    # genomes are lower_snake ("mus_musculus"); taxonomy names are "Mus musculus".
    def _norm_sp(s):
        return (s or "").replace("_", " ").strip().lower()
    ortholog_species = {_norm_sp(o.get("organism")) for o in bundle["orthologs"]}

    homologs = {}
    uni = a.canonical_uniprot
    if uni:
        sims = []
        for ds, src in (("esm2_similarity", "ESM2"), ("diamond_similarity", "Diamond")):
            for t in map_all(uni, f">>uniprot>>{ds}", cap=1):
                try:
                    sc = float(t.get("top_similarity") or 0)
                except (TypeError, ValueError):
                    sc = 0.0
                if t.get("id") and t["id"] != uni:
                    sims.append((sc, t["id"], src))
        sims.sort(key=lambda x: -x[0])
        probed = 0
        for sc, acc, src in sims:
            if len(homologs) >= _HOMOLOG_SPECIES_CAP or probed >= _HOMOLOG_PROBE_CAP:
                break
            probed += 1
            tax = map_all(acc, ">>uniprot>>taxonomy")
            org = (tax[0].get("name") if tax else "") or ""
            if not org or org == "Homo sapiens":      # human paralogs already in `paralogs`
                continue
            if _norm_sp(org) in ortholog_species:     # already in the Compara table above
                continue
            if org not in homologs:
                homologs[org] = {"organism": org, "accession": acc,
                                 "similarity": round(sc, 3), "source": src}
    bundle["cross_species_homologs"] = sorted(homologs.values(),
                                              key=lambda h: -h["similarity"])
    return bundle

SECTION = Section(
    id="5", name="orthologs",
    description="Orthologous genes in model organisms + paralogs (Ensembl Compara)",
    needs=("ensembl_id", "canonical_uniprot"),
    produces=("orthologs", "paralogs", "cross_species_homologs"),
    datasets=DATASETS, chains=CHAINS, collect_fn=collect,
)
