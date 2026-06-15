"""§5 — orthologs: cross-species orthologs (Ensembl Compara) + paralogs, plus
UniProt-wide cross-species homologs (ESM2/Diamond similarity) that reach species
beyond Compara's model-organism set."""
from atlas.biobtree import map_all
from atlas.gene.sections.base import Section

CHAINS = (">>ensembl>>ortholog", ">>ensembl>>paralog",
          ">>uniprot>>diamond_similarity", ">>uniprot>>taxonomy")
DATASETS = ("ensembl", "ortholog", "paralog",
            "uniprot", "diamond_similarity", "taxonomy")

_HOMOLOG_SPECIES_CAP = 40   # distinct non-model species to surface
_HOMOLOG_PROBE_CAP = 100    # bound the per-accession taxonomy resolutions — set to
                            # the Diamond hit-list max (~100) so we can scan the whole
                            # list to reach SPECIES_CAP distinct species when available

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


def _is_virus(org):
    """A viral taxon (e.g. 'Avian leukosis virus' carrying v-erbB, a transduced
    EGFR). Genuine sequence homologs, but not organismal species — de-prioritised
    below the cellular orthologs so they never displace a mammalian ortholog. No
    cellular organism carries 'virus' in its scientific name, so the test is safe."""
    return "virus" in (org or "").lower()


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

    # Cross-species homologs (Diamond similarity) — UniProt-wide, so they reach
    # species the Ensembl Compara orthologs above (model organisms only) miss. We
    # use Diamond ONLY, not ESM2: Diamond is local sequence alignment, so a hit is
    # a genuine sequence homolog with an interpretable % identity (e.g. EGFR's hits
    # span 34–100% identity and are all ERBB-family). ESM2 embedding-cosine
    # saturates near 1.0 even for unrelated proteins (it scored bovine FKBP9 at
    # 0.9991 vs EGFR — a false positive indistinguishable from the real hits), so
    # it's unsafe to assert as "homolog" on a reference page. See BIOBTREE_ISSUES
    # #44. The homolog's gene symbol / Ensembl id aren't in biobtree for non-model
    # species, so we present organism (via >>uniprot>>taxonomy) + UniProt accession
    # + % identity, deduped to the best hit per species. Human hits are skipped
    # (paralogs already cover them); species already in the Compara table above are
    # skipped too, so this block extends BEYOND the model set. Compara genomes are
    # lower_snake ("mus_musculus"); taxonomy names are "Mus musculus".
    def _norm_sp(s):
        return (s or "").replace("_", " ").strip().lower()
    ortholog_species = {_norm_sp(o.get("organism")) for o in bundle["orthologs"]}

    homologs = {}
    uni = a.canonical_uniprot
    if uni:
        sims = []
        for t in map_all(uni, ">>uniprot>>diamond_similarity", cap=1):
            # Diamond schema: id|similarity_count|top_identity|top_bitscore.
            # top_identity is a percentage (0–100); store as a 0–1 fraction so the
            # renderer's `:.1%` shows "91.8%".
            try:
                ident = float(t.get("top_identity") or 0) / 100.0
            except (TypeError, ValueError):
                ident = 0.0
            if t.get("id") and t["id"] != uni and ident > 0:
                sims.append((ident, t["id"]))
        sims.sort(key=lambda x: -x[0])
        probed = n_cellular = 0
        for ident, acc in sims:
            # Stop once the cap can be filled by cellular organisms alone, or the
            # taxonomy-probe budget is spent — viral homologs only fill leftover
            # slots, so there's no need to keep probing past a full cellular set.
            if probed >= _HOMOLOG_PROBE_CAP or n_cellular >= _HOMOLOG_SPECIES_CAP:
                break
            probed += 1
            tax = map_all(acc, ">>uniprot>>taxonomy")
            org = (tax[0].get("name") if tax else "") or ""
            if not org or org == "Homo sapiens":      # human paralogs already in `paralogs`
                continue
            if _norm_sp(org) in ortholog_species:     # already in the Compara table above
                continue
            if org not in homologs:                   # first = highest identity (sorted)
                homologs[org] = {"organism": org, "accession": acc,
                                 "similarity": round(ident, 3), "source": "Diamond"}
                if not _is_virus(org):
                    n_cellular += 1
    # Cellular organisms first (the cross-species orthologs readers expect), then
    # viral homologs (v-erbB-type oncogene captures), each by descending identity;
    # truncate to the species cap so viruses never displace a mammalian ortholog.
    cellular = sorted((h for h in homologs.values() if not _is_virus(h["organism"])),
                      key=lambda h: -h["similarity"])
    viral = sorted((h for h in homologs.values() if _is_virus(h["organism"])),
                   key=lambda h: -h["similarity"])
    bundle["cross_species_homologs"] = (cellular + viral)[:_HOMOLOG_SPECIES_CAP]
    return bundle

SECTION = Section(
    id="5", name="orthologs",
    description="Orthologous genes in model organisms + paralogs (Ensembl Compara)",
    needs=("ensembl_id", "canonical_uniprot"),
    produces=("orthologs", "paralogs", "cross_species_homologs"),
    datasets=DATASETS, chains=CHAINS, collect_fn=collect,
)
