"""§5 — orthologs: cross-species orthologs (Ensembl Compara) + paralogs, plus
UniProt-wide cross-species homologs (ESM2/Diamond similarity) that reach species
beyond Compara's model-organism set."""
from collections import Counter

from atlas.biobtree import map_all
from atlas.gene.sections.base import Section

_MOUSE_MGI_CHAIN = ">>hgnc>>entrez>>orthologentrez>>mgi"
_MOUSE_PHENO_CHAIN = ">>hgnc>>entrez>>orthologentrez>>mgi>>alliance_phenotype"

CHAINS = (">>ensembl>>ortholog", ">>ensembl>>paralog",
          ">>uniprot>>diamond_similarity", ">>uniprot>>taxonomy",
          _MOUSE_MGI_CHAIN, _MOUSE_PHENO_CHAIN)
DATASETS = ("ensembl", "ortholog", "paralog",
            "uniprot", "diamond_similarity", "taxonomy",
            "entrez", "orthologentrez", "mgi", "alliance_phenotype")

# Observed mouse-model phenotypes (MGI, via the Alliance) for the gene's mouse
# ortholog — real knockout/mutant phenotypes, distinct from the HP→uPheno→MP
# TRANSLATION of the human gene's own HPO terms. The human→ortholog hop
# (orthologentrez) was fixed 2026-06-20, so this is now a single chain query
# (human gene → mouse MGI → alliance_phenotype) — no longer the 4-hop step-through
# via the Ensembl ortholog. The >>mgi hop filters orthologentrez's all-species set
# to mouse.
_MOUSE_PHENO_CAP = 25
_PHENO_NOISE = {"no abnormal phenotype detected"}


def _mouse_phenotypes(hgnc_id):
    """[{mp_id, statement, records}] for the gene's mouse ortholog, ranked by MGI
    record count (more annotations = more robustly observed), plus the MGI id(s).
    Drops the 'no abnormal phenotype detected' control rows."""
    mgi_ids = [m["id"] for m in map_all(hgnc_id, _MOUSE_MGI_CHAIN) if m.get("id")]
    counts: Counter = Counter()
    stmts: dict = {}
    # Cap the phenotype pull — top genes carry 600+ records; ~300 (3 pages) is
    # plenty to rank the top _MOUSE_PHENO_CAP by count.
    for r in map_all(hgnc_id, _MOUSE_PHENO_CHAIN, cap=3):
        term = r.get("phenotype_term")
        stmt = (r.get("phenotype_statement") or "").strip()
        if not term or not stmt or stmt.lower() in _PHENO_NOISE:
            continue
        counts[term] += 1
        stmts.setdefault(term, stmt)
    ranked = sorted(counts.items(), key=lambda kv: (-kv[1], kv[0]))
    phenos = [{"mp_id": t, "statement": stmts[t], "records": n}
              for t, n in ranked[:_MOUSE_PHENO_CAP]]
    return phenos, mgi_ids, len(counts)

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

    # Observed mouse-model phenotypes for the mouse ortholog (see _mouse_phenotypes).
    phenos, mgi_ids, distinct = _mouse_phenotypes(a.hgnc_id) if a.hgnc_id else ([], [], 0)
    bundle["mouse_phenotypes"] = phenos
    bundle["mouse_phenotype_total"] = distinct
    bundle["mouse_mgi_ids"] = mgi_ids

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
    produces=("orthologs", "paralogs", "cross_species_homologs",
              "mouse_phenotypes", "mouse_phenotype_total", "mouse_mgi_ids"),
    datasets=DATASETS, chains=CHAINS, collect_fn=collect,
)
