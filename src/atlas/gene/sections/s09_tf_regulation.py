"""§9 — regulation: CollecTRI downstream targets + upstream TFs (direction-
filtered), JASPAR motifs (transcriptional), and miRDB miRNAs targeting the
gene (post-transcriptional). is_tf inferred from CollecTRI/JASPAR presence."""
from atlas.biobtree import map_all
from atlas.gene.sections.base import Section

CHAINS = (
    '>>hgnc>>collectri[tf_gene=="<symbol>"]',
    '>>hgnc>>collectri[target_gene=="<symbol>"]',
    ">>uniprot>>jaspar",
    ">>uniprot>>jaspar>>pubmed",
    ">>hgnc>>refseq>>mirdb",
)
DATASETS = ("collectri", "jaspar", "mirdb", "hgnc", "uniprot", "refseq", "pubmed")

def _f(x):
    try: return float(x)
    except (TypeError, ValueError): return 0.0

def collect(a):
    bundle = {"section": "09_tf_regulation", "symbol": a.symbol}
    # direction-filtered + fully paginated -> high-degree TFs (TP53: 1207 CollecTRI
    # rows) get complete downstream targets, not a 100-row truncation.
    down = map_all(a.hgnc_id, f'>>hgnc>>collectri[tf_gene=="{a.symbol}"]')
    up   = map_all(a.hgnc_id, f'>>hgnc>>collectri[target_gene=="{a.symbol}"]')
    bundle["downstream_targets"] = [{"target": r.get("target_gene"),
                                     "regulation": r.get("regulation")} for r in down]
    bundle["downstream_count"] = len(down)
    bundle["upstream_regulators"] = [{"regulator": r.get("tf_gene"),
                                      "regulation": r.get("regulation")} for r in up]
    bundle["jaspar_motifs"] = [{"id": t["id"], "name": t.get("name"),
                                "class": t.get("class"), "family": t.get("family")}
                               for t in (map_all(a.canonical_uniprot, ">>uniprot>>jaspar")
                                         if a.canonical_uniprot else [])]
    # JASPAR PMIDs — evidence trail for the motifs above. One or two PMIDs
    # per TF gene; light footnote ("see PubMed:NNNN") rather than a table.
    bundle["jaspar_pmids"] = [t["id"] for t in
                              (map_all(a.canonical_uniprot, ">>uniprot>>jaspar>>pubmed")
                               if a.canonical_uniprot else [])
                              if t.get("id")]
    bundle["is_transcription_factor"] = bool(down or bundle["jaspar_motifs"])

    # miRDB — miRNAs that target this gene (post-transcriptional regulators).
    # Sort by max_score (miRDB confidence, 0-100); show top-N. target_count
    # is the miRNA's promiscuity (how many genes total it targets) — lower
    # = more specific, surfaced so the reader can judge.
    mirs = map_all(a.hgnc_id, ">>hgnc>>refseq>>mirdb")
    mirs.sort(key=lambda t: _f(t.get("max_score")), reverse=True)
    bundle["mirna_regulators"] = [{
        "id": t["id"],
        "max_score": t.get("max_score"),
        "avg_score": t.get("avg_score"),
        "target_count": t.get("target_count"),
    } for t in mirs[:30]]
    bundle["mirna_count"] = len(mirs)
    return bundle

SECTION = Section(
    id="9", name="regulation",
    description=("Regulators of the gene — transcriptional (CollecTRI upstream TFs + "
                 "JASPAR motifs) and post-transcriptional (miRDB miRNAs targeting "
                 "the gene), plus the CollecTRI downstream-target fan-out when this "
                 "gene is itself a TF."),
    needs=("hgnc_id", "canonical_uniprot", "symbol"),
    produces=("downstream_targets", "upstream_regulators", "jaspar_motifs",
              "jaspar_pmids", "mirna_regulators", "is_transcription_factor"),
    datasets=DATASETS, chains=CHAINS, collect_fn=collect,
)
