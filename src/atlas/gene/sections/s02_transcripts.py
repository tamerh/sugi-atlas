"""§2 — transcripts: Ensembl + RefSeq + CCDS, MANE-Select, canonical-transcript exons."""
from atlas.biobtree import map_all
from atlas.gene.sections.base import Section

CHAINS = (
    ">>ensembl>>transcript",
    '>>ensembl>>refseq[type=="mRNA"]',
    ">>ensembl>>ccds",
    ">>ensembl>>refseq[is_mane_select==true]",
    ">>refseq>>transcript",
    ">>transcript>>exon",
    ">>hgnc>>entrez>>neighborentrez",   # genomic neighbors (locus context)
)
DATASETS = ("ensembl", "transcript", "refseq", "ccds", "exon",
            "entrez", "neighborentrez")

# Entrez gene-neighborhood `type`s that are real genes worth showing — the raw
# edge is dominated by `biological-region` (recombination/intergenic regions)
# and `pseudo` noise, which we drop.
_NEIGHBOR_GENE_TYPES = {"protein-coding", "ncRNA", "snRNA", "snoRNA",
                        "miRNA", "rRNA", "tRNA", "scRNA", "lncRNA"}
_NEIGHBOR_CAP = 30

def collect(a):
    bundle = {"section": "02_transcripts", "symbol": a.symbol, "ensembl_id": a.ensembl_id}
    if not a.ensembl_id:
        return bundle

    tr = map_all(a.ensembl_id, ">>ensembl>>transcript")
    bundle["ensembl_transcripts"] = [{"id": t["id"], "biotype": t.get("biotype")} for t in tr]
    bundle["ensembl_transcript_count"] = len(tr)

    rs = map_all(a.ensembl_id, '>>ensembl>>refseq[type=="mRNA"]')
    bundle["refseq_mrna"] = [t["id"] for t in rs]
    bundle["refseq_mrna_count"] = len(rs)

    bundle["ccds"] = [t["id"] for t in map_all(a.ensembl_id, ">>ensembl>>ccds")]

    # MANE-Select via dedicated filter (never missed for high-RefSeq genes)
    mane = map_all(a.ensembl_id, ">>ensembl>>refseq[is_mane_select==true]")
    mane_nm = next((t["id"] for t in mane if t.get("type") == "mRNA"), None)
    bundle["mane_select_refseq"] = mane_nm
    bundle["canonical_transcript"] = a.canonical_transcript or (tr[0]["id"] if tr else None)
    if bundle["canonical_transcript"]:
        ex = map_all(bundle["canonical_transcript"], ">>transcript>>exon")
        bundle["canonical_exons"] = [{"id": e["id"], "start": e.get("start"),
                                      "end": e.get("end")} for e in ex]
        bundle["canonical_exon_count"] = len(ex)

    # Genomic neighbors (NCBI/Entrez gene neighborhood) — flanking genes on the
    # chromosome. Filtered to real genes (drop biological-region/pseudo noise and
    # unnamed LOC ids) and self, deduped, capped. Positional context only.
    seen = set()
    neighbors = []
    if a.hgnc_id:
        for t in map_all(a.hgnc_id, ">>hgnc>>entrez>>neighborentrez"):
            sym = (t.get("symbol") or "").strip()
            if (t.get("type") not in _NEIGHBOR_GENE_TYPES or not sym
                    or sym.startswith("LOC") or sym == a.symbol or sym in seen):
                continue
            seen.add(sym)
            neighbors.append({"id": t["id"], "symbol": sym, "type": t.get("type")})
    bundle["genomic_neighbors"] = neighbors[:_NEIGHBOR_CAP]
    bundle["genomic_neighbor_count"] = len(neighbors)
    return bundle

SECTION = Section(
    id="2", name="transcripts",
    description="Ensembl transcripts, RefSeq mRNAs, CCDS, MANE-Select, canonical-transcript exons",
    needs=("ensembl_id", "canonical_transcript", "symbol", "hgnc_id"),
    produces=("ensembl_transcripts", "refseq_mrna", "ccds", "mane_select_refseq",
              "canonical_transcript", "canonical_exons",
              "genomic_neighbors", "genomic_neighbor_count"),
    datasets=DATASETS, chains=CHAINS, collect_fn=collect,
    # refseq_mrna fluctuates with upstream NCBI's REVIEWED-only filter
    # (2026-05-30 refresh: TP53 mRNA 46→25 because NM_/NR_ kept but XM_/XR_
    # PREDICTED dropped). Body_gate demotes those shrinks from regression
    # to drift. ensembl_transcripts can drop similarly when biotypes are
    # recurated.
    shrinkable=("refseq_mrna", "ensembl_transcripts"),
)
