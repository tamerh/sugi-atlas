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
)
DATASETS = ("ensembl", "transcript", "refseq", "ccds", "exon")

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
    return bundle

SECTION = Section(
    id="2", name="transcripts",
    description="Ensembl transcripts, RefSeq mRNAs, CCDS, MANE-Select, canonical-transcript exons",
    needs=("ensembl_id", "canonical_transcript", "symbol"),
    produces=("ensembl_transcripts", "refseq_mrna", "ccds", "mane_select_refseq",
              "canonical_transcript", "canonical_exons"),
    datasets=DATASETS, chains=CHAINS, collect_fn=collect,
)
