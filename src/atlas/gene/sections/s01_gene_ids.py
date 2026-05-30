"""§1 — gene_ids: core HGNC identifiers, cross-refs, ensembl biotype + xref counts.

RNAcentral is added for non-coding RNA genes — it's the closest analogue to
UniProt for ncRNAs (canonical id + classification + length + cross-species
coverage). Protein-coding genes return no RNAcentral row, so the field
elides cleanly."""
from atlas.biobtree import entry, map_all, xref_counts
from atlas.gene.sections.base import Section

CHAINS   = (">>hgnc>>mim", ">>hgnc>>entrez", ">>hgnc>>rnacentral")
DATASETS = ("hgnc", "ensembl", "mim", "entrez", "rnacentral")

def collect(a):
    bundle = {"section": "01_gene_ids", "symbol": a.symbol, "hgnc_id": a.hgnc_id}
    attrs = a.hgnc_entry.get("Attributes", {}).get("Hgnc", {})
    bundle["name"] = (attrs.get("names") or [None])[0]
    bundle["hgnc"] = {
        "symbol": (attrs.get("symbols") or [None])[0],
        "name": (attrs.get("names") or [None])[0],
        "location": attrs.get("location"),
        "locus_type": attrs.get("locus_type"),
        "status": attrs.get("status"),
        "aliases": attrs.get("aliases", []),
    }
    bundle["xref_counts"] = xref_counts(a.hgnc_entry)
    bundle["ensembl_id"] = a.ensembl_id
    bundle["mim"] = [t["id"] for t in map_all(a.hgnc_id, ">>hgnc>>mim")]
    bundle["entrez"] = [t["id"] for t in map_all(a.hgnc_id, ">>hgnc>>entrez")]
    if a.ensembl_id:
        ee = entry(a.ensembl_id, "ensembl")
        ea = ee.get("Attributes", {}).get("Ensembl", {})
        bundle["ensembl"] = {"biotype": ea.get("biotype"), "genome": ea.get("genome")}

    # RNAcentral — canonical id + RNA type + length for ncRNA genes. Returns
    # nothing for protein-coding genes (closes the prior thin coverage of
    # lncRNA / miRNA / etc. gene pages).
    rna_rows = map_all(a.hgnc_id, ">>hgnc>>rnacentral")
    if rna_rows:
        r = rna_rows[0]
        bundle["rnacentral"] = {
            "id": r.get("id"),
            "rna_type": r.get("rna_type"),
            "length": r.get("length"),
            "organism_count": r.get("organism_count"),
        }
    else:
        bundle["rnacentral"] = None
    return bundle

SECTION = Section(
    id="1", name="gene_ids",
    description=("Core gene identifiers (HGNC, Ensembl, OMIM, Entrez) + xref-count "
                 "table. For ncRNA genes also the canonical RNAcentral id (rna_type, "
                 "length, organism count)."),
    needs=("hgnc_id", "hgnc_entry", "ensembl_id"),
    produces=("hgnc", "ensembl_id", "mim", "entrez", "ensembl", "xref_counts",
              "rnacentral"),
    datasets=DATASETS, chains=CHAINS, collect_fn=collect,
)
