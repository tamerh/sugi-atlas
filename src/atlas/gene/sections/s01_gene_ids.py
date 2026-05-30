"""§1 — gene_ids: core HGNC identifiers, cross-refs, ensembl biotype + xref counts."""
from atlas.biobtree import entry, map_all, xref_counts
from atlas.gene.sections.base import Section

CHAINS   = (">>hgnc>>mim", ">>hgnc>>entrez")
DATASETS = ("hgnc", "ensembl", "mim", "entrez")

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
    return bundle

SECTION = Section(
    id="1", name="gene_ids",
    description="Core gene identifiers (HGNC, Ensembl, OMIM, Entrez) + xref-count table",
    needs=("hgnc_id", "hgnc_entry", "ensembl_id"),
    produces=("hgnc", "ensembl_id", "mim", "entrez", "ensembl", "xref_counts"),
    datasets=DATASETS, chains=CHAINS, collect_fn=collect,
)
