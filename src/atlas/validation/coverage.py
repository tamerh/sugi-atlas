#!/usr/bin/env python3
"""Coverage check: does the deterministic collector reproduce the existing
per-section output across the committed genes?

Usage:  python3 coverage.py [section]   (section default "1")

Two axes per gene:
  (A) CALL coverage — every data SOURCE/target the model used in this section,
      vs what the collector's canonical plan reaches. Flags sources the
      collector does NOT cover (potential data gap); pure syntax/path dupes
      that resolve to the same fact are not flagged (target dataset matches).
  (B) FACT coverage — IDs in the historical <prefix>.md that are MISSING from
      the collector's bundle.
"""
import json, re, os, sys, collections
from atlas.gene import collect as C

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
BASE = os.path.join(REPO, "data/validation/gene")

# per section: file prefix, collector fn, datasets that yield a fact, id regex
SECTIONS = {
    "1": dict(prefix="1_gene_ids", fn=C.collect_gene_ids,
              covered={"hgnc", "ensembl", "entrez", "mim", "omim"},
              ids=r"HGNC:\d+|ENSG\d+|ENST\d+|ENSP\d+"),
    "2": dict(prefix="2_transcripts", fn=C.collect_transcripts,
              covered={"ensembl", "transcript", "refseq", "ccds", "exon", "hgnc"},
              ids=r"ENST\d+|ENSE\d+|NM_\d+|NR_\d+|NP_\d+|CCDS\d+|ENSG\d+"),
    "3": dict(prefix="3_protein_ids", fn=C.collect_protein_ids,
              covered={"uniprot", "interpro", "pfam", "refseq", "antibody",
                       "ufeature", "hgnc", "ensembl"},
              ids=r"IPR\d+|PF\d{5}|NP_\d+|[OPQ][0-9][A-Z0-9]{3}[0-9]"),
    "4": dict(prefix="4_structure", fn=C.collect_structure,
              covered={"pdb", "alphafold", "uniprot", "hgnc"},
              ids=r"AF-[A-Z0-9]+-F\d+|\b[1-9][A-Z][A-Z0-9]{2}\b"),
    "5": dict(prefix="5_orthologs", fn=C.collect_orthologs,
              covered={"ortholog", "paralog", "ensembl", "hgnc"},
              ids=r"ENS[A-Z]{0,4}G\d+|FBgn\d+|WBGene\d+"),
    "6": dict(prefix="6_variants", fn=C.collect_variants,
              # mondo intentionally NOT covered (we skip clinvar 'condition')
              covered={"clinvar", "spliceai", "alphamissense", "dbsnp", "entrez",
                       "hgnc", "ensembl", "transcript", "refseq"},
              ids=r"\d+:\d+:[ACGTN-]+:[ACGTN-]+"),
    "7": dict(prefix="7_pathways", fn=C.collect_pathways,
              covered={"reactome", "msigdb", "go", "uniprot", "hgnc", "ensembl"},
              ids=r"R-HSA-\d+|GO:\d+|M\d{4,}"),
    "8": dict(prefix="8_interactions", fn=C.collect_interactions,
              covered={"string", "string_interaction", "intact", "biogrid",
                       "biogrid_interaction", "signor", "esm2_similarity",
                       "diamond_similarity", "uniprot", "hgnc", "ensembl"},
              ids=r"SIGNOR-\d+|EBI-\d+|[OPQ][0-9][A-Z0-9]{3}[0-9]"),
    "9": dict(prefix="9_tf_regulation", fn=C.collect_tf_regulation,
              covered={"collectri", "jaspar", "uniprot", "hgnc", "ensembl"},
              ids=r"MA\d+\.\d+"),
    "10": dict(prefix="10_drugs", fn=C.collect_drugs,
               covered={"chembl_target", "chembl_molecule", "pharmgkb_gene",
                        "pharmgkb", "bindingdb", "clinical_trials", "mondo",
                        "gencc", "clinvar", "uniprot", "hgnc", "ensembl"},
               ids=r"CHEMBL\d+|NCT\d+"),
    "11": dict(prefix="11_expression", fn=C.collect_expression,
               # hpa/gtex/expressionatlas/proteomicsdb/tabula_sapiens are NOT edges
               # (derived xref-ids, Attributes:Empty) — no integrated data to wire.
               covered={"bgee", "bgee_evidence", "scxa", "fantom5_gene",
                        "fantom5_promoter", "uniprot", "hgnc", "ensembl"},
               ids=r"E-\w+-\d+|UBERON:\d+"),
    "12": dict(prefix="12_diseases", fn=C.collect_diseases,
               covered={"mim", "gencc", "mondo", "orphanet", "hpo", "gwas",
                        "clinvar", "hgnc", "entrez"},
               ids=r"MONDO:\d+|HP:\d+|Orphanet:\d+|GCST\d+|MIM:\d+"),
}

def hist_templates(calls):
    out = set()
    for c in calls:
        tool = c.get("tool", "").replace("biobtree_", "")
        a = c.get("args", {})
        if tool == "search":
            out.add(("search", a.get("dataset", "(any)")))
        elif tool == "entry":
            out.add(("entry", a.get("dataset", "?")))
        elif tool == "map":
            chain = a.get("chain", "").replace("&gt;", ">").replace("&lt;", "<")
            tgt = re.split(r">+|::", chain).pop()
            tgt = re.sub(r"\[.*", "", tgt)  # drop filter bracket: clinvar[...] -> clinvar
            out.add(("map->", tgt or "?"))
    return out

def main():
    sec = sys.argv[1] if len(sys.argv) > 1 else "1"
    cfg = SECTIONS[sec]
    genes = sorted(d for d in os.listdir(BASE) if os.path.isdir(os.path.join(BASE, d)))
    all_uncovered = collections.Counter()
    fact_gaps = {}
    ran = 0
    print(f"section {sec} ({cfg['prefix']})\n")
    print(f"{'gene':<10} {'ok':<3} {'callcov':<8} {'factcov':<9} detail")
    print("-" * 90)
    for g in genes:
        cf = os.path.join(BASE, g, f"{cfg['prefix']}_calls.json")
        mf = os.path.join(BASE, g, f"{cfg['prefix']}.md")
        if not os.path.exists(cf):
            continue
        C.CALLS = []
        try:
            bundle = cfg["fn"](g)
        except SystemExit as e:
            print(f"{g:<10} N   collector failed: {e}"); continue
        except Exception as e:
            print(f"{g:<10} N   error: {e}"); continue
        ran += 1

        uncovered = {f"{tool}{ds}" for tool, ds in hist_templates(json.load(open(cf)))
                     if tool != "search" and ds not in cfg["covered"]}
        for u in uncovered:
            all_uncovered[u] += 1

        missing_ids = []
        if os.path.exists(mf):
            md = open(mf).read()
            bstr = json.dumps(bundle)
            ids = set(re.findall(cfg["ids"], md))
            missing_ids = sorted(i for i in ids if i not in bstr)
            if missing_ids:
                fact_gaps[g] = missing_ids

        callcov = "FULL" if not uncovered else f"{len(uncovered)} gap"
        factcov = "FULL" if not missing_ids else f"{len(missing_ids)} miss"
        detail = " | ".join(filter(None, [
            ("src:" + ",".join(sorted(uncovered))) if uncovered else "",
            ("ids:" + ",".join(missing_ids[:6])) if missing_ids else ""]))
        print(f"{g:<10} Y   {callcov:<8} {factcov:<9} {detail}")

    print("\n" + "=" * 60 + f"\nran on {ran} genes (section {sec})")
    print("\nUNCOVERED sources across genes (genes using each):")
    if all_uncovered:
        for src, n in all_uncovered.most_common():
            print(f"  {n:>2}  {src}")
    else:
        print("  (none — plan covers every data source for this section)")
    print(f"\ngenes with missing IDs in bundle: {len(fact_gaps)}")
    for g, ids in fact_gaps.items():
        print(f"  {g}: {ids[:8]}")

if __name__ == "__main__":
    main()
