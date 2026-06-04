"""§4 — Mendelian disease overlap: which cohort genes also cause Mendelian
forms (via OMIM / Orphanet / GenCC). For cancers: which genes also carry
somatic-driver evidence (CIViC / intogen at gene level).

Genes with BOTH GWAS + Mendelian evidence = highest-confidence targets."""
from atlas.section import Section
from atlas.biobtree import map_all
from atlas.disease.cohort import filter_evidence, disease_tokens

CHAINS   = (">>mondo>>orphanet", ">>orphanet>>hgnc", ">>orphanet>>mim",
            ">>mondo>>gencc>>hgnc", ">>hgnc>>gencc", ">>hgnc>>orphanet",
            ">>hgnc>>mim", ">>hgnc>>intogen", ">>hgnc>>civic")
DATASETS = ("mondo", "orphanet", "mim", "gencc", "hgnc", "intogen", "civic")


def collect(a):
    # disease_name lets the render prefer each cohort gene's GenCC record FOR
    # THIS disease over a stronger off-disease one (audit follow-up: BRCA2's
    # Fanconi-D1 record was winning on the medulloblastoma page).
    bundle = {"section": "04_mendelian_overlap", "mondo_id": a.mondo_id,
              "disease_name": a.canonical_name}

    # GenCC: only fan over genes already flagged gencc=True in cohort_evidence.
    # Typical subset is small (<=10 genes); avoids 50-gene fanout.
    gencc_genes = []
    for ga in filter_evidence(a.cohort, a.cohort_evidence, "gencc"):
        for t in map_all(ga.hgnc_id, ">>hgnc>>gencc"):
            gencc_genes.append({
                "symbol": ga.symbol,
                "hgnc_id": ga.hgnc_id,
                "gencc_classification": t.get("classification_title"),
                "mode_of_inheritance": t.get("moi_title"),
                "mondo_disease": t.get("disease_title"),
            })
    bundle["gencc_genes"] = gencc_genes

    # Orphanet: any cohort gene with a >>hgnc>>orphanet hit. Bounded by cohort.
    orphanet_genes = []
    orphanet_hit_symbols = set()
    for ga in a.cohort:
        rows = map_all(ga.hgnc_id, ">>hgnc>>orphanet")
        for t in rows:
            orphanet_genes.append({
                "symbol": ga.symbol,
                "hgnc_id": ga.hgnc_id,
                "orphanet_id": f"Orphanet:{t['id']}" if t.get("id") else None,
                "orphanet_name": t.get("name"),
            })
        if rows:
            orphanet_hit_symbols.add(ga.symbol)
    bundle["orphanet_genes"] = orphanet_genes

    # OMIM Mendelian overlap: intersect per-gene MIM hits with disease's
    # omim_ids — that's the strong-overlap signal (same MIM on both sides).
    disease_mims = set(a.omim_ids or ())
    omim_genes = []
    for ga in a.cohort:
        gene_mims = {t["id"] for t in map_all(ga.hgnc_id, ">>hgnc>>mim") if t.get("id")}
        overlap = gene_mims & disease_mims
        if overlap:
            omim_genes.append({
                "symbol": ga.symbol,
                "hgnc_id": ga.hgnc_id,
                "shared_mim_ids": sorted(f"MIM:{m}" for m in overlap),
            })
    bundle["omim_genes"] = omim_genes

    # Dual evidence: GWAS + ON-DISEASE Mendelian. The GenCC/Orphanet hit must
    # name THIS disease — previously any >>hgnc>>gencc / >>hgnc>>orphanet record
    # counted, so a GWAS-hit gene that also causes an unrelated Mendelian disease
    # (ATXN1 → ataxia, TBX5 → Holt-Oram) was mislabelled a breast-cancer
    # "highest-confidence target". ev.gencc/ev.clinvar already come from the
    # disease's own evidence routes (on-disease by construction); the name-match
    # restricts the Orphanet/GenCC additions to the same disease.
    dn = disease_tokens(a.canonical_name)
    ondisease_gencc = {g["symbol"] for g in gencc_genes
                       if dn and (disease_tokens(g.get("mondo_disease")) & dn)}
    ondisease_orphanet = {g["symbol"] for g in orphanet_genes
                          if dn and (disease_tokens(g.get("orphanet_name")) & dn)}
    dual = []
    for ga in a.cohort:
        ev = a.cohort_evidence.get(ga.hgnc_id) or {}
        if not ev.get("gwas"):
            continue
        if (ev.get("gencc") or ev.get("clinvar")
                or ga.symbol in ondisease_gencc
                or ga.symbol in ondisease_orphanet):
            dual.append(ga.symbol)
    bundle["dual_evidence_genes"] = dual

    # Somatic drivers: cancer only. Fan intogen + civic over full cohort.
    somatic = []
    if a.is_cancer:
        for ga in a.cohort:
            intogen_rows = map_all(ga.hgnc_id, ">>hgnc>>intogen")
            civic_rows = map_all(ga.hgnc_id, ">>hgnc>>civic")
            if not intogen_rows and not civic_rows:
                continue
            intogen = None
            if intogen_rows:
                r = intogen_rows[0]
                intogen = {"role": r.get("role"),
                           "cancer_types": r.get("cancer_types")}
            civic = None
            if civic_rows:
                r = civic_rows[0]
                civic = {"id": r.get("id"), "name": r.get("name")}
            somatic.append({
                "symbol": ga.symbol,
                "hgnc_id": ga.hgnc_id,
                "intogen": intogen,
                "civic": civic,
            })
    bundle["somatic_driver_genes"] = somatic

    return bundle


SECTION = Section(
    id="4", name="mendelian_overlap",
    description=("Cohort genes with Mendelian evidence (GenCC, Orphanet, OMIM); "
                 "for cancers: cohort genes with somatic-driver evidence "
                 "(intOGen, CIViC). High-confidence target subset."),
    needs=("mondo_id", "cohort", "cohort_evidence", "is_cancer", "orphanet_ids", "omim_ids"),
    produces=("disease_name", "gencc_genes", "orphanet_genes", "omim_genes",
              "somatic_driver_genes", "dual_evidence_genes"),
    datasets=DATASETS, chains=CHAINS, collect_fn=collect,
)
