#!/usr/bin/env python3
"""Deterministic markdown renderer for disease bundles — NO model.

Mirrors atlas.gene.render — every fact verbatim from the collector bundle,
zero tokens, zero hallucination. The LLM is reserved for the executive
summary, not section bodies.

§1–§14 use RENDER dict (section_id -> renderer_fn taking the section bundle).
§15 (drug_repurposing), §16 (druggability_pyramid), §17 (undrugged_target_profiles)
are derived views that join multiple §1–§14 bundles, so they live below the
dict and are called explicitly by render_all().
"""
from collections import Counter
from atlas.render_common import table
from atlas.civic import therapy_label

# Shared formatting helpers --------------------------------------------------

def _i(n):
    """Comma-formatted int, falling back to '' for None/empty."""
    if n in (None, ""):
        return ""
    try:
        return f"{int(n):,}"
    except (TypeError, ValueError):
        return str(n)


def _trunc(s, n=80):
    s = (s or "").strip()
    return s if len(s) <= n else s[: n - 1].rstrip() + "…"


# §1 disease_ids ------------------------------------------------------------

def r_disease_ids(b):
    # Skip rows whose value is empty so the ID table only carries facts.
    # Cleaner than rendering empty `| |` cells for diseases that don't have
    # the given cross-reference (most diseases miss at least one of MeSH /
    # OMIM / Orphanet).
    candidate_rows = [
        ("Canonical name", b.get("canonical_name")),
        # Use monarchinitiative.org's stable Mondo redirect — shorter + more
        # human-readable than the OLS4 percent-encoded URL.
        ("Mondo ID",
         f"[{b['mondo_id']}](https://monarchinitiative.org/{b['mondo_id']})"
         if b.get("mondo_id") else None),
        ("EFO", b.get("efo_id")),
        ("MeSH", ", ".join(b.get("mesh_ids") or []) or None),
        ("OMIM", ", ".join(b.get("omim_ids") or []) or None),
        ("Orphanet", ", ".join(b.get("orphanet_ids") or []) or None),
        # Mondo-OBO cross-ontology xrefs — keys present only if data exists.
        ("DOID",     ", ".join((b.get("obo_xrefs") or {}).get("doid")    or []) or None),
        ("ICD-10-CM", ", ".join((b.get("obo_xrefs") or {}).get("icd10cm") or []) or None),
        ("ICD-11",   ", ".join((b.get("obo_xrefs") or {}).get("icd11")   or []) or None),
        ("NCIT",     ", ".join((b.get("obo_xrefs") or {}).get("ncit")    or []) or None),
        ("SNOMED CT", ", ".join((b.get("obo_xrefs") or {}).get("sctid")  or []) or None),
        ("UMLS",     ", ".join((b.get("obo_xrefs") or {}).get("umls")    or []) or None),
        ("MedGen",   ", ".join((b.get("obo_xrefs") or {}).get("medgen")  or []) or None),
        ("GARD",     ", ".join((b.get("obo_xrefs") or {}).get("gard")    or []) or None),
        ("MedDRA",   ", ".join((b.get("obo_xrefs") or {}).get("meddra")  or []) or None),
        ("NORD",     ", ".join((b.get("obo_xrefs") or {}).get("nord")    or []) or None),
        ("Anatomy (UBERON)",
         ", ".join(b.get("anatomy_uberon_ids") or []) or None),
        ("Is cancer (heuristic)", "yes" if b.get("is_cancer") else "no"),
    ]
    rows = [(k, v) for k, v in candidate_rows if v not in (None, "")]
    out = ["## Disease identifiers", "", table(["Field", "Value"], rows)]

    xc = b.get("xref_counts") or {}
    if xc:
        out.append("")
        out.append("**Cross-database coverage (counts from the Mondo entry):**")
        out.append("")
        ordered = sorted(xc.items(), key=lambda kv: -kv[1])
        out.append(table(["Dataset", "Count"], [(k, _i(v)) for k, v in ordered]))

    # Epidemiology — Orphanet curates multi-geography prevalence data per
    # rare disease. Show as a small table when present (validated rows first).
    prevs = b.get("prevalences") or []
    if prevs:
        prevs_sorted = sorted(prevs,
                              key=lambda p: (0 if p.get("validation_status") == "Validated" else 1,
                                             p.get("prevalence_type") or ""))
        out.append("")
        out.append(f"**Epidemiology ({len(prevs)} prevalence records, "
                   "Orphanet):**")
        out.append("")
        out.append(table(["Type", "Class", "Value", "Geography", "Validation"],
                         [(p.get("prevalence_type"), p.get("prevalence_class"),
                           p.get("val_moy"), p.get("geographic"),
                           p.get("validation_status"))
                          for p in prevs_sorted]))

    # Clinical features (HPO phenotypes from Orphanet, frequency-sorted).
    # Show top 30; full list lives in the bundle/sidecar for RAG consumers.
    phs = b.get("phenotypes") or []
    if phs:
        total = b.get("phenotype_count") or len(phs)
        out.append("")
        out.append(f"**Clinical features ({total} HPO phenotypes, Orphanet "
                   f"curated; top {min(30, len(phs))} by frequency):**")
        out.append("")
        out.append(table(["HPO ID", "Term", "Frequency"],
                         [(f"[{p.get('hpo_id')}](https://hpo.jax.org/app/browse/term/{p.get('hpo_id')})",
                           p.get("hpo_term"),
                           p.get("frequency"))
                          for p in phs[:30]]))
    return "\n".join(out)


# §2 gwas_landscape ---------------------------------------------------------

def r_gwas_landscape(b):
    # No common-variant signal (e.g. Mendelian / rare diseases) → say so
    # concisely instead of a "0 across 0 … 0 distinct genes" line that reads
    # like a bug. The real genetic evidence is in §4 (Mendelian) / §3 (ClinVar).
    if not (b.get("assoc_total") or 0):
        return ("## GWAS landscape\n\n*No GWAS associations recorded — common-"
                "variant (GWAS) studies don't cover this disease (typical for "
                "Mendelian / rare diseases). See the curated gene cohort and "
                "Mendelian overlap below.*")
    out = [f"## GWAS landscape", "",
           f"**{_i(b.get('assoc_total'))} GWAS associations across "
           f"{_i(b.get('study_total'))} studies.** "
           f"Cohort touches {_i(b.get('unique_gene_count'))} distinct genes "
           "across the top hits."]
    ta = b.get("top_assocs") or []
    if ta:
        out += ["", "**Top associations by p-value:**", "",
                table(["rsID", "p-value", "Gene", "Risk allele", "Odds ratio"],
                      [(r.get("rsid"), r.get("pvalue"), r.get("gene_symbol"),
                        r.get("risk_allele"), r.get("odds_ratio"))
                       for r in ta[:30]])]
    studies = b.get("studies") or []
    if studies:
        out += ["", "**Top studies (by case count):**", "",
                table(["Study", "Lead author", "Year", "Cases", "Controls", "Title"],
                      [(f"[{s['id']}](https://www.ebi.ac.uk/gwas/studies/{s['id']})"
                        if s.get("id") else "",
                        s.get("lead_author"), s.get("year"),
                        _i(s.get("sample_size_cases")),
                        _i(s.get("sample_size_controls")),
                        _trunc(s.get("title"), 60))
                       for s in studies])]
    return "\n".join(out)


# §3 variant_details --------------------------------------------------------

def r_variant_details(b):
    out = ["## Variant details & genetic-evidence tiers", ""]
    # This section tiers GWAS-derived variants (>>mondo>>gwas>>dbsnp). For
    # diseases with no GWAS (Mendelian / rare) it's empty — render a note, not
    # all-zero tables. Guard each sub-block on real data (non-zero sum), since a
    # dict of zeros is still truthy.
    tv = b.get("top_variants") or []
    tc = b.get("tier_counts") or {}
    md = b.get("maf_distribution") or {}
    cd = b.get("consequence_distribution") or {}
    if not tv and not any((tc or {}).values()):
        out.append("*No GWAS-derived variants to tier here (this block classifies "
                    "common-variant GWAS hits). Rare/coding variants for this "
                    "disease appear under Mendelian overlap and ClinVar.*")
        return "\n".join(out)
    if any(tc.values()):
        out += ["**Tier distribution (top 50 variants):**", ""]
        rows = sorted(tc.items(), key=lambda kv: kv[0])
        out.append(table(["Tier", "Variants"], [(k, _i(v)) for k, v in rows]))
    if any(md.values()):
        out += ["", "**MAF distribution:**", ""]
        out.append(table(["Bucket", "Variants"], [(k, _i(v)) for k, v in md.items()]))
    if any(cd.values()):
        out += ["", "**Functional consequences:**", ""]
        rows = sorted(cd.items(), key=lambda kv: -kv[1])[:15]
        out.append(table(["Consequence", "Count"], [(k, _i(v)) for k, v in rows]))
    if tv:
        out += ["", "**Top variants:**", "",
                table(["rsID", "Chr", "Pos", "Alleles", "MAF",
                       "Consequence", "Gene", "p-value", "Tier"],
                      [(f"[{r['rsid']}](https://www.ncbi.nlm.nih.gov/snp/{r['rsid']})"
                        if r.get("rsid") else "",
                        r.get("chrom"), r.get("pos"), r.get("alleles"),
                        r.get("maf"), r.get("consequence"),
                        r.get("gene_symbol"), r.get("pvalue"), r.get("tier"))
                       for r in tv[:30]])]
    return "\n".join(out)


# §4 mendelian_overlap ------------------------------------------------------

def r_mendelian_overlap(b):
    out = ["## Mendelian disease overlap & somatic drivers", "",
           f"**GenCC: {len(b.get('gencc_genes') or [])} · "
           f"Orphanet: {len(b.get('orphanet_genes') or [])} · "
           f"OMIM-shared: {len(b.get('omim_genes') or [])} · "
           f"Dual-evidence (GWAS+Mendelian): {len(b.get('dual_evidence_genes') or [])}**"]
    dual = b.get("dual_evidence_genes") or []
    if dual:
        # Enrich the dual-evidence table: HGNC link + evidence-route flags so
        # the reader sees *why* each gene qualifies as dual-evidence.
        # cohort_evidence isn't directly on this section's bundle so we
        # surface evidence labels by walking gencc/orphanet/clinvar lists.
        gencc_set = {g.get("symbol") for g in (b.get("gencc_genes") or [])}
        orph_set  = {g.get("symbol") for g in (b.get("orphanet_genes") or [])}
        omim_set  = {g.get("symbol") for g in (b.get("omim_genes") or [])}
        def _routes(sym):
            r = ["GWAS"]
            if sym in gencc_set: r.append("GenCC")
            if sym in orph_set:  r.append("Orphanet")
            if sym in omim_set:  r.append("OMIM")
            return ", ".join(r)
        out += ["", "**Dual-evidence genes (GWAS + Mendelian — highest-confidence targets):**", "",
                table(["Gene", "HGNC", "Evidence routes"],
                      [(f"[{sym}](https://www.genenames.org/tools/search/#!/?query={sym})",
                        sym, _routes(sym))
                       for sym in dual[:50]])]
    sg = b.get("somatic_driver_genes") or []
    if sg:
        # CIViC's `name` field equals the gene symbol — we surface the CIViC
        # ID as a clickable link to the gene's CIViC summary instead, which
        # is the actual value-add.
        out += ["", "**Somatic driver evidence (intOGen + CIViC, cohort fanout):**", "",
                table(["Gene", "intOGen role", "Cancer types", "CIViC"],
                      [(g.get("symbol"),
                        (g.get("intogen") or {}).get("role"),
                        _trunc((g.get("intogen") or {}).get("cancer_types") or "", 50),
                        (f"[CIViC #{(g.get('civic') or {}).get('id')}]"
                         f"(https://civicdb.org/genes/{(g.get('civic') or {}).get('id')}/summary)"
                         if (g.get("civic") or {}).get("id") else ""))
                       for g in sg[:30]])]
    gc = b.get("gencc_genes") or []
    if gc:
        out += ["", "**GenCC Mendelian classification:**", "",
                table(["Gene", "Classification", "Inheritance", "Disease"],
                      [(g.get("symbol"), g.get("gencc_classification"),
                        g.get("mode_of_inheritance"),
                        _trunc(g.get("mondo_disease"), 50)) for g in gc[:30]])]
    og = b.get("orphanet_genes") or []
    if og:
        out += ["", "**Orphanet rare-disease linkage (cohort genes):**", "",
                table(["Gene", "Orphanet ID", "Rare disease"],
                      [(g.get("symbol"),
                        f"[{g['orphanet_id']}](https://www.orpha.net/en/disease/detail/"
                        f"{g['orphanet_id'].split(':')[-1] if g.get('orphanet_id') else ''})"
                        if g.get("orphanet_id") else "",
                        _trunc(g.get("orphanet_name"), 65))
                       for g in og[:50]])]
    omg = b.get("omim_genes") or []
    if omg:
        out += ["", "**OMIM-shared genes (cohort gene's MIM ids overlap the disease's):**", "",
                table(["Gene", "MIM ids"],
                      [(g.get("symbol"), ", ".join(g.get("mim_ids") or []))
                       for g in omg[:30]])]
    return "\n".join(out)


# §5 genes_proteins ---------------------------------------------------------

def r_genes_proteins(b):
    out = ["## Cohort genes → proteins", "",
           f"**{_i(b.get('gene_count'))} cohort genes, "
           f"{_i(b.get('protein_count'))} distinct canonical proteins.**"]
    ev = b.get("evidence_summary") or {}
    if ev:
        out += ["", "**Evidence partition:**", ""]
        out.append(table(["Subset", "Genes"],
                         [(k, _i(v)) for k, v in ev.items()]))
    genes = b.get("genes") or []
    if genes:
        out += ["", "**Cohort genes (full):**", "",
                table(["Symbol", "HGNC", "Ensembl", "UniProt", "Name", "Evidence"],
                      [(g.get("symbol"), g.get("hgnc_id"), g.get("ensembl_id"),
                        g.get("canonical_uniprot"),
                        _trunc(g.get("protein_name") or g.get("hgnc_name"), 50),
                        ",".join(k for k, v in (g.get("evidence") or {}).items() if v))
                       for g in genes[:50]])]

    # Cohort function summary — one-sentence UniProt FUNCTION per gene.
    # Surfaces "what each cohort gene actually does" at a glance, not just
    # HGNC ids. New 2026-05-31, post BIOBTREE_ISSUES #9 resolution.
    fns = b.get("cohort_function_summary") or []
    if fns:
        out += ["", "**Cohort function summary (lead sentence per gene, "
                "UniProt-curated):**", "",
                table(["Symbol", "Protein name", "Function (lead sentence)"],
                      [(f["symbol"], _trunc(f.get("protein_name"), 40),
                        f.get("function_lead", ""))
                       for f in fns[:50]])]
    return "\n".join(out)


# §6 protein_families -------------------------------------------------------

def r_protein_families(b):
    out = ["## Protein-family classification", "",
           f"**Druggable: {_i(b.get('druggable_count'))} · "
           f"Difficult: {_i(b.get('difficult_count'))} · "
           f"Unknown: {_i(b.get('unknown_count'))} · "
           f"Druggable fraction: {b.get('druggable_fraction')}**"]
    fc = b.get("family_counts") or {}
    if fc:
        out += ["", "**Family distribution:**", ""]
        rows = sorted(fc.items(), key=lambda kv: -kv[1])
        out.append(table(["Family", "Genes"], [(k, _i(v)) for k, v in rows]))
    fa = b.get("family_assignments") or []
    if fa:
        out += ["", "**Per-gene assignment:**", "",
                table(["Symbol", "Family", "Druggable?", "EC", "InterPro (top 3)"],
                      [(g.get("symbol"), g.get("assigned_family"),
                        "yes" if g.get("druggable") else "no",
                        g.get("ec") or "",
                        ", ".join((g.get("interpro_names") or [])[:3]))
                       for g in fa[:50]])]
    return "\n".join(out)


# §7 expression_context -----------------------------------------------------

def r_expression_context(b):
    out = ["## Expression context", "",
           f"**Cohort genes with no expression data: "
           f"{_i(b.get('no_expression_count'))}.**"]
    bd = b.get("breadth_distribution") or {}
    if bd:
        out += ["", "**Breadth distribution (Bgee present_calls):**", ""]
        order = ["narrow (1-5 tissues)", "moderate (6-20)", "broad (>20)", "unknown"]
        ordered = [(k, bd[k]) for k in order if k in bd]
        out.append(table(["Bucket", "Genes"], [(k, _i(v)) for k, v in ordered]))
    tt = b.get("cohort_tissue_counts") or {}
    if tt:
        out += ["", "**Top tissues across cohort:**", ""]
        if isinstance(tt, dict):
            items = list(tt.items())
        else:
            items = list(tt)
        items = sorted(items, key=lambda kv: -kv[1])[:20]
        out.append(table(["Tissue", "Cohort genes"],
                         [(k, _i(v)) for k, v in items]))
    pge = b.get("per_gene_expression") or []
    if pge:
        out += ["", "**Per-gene tissue summary (top 30):**", "",
                table(["Symbol", "Bgee breadth", "FANTOM5 breadth", "SCXA", "Top tissues"],
                      [(g.get("symbol"), g.get("bgee_breadth"),
                        g.get("fantom5_breadth"),
                        "yes" if g.get("scxa_present") else "",
                        ", ".join((g.get("top_tissues") or [])[:3]))
                       for g in pge[:30]])]
    return "\n".join(out)


# §8 protein_interactions ---------------------------------------------------

def r_protein_interactions(b):
    out = ["## Protein interactions among cohort", "",
           f"**Intra-cohort edges: {_i(b.get('cohort_edge_count'))}.**"]
    hubs = b.get("hub_genes") or []
    if hubs:
        out += ["", "**Hub genes (top 10 by interactor count):**", "",
                table(["Symbol", "Interactor count"],
                      [(h.get("symbol"), _i(h.get("interactor_count"))) for h in hubs])]
    edges = b.get("cohort_edges") or []
    if edges:
        out += ["", "**Intra-cohort edges:**", "",
                table(["A", "B", "Sources"],
                      [(e.get("a"), e.get("b"), ", ".join(e.get("sources") or []))
                       for e in edges[:50]])]
    return "\n".join(out)


# §9 structural_data --------------------------------------------------------

def r_structural_data(b):
    out = ["## Structural data", "",
           f"**PDB: {_i(b.get('pdb_count'))} · "
           f"AlphaFold-only: {_i(b.get('alphafold_only_count'))} · "
           f"No structure: {_i(b.get('no_structure_count'))}**"]
    pdb = b.get("pdb_genes") or []
    if pdb:
        out += ["", "**Cohort genes with PDB structures (top 30):**", "",
                table(["Symbol", "UniProt", "PDB entries"],
                      [(g.get("symbol"),
                        f"[{g['uniprot']}](https://www.uniprot.org/uniprotkb/{g['uniprot']})"
                        if g.get("uniprot") else "",
                        _i(g.get("pdb_count"))) for g in pdb[:30]])]
    af = b.get("alphafold_only_genes") or []
    if af:
        out += ["", "**AlphaFold-only cohort genes (top 30 by pLDDT):**", "",
                table(["Symbol", "UniProt", "pLDDT"],
                      [(g.get("symbol"),
                        f"[{g['uniprot']}](https://alphafold.ebi.ac.uk/entry/{g['uniprot']})"
                        if g.get("uniprot") else "",
                        g.get("plddt")) for g in af[:30]])]
    return "\n".join(out)


# §10 drug_targets ----------------------------------------------------------

def r_drug_targets(b):
    out = ["## Drug target analysis", "",
           f"**Approved (phase 4): {_i(b.get('approved_count'))} · "
           f"Phase ≥3: {_i(b.get('phase3_count'))} · "
           f"Phased (≥1): {_i(b.get('phased_count'))} · "
           f"Undrugged: {_i(b.get('undrugged_count'))}**"]
    ag = b.get("approved_genes") or []
    if ag:
        out += ["", "**Genes with approved drugs:**", "",
                table(["Symbol", "Lead drug"],
                      [(g.get("symbol"), g.get("drug") or g.get("top_molecule") or "")
                       for g in ag])]
    tt = b.get("top_targets") or []
    if tt:
        out += ["", "**Top cohort targets by molecule count:**", "",
                table(["Symbol", "Molecules", "Max phase"],
                      [(t.get("symbol"), _i(t.get("molecule_count")),
                        t.get("max_phase")) for t in tt])]
    drugs = b.get("drugs") or []
    if drugs:
        out += ["", "**Drugs targeting cohort genes (top 30):**", "",
                table(["Molecule", "Max phase", "Targets in cohort"],
                      [(f"[{d['name'] or d['id']}](https://www.ebi.ac.uk/chembl/"
                        f"target_report_card/{d['id']}/)",
                        d.get("max_phase"),
                        ", ".join((d.get("gene_targets") or [])[:6]))
                       for d in drugs[:30]])]
    return "\n".join(out)


# §11 bioactivity_enzyme ----------------------------------------------------

def r_bioactivity_enzyme(b):
    out = ["## Bioactivity & enzyme data", "",
           f"**Enzyme cohort genes (≥1 EC): {_i(b.get('enzyme_count'))}.**"]
    # Surface every cohort gene with measurable bioactivity (sorted by
    # assay count). Tighter caps lost the migrated page's long-tail signal.
    pgb = b.get("per_gene_bioactivity") or []
    studied_all = sorted(
        [g for g in pgb if (g.get("chembl_assay_total") or 0) > 0],
        key=lambda g: -(g.get("chembl_assay_total") or 0))
    if studied_all:
        out += ["", "**Cohort genes with ChEMBL bioactivity (full, sorted by assay count):**", "",
                table(["Symbol", "Assays", "Type breakdown"],
                      [(g.get("symbol"), _i(g.get("chembl_assay_total")),
                        ", ".join(f"{k}:{v}" for k, v in (g.get("chembl_assay_types") or {}).items()))
                       for g in studied_all[:50]])]
    eg = b.get("enzyme_genes") or []
    if eg:
        out += ["", "**Cohort enzymes (BRENDA EC):**", "",
                table(["Symbol", "EC numbers", "Names"],
                      [(g.get("symbol"),
                        ", ".join(g.get("ec_numbers") or []),
                        _trunc(", ".join(g.get("ec_names") or []), 60))
                       for g in eg[:50]])]
    usp = b.get("undrugged_starting_points") or []
    if usp:
        out += ["", "**Undrugged cohort genes with high screening signal (≥100 ChEMBL assays):**", "",
                table(["Symbol", "ChEMBL assays", "Note"],
                      [(g.get("symbol"), _i(g.get("chembl_assay_total")),
                        g.get("note") or "") for g in usp[:30]])]
    return "\n".join(out)


# §12 pharmacogenomics ------------------------------------------------------

def r_pharmacogenomics(b):
    out = ["## Pharmacogenomics", "",
           f"**Cohort genes with PharmGKB coverage: {_i(b.get('pgx_gene_count'))} "
           f"(VIP: {_i(b.get('vip_count'))}, CPIC: {_i(b.get('cpic_count'))}).**"]
    pg = b.get("pgx_genes") or []
    if pg:
        out += ["", "**PharmGKB coverage by gene:**", "",
                table(["Symbol", "VIP entries", "CPIC entries"],
                      [(g.get("symbol"),
                        _i(g.get("vip_count")), _i(g.get("cpic_count")))
                       for g in pg[:30]])]
    return "\n".join(out)


# §13 clinical_trials -------------------------------------------------------

def r_clinical_trials(b):
    out = ["## Clinical trials", "",
           f"**Total trials: {_i(b.get('trial_count'))}.**"]
    pc = b.get("phase_counts") or {}
    if pc:
        out += ["", "**Phase distribution (across all retrieved trials):**", ""]
        rows = sorted(pc.items(), key=lambda kv: -kv[1])
        out.append(table(["Phase", "Trials"], [(k, _i(v)) for k, v in rows]))
    tt = b.get("top_trials") or []
    if tt:
        out += ["", "**Top trials by phase / activity:**", "",
                table(["NCT", "Phase", "Status", "Title"],
                      [(f"[{t['id']}](https://clinicaltrials.gov/study/{t['id']})"
                        if t.get("id") else "",
                        t.get("phase"), t.get("status"),
                        _trunc(t.get("title"), 65))
                       for t in tt])]
    td = b.get("trial_drugs") or []
    if td:
        out += ["", "**Drugs tested across these trials (top 30):**", "",
                table(["Molecule", "Max phase", "Trials referencing"],
                      [(d.get("name") or d.get("molecule_id"),
                        d.get("max_phase"), _i(d.get("trial_count")))
                       for d in td[:30]])]

    # CIViC precision-subtype map — the drug × molecular subtype × indication
    # triple. Predictive associations only, ranked by CIViC evidence level
    # (A validated → E inferential). The Subtype column is the molecular
    # profile (EGFR T790M, ALK fusion, KRAS G12C, ...); Effect separates
    # Sensitivity/Response from Resistance. Elides for non-cancer diseases.
    ce = b.get("civic_evidence") or []
    if ce:
        etc = b.get("civic_evidence_type_counts") or {}
        extra = ", ".join(f"{n} {k.lower()}" for k, n in etc.items() if k != "Predictive")
        out += ["",
                (f"**Precision-medicine subtype map — drug × molecular subtype "
                 f"(CIViC, {_i(b.get('civic_association_total'))} predictive "
                 f"associations from {_i(b.get('civic_predictive_total'))} curated "
                 f"evidence items"
                 + (f"; also {extra}" if extra else "") + "):**"),
                "",
                table(["Molecular subtype", "Therapy", "Effect", "Level", "CIViC"],
                      [(r["profile"], therapy_label(r["therapy"]), r["significance"],
                        f"CIViC {r['level']}" if r.get("level") else "",
                        f"[EID{r['evidence_id']}](https://civicdb.org/links/evidence_items/{r['evidence_id']})"
                        + (f" +{r['n']-1}" if r.get("n", 1) > 1 else ""))
                       for r in ce])]
        more = (b.get("civic_association_total") or 0) - len(ce)
        if more > 0:
            out += ["", (f"*+{more} more predictive associations (showing top "
                         f"{len(ce)} by evidence level).*")]
    return "\n".join(out)


# §14 pathways --------------------------------------------------------------

def r_pathways(b):
    out = ["## Pathway analysis", "",
           f"**Distinct Reactome pathways touched by cohort: "
           f"{_i(b.get('pathway_count'))}.**"]
    tp = b.get("top_pathways") or []
    if tp:
        out += ["", "**Top pathways by cohort coverage:**", "",
                table(["Pathway", "Genes", "Sample cohort genes"],
                      [(f"[{p.get('name') or p['id']}]"
                        f"(https://reactome.org/PathwayBrowser/#/{p['id']})",
                        _i(p.get("gene_count")),
                        ", ".join((p.get("gene_symbols") or [])[:8]))
                       for p in tp[:30]])]
    return "\n".join(out)


# Derived render-only sections ---------------------------------------------

def r_drug_repurposing(bundles):
    """§15 — compounds with ChEMBL bioactivity against a cohort gene that aren't
    yet in disease-level trials (§10 drugs ∖ §13 trial_drugs).

    The framing is gated on disease class. A bioactivity row is only a SCREENING
    signal — a promiscuous kinase inhibitor assayed against a cohort kinase
    (e.g. sorafenib vs STK36 in primary ciliary dyskinesia) is an off-target
    artifact, not a treatment. So:
      - Cancers (target-inhibition has therapeutic rationale + driver evidence):
        framed as "repurposing candidates", still caveated as a screening signal.
      - Non-cancer / Mendelian / structural diseases: framed honestly as
        "chemical tractability — a research signal, NOT a therapeutic
        recommendation", since target-inhibition repurposing usually has no
        disease mechanism there.
    """
    b10 = bundles.get("10") or {}
    b13 = bundles.get("13") or {}
    is_cancer = bool((bundles.get("1") or {}).get("is_cancer"))
    drugs10 = {d["id"]: d for d in (b10.get("drugs") or []) if d.get("id")}
    trial_mol_ids = {d.get("molecule_id") for d in (b13.get("trial_drugs") or [])}
    repurposable = [d for mid, d in drugs10.items() if mid not in trial_mol_ids]
    repurposable.sort(key=lambda d: -(d.get("max_phase") or 0))

    if is_cancer:
        title = "## Drug repurposing candidates"
        lead = (f"**{len(repurposable)} approved/phased drugs hit cohort targets "
                "but don't yet appear in disease-level clinical trials.** "
                "Target-inhibition rationale is strongest for cancer driver genes; "
                "a bioactivity hit is a screening signal, not a treatment claim.")
    else:
        title = "## Chemical tractability of cohort targets"
        lead = (f"**{len(repurposable)} approved/phased compounds have measured "
                "bioactivity against a cohort gene** (and aren't yet in "
                "disease-level trials). This is a *research / tractability signal, "
                "NOT a therapeutic recommendation* — a bioactivity row often "
                "reflects off-target or screening binding (e.g. promiscuous kinase "
                "inhibitors against a cohort kinase), implying no disease mechanism.")

    out = [title, "", lead]
    if repurposable:
        out += ["", table(["Compound", "Max phase", "Cohort target (bioactivity)"],
                          [(d.get("name") or d.get("id"),
                            d.get("max_phase"),
                            ", ".join((d.get("gene_targets") or [])[:6]))
                           for d in repurposable[:30]])]
    return "\n".join(out)


def r_druggability_pyramid(bundles):
    """§16 — druggability tier counts derived from §6 + §9 + §10.
    Tiers (top → bottom):
      Tier A — approved drug exists
      Tier B — phased (≥1) drug exists, no approved
      Tier C — druggable family + PDB structure, no drug
      Tier D — druggable family, AlphaFold only, no drug
      Tier E — difficult family or no structure, no drug
    """
    b6 = bundles.get("6") or {}
    b9 = bundles.get("9") or {}
    b10 = bundles.get("10") or {}

    fam_by = {g["symbol"]: g for g in (b6.get("family_assignments") or [])}
    struct_by = {g["symbol"]: g for g in (b9.get("per_gene_structure") or [])}
    drug_by = {g["symbol"]: g for g in (b10.get("per_gene_drugs") or [])}

    tiers = Counter()
    tier_members: "dict[str, list[str]]" = {k: [] for k in "ABCDE"}
    for sym, dg in drug_by.items():
        fg = fam_by.get(sym) or {}
        sg = struct_by.get(sym) or {}
        if (dg.get("max_phase") or 0) >= 4:
            t = "A"
        elif (dg.get("max_phase") or 0) >= 1:
            t = "B"
        elif fg.get("druggable") and sg.get("pdb_count"):
            t = "C"
        elif fg.get("druggable") and sg.get("has_alphafold"):
            t = "D"
        else:
            t = "E"
        tiers[t] += 1
        tier_members[t].append(sym)

    labels = {
        "A": "Approved (phase 4 drug)",
        "B": "Phased (≥1) drug, not yet approved",
        "C": "Druggable family + PDB, no drug",
        "D": "Druggable family + AlphaFold only, no drug",
        "E": "Difficult family or no structure, no drug",
    }
    def _members(t):
        m = tier_members[t]
        if not m:
            return ""
        if len(m) <= 10:
            return ", ".join(m)
        return ", ".join(m[:10]) + f" (+{len(m) - 10} more)"

    out = ["## Druggability pyramid", "",
           "Cohort genes binned by druggability tier (high → low):", "",
           table(["Tier", "Definition", "Genes", "Symbols"],
                 [(t, labels[t], _i(tiers.get(t, 0)), _members(t))
                  for t in "ABCDE"])]
    return "\n".join(out)


def r_undrugged_target_profiles(bundles):
    """§17 — undrugged targets sorted by 'how close to a drug start' they are:
    cross-join §8 (interaction with drugged genes) + §10 (drugged flag) +
    §11 (assay depth)."""
    b8 = bundles.get("8") or {}
    b10 = bundles.get("10") or {}
    b11 = bundles.get("11") or {}

    drugged = {g["symbol"] for g in (b10.get("per_gene_drugs") or [])
               if (g.get("max_phase") or 0) >= 1}
    drug_by = {g["symbol"]: g for g in (b10.get("per_gene_drugs") or [])}
    bio_by = {g["symbol"]: g for g in (b11.get("per_gene_bioactivity") or [])}
    inter_by = {g["symbol"]: g for g in (b8.get("per_gene_interactions") or [])}

    rows = []
    for sym, dg in drug_by.items():
        if sym in drugged:
            continue
        interactors = inter_by.get(sym, {}).get("top_interactors") or []
        drugged_partners = [p for p in interactors if p in drugged][:3]
        assay_n = (bio_by.get(sym) or {}).get("chembl_assay_total") or 0
        # Score: prefer high assay count + drugged partner present
        score = (int(assay_n or 0) // 100) + (5 if drugged_partners else 0)
        rows.append({"symbol": sym, "assays": assay_n,
                     "drugged_partners": drugged_partners, "score": score})
    rows.sort(key=lambda r: -r["score"])

    out = ["## Undrugged target profiles", "",
           f"**{len(rows)} cohort genes are undrugged.** "
           "Ranked by 'starting-point quality' (assay depth + drugged-partner adjacency)."]
    if rows:
        out += ["", table(["Symbol", "ChEMBL assays", "Drugged partners (top 3)"],
                          [(r["symbol"], _i(r["assays"]),
                            ", ".join(r["drugged_partners"]) or "—")
                           for r in rows[:30]])]
    return "\n".join(out)


# Registry ------------------------------------------------------------------

RENDER = {
    "1": r_disease_ids,
    "2": r_gwas_landscape,
    "3": r_variant_details,
    "4": r_mendelian_overlap,
    "5": r_genes_proteins,
    "6": r_protein_families,
    "7": r_expression_context,
    "8": r_protein_interactions,
    "9": r_structural_data,
    "10": r_drug_targets,
    "11": r_bioactivity_enzyme,
    "12": r_pharmacogenomics,
    "13": r_clinical_trials,
    "14": r_pathways,
}


def render_all(bundles):
    """Full body = §1..§14 (RENDER dict) + §15..§17 (derived). bundles is
    the full {section_id: section_bundle} dict produced by atlas.disease.collect.collect_all."""
    parts = [RENDER[sid](bundles[sid]) for sid in sorted(RENDER, key=int)]
    parts.append(r_drug_repurposing(bundles))
    parts.append(r_druggability_pyramid(bundles))
    parts.append(r_undrugged_target_profiles(bundles))
    return "\n\n".join(parts)


if __name__ == "__main__":
    import sys
    from atlas.disease.collect import collect_all
    name = sys.argv[1] if len(sys.argv) > 1 else "age-related macular degeneration"
    print(render_all(collect_all(name)))
