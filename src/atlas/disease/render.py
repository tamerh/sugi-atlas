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
import re
from collections import Counter
from atlas.render_common import table, fnum, gencc_rank, phase_label
from atlas.civic import therapy_label
from atlas.page import links

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


_GENCC_STOP = {"disease", "diseases", "syndrome", "cancer", "carcinoma", "tumor",
               "tumour", "neoplasm", "disorder", "disorders", "type", "familial",
               "hereditary", "susceptibility", "predisposition", "complementation",
               "group", "deficiency", "autosomal", "dominant", "recessive",
               "with", "without"}


def _disease_tokens(s):
    """Substantive tokens of a disease name for on-disease matching (drops the
    generic descriptors that would over-match)."""
    return {t for t in re.findall(r"[a-z0-9]+", (s or "").lower())
            if len(t) >= 4 and t not in _GENCC_STOP}


def _dedup_gencc(rows, disease_name=None):
    """Collapse GenCC rows to one per gene. Prefer the record FOR the page's
    disease (so an off-disease but stronger record — BRCA2's Fanconi-D1 on a
    medulloblastoma page — doesn't outrank the on-disease one); among the chosen
    pool keep the strongest classification. Returns
    [(best_row, record_count, on_disease)] with on-disease genes sorted first."""
    dn = _disease_tokens(disease_name)
    by = {}
    for r in rows:
        sym = r.get("symbol")
        if sym:
            by.setdefault(sym, []).append(r)
    out = []
    for rs in by.values():
        on = [r for r in rs if dn and (_disease_tokens(r.get("mondo_disease")) & dn)]
        best = max(on or rs, key=lambda r: gencc_rank(r.get("gencc_classification")))
        out.append((best, len(rs), bool(on)))
    out.sort(key=lambda bn: (not bn[2],
                             -gencc_rank(bn[0].get("gencc_classification")),
                             bn[0].get("symbol") or ""))
    return out


def _data_availability(xc):
    """Compact evidence-volume line from the Mondo xref counts — only the
    datasets that signal *how much data exists* for the disease (true totals).
    Identity/ontology xrefs are excluded (they're in the identifier table)."""
    if not xc:
        return None
    items = []

    def add(key, singular):
        n = xc.get(key) or 0
        if n:
            items.append(f"{_i(n)} {singular}" + ("s" if n != 1 else ""))

    add("clinvar", "ClinVar variant")
    # NB: clinical_trials xref is the RAW (contaminated) count — omitted here;
    # the title-validated count lives in §13.
    add("clingen_variant", "ClinGen variant curation")
    g, gs = xc.get("gwas") or 0, xc.get("gwas_study") or 0
    if g:
        s = f"{_i(g)} GWAS association" + ("s" if g != 1 else "")
        if gs:
            s += f" ({_i(gs)} stud" + ("ies" if gs != 1 else "y") + ")"
        items.append(s)
    add("gencc", "GenCC gene-disease record")
    add("hpo", "HPO phenotype")
    add("cellosaurus", "cell line")
    add("intogen", "intOGen driver record")
    return " · ".join(items) if items else None


# §1 disease_ids ------------------------------------------------------------

def r_disease_ids(b):
    # Skip rows whose value is empty so the ID table only carries facts.
    # Cleaner than rendering empty `| |` cells for diseases that don't have
    # the given cross-reference (most diseases miss at least one of MeSH /
    # OMIM / Orphanet).
    candidate_rows = [
        ("Canonical name", b.get("canonical_name")),
        ("Mondo ID", b.get("mondo_id") or None),
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

    # Synonyms (Mondo) — alternate names people search by ("AT", "Lou Gehrig's
    # disease", brand/abbreviation/historical names). High search/AI value.
    syns = b.get("synonyms") or []
    if syns:
        # " · " separator (not comma) — Mondo synonyms often contain internal
        # commas ("AT, complementation group A"), which a comma-join would blur.
        out += ["", "**Also known as:** " + " · ".join(syns[:20])
                + (f" (+{len(syns) - 20} more)" if len(syns) > 20 else "")]

    # Data availability — a compact "how much evidence exists" line, filtered to
    # the evidence-volume datasets (true totals, beyond the capped per-section
    # displays). The raw Mondo xref-count dump is dropped: its identity-mapping
    # counts of 1 (doid/efo/mesh/…) just duplicate the identifier table above,
    # and the ontology plumbing (mondochild/parent) belongs in disease navigation.
    avail = _data_availability(b.get("xref_counts") or {})
    if avail:
        out += ["", "**Data availability:** " + avail + "."]

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

    # Clinical features (HPO phenotypes from Orphanet, frequency-sorted). These
    # are curated, frequency-annotated, and all clinically meaningful (not noisy
    # like a GWAS tail), so show a generous top 50; the full list is in the
    # bundle/sidecar for RAG consumers.
    phs = b.get("phenotypes") or []
    if phs:
        total = b.get("phenotype_count") or len(phs)
        out.append("")
        out.append(f"**Clinical features ({total} HPO phenotypes, Orphanet "
                   f"curated; top {min(50, len(phs))} by frequency):**")
        out.append("")
        out.append(table(["HPO ID", "Term", "Frequency"],
                         [(p.get("hpo_id"),
                           p.get("hpo_term"),
                           p.get("frequency"))
                          for p in phs[:50]]))
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
           f"Top hits map to {_i(b.get('unique_gene_count'))} distinct genes "
           "(as reported by GWAS)."]
    ta = b.get("top_assocs") or []
    if ta:
        out += ["", "**Top associations by p-value:**", "",
                table(["rsID", "p-value", "Gene", "Risk allele", "Odds ratio"],
                      # fnum on odds_ratio only (audit #7: float32 artifacts like
                      # 1.4348061); p-value is exponent-form (3e-67) and must NOT
                      # be rounded — fnum(3e-67, 2) would collapse it to 0.
                      [(r.get("rsid"), r.get("pvalue"), r.get("gene_symbol"),
                        r.get("risk_allele"), fnum(r.get("odds_ratio")))
                       for r in ta[:30]])]
    studies = b.get("studies") or []
    if studies:
        out += ["", "**Top studies (by case count):**", "",
                table(["Study", "Lead author", "Year", "Cases", "Controls", "Title"],
                      [(s.get("id") or "",
                        s.get("lead_author"), s.get("year"),
                        _i(s.get("sample_size_cases")),
                        _i(s.get("sample_size_controls")),
                        _trunc(s.get("title"), 60))
                       for s in studies])]
    return "\n".join(out)


# §3 variant_details --------------------------------------------------------

def r_variant_details(b):
    out = ["## Variant details and genetic-evidence tiers", ""]
    # This section tiers GWAS-derived variants (>>mondo>>gwas>>dbsnp). For
    # diseases with no GWAS (Mendelian / rare) it's empty — render a note, not
    # all-zero tables. Guard each sub-block on real data (non-zero sum), since a
    # dict of zeros is still truthy.
    tv = b.get("top_variants") or []
    tc = b.get("tier_counts") or {}
    md = b.get("maf_distribution") or {}
    cd = b.get("consequence_distribution") or {}
    cv = b.get("clinvar_variants") or []
    gwas_present = bool(tv) or any((tc or {}).values())
    if not gwas_present and not cv:
        out.append("*No tiered GWAS variants or ClinVar records for this disease.*")
        return "\n".join(out)

    # ClinVar germline variants — the rare/coding genetic evidence (primary for
    # Mendelian diseases, which have no GWAS). Rendered first when there's no
    # GWAS signal; otherwise after the GWAS tiers.
    def _clinvar_block():
        cc = b.get("clinvar_class_counts") or {}
        # clinvar_total here is a paginated fetch count (caps at ~600), not the
        # disease's full ClinVar set — present it as a retrieved sample / floor
        # so it doesn't contradict the accurate xref total in "At a glance".
        bl = [f"**ClinVar germline variants ({_i(b.get('clinvar_total'))} retrieved; "
              f"paginated sample, class counts are floors):**"]
        if cc:
            order = sorted(cc.items(), key=lambda kv: -kv[1])
            bl.append("\n" + ", ".join(f"{_i(v)} {k.lower()}" for k, v in order))
        bl += ["", table(["ClinVar", "Variant (HGVS)", "Gene", "Classification", "Review"],
                         [(v.get("id") or "",
                           v.get("hgvs"), v.get("gene"), v.get("classification"),
                           v.get("review_status")) for v in cv])]
        return bl

    if cv and not gwas_present:
        out += _clinvar_block()
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
                      [(r.get("rsid") or "",
                        r.get("chrom"), r.get("pos"), r.get("alleles"),
                        fnum(r.get("maf"), 3), r.get("consequence"),
                        r.get("gene_symbol"), r.get("pvalue"), r.get("tier"))
                       for r in tv[:30]])]
    if cv:  # both GWAS and ClinVar present → ClinVar after the GWAS tiers
        out += ["", *_clinvar_block()]
    return "\n".join(out)


# §4 mendelian_overlap ------------------------------------------------------

def r_mendelian_overlap(b):
    out = ["## Mendelian disease overlap and somatic drivers", "",
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
                      [(sym, sym, _routes(sym))
                       for sym in dual[:50]])]
    sg = b.get("somatic_driver_genes") or []
    if sg:
        # CIViC's `name` field equals the gene symbol — we surface the CIViC
        # gene ID (plain text) as the actual value-add.
        out += ["", "**Somatic driver evidence (intOGen + CIViC, cohort fanout):**", "",
                table(["Gene", "intOGen role", "Cancer types", "CIViC"],
                      [(g.get("symbol"),
                        (g.get("intogen") or {}).get("role"),
                        _trunc((g.get("intogen") or {}).get("cancer_types") or "", 50),
                        (f"CIViC #{(g.get('civic') or {}).get('id')}"
                         if (g.get("civic") or {}).get("id") else ""))
                       for g in sg[:30]])]
    gc = b.get("gencc_genes") or []
    if gc:
        # Collapse to one row per gene (audit #13: cohort fan-out pulls every
        # GenCC submission — Lynch syndrome had MLH1 ×19); prefer each gene's
        # record FOR this disease, else its strongest, and surface the count.
        ded = _dedup_gencc(gc, b.get("disease_name"))
        out += ["", "**GenCC gene–disease validity (cohort genes):** *the Disease "
                "column is the GenCC-asserted condition — a cohort gene's "
                "strongest validity may be for a related predisposition syndrome.*", "",
                table(["Gene", "Classification", "Inheritance", "Disease", "Records"],
                      [(g.get("symbol"), g.get("gencc_classification"),
                        g.get("mode_of_inheritance"),
                        _trunc(g.get("mondo_disease"), 50),
                        str(n) if n > 1 else "") for g, n, _on in ded[:30]])]
    og = b.get("orphanet_genes") or []
    if og:
        out += ["", "**Orphanet rare-disease linkage (cohort genes):**", "",
                table(["Gene", "Orphanet ID", "Rare disease"],
                      [(g.get("symbol"),
                        g.get("orphanet_id") or "",
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
                      [(links.maybe_link(g.get("symbol"), links.gene_url(symbol=g.get("symbol"), hgnc_id=g.get("hgnc_id"))),
                        g.get("hgnc_id"), g.get("ensembl_id"),
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
        items = sorted(items, key=lambda kv: -kv[1])[:40]
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
                        g.get("uniprot") or "",
                        _i(g.get("pdb_count"))) for g in pdb[:30]])]
    af = b.get("alphafold_only_genes") or []
    if af:
        out += ["", "**AlphaFold-only cohort genes (top 30 by pLDDT):**", "",
                table(["Symbol", "UniProt", "pLDDT"],
                      [(g.get("symbol"),
                        g.get("uniprot") or "",
                        g.get("plddt")) for g in af[:30]])]
    return "\n".join(out)


# §10 drug_targets ----------------------------------------------------------

def r_drug_targets(b):
    out = ["## Drug target analysis", "",
           f"**Approved (phase 4): {_i(b.get('approved_count'))} · "
           f"Phase ≥3: {_i(b.get('phase3_count'))} · "
           f"Phased (≥1): {_i(b.get('phased_count'))} · "
           f"Undrugged: {_i(b.get('undrugged_count'))}**"]
    es, dr = b.get("enrichment_size") or 0, b.get("enrichment_druggable") or 0
    if es:
        pct = round(100 * dr / es) if es else 0
        out += ["", f"**Druggability breadth:** {_i(dr)} of {_i(es)} "
                f"evidence-associated genes ({pct}%) have a ChEMBL target "
                f"(buckets above are over the deeply-mined display cohort)."]
    ag = b.get("approved_genes") or []
    if ag:
        out += ["", "**Genes with approved drugs:**", "",
                table(["Symbol", "Lead drug"],
                      [(links.maybe_link(g.get("symbol"), links.gene_url(symbol=g.get("symbol"))),
                        g.get("drug") or g.get("top_molecule") or "")
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
                      [(links.maybe_link(d.get("name") or d.get("id") or "",
                                         links.drug_url(chembl_id=d.get("id"), name=d.get("name"))),
                        d.get("max_phase"),
                        ", ".join((d.get("gene_targets") or [])[:6]))
                       for d in drugs[:30]])]
    return "\n".join(out)


# §11 bioactivity_enzyme ----------------------------------------------------

def r_bioactivity_enzyme(b):
    out = ["## Bioactivity and enzyme data", "",
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
    # NB: the upstream PharmGKB `is_vip` flag is degenerate here — it returns
    # true for every gene record (so vip_count == "has a PharmGKB record", not a
    # real Very-Important-Pharmacogene designation). Surfacing it as "VIP" said
    # all 50 cohort genes were VIPs, which is wrong — so we drop it and report
    # only PharmGKB coverage + CPIC dosing-guideline counts (the real signal).
    cpic = b.get("cpic_count") or 0
    out = ["## Pharmacogenomics", "",
           f"**Cohort genes with a PharmGKB record: {_i(b.get('pgx_gene_count'))}; "
           f"with CPIC/DPWG dosing guidelines: {_i(cpic)}.**"]
    pg = [g for g in (b.get("pgx_genes") or []) if (g.get("cpic_count") or 0)]
    if pg:
        out += ["", "**Cohort genes with a CPIC/DPWG dosing guideline:**", "",
                table(["Symbol", "CPIC guidelines"],
                      [(g.get("symbol"), _i(g.get("cpic_count"))) for g in pg[:30]])]
    elif not cpic:
        out += ["", "*No cohort gene has a CPIC/DPWG genotype-guided dosing "
                "guideline (PharmGKB).*"]
    return "\n".join(out)


# §13 clinical_trials -------------------------------------------------------

def r_clinical_trials(b):
    n = b.get("trial_count") or 0
    out = ["## Clinical trials", "",
           f"**Clinical trials: {_i(n)}.**"]
    pc = b.get("phase_counts") or {}
    if pc:
        out += ["", "**Phase distribution (across all retrieved trials):**", ""]
        rows = sorted(pc.items(), key=lambda kv: -kv[1])
        out.append(table(["Phase", "Trials"], [(k, _i(v)) for k, v in rows]))
    tt = b.get("top_trials") or []
    if tt:
        out += ["", "**Top trials by phase / activity:**", "",
                table(["NCT", "Phase", "Status", "Title"],
                      [(t.get("id") or "",
                        phase_label(t.get("phase")), t.get("status"),
                        _trunc(t.get("title"), 65))
                       for t in tt])]
    td = b.get("trial_drugs") or []
    if td:
        out += ["", "**Drugs tested across these trials (top 30):**", "",
                table(["Molecule", "Max phase", "Trials referencing"],
                      [(links.maybe_link(d.get("name") or d.get("molecule_id"),
                                         links.drug_url(chembl_id=d.get("molecule_id"), name=d.get("name"))),
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
                      [(r["profile"],
                        links.maybe_link(therapy_label(r["therapy"]), links.drug_url(name=therapy_label(r["therapy"]))),
                        r["significance"],
                        f"CIViC {r['level']}" if r.get("level") else "",
                        f"EID{r['evidence_id']}"
                        + (f" +{r['n']-1}" if r.get("n", 1) > 1 else ""))
                       for r in ce])]
        more = (b.get("civic_association_total") or 0) - len(ce)
        if more > 0:
            out += ["", (f"*+{more} more predictive associations (showing top "
                         f"{len(ce)} by evidence level).*")]
    return "\n".join(out)


# §14 pathways --------------------------------------------------------------

def r_pathways(b):
    es, gp = b.get("enrichment_size") or 0, b.get("genes_with_pathways") or 0
    over = (f" Enrichment computed across {_i(es)} evidence-associated genes "
            f"({_i(gp)} with Reactome annotation)." if es else "")
    out = ["## Pathway analysis", "",
           f"**Distinct Reactome pathways touched by cohort: "
           f"{_i(b.get('pathway_count'))}.**{over}"]
    tp = b.get("top_pathways") or []
    if tp:
        out += ["", "**Top pathways by cohort coverage:**", "",
                table(["Pathway", "Genes", "Sample cohort genes"],  # table() dedups
                      [(p.get("name") or p.get("id") or "",
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
    """Disease page body in the FROZEN canonical H2 order (docs/PAGE_CONTRACT.md):
    Identifiers → Genetics & variants → Genes & proteins → Function → Therapeutics
    → Clinical trials & evidence. (Summary wrapped by assemble_page; Related
    appended after.) The current 18 section-renderers fold into these as H3."""
    from atlas.render_common import demote, emit_canonical, with_heading_id

    def S(s, anchor):
        return with_heading_id(demote(RENDER[s](bundles[s])), anchor)

    def D(md, anchor):
        return with_heading_id(demote(md), anchor)

    def join(*parts):
        return "\n\n".join(p for p in parts if p and p.strip())

    spec = [
        ("Identifiers", "identifiers", S("1", "disease-ids"), None),
        ("Genetics & variants", "genetics",
         join(S("2", "gwas"), S("3", "variant-tiers")),
         "No common-variant (GWAS) or curated variant data for this disease."),
        ("Genes & proteins", "genes",
         join(S("4", "mendelian"), S("5", "cohort-genes"), S("6", "protein-families"),
              S("7", "expression"), S("8", "interactions"), S("9", "structural")),
         "No associated genes curated for this disease."),
        ("Function", "function", S("14", "pathways"),
         "No pathway enrichment — requires an associated-gene cohort."),
        ("Therapeutics", "drugs",
         join(S("10", "drug-targets"), S("11", "bioactivity"), S("12", "pharmacogenomics"),
              D(r_drug_repurposing(bundles), "tractability"),
              D(r_druggability_pyramid(bundles), "druggability"),
              D(r_undrugged_target_profiles(bundles), "undrugged")),
         "No druggable-target or therapeutic data for this disease's cohort."),
        ("Clinical trials & evidence", "trials", S("13", "clinical-trials"),
         "No clinical trials or CIViC evidence naming this disease."),
    ]
    return emit_canonical(spec)


if __name__ == "__main__":
    import sys
    from atlas.disease.collect import collect_all
    name = sys.argv[1] if len(sys.argv) > 1 else "age-related macular degeneration"
    print(render_all(collect_all(name)))
