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
    # Display truncation is the frontend's job (CSS ellipsis); the source markdown
    # must carry the full value. Web-team report: titles/names were arriving
    # pre-clipped with a baked-in "…". Signature kept so the call sites are
    # unchanged; `n` is now ignored.
    return (s or "").strip()


# Shared with the dual-evidence on-disease filter (atlas.disease.cohort) so both
# the GenCC dedup here and §4 use identical disease-name matching.
from atlas.disease.cohort import disease_tokens as _disease_tokens


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
    return "\n".join(out)


def r_epidemiology(b):
    """Epidemiology — Orphanet multi-geography prevalence (validated first).
    Moved from §1 into the Clinical features zone (it's clinical profile, not an
    identifier). Returns '' when no prevalence data so the zone elides."""
    prevs = b.get("prevalences") or []

    def _informative(p):
        # Drop rows with no quantitative content — blank value AND a class that's
        # blank or "Unknown" (e.g. "Prevalence at birth | Unknown | | Hungary").
        cls = str(p.get("prevalence_class") or "").strip().lower()
        val = str(p.get("val_moy") or "").strip()
        return bool(val) or (cls and cls != "unknown")

    prevs = [p for p in prevs if _informative(p)]
    if not prevs:
        return ""

    def _geo_rank(g):                       # broadest geography first
        g = (g or "").lower()
        return 0 if "worldwide" in g else (1 if "europe" in g else 2)

    prevs_sorted = sorted(
        prevs,
        key=lambda p: (0 if p.get("validation_status") == "Validated" else 1,
                       _geo_rank(p.get("geographic")), p.get("prevalence_type") or ""))
    # Cap — rare diseases (cystic fibrosis: 62 records, every EU country ×
    # point/birth prevalence) would otherwise dwarf the symptoms table below.
    shown = prevs_sorted[:20]
    note = (f", top {len(shown)} (validated / broadest geography first)"
            if len(prevs_sorted) > len(shown) else "")
    return "\n".join(
        ["## Epidemiology", "",
         "### Prevalence records {#prevalence}", "",
         f"{len(prevs)} prevalence record(s), Orphanet{note}:", "",
         table(["Type", "Class", "Value", "Geography", "Validation"],
               [(p.get("prevalence_type"), p.get("prevalence_class"),
                 p.get("val_moy"), p.get("geographic"), p.get("validation_status"))
                for p in shown])])


def r_symptoms(b):
    """Signs & symptoms — Orphanet-curated HPO clinical features, frequency-sorted
    (top 50; full list in the bundle/JSON-LD sidecar). The headline clinical
    presentation of a disease, so it gets a first-class section right after the
    summary rather than being buried under identifiers."""
    phs = b.get("phenotypes") or []
    if not phs:
        return ""
    total = b.get("phenotype_count") or len(phs)
    return "\n".join(
        ["## Signs & symptoms", "",
         "### Clinical features (HPO) {#hpo-features}", "",
         f"{total} HPO clinical features (Orphanet curated; top "
         f"{min(50, len(phs))} by frequency):", "",
         table(["HPO ID", "Term", "Frequency"],
               [(p.get("hpo_id"), p.get("hpo_term"), p.get("frequency"))
                for p in phs[:50]])])


def r_disease_family(b1, b5):
    """Mondo ontology family — broader term (parent) + subtypes (children), with
    a pointer to the parent when this page's own cohort is sparse. Granular
    subtype terms (IDH-wildtype glioblastoma, MONDO:0850335) often carry little
    direct evidence while the broader term (glioblastoma) holds the cohort,
    trials, and CIViC — so we route the reader there. All links are
    manifest-gated (render only when the target is a built page). Body only — the
    "## Disease family {#family}" heading comes from the canonical zone."""
    b1 = b1 or {}
    parent = b1.get("parent") or None
    ancestors = b1.get("ancestors") or []
    children = b1.get("children") or []
    siblings = b1.get("siblings") or []
    if not (parent or ancestors or children or siblings):
        return ""

    def _dl(t):  # manifest-gated disease link by id (label = name)
        return links.maybe_link(t.get("name"),
                                links.disease_url(mondo_id=t.get("id"), name=t.get("name")))

    out = []
    parent_url = (links.disease_url(mondo_id=parent.get("id"), name=parent.get("name"))
                  if parent else None)
    gene_count = (b5 or {}).get("gene_count") or 0
    # Lead: route a sparse subtype to its parent; else flag an umbrella term.
    if parent and parent_url and gene_count == 0:
        out += [f"This is a subtype of **{links.maybe_link(parent.get('name'), parent_url)}**. "
                "Genetic, therapeutic, and trial evidence is largely curated at the "
                "broader-term level — see the parent page for the associated-gene "
                "cohort and molecular evidence.", ""]
    elif b1.get("child_count"):
        n = b1["child_count"]
        out += [f"An umbrella term covering {n} Mondo subtype"
                f"{'s' if n != 1 else ''}.", ""]
    # Breadcrumb: root → … → parent → THIS (current term bold, not linked).
    if ancestors:
        crumb = " › ".join(_dl(t) for t in reversed(ancestors))
        cur = b1.get("canonical_name") or b1.get("name") or "this disease"
        out.append(f"**Classification path:** {crumb} › **{cur}**")
    # Siblings (co-subtypes under the same parent) — lateral navigation.
    if siblings:
        out.append(f"\n**Related subtypes ({len(siblings)}):** "
                   + ", ".join(_dl(t) for t in siblings))
    # Children (this term's own subtypes).
    if children:
        out.append(f"\n**Subtypes ({len(children)}):** "
                   + ", ".join(_dl(t) for t in children))
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
        out += ["", "### Top associations by p-value {#gwas-associations}", "",
                table(["rsID", "p-value", "Gene", "Risk allele", "Odds ratio"],
                      # fnum on odds_ratio only (audit #7: float32 artifacts like
                      # 1.4348061); p-value is exponent-form (3e-67) and must NOT
                      # be rounded — fnum(3e-67, 2) would collapse it to 0.
                      [(r.get("rsid"), r.get("pvalue"), r.get("gene_symbol"),
                        r.get("risk_allele"), fnum(r.get("odds_ratio")))
                       for r in ta[:30]])]
    studies = b.get("studies") or []
    if studies:
        out += ["", "### Top studies (by case count) {#gwas-studies}", "",
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
        bl = ["### ClinVar germline variants {#clinvar-variants}", "",
              f"{_i(b.get('clinvar_total'))} retrieved; paginated sample, class counts "
              f"are floors:"]
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
        out += ["### Tier distribution (top 50 variants) {#tier-distribution}", ""]
        rows = sorted(tc.items(), key=lambda kv: kv[0])
        out.append(table(["Tier", "Variants"], [(k, _i(v)) for k, v in rows]))
    if any(md.values()):
        out += ["", "### MAF distribution {#maf-distribution}", ""]
        out.append(table(["Bucket", "Variants"], [(k, _i(v)) for k, v in md.items()]))
    if any(cd.values()):
        out += ["", "### Functional consequences {#consequences}", ""]
        rows = sorted(cd.items(), key=lambda kv: -kv[1])[:15]
        out.append(table(["Consequence", "Count"], [(k, _i(v)) for k, v in rows]))
    if tv:
        out += ["", "### Top variants {#top-variants}", "",
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
        out += ["", "### Dual-evidence genes (GWAS + Mendelian — highest-confidence targets) {#dual-evidence}", "",
                table(["Gene", "HGNC", "Evidence routes"],
                      [(sym, sym, _routes(sym))
                       for sym in dual[:50]])]
    sg = b.get("somatic_driver_genes") or []
    if sg:
        # CIViC's `name` field equals the gene symbol — we surface the CIViC
        # gene ID (plain text) as the actual value-add.
        out += ["", "### Somatic driver evidence (intOGen + CIViC, cohort fanout) {#somatic-drivers}", "",
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
        out += ["", "### GenCC gene–disease validity (cohort genes) {#gencc-validity}", "",
                "*the Disease column is the GenCC-asserted condition — a cohort gene's "
                "strongest validity may be for a related predisposition syndrome.*", "",
                table(["Gene", "Classification", "Inheritance", "Disease", "Records"],
                      [(g.get("symbol"), g.get("gencc_classification"),
                        g.get("mode_of_inheritance"),
                        _trunc(g.get("mondo_disease"), 50),
                        str(n) if n > 1 else "") for g, n, _on in ded[:30]])]
    og = b.get("orphanet_genes") or []
    if og:
        out += ["", "### Orphanet rare-disease linkage (cohort genes) {#orphanet-linkage}", "",
                table(["Gene", "Orphanet ID", "Rare disease"],
                      [(g.get("symbol"),
                        g.get("orphanet_id") or "",
                        _trunc(g.get("orphanet_name"), 65))
                       for g in og[:50]])]
    omg = b.get("omim_genes") or []
    if omg:
        out += ["", "### OMIM-shared genes (cohort gene's MIM ids overlap the disease's) {#omim-shared}", "",
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
    if ev and any(ev.values()):              # skip the all-zero partition (no cohort)
        out += ["", "### Evidence partition {#evidence-partition}", ""]
        out.append(table(["Subset", "Genes"],
                         [(k, _i(v)) for k, v in ev.items() if v]))
    genes = b.get("genes") or []
    if genes:
        out += ["", "### Cohort genes (full) {#cohort-genes-full}", "",
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
        out += ["", "### Cohort function summary {#cohort-function}", "",
                "Lead sentence per gene, UniProt-curated.", "",
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
        fe = b.get("family_enrichment") or {}

        def _fdr(q):
            return "—" if q is None else (f"{q:.0e}" if q < 1e-3 else f"{q:.3f}")

        if fe and any((fe.get(k) or {}).get("fdr") is not None for k in fc):
            def _fkey(kv):                       # enriched (lowest FDR) first
                e = fe.get(kv[0]) or {}
                tested = e.get("fdr") is not None
                return (0 if tested else 1, e.get("fdr") if tested else 0.0,
                        -(e.get("fold") or 0.0), -kv[1], kv[0])
            out += ["", "### Family distribution {#family-distribution}", "",
                    "Cohort families vs a genome-wide background (hypergeometric, "
                    "BH-FDR; fold = observed/expected). Counts kept; sorted by "
                    "enrichment, so the catch-all Other/Unknown bucket no longer leads.",
                    "",
                    table(["Family", "Genes", "Fold", "FDR"],
                          [(k, _i(v),
                            f"{(fe.get(k) or {}).get('fold'):.1f}×" if (fe.get(k) or {}).get("fold") else "—",
                            _fdr((fe.get(k) or {}).get("fdr")))
                           for k, v in sorted(fc.items(), key=_fkey)])]
        else:
            out += ["", "### Family distribution {#family-distribution}", ""]
            rows = sorted(fc.items(), key=lambda kv: -kv[1])
            out.append(table(["Family", "Genes"], [(k, _i(v)) for k, v in rows]))
    fa = b.get("family_assignments") or []
    if fa:
        out += ["", "### Per-gene assignment {#family-assignment}", "",
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
    scm = b.get("sc_marker_gene_count") or 0
    if scm:
        out.append(f"\n**{_i(scm)} cohort gene{'s' if scm != 1 else ''} are a "
                   f"single-cell marker in ≥1 SCXA experiment.**")
    bd = b.get("breadth_distribution") or {}
    if bd and any(bd.values()):              # skip the all-zero distribution (no cohort)
        out += ["", "### Breadth distribution (Bgee present_calls) {#breadth-distribution}", ""]
        order = ["narrow (1-5 tissues)", "moderate (6-20)", "broad (>20)", "unknown"]
        ordered = [(k, bd[k]) for k in order if k in bd]
        out.append(table(["Bucket", "Genes"], [(k, _i(v)) for k, v in ordered]))
    tt = b.get("cohort_tissue_counts") or {}
    if tt:
        out += ["", "### Top tissues across cohort {#cohort-tissues}", ""]
        if isinstance(tt, dict):
            items = list(tt.items())
        else:
            items = list(tt)
        items = sorted(items, key=lambda kv: -kv[1])[:40]
        out.append(table(["Tissue", "Cohort genes"],
                         [(k, _i(v)) for k, v in items]))
    pge = b.get("per_gene_expression") or []
    if pge:
        out += ["", "### Per-gene tissue summary (top 30) {#per-gene-tissue}", "",
                table(["Symbol", "Bgee breadth", "FANTOM5 breadth", "SCXA", "Top tissues"],
                      [(g.get("symbol"), g.get("bgee_breadth"),
                        g.get("fantom5_breadth"),
                        ("marker" if g.get("scxa_marker")
                         else "yes" if g.get("scxa_present") else ""),
                        ", ".join((g.get("top_tissues") or [])[:3]))
                       for g in pge[:30]])]
    return "\n".join(out)


# §8 protein_interactions ---------------------------------------------------

def r_protein_interactions(b):
    out = ["## Protein interactions among cohort", "",
           f"**Intra-cohort edges: {_i(b.get('cohort_edge_count'))}.**"]
    hubs = b.get("hub_genes") or []
    if hubs:
        out += ["", "### Hub genes (top 10 by interactor count) {#hub-genes}", "",
                table(["Symbol", "Interactor count"],
                      [(h.get("symbol"), _i(h.get("interactor_count"))) for h in hubs])]
    edges = b.get("cohort_edges") or []
    if edges:
        out += ["", "### Intra-cohort edges {#intra-cohort-edges}", "",
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
        out += ["", "### Cohort genes with PDB structures (top 30) {#cohort-pdb}", "",
                table(["Symbol", "UniProt", "PDB entries"],
                      [(g.get("symbol"),
                        g.get("uniprot") or "",
                        _i(g.get("pdb_count"))) for g in pdb[:30]])]
    af = b.get("alphafold_only_genes") or []
    if af:
        out += ["", "### AlphaFold-only cohort genes (top 30 by pLDDT) {#cohort-alphafold}", "",
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
        out += ["", "### Genes with an approved drug {#approved-drug-genes}", "",
                "The molecule shown is one approved compound that hits the gene — not "
                "necessarily a drug of choice or one indicated for this disease.", "",
                table(["Symbol", "Example approved molecule"],
                      [(links.maybe_link(g.get("symbol"), links.gene_url(symbol=g.get("symbol"))),
                        g.get("drug") or g.get("top_molecule") or "")
                       for g in ag])]
    tt = b.get("top_targets") or []
    if tt:
        out += ["", "### Top cohort targets by molecule count {#top-targets}", "",
                table(["Symbol", "Molecules", "Max phase"],
                      [(t.get("symbol"), _i(t.get("molecule_count")),
                        t.get("max_phase")) for t in tt])]
    drugs = b.get("drugs") or []
    if drugs:
        out += ["", "### Drugs targeting cohort genes (top 30) {#cohort-drugs}", "",
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
        out += ["", "### Cohort genes with ChEMBL bioactivity (full, sorted by assay count) {#cohort-bioactivity}", "",
                table(["Symbol", "Assays", "Type breakdown"],
                      [(g.get("symbol"), _i(g.get("chembl_assay_total")),
                        ", ".join(f"{k}:{v}" for k, v in (g.get("chembl_assay_types") or {}).items()))
                       for g in studied_all[:50]])]
    eg = b.get("enzyme_genes") or []
    if eg:
        out += ["", "### Cohort enzymes (BRENDA EC) {#cohort-enzymes}", "",
                table(["Symbol", "EC numbers", "Names"],
                      [(g.get("symbol"),
                        ", ".join(g.get("ec_numbers") or []),
                        _trunc(", ".join(g.get("ec_names") or []), 60))
                       for g in eg[:50]])]
    usp = b.get("undrugged_starting_points") or []
    if usp:
        # NB: this is a studied-ness signal, NOT a drugged/undrugged claim — the
        # list includes well-drugged targets (ESR1, EGFR…). Real drugged-status
        # resolution is done in §17 against the approved-drug set; here we only
        # report screening volume and point readers to Therapeutics for status.
        out += ["", "### Cohort genes with high screening signal {#screening-signal}", "",
                "≥100 ChEMBL assays — a studied-ness signal; see Therapeutics for "
                "approved-drug status.", "",
                table(["Symbol", "ChEMBL assays"],
                      [(g.get("symbol"), _i(g.get("chembl_assay_total")))
                       for g in usp[:30]])]
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
        out += ["", "### Cohort genes with a CPIC/DPWG dosing guideline {#cohort-pgx}", "",
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
        out += ["", "### Phase distribution (across all retrieved trials) {#trial-phases}", ""]
        rows = sorted(pc.items(), key=lambda kv: -kv[1])
        out.append(table(["Phase", "Trials"], [(k, _i(v)) for k, v in rows]))
    tt = b.get("top_trials") or []
    if tt:
        out += ["", "### Top trials by phase / activity {#top-trials}", "",
                table(["NCT", "Phase", "Status", "Title"],
                      [(t.get("id") or "",
                        phase_label(t.get("phase")), t.get("status"),
                        _trunc(t.get("title"), 65))
                       for t in tt])]
    td = b.get("trial_drugs") or []
    if td:
        out += ["", "### Drugs tested across these trials (top 30) {#trial-drugs}", "",
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
                "### Precision-medicine subtype map (CIViC) {#civic}",
                "",
                (f"Drug × molecular subtype: {_i(b.get('civic_association_total'))} "
                 f"predictive associations from {_i(b.get('civic_predictive_total'))} "
                 f"curated evidence items"
                 + (f"; also {extra}" if extra else "") + "."),
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

    def _samp(p):  # reconcile the sample with the Genes count (audit: 8 vs 9)
        syms = p.get("gene_symbols") or []
        gc = p.get("gene_count") or len(syms)
        shown = syms[:8]
        extra = gc - len(shown)
        return ", ".join(shown) + (f" (+{extra} more)" if extra > 0 else "")

    def _fdr(q):
        return "—" if q is None else (f"{q:.0e}" if q < 1e-3 else f"{q:.3f}")

    def _ora_block(rows, anchor, heading, col0, n_annot):
        """ORA-ranked table (fold + FDR, sorted by enrichment) when the background
        is present; count-only fallback otherwise. Counts + members kept either
        way (ground-truth for human/agent consumers)."""
        if not rows:
            return []
        if any(p.get("fdr") is not None for p in rows):
            return ["", f"### {heading} {{#{anchor}}}", "",
                    "Over-representation of cohort genes vs the genome-wide "
                    "background (hypergeometric test, Benjamini-Hochberg FDR; fold = "
                    f"observed/expected over {_i(n_annot)} annotated cohort genes). "
                    "Counts and members are kept as ground-truth; sorted by enrichment.",
                    "",
                    table([col0, "Cohort genes", "Fold", "FDR", "Sample cohort genes"],
                          [(p.get("name") or p.get("id") or "", _i(p.get("gene_count")),
                            f"{p['fold']:.1f}×" if p.get("fold") else "—",
                            _fdr(p.get("fdr")), _samp(p)) for p in rows[:30]])]
        return ["", f"### {heading} {{#{anchor}}}", "",
                table([col0, "Genes", "Sample cohort genes"],
                      [(p.get("name") or p.get("id") or "", _i(p.get("gene_count")), _samp(p))
                       for p in rows[:30]])]

    out += _ora_block(b.get("top_pathways") or [], "cohort-pathways",
                      "Pathways by enrichment", "Pathway", gp)
    out += _ora_block(b.get("top_go") or [], "go-enrichment",
                      "GO biological processes by enrichment", "GO term",
                      b.get("genes_with_go") or 0)
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

    # Cohort-empty diseases (no gene cohort) would otherwise render an all-zero
    # 5-row pyramid — skip it; the zone placeholder + Disease family parent
    # pointer cover the "no data here" case.
    if not tiers:
        return ""

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


def _cohort_empty_note(bundles):
    """When no gene cohort resolves, explain (at the top of the Genes zone) that
    the molecular/therapeutic sections are gene-derived and therefore empty —
    so an antibody-mediated / autoimmune / non-gene-defined disease doesn't read
    as broken. Parent pointer is intentionally omitted here (Disease family
    already carries it). '' when the cohort has genes."""
    if (bundles.get("5") or {}).get("gene_count"):
        return ""
    return ("*No associated-gene cohort resolved for this disease. Atlas builds "
            "the molecular and therapeutic sections — associated genes, protein "
            "families, druggability, pathways, interactions, and drug associations "
            "— by aggregating over a disease's associated genes (resolved via "
            "GWAS / GenCC / ClinVar / CIViC), and none resolved here. This is "
            "expected for antibody-mediated, autoimmune, or otherwise "
            "non-gene-defined conditions; the curated evidence for this disease is "
            "its clinical features, GWAS susceptibility, and clinical trials "
            "(above).*")


def r_drugs_indicated(bundles):
    """Drugs DIRECTLY indicated for this disease (ChEMBL drug_indication), from
    the merge-phase indication index (links.indicated_drugs → injected as
    bundles['_indicated_drugs']). Disease-direct therapeutic edge: unlike the
    cohort-target tables below it needs no associated-gene cohort, so autoimmune
    / clinical-only diseases (anti-NMDA aside — off-label drugs aren't ChEMBL
    indications) can still surface real registered drugs.

    TIERED so "indicated" never overstates: phase ≥3 (approved + late-stage) in
    the main table, approved first; phase 2 listed separately as investigational
    candidates (efficacy unproven, ~70% never reach approval). Phase ≤1 is
    excluded upstream. Empty → '' so the Therapeutics zone elides cleanly."""
    rows = bundles.get("_indicated_drugs") or []
    if not rows:
        return ""
    # Tier: only FDA-approved (phase 4) is an INDICATION; phase 2-3 are
    # investigational trials, shown but clearly disclaimed (the aspirin-vs-cancer
    # class). Phase ≤1 isn't in the index.
    approved = [r for r in rows if (r.get("max_phase") or 0) >= 4]
    trials = [r for r in rows if (r.get("max_phase") or 0) in (2, 3)]

    out = ["## Drugs indicated or in trials for this disease", ""]
    if approved:
        out += [f"**{_i(len(approved))} approved drug{'s' if len(approved) != 1 else ''}** "
                "(FDA phase 4) — disease-direct ChEMBL indications, not inferred from the "
                "associated-gene cohort below.", "",
                table(["Drug", "Status"],
                      [(links.maybe_link(r.get("name"), r.get("url")), "Approved (phase 4)")
                       for r in approved])]
    else:
        out += ["No drug is approved (FDA phase 4) with a disease-direct ChEMBL "
                "indication for this disease.", ""]
    if trials:
        out += ["",
                f"**{_i(len(trials))} drug{'s' if len(trials) != 1 else ''} in clinical "
                "trials for this disease (phase 2–3, investigational):** efficacy not "
                "established — a trial record, *not* an indication.", "",
                table(["Drug", "Highest phase"],
                      [(links.maybe_link(r.get("name"), r.get("url")),
                        f"Phase {r.get('max_phase')}")
                       for r in sorted(trials, key=lambda r: -(r.get("max_phase") or 0))])]
    return "\n".join(out)


def render_all(bundles):
    """Disease page body in the FROZEN canonical H2 order (docs/PAGE_CONTRACT.md):
    Clinical features → Identifiers → Genetics & variants → Genes & proteins →
    Function → Therapeutics → Clinical trials & evidence. (Summary wrapped by
    assemble_page; Related appended after.) Clinical features (epidemiology +
    HPO signs/symptoms) leads — it's the headline content for a disease."""
    from atlas.render_common import demote, emit_canonical, with_heading_id

    def S(s, anchor):
        return with_heading_id(demote(RENDER[s](bundles[s])), anchor)

    def D(md, anchor):
        return with_heading_id(demote(md), anchor)

    def join(*parts):
        return "\n\n".join(p for p in parts if p and p.strip())

    # Cohort-derived sections fan out over the disease's associated genes; with no
    # cohort they would render only all-zero "0 cohort genes" leads. Skip them and
    # let the zone fall back to its placeholder (Genes & proteins keeps the
    # explanatory cohort-empty note). Disease-direct content (GWAS, variants,
    # indicated drugs, trials, CIViC) is cohort-independent and always renders.
    has_cohort = bool((bundles.get("5") or {}).get("gene_count"))
    cohort_genes = (join(S("4", "mendelian"), S("5", "cohort-genes"), S("6", "protein-families"),
                         S("7", "expression"), S("8", "interactions"), S("9", "structural"))
                    if has_cohort else "")
    cohort_drugs = (join(S("10", "drug-targets"), S("11", "bioactivity"), S("12", "pharmacogenomics"),
                         D(r_drug_repurposing(bundles), "tractability"),
                         D(r_druggability_pyramid(bundles), "druggability"),
                         D(r_undrugged_target_profiles(bundles), "undrugged"))
                    if has_cohort else "")

    spec = [
        ("Clinical features", "clinical",
         join(D(r_epidemiology(bundles["1"]), "epidemiology"),
              D(r_symptoms(bundles["1"]), "symptoms")),
         "No curated clinical features (Orphanet) for this disease."),
        ("Identifiers", "identifiers", S("1", "disease-ids"), None),
        ("Disease family", "family",
         r_disease_family(bundles.get("1"), bundles.get("5")),
         "No broader Mondo term or subtypes recorded for this disease."),
        ("Genetics & variants", "genetics",
         join(S("2", "gwas"), S("3", "variant-tiers")),
         "No common-variant (GWAS) or curated variant data for this disease."),
        ("Genes & proteins", "genes",
         join(_cohort_empty_note(bundles), cohort_genes),
         "No associated genes curated for this disease."),
        ("Function", "function", S("14", "pathways") if has_cohort else "",
         "No pathway enrichment — requires an associated-gene cohort."),
        ("Therapeutics", "drugs",
         join(D(r_drugs_indicated(bundles), "indicated"), cohort_drugs),
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
