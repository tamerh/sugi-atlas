#!/usr/bin/env python3
"""Deterministic markdown renderer for collector bundles — NO model.

Renders the structured §1-§6 bundles into the same markdown the agentic
pipeline produced, with templates only. Every fact comes verbatim from the
bundle: zero tokens, zero hallucination, perfectly uniform output. The LLM is
reserved for the synthesis/executive-summary layer, not this.

  python3 render.py TP53 2        # collect §2 for TP53 and render it
  python3 render.py TP53 all      # render §1-§6
"""
import sys, os, html
from atlas.gene import collect as C
from atlas.render_common import table


def _cap(n):
    return ""  # pagination works (p= param) — counts are real, not capped at 100


def r_gene_ids(b):
    h = b.get("hgnc", {})
    rows = [("HGNC ID", b.get("hgnc_id")), ("Approved symbol", h.get("symbol")),
            ("Name", h.get("name")), ("Location", h.get("location")),
            ("Locus type", h.get("locus_type")), ("Status", h.get("status")),
            ("Aliases", ", ".join(h.get("aliases", []))),
            ("Ensembl gene", b.get("ensembl_id")),
            ("Ensembl biotype", (b.get("ensembl") or {}).get("biotype")),
            ("OMIM", ", ".join(b.get("mim", []))),
            ("Entrez", ", ".join(b.get("entrez", [])))]
    # RNAcentral row — ncRNA genes only; protein-coding genes have None
    # and the row elides. Builds a clickable RNAcentral URL.
    rc = b.get("rnacentral")
    if rc:
        rid = rc.get("id")
        rt = rc.get("rna_type") or ""
        ln = rc.get("length")
        oc = rc.get("organism_count")
        rows.append(("RNAcentral",
                     f"[{rid}](https://rnacentral.org/rna/{rid}) — {rt}"
                     + (f", {ln} nt" if ln else "")
                     + (f", {oc} organism(s)" if oc else "")))
    return "## Gene identifiers\n\n" + table(["Field", "Value"], rows)


def r_transcripts(b):
    L = ["## Transcript identifiers", ""]
    L.append(f"**Ensembl transcripts: {b.get('ensembl_transcript_count', 0)}**\n")
    L.append(table(["Transcript", "Biotype"],
                   [(t["id"], t.get("biotype")) for t in b.get("ensembl_transcripts", [])]))
    n = b.get("refseq_mrna_count", 0)
    L.append(f"\n**RefSeq mRNA: {n}{_cap(n)}** — MANE Select: `{b.get('mane_select_refseq')}`")
    L.append(", ".join(f"`{x}`" for x in b.get("refseq_mrna", [])))
    L.append(f"\n**CCDS:** " + ", ".join(f"`{x}`" for x in b.get("ccds", [])))
    L.append(f"\n**Canonical transcript `{b.get('canonical_transcript')}` — "
             f"{b.get('canonical_exon_count', 0)} exons**\n")
    L.append(table(["Exon", "Start", "End"],
                   [(e["id"], e.get("start"), e.get("end")) for e in b.get("canonical_exons", [])]))
    return "\n".join(L)


_CC_ORDER = [
    ("function", "Function"),
    ("subunit", "Subunit / interactions"),
    ("subcellular_location", "Subcellular location"),
    ("tissue_specificity", "Tissue specificity"),
    ("ptm", "Post-translational modifications"),
    ("disease", "Disease relevance"),
    ("activity_regulation", "Activity regulation"),
    ("cofactor", "Cofactor"),
    ("domain", "Domain organisation"),
    ("induction", "Induction"),
    ("pathway", "Pathway"),
    ("polymorphism", "Polymorphism"),
    ("miscellaneous", "Miscellaneous"),
    ("similarity", "Similarity"),
]


def r_protein_ids(b):
    L = ["## Protein identifiers", ""]
    pn = b.get("protein_name")
    if pn:
        L.append(f"**{pn}** — `{b.get('canonical_uniprot')}` "
                 f"(reviewed: {', '.join(b.get('reviewed_uniprot', []))})")
    else:
        L.append(f"**Canonical reviewed UniProt:** `{b.get('canonical_uniprot')}`"
                 f" (reviewed: {', '.join(b.get('reviewed_uniprot', []))})")
    alt = b.get("alternative_names") or []
    if alt:
        L.append(f"\n**Alternative names:** " + ", ".join(alt))
    L.append(f"\n**All UniProt accessions ({b.get('uniprot_count', 0)}):** "
             + ", ".join(f"`{x}`" for x in b.get("uniprot_all", [])))

    # UniProt CC narratives (curated, evidence-codes stripped at collector layer).
    # The audit's #1 content gap — biobtree's 2026-05-31 refresh shipped these.
    cc = b.get("cc") or {}
    if cc:
        L.append("\n### UniProt curated annotations\n")
        for key, label in _CC_ORDER:
            text = cc.get(key)
            if not text:
                continue
            L.append(f"**{label}.** {text}")
            L.append("")  # blank between paragraphs

    # Named isoforms — surfaces p53α/β/γ, K-Ras4A/4B, p16INK4a/p14ARF, etc.
    isoforms = b.get("isoforms") or []
    if isoforms:
        L.append(f"\n**Isoforms ({len(isoforms)}):**\n")
        L.append(table(["UniProt ID", "Names", "Canonical?"],
                       [(f"[{iso['id']}](https://www.uniprot.org/uniprotkb/"
                         f"{iso['id'].split('-')[0]}#sequences)" if iso.get("id") else "",
                         ", ".join(iso.get("names") or []),
                         "yes" if iso.get("is_canonical") else "")
                        for iso in isoforms]))

    rp = b.get("refseq_protein", [])
    L.append(f"\n**RefSeq proteins ({b.get('refseq_protein_count', 0)}):** "
             + ", ".join(f"`{p['id']}`" + ("*" if p.get("mane") else "") for p in rp) + "  (*=MANE)")
    L.append("\n**Domains & families (InterPro):**\n")
    L.append(table(["ID", "Name", "Type"],
                   [(d["id"], d.get("name"), d.get("type")) for d in b.get("interpro", [])]))
    L.append(f"\n**Pfam:** " + ", ".join(f"`{x}`" for x in b.get("pfam", [])))
    L.append(f"\n**Antibody resources:** {b.get('antibody_count', 0)}")
    # BRENDA enzyme classification — EC number + name + summary stats.
    # Non-enzyme proteins (TFs, inhibitors) have empty brenda_ec; block elides.
    brenda = b.get("brenda_ec") or []
    if brenda:
        L.append(f"\n**Enzyme classification (BRENDA):**\n")
        for e in brenda:
            ec = e.get("ec") or ""
            L.append(f"- [EC {ec}](https://www.brenda-enzymes.org/enzyme.php?ecno={ec}) — "
                     f"{e.get('name') or ''} *(BRENDA: {e.get('organism_count', 0)} organisms, "
                     f"{e.get('substrate_count', 0)} substrates, {e.get('inhibitor_count', 0)} inhibitors, "
                     f"{e.get('km_count', 0)} Km, {e.get('kcat_count', 0)} kcat entries)*")
    # UniProt sequence features (annotated functional + structural sites)
    uc = b.get("ufeature_counts") or {}
    if uc:
        total = sum(uc.values())
        L.append(f"\n**UniProt features ({total} total):** "
                 + ", ".join(f"{t} {n}" for t, n in sorted(uc.items(), key=lambda x: -x[1])))
        # functional features worth listing with locations (skip bulk sequence-variant
        # — covered in §6 ClinVar — and secondary-structure rows — §4 territory).
        skip = {"sequence variant", "strand", "helix", "turn"}
        feats = [f for f in b.get("ufeatures", []) if f.get("type") not in skip]
        if feats:
            L.append("\n**Annotated functional features (top 25):**\n")
            L.append(table(["Type", "Location", "Description"],
                           [(f["type"], f"{f.get('begin')}–{f.get('end')}",
                             (f.get("description") or "")[:80]) for f in feats[:25]]))
    return "\n".join(L)


def r_structure(b):
    L = ["## Structure", ""]
    n = b.get("pdb_count", 0)
    L.append(f"**Experimental structures (PDB): {n}{_cap(n)}**\n")
    L.append(table(["PDB", "Method", "Resolution (Å)"],
                   [(p["id"], p.get("method"), p.get("resolution")) for p in b.get("pdb", [])]))
    L.append("\n**Predicted structure (AlphaFold):**\n")
    afs = b.get("alphafold", []) or []
    present_rows = [a for a in afs if a.get("present", True) and a.get("id")]
    missing_rows = [a for a in afs if not a.get("present", True)]
    if present_rows:
        L.append(table(["Model", "pLDDT", "Fraction very-high"],
                       [(f"[{a['id']}](https://alphafold.ebi.ac.uk/entry/{a['uniprot']})",
                         a.get("plddt"), a.get("fraction_plddt_very_high"))
                        for a in present_rows]))
    for a in missing_rows:
        L.append(f"\n*No AlphaFold model available for {a['uniprot']} — "
                 "AlphaFold DB does not currently provide models for proteins "
                 "above ~3000 aa.*")
    return "\n".join(L)


def r_orthologs(b):
    L = ["## Cross-species orthologs", "", f"**{b.get('ortholog_count', 0)} orthologs**\n"]
    L.append(table(["Organism", "Symbol", "Gene ID"],
                   [(o.get("organism"), o.get("symbol"), o["id"]) for o in b.get("orthologs", [])]))
    para = b.get("paralogs", [])
    if para:
        L.append(f"\n**Paralogs ({b.get('paralog_count', 0)}):** "
                 + ", ".join(f"{p.get('symbol')} ({p['id']})" for p in para[:20]))
    return "\n".join(L)


def r_variants(b):
    L = ["## Clinical variants & AI predictions", ""]
    bd = b.get("clinvar_breakdown", {})
    L.append(f"**ClinVar: {b.get('clinvar_total', 0)} variants total.** "
             f"Per-class counts are floors (≥ shown; pagination cap):\n")
    L.append(table(["Classification", "Count (floor)"], list(bd.items())))
    L.append(f"\n**Top pathogenic / likely-pathogenic ({len(b.get('top_pathogenic', []))}):**\n")
    L.append(table(["Variant ID", "HGVS", "Classification"],
                   [(v["id"], v.get("hgvs"), v.get("classification")) for v in b.get("top_pathogenic", [])]))
    L.append(f"\n**SpliceAI: {b.get('spliceai_total', 0)} predictions.** Top by Δscore:\n")
    L.append(table(["Variant", "Effect", "Δscore"],
                   [(v["id"], v.get("effect"), v.get("score")) for v in b.get("top_spliceai", [])]))
    L.append(f"\n**AlphaMissense: {b.get('alphamissense_total', 0)} scored.** "
             f"Top likely-pathogenic:\n")
    L.append(table(["Variant", "Protein change", "am_pathogenicity"],
                   [(v["id"], v.get("variant"), v.get("am_pathogenicity")) for v in b.get("top_alphamissense", [])]))
    ds = b.get("dbsnp_sample", [])
    if ds:
        L.append(f"\n**dbSNP variants (sampled {b.get('dbsnp_sampled', 0)} via entrez):** "
                 + ", ".join(f"{d['id']} ({d['pos']} {d['change']})" for d in ds[:15]))
    return "\n".join(L)


def r_pathways(b):
    L = ["## Pathways & Gene Ontology", "",
         f"**Reactome pathways: {b.get('reactome_count', 0)}**\n"]
    L.append(table(["ID", "Pathway"], [(p["id"], p.get("name")) for p in b.get("reactome", [])[:30]]))
    L.append(f"\n**MSigDB gene sets: {b.get('msigdb_total', 0)}** (showing top):")
    L.append(", ".join(f"`{m['name']}`" for m in b.get("msigdb", [])[:15]))
    go = b.get("go", {})
    for cat in ("biological_process", "molecular_function", "cellular_component"):
        terms = go.get(cat, [])
        L.append(f"\n**GO {cat.replace('_', ' ').title()} ({len(terms)}):**")
        L.append(", ".join(f"{t['name']} ({t['id']})" for t in terms[:20]))

    # Top-level parent rollups — give the page a hierarchical-navigation view.
    # Reactome's hierarchy is tight (1-2 parents per pathway); GO's is broader.
    rp = b.get("reactome_parent_rollup") or []
    if rp:
        L.append(f"\n**Reactome top-level categories (rollup of top-{len(rp)} pathways):**\n")
        L.append(table(["Category", "Pathways"],
                       [(p.get("name") or p["id"], p.get("pathway_count")) for p in rp[:15]]))
    gp = b.get("go_parent_rollup") or []
    if gp:
        L.append(f"\n**GO top-level categories (rollup of top GO terms by namespace):**\n")
        L.append(table(["Category", "Terms"],
                       [(p.get("name") or p["id"], p.get("term_count")) for p in gp[:20]]))
    return "\n".join(L)


def r_interactions(b):
    L = ["## Protein interactions & networks", ""]
    L.append(f"**STRING ({b.get('string_count', 0)}), top by confidence (×1000):**\n")
    L.append(table(["Partner", "Score"], [(s.get("partner"), s.get("score")) for s in b.get("string", [])[:20]]))
    L.append(f"\n**IntAct ({b.get('intact_count', 0)}), top by confidence:**\n")
    L.append(table(["A", "B", "Type", "Score"],
                   [(i.get("a"), i.get("b"), i.get("type"), i.get("score")) for i in b.get("intact", [])[:20]]))
    L.append(f"\n**BioGRID ({b.get('biogrid_count', 0)}):** "
             + ", ".join(f"{x.get('partner')} ({x.get('method')})" for x in b.get("biogrid", [])[:15]))
    L.append(f"\n**ESM2 similar proteins:** " + ", ".join(f"`{p}`" for p in b.get("esm2_similar", [])[:20]))
    L.append(f"\n**Diamond homologs:** " + ", ".join(f"`{p}`" for p in b.get("diamond_similar", [])[:20]))
    L.append(f"\n**SIGNOR signaling ({b.get('signor_count', 0)}):**\n")
    L.append(table(["A", "Effect", "B", "Mechanism"],
                   [(s.get("a"), s.get("effect"), s.get("b"), s.get("mechanism")) for s in b.get("signor", [])[:30]]))
    return "\n".join(L)


def r_tf_regulation(b):
    L = ["## Regulation", "",
         f"**Is transcription factor: {b.get('is_transcription_factor')}**\n"]
    L.append(f"**Downstream targets (CollecTRI): {b.get('downstream_count', 0)}**\n")
    L.append(table(["Target", "Regulation"],
                   [(t.get("target"), t.get("regulation")) for t in b.get("downstream_targets", [])[:30]]))
    L.append("\n**JASPAR motifs:**\n")
    L.append(table(["Motif", "Name", "Family"],
                   [(m["id"], m.get("name"), m.get("family")) for m in b.get("jaspar_motifs", [])]))
    # JASPAR PMIDs — evidence trail for the motifs above.
    pmids = b.get("jaspar_pmids") or []
    if pmids:
        links = ", ".join(f"[PubMed:{p}](https://pubmed.ncbi.nlm.nih.gov/{p}/)"
                          for p in pmids[:10])
        L.append(f"\n*JASPAR matrix evidence (PMIDs):* {links}")
    L.append(f"\n**Upstream regulators (CollecTRI, top):** "
             + ", ".join(r.get("regulator") for r in b.get("upstream_regulators", [])[:20]))
    # miRDB miRNAs targeting this gene — post-transcriptional regulators.
    # Sorted by max_score (miRDB confidence). target_count is the miRNA's
    # promiscuity across all genes (lower = more specific target relationship).
    n_mir = b.get("mirna_count", 0)
    if n_mir:
        L.append(f"\n**miRNA regulators (miRDB): {n_mir}** targeting "
                 f"{b.get('symbol')}, top 30 by miRDB confidence (max_score; "
                 f"target_count = how many genes the miRNA targets in total — "
                 f"lower means more specific):\n")
        L.append(table(["miRNA", "Max score", "Avg score", "miRNA target_count"],
                       [(m["id"], m.get("max_score"), m.get("avg_score"),
                         m.get("target_count")) for m in b.get("mirna_regulators", [])]))
    return "\n".join(L)


def r_drugs(b):
    L = ["## Drug & pharmacology data", "",
         f"**Is drug target: {b.get('is_drug_target')}**\n"]
    L.append(f"**ChEMBL targets ({len(b.get('chembl_targets', []))}):** "
             + ", ".join(f"{t['id']} ({t.get('type')})" for t in b.get("chembl_targets", [])[:10]))
    pt = b.get("patent_total", 0)
    patent_note = (f" Patent mentions across the top 20 by phase: **{pt:,}** "
                   f"(via chembl_molecule>>patent_compound — counts attach to the compound, "
                   f"not the gene–compound relationship, so off-target/promiscuous "
                   f"molecules can dominate). " if pt else "")
    L.append(f"\n**Molecules with ChEMBL bioactivity (phase ≥1): {b.get('molecule_count', 0)}**, "
             f"by development phase (incl. off-target/promiscuous compounds).{patent_note}\n")
    L.append(table(["Molecule", "Name", "Phase", "Patents"],
                   [(m["id"], m.get("name"), m.get("phase"),
                     f"{m['patent_count']:,}" if m.get("patent_count") else "")
                    for m in b.get("molecules", [])[:30]]))
    pg = b.get("pharmgkb", [])
    L.append(f"\n**PharmGKB:** {len(pg)} entr{'y' if len(pg)==1 else 'ies'}"
             + (f" (VIP={pg[0].get('vip')}, CPIC={pg[0].get('cpic_guideline')})" if pg else ""))
    bd = b.get("bindingdb_sample", [])
    if bd:
        L.append(f"\n**Binding affinities (BindingDB, sampled {b.get('bindingdb_sampled', 0)}):**\n")
        L.append(table(["Ligand", "Ki", "IC50"],
                       [(x.get("ligand"), x.get("ki"), x.get("ic50")) for x in bd[:15]]))

    # ChEMBL bioactivities (pchembl-ranked). pchembl is the gold potency
    # metric — directly comparable across assay types. Renders even when
    # PubChem activity is empty (e.g. KRAS — see BIOBTREE_ISSUES.md #12).
    ca = b.get("chembl_activities") or []
    if ca:
        L.append(f"\n**ChEMBL bioactivities ({b.get('chembl_activity_potent_count', 0)} "
                 f"potent at pChembl≥5 of {b.get('chembl_activity_total', 0)} total), "
                 f"top 30 by pChembl (potency: 10 = 0.1 nM, 6 = 1 µM):**\n")
        L.append(table(["pChembl", "Type", "Value", "Unit", "Activity ID"],
                       [(r.get("pchembl"), r.get("type"), r.get("value"),
                         r.get("unit"), r["id"]) for r in ca]))

    # PubChem BioAssay actives — sorted by potency. CID/AID get clickable
    # PubChem URLs so an AI agent (or human) can drill into the assay record
    # directly without needing a separate entry fetch.
    pba = b.get("pubchem_bioassay") or []
    if pba:
        L.append(f"\n**PubChem BioAssay actives ({b.get('pubchem_bioassay_active_count', 0)} "
                 f"Active w/ measured affinity, of {b.get('pubchem_bioassay_total', 0)} total "
                 f"PubChem activities), top 30 by potency:**\n")
        L.append(table(["CID", "AID", "Type", "Value", "Unit"],
                       [(f"[{p['cid']}](https://pubchem.ncbi.nlm.nih.gov/compound/{p['cid']})" if p.get('cid') else '',
                         f"[{p['aid']}](https://pubchem.ncbi.nlm.nih.gov/bioassay/{p['aid']})" if p.get('aid') else '',
                         p.get("activity_type"), p.get("value"), p.get("unit"))
                        for p in pba]))

    # CTD literature-mined chemical-gene interactions — Comparative
    # Toxicogenomics Database. Each row: a chemical (drug, toxin,
    # environmental compound) + CV-coded action verbs + PubMed-count support.
    # High AI value: every claim is anchored by literature counts.
    ctd = b.get("ctd_interactions") or []
    if ctd:
        L.append(f"\n**CTD chemical–gene interactions (human, "
                 f"{b.get('ctd_interaction_total', 0)} total), top 30 by PubMed support:**\n")
        L.append(table(["Chemical", "Actions (CV verbs)", "PubMed papers"],
                       [(f"[{r['chemical']}](http://ctdbase.org/detail.go?type=chem&acc={r['chemical_id']})",
                         ", ".join(r["actions"]),
                         r["pmids"]) for r in ctd]))

    # ChEMBL screening-assay depth — how heavily this target has been
    # profiled, and via what experimental style (Binding / Functional / ADMET
    # / Toxicity). Aggregated across all ChEMBL targets for this protein.
    cat = b.get("chembl_assay_total") or 0
    if cat:
        type_counts = b.get("chembl_assay_type_counts") or {}
        breakdown = ", ".join(f"{n} {k.lower()}" for k, n in type_counts.items())
        L.append(f"\n**ChEMBL screening assays ({cat} unique, capped per target):** "
                 f"{breakdown}")
        samples = b.get("chembl_assay_samples") or []
        if samples:
            L.append("\nRepresentative assays (with source publication via chembl_document):\n")
            rows = []
            for s in samples:
                doc = ""
                if s.get("doc_id"):
                    title = (s.get("doc_title") or s["doc_id"]).strip()
                    if len(title) > 90:
                        title = title[:87] + "…"
                    doc = (f"[{title}](https://www.ebi.ac.uk/chembl/document_report_card/{s['doc_id']}/)"
                           + (f" — *{s['doc_journal']}*" if s.get("doc_journal") else ""))
                rows.append((
                    f"[{s['id']}](https://www.ebi.ac.uk/chembl/assay_report_card/{s['id']}/)",
                    s["type"], s["desc"], doc))
            L.append(table(["Assay ID", "Type", "Description", "Source paper"], rows))

    # Cellosaurus — cell lines associated with this gene (mutated, deficient,
    # expressed-in, model-of). Counts can run into thousands for heavily-
    # studied tumor suppressors / oncogenes. Sample list is the first 10
    # returned by biobtree (id-ordered, not curated) — meant as a starting
    # point for downstream lookups, not a "top 10".
    cct = b.get("cellosaurus_total") or 0
    if cct:
        cats = b.get("cellosaurus_category_counts") or {}
        breakdown = ", ".join(f"{n:,} {k.lower()}" for k, n in cats.items())
        L.append(f"\n**Cellosaurus cell lines ({cct:,}):** {breakdown}")
        cs = b.get("cellosaurus_samples") or []
        if cs:
            L.append("\nFirst 10 cell lines (id-ordered, not curated):\n")
            L.append(table(["Cellosaurus", "Name", "Category", "Sex"],
                           [(f"[{c['id']}](https://www.cellosaurus.org/{c['id']})",
                             c.get("name"), c.get("category"), c.get("sex")) for c in cs]))

    ct = b.get("disease_trials", [])
    L.append(f"\n**Clinical trials for the gene's associated diseases "
             f"({b.get('disease_trial_count', 0)}, via MONDO — disease-level, not drug-specific):**\n")
    L.append(table(["Trial", "Phase", "Status", "Title"],
                   [(t["id"], t.get("phase"), t.get("status"), (t.get("title") or "")[:55]) for t in ct[:20]]))
    return "\n".join(L)


def r_expression(b):
    L = ["## Expression profiles", ""]
    bg = b.get("bgee")
    if bg:
        L.append(f"**Bgee:** expression breadth **{bg.get('breadth')}**, "
                 f"{bg.get('present_calls')} present calls, max score {bg.get('max_expression_score')}.")
    f5 = b.get("fantom5")
    if f5:
        L.append(f"\n**FANTOM5 (CAGE):** breadth **{f5.get('breadth')}**, "
                 f"TPM avg {f5.get('tpm_average')} / max {f5.get('tpm_max')}, "
                 f"expressed in {f5.get('samples_expressed')} samples.")
    fp = b.get("fantom5_promoters") or []
    if fp:
        L.append(f"\n**FANTOM5 promoters ({len(fp)} alternative TSS):**\n")
        L.append(table(["Promoter ID", "TPM avg", "Samples expressed"],
                       [(p["id"], p.get("tpm_average"), p.get("samples_expressed")) for p in fp[:10]]))
    # Tissue/cell name -> federated link via bioregistry (resolves UBERON or
    # CL ids to their authority page). Bare names without an anatomy_id
    # (rare) render unwrapped.
    def _t_link(t):
        name = t.get("tissue") or ""
        aid = t.get("anatomy_id")
        return f"[{name}](https://bioregistry.io/{aid})" if aid and name else name
    L.append(f"\n**Top tissues by expression ({b.get('tissue_count', 0)} total):**\n")
    L.append(table(["Tissue", "Anatomy ID", "Score", "Rank", "Quality"],
                   [(_t_link(t), t.get("anatomy_id") or "",
                     t.get("score"), t.get("rank"), t.get("quality"))
                    for t in b.get("top_tissues", [])[:30]]))
    sc = b.get("single_cell_datasets", [])
    L.append(f"\n**Single-cell datasets ({len(sc)}):**\n")
    L.append(table(["Dataset", "Description", "Cells"],
                   [(d["id"], (d.get("description") or "")[:60], d.get("cells")) for d in sc[:15]]))
    return "\n".join(L)


def r_cancer_overview(bundle):
    """Cancer-significance block — emitted between the page lead and the
    LLM exec summary IFF the gene has CIViC and/or intOGen data. Folds the
    two most-grounded cancer signals into a single early narrative block,
    which is what AI agents extract as the headline interpretation.

    Returns "" for non-cancer genes so the block elides cleanly.
    """
    b12 = bundle.get("12") or {}
    civic = b12.get("civic")
    intogen = b12.get("intogen")
    if not civic and not intogen:
        return ""
    L = ["## Cancer significance"]
    if civic:
        desc = (civic.get("description") or "").strip()
        if desc:
            L.append(f"*From [CIViC](https://civicdb.org/genes/{civic.get('id')}/summary) — curated cancer-variant interpretation:*")
            L.append(desc)
    if intogen:
        role = intogen.get("role") or ""
        role_human = {
            "Act": "**activating** (oncogene-like)",
            "LoF": "**loss-of-function** (tumor-suppressor-like)",
            "ambiguous": "**ambiguous** (mixed evidence)",
        }.get(role, f"**{role}**")
        cancers = [c for c in (intogen.get("cancer_types") or "").split(",") if c]
        head = (f"*From [intOGen](https://www.intogen.org/search?gene={intogen.get('symbol')}) "
                f"— cancer-driver classification:* {role_human} "
                f"across **{len(cancers)} cancer types**")
        if cancers:
            shown = ", ".join(cancers[:12]) + (f"…(+{len(cancers)-12} more)" if len(cancers) > 12 else "")
            head += f" — {shown}"
        L.append(head + ".")
    return "\n\n".join(L)


def r_diseases(b):
    L = ["## Disease associations", ""]
    L.append(f"**OMIM:** gene `{', '.join(b.get('gene_omim', []))}` | "
             f"disease phenotypes: {', '.join(b.get('disease_omim', [])[:20])}")
    L.append("\n**GenCC curated gene-disease:**\n")
    L.append(table(["Disease", "Classification", "Inheritance"],
                   [(g.get("disease"), g.get("classification"), g.get("inheritance")) for g in b.get("gencc", [])[:20]]))
    L.append(f"\n**Mondo ({len(b.get('mondo', []))}):** "
             + ", ".join(f"{m.get('name')} ({m['id']})" for m in b.get("mondo", [])[:15]))
    L.append(f"\n**Orphanet ({len(b.get('orphanet', []))}):** "
             + ", ".join(f"{o.get('name') or ''} ({o['id']})" for o in b.get("orphanet", [])[:15]))
    L.append(f"\n**HPO phenotypes: {b.get('hpo_total', 0)}** (top):\n")
    L.append(table(["HPO", "Term"], [(h["id"], h.get("name")) for h in b.get("hpo", [])[:30]]))
    L.append(f"\n**GWAS associations: {b.get('gwas_total', 0)}** (top):\n")
    L.append(table(["Study", "Trait", "p-value"],
                   [(g["id"], g.get("trait"), g.get("p_value")) for g in b.get("gwas", [])[:30]]))

    # EFO canonical trait names — normalized vocabulary for the free-text
    # GWAS-catalog traits above. The same disease can appear under multiple
    # synonymous GWAS-catalog spellings; EFO collapses them. Useful for
    # AI agents doing cross-resource joins on trait IDs.
    efo = b.get("efo_traits") or []
    if efo:
        L.append(f"\n**EFO canonical traits ({len(efo)}, from GWAS):**\n")
        L.append(table(["EFO ID", "Trait name"],
                       [(f"[{e['id']}](https://www.ebi.ac.uk/efo/{e['id'].replace(':','_')})",
                         e.get("name")) for e in efo[:30]]))

    # MeSH disease descriptors — NLM's controlled disease vocabulary.
    # Tree numbers (e.g. C04.700.600) classify into MeSH categories
    # (C04=Neoplasms, C16=Congenital, C18=Nutritional/Metabolic, etc.) —
    # useful for grouping diseases at a coarser level.
    mesh = b.get("mesh_descriptors") or []
    if mesh:
        L.append(f"\n**MeSH disease descriptors ({len(mesh)}):**\n")
        L.append(table(["Descriptor", "Name", "Tree numbers"],
                       [(f"[{m['id']}](https://www.ncbi.nlm.nih.gov/mesh/?term={m['id']})",
                         m["name"] + (" *(supp.)*" if m.get("is_supplementary") else ""),
                         "; ".join(m.get("tree_numbers") or []))
                        for m in mesh[:30]]))
    return "\n".join(L)


RENDER = {"1": r_gene_ids, "2": r_transcripts, "3": r_protein_ids,
          "4": r_structure, "5": r_orthologs, "6": r_variants,
          "7": r_pathways, "8": r_interactions, "9": r_tf_regulation,
          "10": r_drugs, "11": r_expression, "12": r_diseases}

if __name__ == "__main__":
    sym = sys.argv[1] if len(sys.argv) > 1 else "TP53"
    sec = sys.argv[2] if len(sys.argv) > 2 else "1"
    secs = list(RENDER) if sec == "all" else [sec]
    for s in secs:
        print(RENDER[s](C.SECTIONS[s](sym)))
        print()
