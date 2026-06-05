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
from atlas.render_common import table, phase_label, fnum
from atlas.civic import therapy_label
from atlas.page import links


def _cap(n):
    return ""  # pagination works (p= param) — counts are real, not capped at 100


def _labeled(label, items):
    """'**label:** a, b, c' — or '' when items is empty, so an empty list doesn't
    leave a dangling '**label:** ' with no content (audit: lncRNA gene pages)."""
    items = [x for x in items]
    return (f"\n**{label}:** " + ", ".join(items)) if items else ""


def _ctd_actions(actions, cap=5):
    """CTD encodes each action as `verb^object` ("increases^expression"); a
    chemical can carry 20+, producing an unreadable wall. Humanize the caret and
    cap the chain (the full set is in the bundle/JSON sidecar)."""
    shown = [str(a).replace("^", " ") for a in (actions or [])[:cap]]
    extra = len(actions or []) - len(shown)
    return ", ".join(shown) + (f" (+{extra} more)" if extra > 0 else "")


def _dedup_sentences(text):
    """Collapse duplicate period-separated segments (UniProt subcellular-location
    CC concatenates per-isoform notes, which repeat the same location many times)."""
    parts = [p.strip().rstrip(".") for p in (text or "").split(". ")]
    seen, out = set(), []
    for p in parts:
        k = p.strip().lower()
        if not k or k in seen:
            continue
        seen.add(k)
        out.append(p)
    return (". ".join(out) + ".") if out else (text or "")


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
                     f"{rid} — {rt}"
                     + (f", {ln} nt" if ln else "")
                     + (f", {oc} organism(s)" if oc else "")))
    # Skip empty-value rows (e.g. "| Aliases |  |", "| OMIM |  |" on lncRNA genes)
    # — the disease/drug identifier tables already do this.
    rows = [(k, v) for k, v in rows if v not in (None, "")]
    return "## Gene identifiers\n\n" + table(["Field", "Value"], rows)


def r_transcripts(b):
    L = ["## Transcript identifiers", ""]
    # Inline ID list + biotype breakdown (consistent with the RefSeq lists
    # below — an ID list with one mostly-uniform attribute doesn't need a table).
    ts = b.get("ensembl_transcripts", [])
    bc = {}
    for t in ts:
        bt = t.get("biotype") or "?"
        bc[bt] = bc.get(bt, 0) + 1
    breakdown = ", ".join(f"{c} {bt}" for bt, c in sorted(bc.items(), key=lambda kv: -kv[1]))
    L.append(f"**Ensembl transcripts: {b.get('ensembl_transcript_count', 0)}**"
             + (f" — {breakdown}" if breakdown else ""))
    if ts:
        L.append("\n" + ", ".join(f"`{t['id']}`" for t in ts))
    n = b.get("refseq_mrna_count", 0)
    L.append(f"\n**RefSeq mRNA: {n}{_cap(n)}** — MANE Select: `{b.get('mane_select_refseq')}`")
    L.append(", ".join(f"`{x}`" for x in b.get("refseq_mrna", [])))
    L.append(_labeled("CCDS", (f"`{x}`" for x in b.get("ccds", []))))
    L.append("\n### Canonical transcript exons {#canonical-exons}\n")
    L.append(f"`{b.get('canonical_transcript')}` — {b.get('canonical_exon_count', 0)} exons\n")
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
    # Emitted in full: the markdown carries the complete narrative and the frontend
    # clamps for display (web-team request). A single section-level deep link to
    # the full UniProt entry remains as the source pointer.
    cc = b.get("cc") or {}
    if cc:
        canon = b.get("canonical_uniprot")
        uurl = f"https://www.uniprot.org/uniprotkb/{canon}/entry" if canon else None
        # Single deep link in the subtitle (not one per truncated block).
        head = "### UniProt curated annotations"
        if uurl:
            head += f" — [full annotation on UniProt →]({uurl})"
        head += " {#uniprot-cc}"
        L.append("\n" + head + "\n")
        for key, label in _CC_ORDER:
            text = cc.get(key)
            if not text:
                continue
            # Subcellular-location CC concatenates per-isoform notes, producing
            # run-ons ("Cytoplasm. Nucleus. Nucleus. … Cytoplasm Nucleus. ×6").
            # Collapse duplicate period-separated segments.
            if key == "subcellular_location":
                text = _dedup_sentences(text)
            # Full CC narrative: the section-level UniProt link above is the
            # source pointer, and display clamping is the frontend's job, so the
            # markdown carries the complete text (no generator-side truncation).
            L.append(f"**{label}.** {text.strip()}")
            L.append("")  # blank between paragraphs

    # (NCBI-curated gene summary moved to the top-of-page Overview, and the
    #  gene-level fitness/dosage + GeneRIF blocks moved out to their own
    #  sections — see r_functional_genomics / r_generifs. §3 is protein
    #  identifiers + protein-scoped annotation only.)

    # Named isoforms — surfaces p53α/β/γ, K-Ras4A/4B, p16INK4a/p14ARF, etc.
    isoforms = b.get("isoforms") or []
    if isoforms:
        L.append(f"\n### Isoforms ({len(isoforms)}) {{#isoforms}}\n")
        L.append(table(["UniProt ID", "Names", "Canonical?"],
                       [(iso.get("id") or "",
                         ", ".join(iso.get("names") or []),
                         "yes" if iso.get("is_canonical") else "")
                        for iso in isoforms]))

    rp = b.get("refseq_protein", [])
    L.append(f"\n**RefSeq proteins ({b.get('refseq_protein_count', 0)}):** "
             + ", ".join(f"`{p['id']}`" + ("*" if p.get("mane") else "") for p in rp) + "  (*=MANE)")
    L.append("\n### Domains & families (InterPro) {#domains}\n")
    L.append(table(["ID", "Name", "Type"],
                   [(d["id"], d.get("name"), d.get("type")) for d in b.get("interpro", [])]))
    L.append(_labeled("Pfam", (f"`{x}`" for x in b.get("pfam", []))))
    # Therapeutic antibodies targeting this protein (TheraSAbDab/SAbDab via the
    # antibody dataset). Only shown when present — the reverse edge is sparse, so
    # a universal "0" otherwise reads as a research-reagent count, which it isn't.
    if b.get("antibody_count"):
        L.append(f"\n**Therapeutic antibodies (targeting this protein):** "
                 f"{b.get('antibody_count')}")
    # BRENDA enzyme classification — EC number + name + summary stats.
    # Non-enzyme proteins (TFs, inhibitors) have empty brenda_ec; block elides.
    brenda = b.get("brenda_ec") or []
    if brenda:
        L.append(f"\n**Enzyme classification (BRENDA):**\n")
        for e in brenda:
            ec = e.get("ec") or ""
            L.append(f"- EC {ec} — "
                     f"{e.get('name') or ''} *(BRENDA: {e.get('organism_count', 0)} organisms, "
                     f"{e.get('substrate_count', 0)} substrates, {e.get('inhibitor_count', 0)} inhibitors, "
                     f"{e.get('km_count', 0)} Km, {e.get('kcat_count', 0)} kcat entries)*")
        # Per-substrate kinetics (Tier 2: Km spans organisms/conditions).
        kin = b.get("brenda_kinetics") or []
        if kin:
            L.append("\n### Substrate kinetics (BRENDA) {#brenda}\n")
            L.append(f"{b.get('brenda_kinetics_total', 0)} substrates with measured Km, "
                     f"best-characterized {len(kin)}. *Km ranges are aggregated across "
                     f"organisms/conditions.*\n")

            def _km(r):
                def f(v):
                    try:
                        x = float(v)
                        return x if x > 0 else None     # 0 = missing-value floor
                    except (TypeError, ValueError):
                        return None
                lo, hi = f(r.get("min_km")), f(r.get("max_km"))
                # Show a value only with a clean lower bound — a max-only figure
                # from a 0-floored aggregate (DHF 21.7 mM) would misrepresent a
                # µM-affinity substrate. The measurement count still conveys depth.
                if not lo:
                    return "—"
                return f"{lo:g}–{hi:g}" if (hi and hi != lo) else f"{lo:g}"
            L.append(table(["Substrate", "Km (mM)", "Measurements"],
                           [((k.get("substrate") or "")[:48], _km(k), k.get("km_count"))
                            for k in kin]))
    # Catalyzed reactions (Rhea) — the enzyme mechanism as substrate → product.
    # Independent of BRENDA (a protein can carry Rhea reactions with no BRENDA
    # EC row). biobtree now projects the human-readable equation (#29).
    rhea = b.get("rhea_reactions") or []
    if rhea:
        L.append(f"\n**Catalyzed reactions (Rhea), {len(rhea)} shown:**\n")
        for r in rhea:
            rid = (r.get("id") or "").replace("RHEA:", "")
            L.append(f"- {r.get('equation')}" + (f" *(RHEA:{rid})*" if rid else ""))
    # UniProt sequence features — one-line census here; the structured
    # druggable-residue breakdown lives in r_residue_map (the protein zone),
    # so the former flat "top 40 features" table is no longer duplicated here.
    uc = b.get("ufeature_counts") or {}
    if uc:
        total = sum(uc.values())
        L.append(f"\n**UniProt features ({total} total):** "
                 + ", ".join(f"{t} {n}" for t, n in sorted(uc.items(), key=lambda x: -x[1])))
    return "\n".join(L)


def r_functional_genomics(b):
    """ClinGen dosage sensitivity + DepMap CRISPR fitness — gene-level signals
    for variant interpretation (dosage) and drug-target prioritization (depmap).
    Carved out of §3 (these are gene-level, not protein identifiers). The
    one-line verdicts also surface in the top-of-page "At a glance" block;
    this section carries the full scale + provenance links. Returns "" when
    neither signal is present."""
    cd = b.get("clingen_dosage") or {}
    dm = b.get("depmap") or {}
    # DepMap only counts as content when notable — a low pct_dependent (e.g.
    # 0.2%) just means "not a dependency" and would leave a lone, low-value
    # one-line section. Show it only if common-essential / strongly-selective
    # / ≥10% (mirrors the "At a glance" gate).
    pct = dm.get("pct_dependent", "")
    try:
        dm_notable = float(pct) >= 10
    except (TypeError, ValueError):
        dm_notable = False
    dm_notable = (dm_notable or dm.get("strongly_selective") == "true"
                  or dm.get("common_essential") == "true")
    if not cd and not dm_notable:
        return ""
    L = ["## Functional genomics", ""]
    if cd:
        haplo = cd.get("haplo_score", "")
        triplo = cd.get("triplo_score", "")
        # ClinGen scale: 3=sufficient evidence, 2=emerging, 1=little, 0=no evidence,
        # 30=AR (recessive — loss of one copy doesn't cause disease),
        # 40=dosage sensitivity unlikely.
        scale = {"3": "sufficient evidence",
                 "2": "emerging evidence",
                 "1": "little evidence",
                 "0": "no evidence",
                 "30": "autosomal recessive",
                 "40": "dosage sensitivity unlikely"}
        L.append(f"**ClinGen dosage:** haploinsufficiency `{haplo}` "
                 f"({scale.get(haplo, 'unscored')}), "
                 f"triplosensitivity `{triplo}` ({scale.get(triplo, 'unscored')}). "
                 "[ClinGen Gene Dosage Map](https://search.clinicalgenome.org/kb/gene-dosage)")
    if dm_notable:
        sel = dm.get("strongly_selective", "")
        ce = dm.get("common_essential", "")
        L.append(f"**DepMap (CRISPR cell-line fitness):** "
                 f"dependent in {pct}% of screened cell lines"
                 + (", strongly selective" if sel == "true" else "")
                 + (", common-essential" if ce == "true" else "")
                 + ".")
    return "\n".join(L)


def r_generifs(b):
    """GeneRIFs — NCBI per-gene PMID-anchored claims (top 30 of N).
    Citation-grounded gene knowledge, dense with disease + mechanism + clinical
    context. Carved out of §3 (gene-level literature, not protein IDs).
    Returns "" when none present."""
    rifs = b.get("generifs") or []
    if not rifs:
        return ""
    L = [f"## Literature-anchored findings (GeneRIF, showing {len(rifs)})", ""]
    for r in rifs:
        # Escape asterisks: GeneRIF prose occasionally has literal '**' (e.g.
        # "-Leu**Pro-"), which Markdown would read as an (unbalanced) bold marker.
        text = (r.get("text") or "").strip().replace("*", r"\*")
        pmid = r.get("pmid")
        cite = f" (PMID:{pmid})" if pmid else ""
        L.append(f"- {text}{cite}")
    return "\n".join(L)


# Functional-residue map (audit-enrichment #1 / MOLECULAR_ENRICHMENT layer B):
# regroup UniProt sequence features into a drug-discovery view. (group label,
# feature types, show per-residue description?)
_RESIDUE_GROUPS = [
    ("Catalytic / active sites", ("active site", "site"), True),
    ("Ligand- & substrate-binding residues", ("binding site",), True),
    ("Post-translational modifications", ("modified residue",
     "lipid moiety-binding region", "cross-link"), False),
    ("Disulfide bonds", ("disulfide bond",), False),
    ("Glycosylation sites", ("glycosylation site",), False),
]


def _res_loc(f):
    b, e = f.get("begin"), f.get("end")
    return str(b) if (e in (None, "", b)) else f"{b}–{e}"


def r_residue_map(b):
    """Structure-function residue view over UniProt features — catalytic,
    ligand-binding, PTM, disulfide, glycosylation, and mutagenesis-validated
    positions, grouped per reviewed product. Pure restructure of data §3 already
    collects (each feature is accession-stamped); '' when no features."""
    feats = b.get("ufeatures") or []
    canon = b.get("canonical_uniprot")

    def _product_block(ufe):
        """Residue-group lines for one product; [] when it has no mappable
        residues (avoids an orphan product header)."""
        lines = []
        for label, types, show_desc in _RESIDUE_GROUPS:
            rows = [f for f in ufe if f.get("type") in types]
            if not rows:
                continue
            if show_desc:
                items = "; ".join(
                    f"**{_res_loc(f)}**" + (f" ({f['description']})" if f.get("description") else "")
                    for f in rows[:12])
                more = " …" if len(rows) > 12 else ""
                lines.append(f"\n**{label} ({len(rows)}):** {items}{more}")
            else:
                locs = ", ".join(_res_loc(f) for f in rows[:20])
                more = " …" if len(rows) > 20 else ""
                lines.append(f"\n**{label} ({len(rows)}):** {locs}{more}")
        mut = [f for f in ufe if f.get("type") == "mutagenesis site"]
        if mut:
            lines.append(f"\n**Mutagenesis-validated functional residues ({len(mut)}):**\n")
            lines.append(table(["Position", "Phenotype"],
                               [(_res_loc(f), (f.get("description") or "")[:120]) for f in mut[:25]]))
        return lines

    blocks = []
    for u in (b.get("reviewed_uniprot") or []):
        ufe = [f for f in feats if f.get("uniprot") == u]
        body = _product_block(ufe)
        if body:
            blocks.append((u, body))
    if not blocks:
        return ""
    L = ["## Functional residue map", "",
         "Curated UniProt residues grouped by drug-discovery relevance — "
         "catalytic, ligand-binding, modification, and mutation-validated "
         "positions. *Source: UniProtKB sequence features.*"]
    multi = len(blocks) > 1
    for u, body in blocks:
        if multi:
            L.append(f"\n**{u}{' (canonical)' if u == canon else ''}**")
        L.extend(body)
    return "\n".join(L)


def r_structure(b):
    L = ["## Structure", ""]
    n = b.get("pdb_count", 0)
    pdb = b.get("pdb", []) or []
    # Cap the dump — TP53 carries 313 PDB entries, EGFR 351; an uncapped table
    # buries the rest of the Protein section. Best resolution first (None last),
    # like the other capped tables. Full set is paginated in biobtree (p= param).
    def _res(p):                            # resolution arrives as a string;
        try:                                # coerce for sorting, push N/A last
            return float(p.get("resolution"))
        except (TypeError, ValueError):
            return float("inf")
    pdb_sorted = sorted(pdb, key=_res)
    shown = pdb_sorted[:30]
    note = f", top {len(shown)} by resolution" if len(pdb) > len(shown) else ""
    L.append("\n### Experimental structures (PDB) {#pdb}\n")
    L.append(f"{n} structures{note}.\n")
    L.append(table(["PDB", "Method", "Resolution (Å)"],
                   [(p["id"], p.get("method"), fnum(p.get("resolution")))  # 2-dp; raw is 1.83549
                    for p in shown]))
    L.append("\n### Predicted structure (AlphaFold) {#alphafold}\n")
    afs = b.get("alphafold", []) or []
    present_rows = [a for a in afs if a.get("present", True) and a.get("id")]
    missing_rows = [a for a in afs if not a.get("present", True)]
    if present_rows:
        L.append(table(["Model", "pLDDT", "Fraction very-high"],
                       [(a.get("id") or "",
                         a.get("plddt"), a.get("fraction_plddt_very_high"))
                        for a in present_rows]))
    for a in missing_rows:
        L.append(f"\n*No AlphaFold model available for {a['uniprot']} — "
                 "AlphaFold DB does not currently provide models for proteins "
                 "above ~3000 aa.*")
    # Antibody-complex structures (SAbDab) — antibodies co-crystallized with this
    # protein; a high count marks a validated antibody target. Structural (not
    # necessarily therapeutic) antibodies.
    ab = b.get("antibody_structures") or []
    if ab:
        L.append(f"\n**Antibody-complex structures (SAbDab): {len(ab)}** — "
                 + ", ".join(f"`{p}`" for p in ab[:25])
                 + (f" (+{len(ab) - 25} more)" if len(ab) > 25 else ""))
    return "\n".join(L)


def r_orthologs(b):
    L = ["## Cross-species orthologs", "", f"**{b.get('ortholog_count', 0)} orthologs**\n"]
    L.append(table(["Organism", "Symbol", "Gene ID"],
                   [(o.get("organism"), o.get("symbol"), o["id"]) for o in b.get("orthologs", [])]))
    para = b.get("paralogs", [])
    if para:
        L.append(f"\n**Paralogs ({b.get('paralog_count', 0)}):** "
                 + ", ".join(f"{p.get('symbol')} ({p['id']})" for p in para[:40]))
    return "\n".join(L)


def r_variants(b):
    L = ["## Clinical variants and AI predictions", ""]
    bd = b.get("clinvar_breakdown", {})
    L.append("\n### ClinVar {#clinvar}\n")
    L.append(f"{b.get('clinvar_total', 0)} variants total. Per-class counts are floors "
             f"(≥ shown; pagination cap):\n")
    L.append(table(["Classification", "Count (floor)"], list(bd.items())))
    L.append(f"\n### Top pathogenic / likely-pathogenic ({len(b.get('top_pathogenic', []))}) {{#top-pathogenic}}\n")
    L.append(table(["Variant ID", "HGVS", "Classification"],
                   [(v["id"], v.get("hgvs"), v.get("classification")) for v in b.get("top_pathogenic", [])]))
    L.append("\n### SpliceAI {#spliceai}\n")
    L.append(f"{b.get('spliceai_total', 0)} predictions. Top by Δscore:\n")
    L.append(table(["Variant", "Effect", "Δscore"],
                   [(v["id"], v.get("effect"), v.get("score")) for v in b.get("top_spliceai", [])]))
    L.append("\n### AlphaMissense {#alphamissense}\n")
    L.append(f"{b.get('alphamissense_total', 0)} scored. Top likely-pathogenic:\n")
    L.append(table(["Variant", "Protein change", "am_pathogenicity"],
                   [(v["id"], v.get("variant"), v.get("am_pathogenicity")) for v in b.get("top_alphamissense", [])]))
    ds = b.get("dbsnp_sample", [])
    if ds:
        L.append(f"\n**dbSNP variants (sampled {b.get('dbsnp_sampled', 0)} via entrez):** "
                 + ", ".join(f"{d['id']} ({d['pos']} {d['change']})" for d in ds[:15]))
    return "\n".join(L)


def r_pathways(b):
    L = ["## Pathways and Gene Ontology", "",
         "### Reactome pathways {#reactome}", "",
         f"{b.get('reactome_count', 0)} pathways\n"]
    L.append(table(["ID", "Pathway"], [(p["id"], p.get("name")) for p in b.get("reactome", [])[:30]]))
    L.append(f"\n**MSigDB gene sets: {b.get('msigdb_total', 0)}** (showing top):")
    L.append(", ".join(f"`{m['name']}`" for m in b.get("msigdb", [])[:15]))
    go = b.get("go", {})
    for cat in ("biological_process", "molecular_function", "cellular_component"):
        terms = go.get(cat, [])
        L.append(f"\n**GO {cat.replace('_', ' ').title()} ({len(terms)}):**")
        L.append(", ".join(f"{t['name']} ({t['id']})" for t in terms[:40]))

    # Top-level parent rollups — give the page a hierarchical-navigation view.
    # Reactome's hierarchy is tight (1-2 parents per pathway); GO's is broader.
    rp = b.get("reactome_parent_rollup") or []
    if rp:
        L.append("\n### Reactome top-level categories {#reactome-categories}\n")
        L.append(f"Rollup of top-{len(rp)} pathways:\n")
        L.append(table(["Category", "Pathways"],
                       [(p.get("name") or p["id"], p.get("pathway_count")) for p in rp[:15]]))
    gp = b.get("go_parent_rollup") or []
    if gp:
        L.append("\n### GO top-level categories {#go-categories}\n")
        L.append("Rollup of top GO terms by namespace:\n")
        L.append(table(["Category", "Terms"],
                       [(p.get("name") or p["id"], p.get("term_count")) for p in gp[:40]]))
    return "\n".join(L)


def r_interactions(b):
    L = ["## Protein interactions and networks", ""]
    L.append("\n### STRING {#string}\n")
    L.append(f"{b.get('string_count', 0)} interactions, top by confidence (×1000):\n")
    # Show both sides of each interaction (this protein ↔ partner) + the partner's
    # UniProt accession. Partner is the non-query side (biobtree #34 workaround).
    self_sym = b.get("symbol") or "—"
    L.append(table(["Protein A", "Protein B", "Partner UniProt", "Score"],
                   [(self_sym, s.get("partner_symbol") or s.get("partner"),
                     s.get("partner"), s.get("score")) for s in b.get("string", [])[:40]]))
    L.append("\n### IntAct {#intact}\n")
    L.append(f"{b.get('intact_count', 0)} interactions, top by confidence:\n")
    L.append(table(["A", "B", "Type", "Score"],
                   [(i.get("a"), i.get("b"), i.get("type"), i.get("score")) for i in b.get("intact", [])[:40]]))
    L.append(_labeled(f"BioGRID ({b.get('biogrid_count', 0)})",
                      (f"{x.get('partner')} ({x.get('method')})" for x in b.get("biogrid", [])[:15])))
    L.append(_labeled("ESM2 similar proteins", (f"`{p}`" for p in b.get("esm2_similar", [])[:40])))
    L.append(_labeled("Diamond homologs", (f"`{p}`" for p in b.get("diamond_similar", [])[:40])))
    L.append("\n### SIGNOR signaling {#signor}\n")
    L.append(f"{b.get('signor_count', 0)} interactions.\n")
    L.append(table(["A", "Effect", "B", "Mechanism"],
                   [(s.get("a"), s.get("effect"), s.get("b"), s.get("mechanism")) for s in b.get("signor", [])[:30]]))
    return "\n".join(L)


def r_tf_regulation(b):
    L = ["## Regulation", "",
         f"**Is transcription factor: {'yes' if b.get('is_transcription_factor') else 'no'}**\n"]
    L.append("\n### Downstream targets (CollecTRI) {#collectri}\n")
    L.append(f"{b.get('downstream_count', 0)} targets.\n")
    L.append(table(["Target", "Regulation"],
                   [(t.get("target"), t.get("regulation")) for t in b.get("downstream_targets", [])[:30]]))
    # JASPAR motifs — only render when present (most non-TF genes have none,
    # which would otherwise leave an empty header-only table).
    motifs = b.get("jaspar_motifs") or []
    if motifs:
        L.append("\n### JASPAR motifs {#jaspar}\n")
        L.append(table(["Motif", "Name", "Family"],
                       [(m["id"], m.get("name"), m.get("family")) for m in motifs]))
    # JASPAR PMIDs — evidence trail for the motifs above.
    pmids = b.get("jaspar_pmids") or []
    if pmids:
        links = ", ".join(f"PMID:{p}" for p in pmids[:10])
        L.append(f"\n*JASPAR matrix evidence (PMIDs):* {links}")
    L.append(_labeled("Upstream regulators (CollecTRI, top)",
                      (r.get("regulator") for r in b.get("upstream_regulators", [])[:40])))
    # miRDB miRNAs targeting this gene — post-transcriptional regulators.
    # Sorted by max_score (miRDB confidence). target_count is the miRNA's
    # promiscuity across all genes (lower = more specific target relationship).
    n_mir = b.get("mirna_count", 0)
    if n_mir:
        L.append("\n### miRNA regulators (miRDB) {#mirdb}\n")
        L.append(f"{n_mir} targeting {b.get('symbol')}, top 30 by miRDB confidence "
                 f"(max_score; target_count = how many genes the miRNA targets in total "
                 f"— lower means more specific):\n")
        L.append(table(["miRNA", "Max score", "Avg score", "miRNA target_count"],
                       [(m["id"], m.get("max_score"), m.get("avg_score"),
                         m.get("target_count")) for m in b.get("mirna_regulators", [])]))
    return "\n".join(L)


def r_drugs(b):
    L = ["## Drug and pharmacology data", "",
         f"**Is drug target: {'yes' if b.get('is_drug_target') else 'no'}**\n"]
    # ChEMBL target/bioactivity blocks only render when populated — non-target
    # genes (e.g. APOE) would otherwise show "(0)" lines + empty tables. CIViC /
    # PharmGKB / CTD blocks below still carry the section for such genes.
    tgts = b.get("chembl_targets", [])
    if tgts:
        L.append(f"**ChEMBL targets ({len(tgts)}):** "
                 + ", ".join(f"{t['id']} ({t.get('type')})" for t in tgts[:10]))
    mols = b.get("molecules", [])
    mc = b.get("molecule_count", 0)
    if mols or mc:
        pt = b.get("patent_total", 0)
        patent_note = (f" Patent mentions across the top 20 by phase: **{pt:,}** "
                       f"(via chembl_molecule>>patent_compound — counts attach to the compound, "
                       f"not the gene–compound relationship, so off-target/promiscuous "
                       f"molecules can dominate). " if pt else "")
        L.append("\n### Molecules with ChEMBL bioactivity {#chembl-molecules}\n")
        L.append(f"{mc} molecules (phase ≥1), by development phase (incl. off-target/"
                 f"promiscuous compounds).{patent_note}\n")
        L.append(table(["Molecule", "Name", "Phase", "Patents"],
                       [(m["id"],
                         links.maybe_link(m.get("name"), links.drug_url(chembl_id=m["id"], name=m.get("name"))),
                         m.get("phase"),
                         f"{m['patent_count']:,}" if m.get("patent_count") else "")
                        for m in mols[:30]]))
    # CIViC clinical evidence — drug × variant × indication (the precision-
    # medicine triple). Predictive associations only, deduped + ranked by CIViC
    # evidence level (A validated → E inferential). The Effect column separates
    # Sensitivity/Response from Resistance (opposite clinical meaning, both
    # actionable). Each row links to a representative CIViC evidence item;
    # "+N" flags additional supporting items for the same association. Empty
    # for non-cancer genes so the block elides cleanly.
    ce = b.get("civic_evidence") or []
    if ce:
        etc = b.get("civic_evidence_type_counts") or {}
        extra = ", ".join(f"{n} {k.lower()}" for k, n in etc.items() if k != "Predictive")
        L.append("\n### Clinical evidence (CIViC) {#civic}\n")
        L.append(f"Drug × variant × indication: {b.get('civic_association_total', 0)} "
                 f"predictive associations from {b.get('civic_predictive_total', 0)} curated "
                 f"evidence items" + (f"; also {extra}" if extra else "") + ".\n")
        L.append(table(["Variant", "Therapy", "Indication", "Effect", "Level", "CIViC"],
                       [(r["profile"],
                         links.maybe_link(therapy_label(r["therapy"]), links.drug_url(name=therapy_label(r["therapy"]))),
                         links.maybe_link(r["disease"], links.disease_url(name=r["disease"])),
                         r["significance"],
                         f"CIViC {r['level']}" if r.get("level") else "",
                         f"EID{r['evidence_id']}"
                         + (f" +{r['n']-1}" if r.get("n", 1) > 1 else ""))
                        for r in ce]))
        more = b.get("civic_association_total", 0) - len(ce)
        if more > 0:
            L.append(f"\n*+{more} more predictive associations (showing top "
                     f"{len(ce)} by evidence level).*")

    pg = b.get("pharmgkb", [])
    L.append(f"\n**PharmGKB:** {len(pg)} entr{'y' if len(pg)==1 else 'ies'}"
             + (f" (VIP={pg[0].get('vip')}, CPIC={pg[0].get('cpic_guideline')})" if pg else ""))

    # PharmGKB clinical annotations — variant + drug + phenotype tuples
    # with PharmGKB's evidence-level rating (1A strongest → 4 weakest).
    pgc = b.get("pharmgkb_clinical") or []
    if pgc:
        L.append("\n### PharmGKB clinical annotations {#pharmgkb-clinical}\n")
        L.append(f"{len(pgc)} annotations.\n")
        L.append(table(
            ["Variant", "Type", "Level", "Drugs", "Phenotypes"],
            [(c.get("variant"), c.get("type"), c.get("level_of_evidence"),
              links.link_csv(c.get("chemicals"), lambda s: links.drug_url(name=s)),
              c.get("phenotypes"))
             for c in pgc[:40]]))

    # PharmGKB variant pages — variant-level aggregations with PharmGKB's
    # composite score + count of clinical annotations.
    pgv = b.get("pharmgkb_variant") or []
    if pgv:
        L.append("\n### PharmGKB variants {#pharmgkb-variants}\n")
        L.append(f"{len(pgv)} variants.\n")
        L.append(table(
            ["Variant", "Genes", "Level", "Score", "#Clin annots", "Drugs"],
            [(v.get("name"),
              links.link_csv(v.get("gene_symbols"), lambda s: links.gene_url(symbol=s)),
              v.get("level_of_evidence"),
              v.get("score"), v.get("clinical_annotation_count"),
              links.link_csv(v.get("associated_drugs"), lambda s: links.drug_url(name=s)))
             for v in pgv[:40]]))

    # PharmGKB guidelines — CPIC / DPWG / CPNDS dosing guidance per
    # gene+drug pair. Canonical pharmacogenes have tens of guidelines
    # (CYP2D6: 69, CYP2C19: 37); non-pharmacogenes have none.
    pgg = b.get("pharmgkb_guideline") or []
    if pgg:
        # CPIC > DPWG > others as canonical priority.
        order = {"CPIC": 0, "DPWG": 1, "CPNDS": 2}
        pgg_sorted = sorted(pgg, key=lambda g: (order.get(g.get("source"), 99),
                                                 g.get("chemical_names") or ""))
        L.append("\n### PharmGKB dosing guidelines {#pharmgkb-guidelines}\n")
        L.append(f"{len(pgg)} guidelines.\n")
        L.append(table(
            ["Source", "Drug", "Guideline", "Dosing?", "Recommendation?"],
            [(g.get("source"),
              links.link_csv(g.get("chemical_names"), lambda s: links.drug_url(name=s)),
              g.get("name"),
              "yes" if g.get("has_dosing_info") else "",
              "yes" if g.get("has_recommendation") else "")
             for g in pgg_sorted[:30]]))
    # GtoPdb / IUPHAR — hand-curated pharmacology (Tier 1: leads, plainly).
    gt = b.get("gtopdb_target")
    gi = b.get("gtopdb_interactions") or []
    if gt or gi:
        L.append("\n### GtoPdb / IUPHAR curated pharmacology {#gtopdb}\n")
        L.append("*(IUPHAR/BPS Guide to Pharmacology — expert-curated)*")
        if gt:
            cls = (gt.get("type") or "").replace("_", " ")
            fam = gt.get("family")
            L.append(f"\n**Target class:** {cls}" + (f" — *{fam}*" if fam else ""))
        if gi:
            L.append(f"\n**Most potent curated ligand interactions "
                     f"({b.get('gtopdb_interaction_count', 0)} total), top {len(gi)}:**\n")
            def _gact(x):
                # GtoPdb action AND type can both be the literal 'None'/'Unknown'
                # — take the first meaningful one, else blank (no null-looking cell).
                for v in (x.get("action"), x.get("type")):
                    v = (v or "").strip()
                    if v and v.lower() not in ("none", "unknown"):
                        return v
                return ""
            L.append(table(["Ligand", "Action", "Affinity", "Parameter"],
                           [(links.maybe_link(x.get("ligand"), links.drug_url(name=x.get("ligand"))),
                             _gact(x), x.get("affinity"), x.get("parameter")) for x in gi]))

    # BindingDB — measured affinities (Tier 2: heterogeneous assays; the note
    # carries the provenance so the ranking isn't read as one comparable scale).
    br = b.get("bindingdb_ranked") or []
    if br:
        L.append("\n### Binding affinities (BindingDB) {#bindingdb}\n")
        L.append(f"{b.get('bindingdb_measured', 0)} measured of {b.get('bindingdb_human', 0)} "
                 f"human assays ({b.get('bindingdb_total', 0)} total across all organisms); "
                 f"most potent {len(br)} below. *Values come from heterogeneous assays and "
                 f"are not directly comparable.*\n")
        # Patent column when a displayed compound carries a source patent — the
        # hyphenated number (US-8524722) linked to Google Patents, plus the
        # invention title (more readable than the IUPAC). Blank for the rest.
        if any(x.get("patent") for x in br):
            def _pat(x):
                pn = x.get("patent")
                if not pn:
                    return ""
                disp = pn[:2] + "-" + pn[2:] if len(pn) > 2 and pn[:2].isalpha() else pn
                cell = links.maybe_link(disp, f"https://patents.google.com/patent/{pn}")
                return f"{cell}: {x['patent_title']}" if x.get("patent_title") else cell
            L.append(table(["Ligand", "Measure", "Value", "Patent"],
                           [(x.get("ligand"), x.get("measure"), x.get("value"), _pat(x))
                            for x in br]))
        else:
            L.append(table(["Ligand", "Measure", "Value"],
                           [(x.get("ligand"), x.get("measure"), x.get("value")) for x in br]))

    # ChEMBL bioactivities (pchembl-ranked). pchembl is the gold potency
    # metric — directly comparable across assay types.
    ca = b.get("chembl_activities") or []
    if ca:
        L.append("\n### ChEMBL bioactivities {#chembl-bioactivity}\n")
        L.append(f"{b.get('chembl_activity_potent_count', 0)} potent at pChembl≥5 of "
                 f"{b.get('chembl_activity_total', 0)} total, top {len(ca)} by pChembl "
                 f"(potency: 10 = 0.1 nM, 6 = 1 µM).\n")
        L.append(table(["pChembl", "Type", "Value", "Unit", "Molecule"],
                       [(r.get("pchembl"), r.get("type"), r.get("value"), r.get("unit"),
                         links.maybe_link(r.get("molecule_name") or r.get("molecule_id") or "",
                                          f"https://www.ebi.ac.uk/chembl/compound_report_card/{r['molecule_id']}/")
                         if r.get("molecule_id") else "")
                        for r in ca]))

    # PubChem BioAssay actives — sorted by potency. CID/AID get clickable
    # PubChem URLs so an AI agent (or human) can drill into the assay record
    # directly without needing a separate entry fetch.
    pba = b.get("pubchem_bioassay") or []
    if pba:
        L.append("\n### PubChem BioAssay actives {#pubchem-bioassay}\n")
        L.append(f"{b.get('pubchem_bioassay_active_count', 0)} with measured affinity, of "
                 f"{b.get('pubchem_bioassay_total', 0)} total; {len(pba)} most potent distinct "
                 f"compounds. *Largely complementary to BindingDB; screening values are coarse "
                 f"(µM, 4 dp), so sub-nM hits tie at the floor.*\n")
        def _cmpd(p):
            return p.get("name") or (
                links.maybe_link(p["cid"], f"https://pubchem.ncbi.nlm.nih.gov/compound/{p['cid']}/")
                if p.get("cid") else "")

        def _assay(p):
            if not p.get("aid"):
                return ""
            link = links.maybe_link(p["aid"], f"https://pubchem.ncbi.nlm.nih.gov/bioassay/{p['aid']}/")
            nm = (p.get("assay_name") or "").replace("*", r"\*")   # what the assay measured
            return f"{link}: {nm}" if nm else link
        L.append(table(["Compound", "Assay", "Type", "Value", "Unit"],
                       [(_cmpd(p), _assay(p), p.get("activity_type"), p.get("value"), p.get("unit"))
                        for p in pba]))

    # CTD literature-mined chemical-gene interactions — Comparative
    # Toxicogenomics Database. Each row: a chemical (drug, toxin,
    # environmental compound) + CV-coded action verbs + PubMed-count support.
    # High AI value: every claim is anchored by literature counts.
    ctd = b.get("ctd_interactions") or []
    if ctd:
        L.append("\n### CTD chemical–gene interactions {#ctd}\n")
        L.append(f"{b.get('ctd_interaction_total', 0)} total (human), top {len(ctd)} by "
                 f"PubMed support.\n")
        L.append(table(["Chemical", "Actions (top 5)", "PubMed papers"],
                       [(r["chemical"], _ctd_actions(r["actions"]),
                         r["pmids"]) for r in ctd]))

    # ChEMBL screening-assay depth — how heavily this target has been
    # profiled, and via what experimental style (Binding / Functional / ADMET
    # / Toxicity). Aggregated across all ChEMBL targets for this protein.
    cat = b.get("chembl_assay_total") or 0
    if cat:
        type_counts = b.get("chembl_assay_type_counts") or {}
        breakdown = ", ".join(f"{n} {k.lower()}" for k, n in type_counts.items())
        L.append("\n### ChEMBL screening assays {#chembl-assays}\n")
        L.append(f"{cat} unique, capped per target: {breakdown}")
        samples = b.get("chembl_assay_samples") or []
        if samples:
            L.append("\nRepresentative assays (with source publication via chembl_document):\n")
            rows = []
            for s in samples:
                doc = ""
                if s.get("doc_id"):
                    title = (s.get("doc_title") or s["doc_id"]).strip()
                    doc = (title
                           + (f" — *{s['doc_journal']}*" if s.get("doc_journal") else ""))
                rows.append((s["id"], s["type"], s["desc"], doc))
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
        L.append("\n### Cellosaurus cell lines {#cellosaurus}\n")
        L.append(f"{cct:,} cell lines: {breakdown}")
        cs = b.get("cellosaurus_samples") or []
        if cs:
            L.append("\nFirst 10 cell lines (id-ordered, not curated):\n")
            L.append(table(["Cellosaurus", "Name", "Category", "Sex"],
                           [(c["id"], c.get("name"), c.get("category"), c.get("sex"))
                            for c in cs]))

    ct = b.get("disease_trials", [])
    L.append("\n### Clinical trials (associated diseases) {#gene-trials}\n")
    L.append(f"{b.get('disease_trial_count', 0)} trials via MONDO — disease-level, not "
             f"drug-specific.\n")
    L.append(table(["Trial", "Phase", "Status", "Title"],
                   # phase_label (audit #8): biobtree emits 'NaN' for trials with
                   # no interventional phase; render it as 'Not specified', not
                   # a leaked 'nan'/'NAN'.
                   [(t["id"], phase_label(t.get("phase")), t.get("status"),
                     (t.get("title") or "").strip()) for t in ct[:40]]))
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
        L.append(f"\n### FANTOM5 promoters ({len(fp)} alternative TSS) {{#fantom5-promoters}}\n")
        L.append(table(["Promoter ID", "TPM avg", "Samples expressed"],
                       [(p["id"], p.get("tpm_average"), p.get("samples_expressed")) for p in fp[:10]]))
    # Tissue name as plain text — the UBERON/CL id is shown in its own column;
    # we don't link every id (consistency: links are reserved for primary
    # targets, not decorative on every row).
    # Bgee expression score (0-100, higher = more highly expressed; derived from
    # the per-condition rank). We drop the raw Rank column — it's inversely
    # redundant with the score and the pair reads as a contradiction (high score
    # next to a mid-looking rank) when it isn't.
    L.append("\n### Top tissues by expression {#tissue-expression}\n")
    L.append(f"{b.get('tissue_count', 0)} total, by Bgee expression score (0-100, higher "
             f"= more expressed):\n")
    L.append(table(["Tissue", "Anatomy ID", "Expression score", "Quality"],
                   [(t.get("tissue") or "", t.get("anatomy_id") or "",
                     t.get("score"), t.get("quality"))
                    for t in b.get("top_tissues", [])[:30]]))
    # Single-cell (SCXA) — per-gene marker status + max expression across
    # single-cell experiments (biobtree #31: via the scxa_expression node).
    sc = b.get("single_cell") or {}
    if sc.get("total_experiments"):
        L.append("\n### Single-cell (SCXA) {#scxa}\n")
        L.append(f"Detected in {sc['total_experiments']} experiment(s), a significant "
                 f"marker in {sc.get('marker_experiments', 0)}.\n")
        exps = sc.get("experiments") or []
        if exps:
            L.append(table(["Experiment", "Marker?", "Max mean expression"],
                           [(_scxa_link(e.get("experiment_id")),
                             "yes" if e.get("is_marker") else "",
                             e.get("max_expression")) for e in exps[:15]]))
    return "\n".join(L)


def _scxa_link(eid):
    """Link an SCXA experiment accession to its EBI Single Cell Expression Atlas
    page; plain text when missing."""
    return (f"[{eid}](https://www.ebi.ac.uk/gxa/sc/experiments/{eid})"
            if eid else "")


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
             f"disease phenotypes: {', '.join(b.get('disease_omim', [])[:40])}")
    L.append("\n### GenCC curated gene-disease {#gencc}\n")
    # Collapse to one row per disease (audit #13: GenCC has many submissions per
    # gene–disease pair — FANCB ×4); keep the strongest classification.
    from atlas.render_common import gencc_rank
    _gc = {}
    for g in b.get("gencc", []):
        d = g.get("disease") or ""
        if d and (d not in _gc or gencc_rank(g.get("classification")) > gencc_rank(_gc[d].get("classification"))):
            _gc[d] = g
    _gc_rows = sorted(_gc.values(),
                      key=lambda g: -gencc_rank(g.get("classification")))
    L.append(table(["Disease", "Classification", "Inheritance"],
                   [(links.maybe_link(g.get("disease"), links.disease_url(name=g.get("disease"))),
                     g.get("classification"), g.get("inheritance")) for g in _gc_rows[:40]]))

    # ClinGen Gene-Disease Validity — expert-panel curated relationship strength.
    # Distinct from GenCC: ClinGen is the canonical authority for gene-disease
    # validity classifications used by clinical variant labs (ACMG/AMP guidelines).
    cgv = b.get("clingen_validity") or []
    if cgv:
        L.append(f"\n### ClinGen Gene-Disease Validity ({len(cgv)}) {{#clingen}}\n")
        L.append("Expert-panel classifications — Definitive > Strong > Moderate > "
                 "Limited > Disputed > Refuted.\n")
        L.append(table(["Disease", "Classification", "Inheritance"],
                       [(links.maybe_link(c.get("disease"), links.disease_url(name=c.get("disease"))),
                         c.get("classification"), c.get("moi"))
                        for c in cgv[:40]]))
    L.append(f"\n**Mondo ({len(b.get('mondo', []))}):** "
             + ", ".join(f"{links.maybe_link(m.get('name'), links.disease_url(mondo_id=m['id'], name=m.get('name')))} ({m['id']})"
                         for m in b.get("mondo", [])[:15]))
    L.append(f"\n**Orphanet ({len(b.get('orphanet', []))}):** "
             + ", ".join(f"{o.get('name') or ''} ({o['id']})" for o in b.get("orphanet", [])[:15]))
    # No relevance/frequency ordering on the gene→HPO route, so don't claim
    # "(top)" — these are the first 30 by HPO id. (Frequency-ranked clinical
    # features live on the disease pages.)
    _hpo = b.get("hpo", [])
    L.append("\n### HPO phenotypes {#hpo}\n")
    L.append(f"{b.get('hpo_total', 0)} total ({min(30, len(_hpo))} of "
             f"{b.get('hpo_total', 0)} shown, HPO-id order):\n")
    L.append(table(["HPO", "Term"], [(h["id"], h.get("name")) for h in _hpo[:30]]))
    L.append("\n### GWAS associations {#gwas-assoc}\n")
    L.append(f"{b.get('gwas_total', 0)} associations (top):\n")
    L.append(table(["Study", "Trait", "p-value"],
                   [(g["id"], g.get("trait"), g.get("p_value")) for g in b.get("gwas", [])[:30]]))

    # EFO canonical trait names — normalized vocabulary for the free-text
    # GWAS-catalog traits above. The same disease can appear under multiple
    # synonymous GWAS-catalog spellings; EFO collapses them. Useful for
    # AI agents doing cross-resource joins on trait IDs.
    efo = b.get("efo_traits") or []
    if efo:
        L.append(f"\n### EFO canonical traits ({len(efo)}, from GWAS) {{#efo}}\n")
        L.append(table(["EFO ID", "Trait name"],
                       [(e["id"], e.get("name")) for e in efo[:30]]))

    # MeSH disease descriptors — NLM's controlled disease vocabulary.
    # Tree numbers (e.g. C04.700.600) classify into MeSH categories
    # (C04=Neoplasms, C16=Congenital, C18=Nutritional/Metabolic, etc.) —
    # useful for grouping diseases at a coarser level.
    # Drop nameless descriptors — biobtree occasionally returns a MeSH id with
    # no name (e.g. `C538339`); a row with a bare id misleads LLMs.
    mesh = [m for m in (b.get("mesh_descriptors") or []) if (m.get("name") or "").strip()]
    if mesh:
        L.append(f"\n### MeSH disease descriptors ({len(mesh)}) {{#mesh}}\n")
        L.append(table(["Descriptor", "Name", "Tree numbers"],
                       [(m["id"],
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
