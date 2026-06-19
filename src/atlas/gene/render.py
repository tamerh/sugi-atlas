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
from atlas.render_common import table, phase_label, fnum, more_line, capped_table, pval
from atlas.civic import therapy_label, LEGEND as CIVIC_LEGEND
from atlas.page import links

# Display cap for "curated, naturally bounded" tables (cohort/pathway/GenCC/etc.):
# most pages show everything; the cap only bites rare outliers, and the "+N more"
# line discloses when it does. NOT used for score-ranked firehose tables (STRING,
# IntAct, GWAS, SpliceAI, AlphaMissense, CollecTRI targets, miRNA, tissues) — those
# keep their own literal caps pending a confidence/significance-threshold redesign.
ROW_CAP = 60


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


def r_genomic_neighbors(b):
    nb = b.get("genomic_neighbors") or []
    if not nb:
        return ""
    total = b.get("genomic_neighbor_count") or len(nb)
    L = ["## Genomic neighbors (locus context)", "",
         "Genes flanking this one on the chromosome (NCBI/Entrez gene "
         "neighborhood). Positional context only — some neighbors share "
         "regulatory elements (bidirectional promoters, antisense/overlapping "
         "transcripts), but chromosomal adjacency is not itself a functional "
         "relationship (see Interactions for that).", "",
         table(["Symbol", "Type", "Entrez"],
               [(links.maybe_link(n.get("symbol"), links.gene_url(symbol=n.get("symbol"))),
                 n.get("type"), n.get("id")) for n in nb])]
    L.append(more_line(total, len(nb)))
    return "\n".join(L)


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
        L.append("\n### Enzyme classification (BRENDA) {#enzyme-ec}\n")
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
        L.append("\n### Catalyzed reactions (Rhea) {#rhea}\n")
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
    gc = b.get("gnomad_constraint") or {}
    has_gc = gc.get("pli") not in (None, "") or gc.get("loeuf") not in (None, "")
    if not cd and not dm_notable and not has_gc:
        return ""
    L = ["## Functional genomics", ""]
    if has_gc:
        bits = []
        for label, key in (("pLI", "pli"), ("LOEUF", "loeuf"), ("missense Z", "mis_z")):
            v = gc.get(key)
            if v not in (None, ""):
                try:
                    bits.append(f"{label} {float(v):.3g}")
                except (TypeError, ValueError):
                    pass
        if bits:
            L.append(f"**gnomAD constraint:** {', '.join(bits)} — population "
                     "loss-of-function intolerance (pLI ≥ 0.9 = LoF-intolerant; lower "
                     "LOEUF = more constrained). "
                     "[gnomAD](https://gnomad.broadinstitute.org/)")
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
        # Which models depend on it most (gene-effect, lower = stronger dependency).
        lines = b.get("depmap_lines") or []
        if lines:
            names = ", ".join(f"{d['cell_line']} ({d['gene_effect']})" for d in lines[:15])
            L.append(f"\n**Most-dependent cell lines** *(DepMap gene-effect, lower = "
                     f"stronger; ≤ −0.5 shown)*: {names}.")
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
# feature types, display mode, anchor). Display modes:
#   "perres" — one entry per residue with its inline description (catalytic /
#              ligand-binding: each position carries a distinct description).
#   "bytype" — group positions by modification identity (PTM / disulfide /
#              glycosylation: many residues share a few modification types, so
#              "Phosphoserine — 229, 695" beats a bare, meaningless number list).
_RESIDUE_GROUPS = [
    ("Catalytic / active sites", ("active site", "site"), "perres", "residue-catalytic"),
    ("Ligand- & substrate-binding residues", ("binding site",), "perres", "residue-binding"),
    ("Post-translational modifications", ("modified residue",
     "lipid moiety-binding region", "cross-link"), "bytype", "residue-ptm"),
    ("Disulfide bonds", ("disulfide bond",), "bytype", "residue-disulfide"),
    ("Glycosylation sites", ("glycosylation site",), "bytype", "residue-glycosylation"),
]


def _res_loc(f):
    b, e = f.get("begin"), f.get("end")
    return str(b) if (e in (None, "", b)) else f"{b}–{e}"


def _clip(text, n=140):
    """Truncate on a word boundary with an ellipsis — never mid-word (a naive
    [:n] left mutagenesis phenotypes ending '…downstream kinases. reduced tran')."""
    text = (text or "").strip()
    if len(text) <= n:
        return text
    cut = text[:n].rsplit(" ", 1)[0].rstrip(",;.") or text[:n]
    return cut + " …"


def _mod_label(desc):
    """The modification identity from a UniProt feature description — the part
    before the first ';' (drops the '; by <kinase>' / evidence tail), with the
    leading letter capitalized to UniProt's display convention. '' when the
    feature carries no description (most intrachain disulfide bonds)."""
    head = (desc or "").split(";")[0].strip()
    return head[:1].upper() + head[1:] if head else ""


def _join_locs(locs, cap=40):
    """Comma-join residue positions for a table cell, capped so one modification
    type can't blow out the cell; the H4 heading still carries the true count."""
    return ", ".join(locs[:cap]) + (f" … (+{len(locs) - cap})" if len(locs) > cap else "")


def _residue_type_table(rows, fallback):
    """Positions grouped by modification identity → a `| Type | Positions |`
    table. e.g. 'Phosphoserine | 229, 695'. Features with no description
    (unlabeled intrachain disulfide bonds) collapse onto one `fallback`-labelled
    row, so the table never has a blank Type cell."""
    groups = {}                              # label -> [locs], insertion order
    for f in rows:
        groups.setdefault(_mod_label(f.get("description")) or fallback, []).append(_res_loc(f))
    return table(["Type", "Positions"],
                 [(lab, _join_locs(locs)) for lab, locs in groups.items()])


def r_residue_map(b):
    """Structure-function residue view over UniProt features — catalytic,
    ligand-binding, PTM, disulfide, glycosylation, and mutagenesis-validated
    positions, grouped per reviewed product. Pure restructure of data §3 already
    collects (each feature is accession-stamped); '' when no features."""
    feats = b.get("ufeatures") or []
    canon = b.get("canonical_uniprot")

    def _product_block(ufe, u="", multi=False):
        """Residue-group sub-sub-sections (H4) for one product; [] when it has no
        mappable residues. Each category is its own H4 heading (the web renderer
        wraps each into a collapsible sub-block) with the body underneath; the
        count stays in the heading. For dual-product genes the visible suffix and
        the {#id} both carry the product (' — P04637' / '-p04637') so the repeated
        category headings stay unique. The residue-* anchor prefix is exempt from
        the frozen H4 contract (the set varies per gene/product)."""
        suffix = f" — {u}{' (canonical)' if u == canon else ''}" if multi else ""
        asuf = f"-{u.lower()}" if multi else ""
        lines = []
        for label, types, mode, anchor in _RESIDUE_GROUPS:
            rows = [f for f in ufe if f.get("type") in types]
            if not rows:
                continue
            if mode == "perres":             # one row per residue, with its role
                # binding sites carry the ligand in a clean field (ATP, Mg(2+)…);
                # active sites carry their role in description (proton acceptor).
                body = table(["Position", "Role"],
                             [(_res_loc(f), f.get("ligand") or f.get("description") or "")
                              for f in rows[:ROW_CAP]])
                if len(rows) > ROW_CAP:      # heading carries the true count; flag the hidden tail
                    body += f"\n\n*…and {len(rows) - ROW_CAP} more — see UniProt.*"
            else:                            # bytype — positions grouped by modification identity
                singular = label[:-1] if label.endswith("s") else label
                body = _residue_type_table(rows, singular)
            lines.append(f"\n### {label} ({len(rows)}){suffix} {{#{anchor}{asuf}}}\n\n{body}")
        mut = [f for f in ufe if f.get("type") == "mutagenesis site"]
        if mut:
            lines.append(f"\n### Mutagenesis-validated functional residues ({len(mut)}){suffix} "
                         f"{{#residue-mutagenesis{asuf}}}\n")
            lines.append(table(["Position", "Phenotype"],
                               [(_res_loc(f), _clip(f.get("description"))) for f in mut[:25]]))
        return lines

    prods = [(u, [f for f in feats if f.get("uniprot") == u])
             for u in (b.get("reviewed_uniprot") or [])]
    prods = [(u, ufe) for u, ufe in prods if _product_block(ufe)]   # drop residue-less products
    if not prods:
        return ""
    L = ["## Functional residue map", "",
         "Curated UniProt residues grouped by drug-discovery relevance — "
         "catalytic, ligand-binding, modification, and mutation-validated "
         "positions. *Source: UniProtKB sequence features.*"]
    multi = len(prods) > 1
    for u, ufe in prods:
        L.extend(_product_block(ufe, u, multi))
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
    L.append("\n### Experimental structures (PDB) {#pdb}\n")
    L.append(capped_table(["PDB", "Title", "Method", "Resolution (Å)"],
                          [(p["id"], _clip(p.get("title"), 90), p.get("method"),
                            fnum(p.get("resolution")))  # 2-dp; raw is 1.83549
                           for p in pdb_sorted],
                          30, total=n, noun="structures by resolution"))
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
        L.append("\n### Antibody-complex structures (SAbDab) {#sabdab}\n")
        L.append(f"{len(ab)} antibodies co-crystallized with this protein: "
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
                 + ", ".join(f"{p.get('symbol')} ({p['id']})" for p in para[:ROW_CAP]))
        L.append(more_line(b.get("paralog_count"), len(para[:ROW_CAP])))
    # Cross-species homologs (UniProt-wide Diamond local-alignment similarity) —
    # reach species beyond Ensembl Compara's model set (the table above). Diamond
    # only (not ESM2 embedding similarity, which over-calls — see s05_orthologs);
    # one row per species, best hit. The non-model gene symbol/Ensembl id aren't in
    # biobtree, so we show the UniProt accession + % sequence identity.
    # Observed mouse-model phenotypes (MGI knockout/mutant), via the gene's mouse
    # ortholog — real model-organism biology, NOT a translation of the human
    # gene's HPO terms. Ranked by MGI annotation count (robustness).
    mp = b.get("mouse_phenotypes") or []
    if mp:
        total = b.get("mouse_phenotype_total") or len(mp)
        mgi = (b.get("mouse_mgi_ids") or [None])[0]
        mgi_link = (f" ([{mgi}](https://www.informatics.jax.org/marker/{mgi}))"
                    if mgi else "")
        L.append("\n### Mouse-model phenotypes {#mouse-phenotypes}\n")
        L.append(f"Observed phenotypes of the mouse ortholog{mgi_link} from MGI "
                 "(Alliance of Genome Resources) — knockout/mutant phenotypes, "
                 f"ranked by annotation count. {total} distinct phenotype terms.\n")
        L.append(table(["Mouse phenotype (MP)", "Records"],
                       [(p.get("statement"), p.get("records")) for p in mp]))
        L.append(more_line(total, len(mp)))
    hom = b.get("cross_species_homologs") or []
    if hom:
        L.append("\n### Additional cross-species homologs {#cross-species-homologs}\n")
        L.append("Sequence homologs in species beyond Ensembl Compara's model "
                 "organisms, by Diamond local protein alignment (% sequence identity; "
                 "best hit per species). Lower-identity rows are more distant homologs. "
                 "The matched accessions' gene symbols aren't in biobtree yet, so each "
                 "row shows the UniProt accession.\n")
        L.append(table(["Organism", "UniProt", "Identity", "Method"],
                       [(h.get("organism"),
                         links.maybe_link(h.get("accession"), links.uniprot_url(h.get("accession"))),
                         f"{h.get('similarity'):.1%}" if isinstance(h.get("similarity"), float) else h.get("similarity"),
                         h.get("source")) for h in hom]))
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
    # ClinGen VCEP expert-panel interpretations — ACMG calls reviewed by a Variant
    # Curation Expert Panel; a higher-authority tier than individual ClinVar
    # submissions. Summary breakdown (not a per-variant dump).
    cg_total = b.get("clingen_variant_total")
    if cg_total:
        L.append(f"\n### ClinGen expert-panel interpretations ({cg_total}) {{#clingen-variants}}\n")
        vceps = b.get("clingen_variant_vceps") or []
        dis = b.get("clingen_variant_diseases") or []
        ctx = "Expert ACMG interpretations"
        if vceps:
            ctx += " from " + ", ".join(vceps)
        if dis:
            ctx += " for " + "; ".join(dis[:3]) + (" …" if len(dis) > 3 else "")
        L.append(ctx + " — a higher-authority tier than individual ClinVar "
                 "submissions.\n")
        L.append(table(["Assertion", "Variants"], b.get("clingen_variant_breakdown") or []))
    L.append("\n### SpliceAI {#spliceai}\n")
    L.append(capped_table(["Variant", "Effect", "Δscore"],
                          [(v["id"], v.get("effect"), v.get("score")) for v in b.get("top_spliceai", [])],
                          None, total=b.get("spliceai_total"), noun="predictions by Δscore"))
    L.append("\n### AlphaMissense {#alphamissense}\n")
    L.append(capped_table(["Variant", "Protein change", "am_pathogenicity"],
                          [(v["id"], v.get("variant"), v.get("am_pathogenicity")) for v in b.get("top_alphamissense", [])],
                          None, total=b.get("alphamissense_total"), noun="scored, likely-pathogenic"))
    ds = b.get("dbsnp_sample", [])
    if ds:
        shown = min(50, len(ds))
        sampled = b.get("dbsnp_sampled", 0)
        of = f" of ~{sampled:,} sampled via entrez" if sampled > shown else " via entrez"
        L.append(f"\n### dbSNP variants (showing {shown}{of}) {{#dbsnp}}\n")
        L.append("Population variants with gnomAD minor-allele frequency where "
                 "available (frequency-bearing variants shown first); blank MAF = "
                 "not reported in gnomAD (typically very rare).\n")

        def _maf(d):
            f = d.get("gnomad")
            if f is None:
                return ""
            return f"{f:.2%}" if f >= 0.0001 else f"{f:.1e}"
        def _rs(i):                       # biobtree returns "RS123"; dbSNP convention is "rs123"
            return ("rs" + i[2:]) if i[:2].upper() == "RS" else i
        L.append(table(["Variant", "Position", "Change", "gnomAD MAF", "Class"],
                       [(links.maybe_link(_rs(d["id"]), links.variant_link(_rs(d["id"]))),
                         d["pos"], d["change"], _maf(d),
                         d.get("variant_class") or "") for d in ds[:50]]))
    return "\n".join(L)


def r_noncoding_genesets(b7):
    """Function-zone content for NON-CODING genes: MSigDB gene-set membership.
    MSigDB sets are keyed by gene id, so membership is intrinsic to this gene
    (unlike the positional variant/disease data we scrub) — it's the one
    functional signal valid for a lncRNA with no protein pathway/GO/interaction
    data. '' when the gene is in no sets."""
    sets = b7.get("msigdb") or []
    names = [s.get("name") for s in sets if s.get("name")]
    if not names:
        return ""
    total = b7.get("msigdb_total") or len(names)
    return "\n".join(
        ["## Gene-set membership (MSigDB)", "",
         f"This non-coding RNA is a member of {total:,} curated MSigDB gene set"
         f"{'s' if total != 1 else ''} — regulatory-target and signature collections "
         "that include it (it has no protein-based pathway / GO / interaction data).",
         "", ", ".join(f"`{n}`" for n in names[:ROW_CAP])])


def r_ncrna_function(b14):
    """Non-coding RNA function — Rfam family + Rfam-derived GO (RNAcentral).
    Structured RNAs (rRNA/tRNA/snoRNA/miRNA) carry an Rfam family + GO; bare
    lncRNAs usually carry neither (use the disease/interaction layers instead).
    '' when the gene has no Rfam and no GO."""
    rfam = b14.get("rfam") or []
    go = b14.get("go") or []
    if not rfam and not go:
        return ""
    L = ["## Non-coding RNA function {#ncrna-function}", "",
         "Functional annotation for the RNA itself (RNAcentral / Rfam), "
         "independent of any protein product."]
    if rfam:
        L.append("\n**Rfam family:** " + ", ".join(
            f"{r.get('rfam_description') or r.get('rfam_id')} (`{r.get('rfam_id')}`)"
            for r in rfam))
    if go:
        L.append("\n### Gene Ontology (Rfam-derived) {#ncrna-go}\n")
        L.append(table(["GO ID", "Aspect", "Term"],
                       [(g.get("id"), (g.get("type") or "").replace("_", " "), g.get("name"))
                        for g in go[:ROW_CAP]]))
    return "\n".join(L)


def r_ncrna_disease(b14):
    """ncRNA -> disease associations (LncRNADisease v3.0 + HMDD). Curated,
    symbol-keyed (not positional), so valid for a non-coding gene. '' when none."""
    dis = b14.get("diseases") or []
    if not dis:
        return ""
    total = b14.get("disease_total") or len(dis)
    L = ["## ncRNA disease associations {#ncrna-disease}", "",
         "Curated non-coding-RNA → disease associations (LncRNADisease, HMDD) — "
         "these are RNA-level associations for this gene, not the positional "
         "variant data omitted above.", ""]
    L.append(capped_table(["Disease", "Causal?", "Validated by"],
                          [(d.get("disease_name"), d.get("causality"),
                            (d.get("validated_method") or "").replace("//", ", "))
                           for d in dis],
                          ROW_CAP, total=total, noun="ncRNA–disease associations"))
    return "\n".join(L)


def r_ncrna_interactions(b14):
    """ncRNA -> partner interactions (NPInter v5): RNA-RNA / RNA-protein. '' when none."""
    inter = b14.get("interactions") or []
    if not inter:
        return ""
    total = b14.get("interaction_total") or len(inter)
    L = ["## ncRNA interactions {#ncrna-interactions}", "",
         "Experimentally supported non-coding-RNA interaction partners (NPInter).", ""]
    L.append(capped_table(["Partner", "Partner type", "Level", "Source"],
                          [(links.maybe_link(i.get("partner_name"),
                                             links.gene_url(symbol=i.get("partner_name"))),
                            i.get("partner_type"), i.get("level"), i.get("datasource"))
                           for i in inter],
                          ROW_CAP, total=total, noun="interaction partners"))
    return "\n".join(L)


def r_ncrna_drugs(b14):
    """ncRNA -> drug associations (resistance / target). '' when none."""
    drg = b14.get("drugs") or []
    if not drg:
        return ""
    total = b14.get("drug_total") or len(drg)
    L = ["## ncRNA–drug associations {#ncrna-drugs}", "",
         "Curated non-coding-RNA → drug associations (drug response / resistance / "
         "target).", ""]
    L.append(capped_table(["Drug", "Relation", "Effect", "Condition"],
                          [(links.maybe_link(d.get("drug_name"),
                                             links.drug_url(name=d.get("drug_name"))),
                            d.get("relation"), d.get("effect"), d.get("condition"))
                           for d in drg],
                          ROW_CAP, total=total, noun="ncRNA–drug associations"))
    return "\n".join(L)


def r_pathways(b):
    L = ["## Pathways and Gene Ontology", "",
         "### Reactome pathways {#reactome}", ""]
    L.append(capped_table(["ID", "Pathway"], [(p["id"], p.get("name")) for p in b.get("reactome", [])],
                          ROW_CAP, total=b.get("reactome_count"), noun="pathways"))
    L.append(f"\n**MSigDB gene sets: {b.get('msigdb_total', 0)}** (showing top):")
    L.append(", ".join(f"`{m['name']}`" for m in b.get("msigdb", [])[:15]))
    go = b.get("go", {})
    for cat in ("biological_process", "molecular_function", "cellular_component"):
        terms = go.get(cat, [])
        L.append(f"\n**GO {cat.replace('_', ' ').title()} ({len(terms)}):**")
        L.append(", ".join(f"{t['name']} ({t['id']})" for t in terms[:ROW_CAP]))
        L.append(more_line(len(terms), ROW_CAP))

    # Top-level parent rollups — give the page a hierarchical-navigation view.
    # Reactome's hierarchy is tight (1-2 parents per pathway); GO's is broader.
    rp = b.get("reactome_parent_rollup") or []
    if rp:
        L.append("\n### Reactome top-level categories {#reactome-categories}\n")
        L.append(capped_table(["Category", "Pathways"],
                              [(p.get("name") or p["id"], p.get("pathway_count")) for p in rp],
                              ROW_CAP, noun="top-level categories"))
    gp = b.get("go_parent_rollup") or []
    if gp:
        L.append("\n### GO top-level categories {#go-categories}\n")
        L.append(capped_table(["Category", "Terms"],
                              [(p.get("name") or p["id"], p.get("term_count")) for p in gp],
                              ROW_CAP, noun="GO namespaces"))
    return "\n".join(L)


def r_interactions(b):
    L = ["## Protein interactions and networks", ""]
    L.append("\n### STRING {#string}\n")
    # Show both sides of each interaction (this protein ↔ partner) + the partner's
    # UniProt accession. Partner is the non-query side (biobtree #34 workaround).
    self_sym = b.get("symbol") or "—"
    L.append(capped_table(["Protein A", "Protein B", "Partner UniProt", "Score"],
                          [(self_sym, s.get("partner_symbol") or s.get("partner"),
                            s.get("partner"), s.get("score")) for s in b.get("string", [])],
                          40, total=b.get("string_count"), noun="interactions by confidence (×1000)"))
    L.append("\n### IntAct {#intact}\n")
    L.append(capped_table(["A", "B", "Type", "Score"],
                          [(i.get("a"), i.get("b"), i.get("type"), i.get("score")) for i in b.get("intact", [])],
                          40, total=b.get("intact_count"), noun="interactions by confidence"))
    L.append(_labeled(f"BioGRID ({b.get('biogrid_count', 0)})",
                      (f"{x.get('partner')} ({x.get('method')})" for x in b.get("biogrid", [])[:15])))
    L.append("\n### SIGNOR signaling {#signor}\n")
    L.append(capped_table(["A", "Effect", "B", "Mechanism"],
                          [(s.get("a"), s.get("effect"), s.get("b"), s.get("mechanism")) for s in b.get("signor", [])],
                          ROW_CAP, total=b.get("signor_count"), noun="signaling interactions"))
    # CORUM — named, curated protein complexes this protein is a subunit of. A
    # discrete biological complement to the score-ranked PPI firehose above.
    cor = b.get("corum") or []
    if cor:
        L.append(f"\n### Protein complexes (CORUM) ({b.get('corum_count', len(cor))}) {{#corum}}\n")
        L.append("Experimentally-characterized complexes this protein is a subunit "
                 "of (★ = the complex contains a known drug target).\n")
        L.append(capped_table(["Complex", "Subunits", "Members"],
                              [("★ " + c["name"] if c.get("has_drug_targets") else c["name"],
                                c.get("subunit_count"),
                                (c.get("subunits") or "").replace(";", ", "))
                               for c in cor],
                              ROW_CAP, total=b.get("corum_count"), noun="complexes by size"))
    # CellPhoneDB — curated ligand–receptor pairs (cell-cell communication), a
    # class the PPI firehose above doesn't frame. Partner gene linked; role is
    # this gene's side of the pair (ligand / receptor).
    lr = b.get("cellphonedb") or []
    if lr:
        L.append(f"\n### Ligand–receptor pairs (CellPhoneDB) ({b.get('cellphonedb_count', len(lr))}) {{#cellphonedb}}\n")
        L.append("Curated ligand–receptor partners (cell-cell communication), "
                 "distinct from the protein–protein interactions above.\n")
        L.append(capped_table(["Partner", "This gene's role", "Signaling class"],
                              [(links.maybe_link(p.get("partner"), links.gene_url(symbol=p.get("partner"))),
                                p.get("role") or "", p.get("classification") or "")
                               for p in lr],
                              ROW_CAP, total=b.get("cellphonedb_count"), noun="ligand–receptor pairs"))
    L.extend(_interactome_enrichment(b))
    return "\n".join(L)


def _interactome_enrichment(b):
    """Reactome / GO-BP over-representation among this gene's IntAct interaction
    partners — a functional readout of the neighbourhood, distinct from the gene's
    own pathway/GO memberships. Standard gene-set-size band + FDR + fold via
    atlas.ora; counts/members kept. [] when too few partners (degrade gracefully)."""
    from atlas import ora
    partners = b.get("interaction_partners") or []
    rc = ora.interactome_enrichment(partners, "reactome")
    go = ora.interactome_enrichment(partners, "go")
    if not rc and not go:
        return []
    L = ["\n### Enriched among interaction partners {#interactome-enrichment}\n",
         f"Reactome pathways and GO biological processes over-represented among this "
         f"gene's {len(partners)} IntAct physical interaction partners (hypergeometric "
         f"vs the genome-wide background, BH-FDR, gene-set size 15–500, ranked by fold). "
         f"A functional readout of the neighbourhood — *distinct from this gene's own "
         f"memberships above, and biased toward well-studied / hub proteins, so read it "
         f"as themes rather than proof.*"]

    def _tbl(rows, col0):
        return table([col0, "Partners", "Fold", "FDR"],
                     [(it.get("name") or it["id"], it["k"],
                       f"{it['fold']:.1f}×", f"{it['fdr']:.0e}") for it in rows])
    if rc:
        L += ["\n**Reactome pathways:**\n", _tbl(rc, "Pathway")]
    if go:
        L += ["\n**GO biological processes:**\n", _tbl(go, "GO term")]
    return L


def r_tf_regulation(b):
    L = ["## Regulation", "",
         f"**Is transcription factor: {'yes' if b.get('is_transcription_factor') else 'no'}**\n"]
    dt = b.get("downstream_targets") or []
    if dt:                                          # elide when no CollecTRI targets
        L.append("\n### Downstream targets (CollecTRI) {#collectri}\n")
        L.append(capped_table(["Target", "Regulation"],
                              [(t.get("target"), t.get("regulation")) for t in dt],
                              30, total=b.get("downstream_count"), noun="targets"))
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
        L.append(f"Targeting {b.get('symbol')} by miRDB confidence (max_score; target_count "
                 f"= how many genes the miRNA targets in total — lower means more specific).\n")
        L.append(capped_table(["miRNA", "Max score", "Avg score", "miRNA target_count"],
                              [(m["id"], m.get("max_score"), m.get("avg_score"),
                                m.get("target_count")) for m in b.get("mirna_regulators", [])],
                              None, total=n_mir, noun="miRNAs"))
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
        L.append(f"Phase ≥1, by development phase (incl. off-target/promiscuous "
                 f"compounds).{patent_note}\n")
        L.append(capped_table(["Molecule", "Name", "Phase", "Patents"],
                              [(m["id"],
                                links.maybe_link(m.get("name"), links.drug_url(chembl_id=m["id"], name=m.get("name"))),
                                m.get("phase"),
                                f"{m['patent_count']:,}" if m.get("patent_count") else "")
                               for m in mols],
                              ROW_CAP, total=mc, noun="molecules by phase"))
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
        L.append(f"Drug × variant × indication, from {b.get('civic_predictive_total', 0)} "
                 f"curated evidence items" + (f"; also {extra}" if extra else "") + ". "
                 + CIVIC_LEGEND + "\n")
        L.append(capped_table(["Variant", "Therapy", "Indication", "Effect", "Level", "CIViC"],
                              [(r["profile"],
                                links.maybe_link(therapy_label(r["therapy"]), links.drug_url(name=therapy_label(r["therapy"]))),
                                links.maybe_link(r["disease"], links.disease_url(name=r["disease"])),
                                r["significance"],
                                f"CIViC {r['level']}" if r.get("level") else "",
                                f"EID{r['evidence_id']}"
                                + (f" +{r['n']-1}" if r.get("n", 1) > 1 else ""))
                               for r in ce],
                              None, total=b.get("civic_association_total", 0),
                              noun="predictive associations by evidence level"))
        if any(r.get("n", 1) > 1 for r in ce):
            L.append("\n*CIViC column: a representative evidence item (EID); "
                     "+N = additional CIViC items supporting the same association.*")

    # CIViC curated clinical variants — the named-variant catalogue (with variant
    # type) beneath the predictive evidence above. Compact name list.
    cvs = b.get("civic_variants") or []
    if cvs:
        total = b.get("civic_variant_total", len(cvs))
        names = ", ".join(f"{v['name']}" + (f" ({v['type']})" if v.get("type") else "")
                          for v in cvs[:ROW_CAP])
        L.append(f"\n### CIViC curated variants ({total}) {{#civic-variants}}\n")
        L.append("Named clinical variants curated in CIViC for this gene "
                 "(the catalogue beneath the predictive evidence above):\n")
        L.append(names + (f" …(+{total - ROW_CAP} more)" if total > ROW_CAP else "") + ".")

    pg = b.get("pharmgkb", [])
    # Don't surface is_vip: it's broken upstream (always true — ACTB/GAPDH/TTN
    # all report VIP=true alongside CYP2D6), so it carries no signal. CPIC is
    # real and discriminates, so keep it.
    L.append(f"\n**PharmGKB:** {len(pg)} entr{'y' if len(pg)==1 else 'ies'}"
             + (f" (CPIC guideline: {'yes' if pg[0].get('cpic_guideline') == 'true' else 'no'})"
                if pg else ""))

    _variant_link = links.variant_link

    # PharmGKB clinical annotations — variant × association × phenotype × drug,
    # ranked by PharmGKB evidence level (1 strongest → 4). Same dataset and
    # presentation as the drug page's clinical-annotations table, mirrored from
    # the gene's side (here the gene is fixed and we list the drugs).
    pgc = b.get("pharmgkb_clinical") or []
    if pgc:
        pgc = sorted(pgc, key=lambda c: str(c.get("level_of_evidence") or "9"))
        L.append("\n### PharmGKB clinical annotations {#pharmgkb-clinical}\n")
        L.append(capped_table(
            ["Variant", "Association", "Level", "Drugs", "Phenotypes"],
            [(links.maybe_link(c.get("variant") or "", _variant_link(c.get("variant"))),
              c.get("type"), c.get("level_of_evidence"),
              links.link_csv(c.get("chemicals"), lambda s: links.drug_url(name=s)),
              c.get("phenotypes"))
             for c in pgc],
            ROW_CAP, noun="annotations by evidence level (1 strongest → 4)"))

    # PharmGKB variant pages — variant-level aggregations with PharmGKB's
    # composite score + count of clinical annotations.
    pgv = b.get("pharmgkb_variant") or []
    if pgv:
        L.append("\n### PharmGKB variants {#pharmgkb-variants}\n")
        L.append(capped_table(
            ["Variant", "Genes", "Level", "Score", "#Clin annots", "Drugs"],
            [(links.maybe_link(v.get("name") or "", _variant_link(v.get("name"))),
              links.link_csv(v.get("gene_symbols"), lambda s: links.gene_url(symbol=s)),
              v.get("level_of_evidence"),
              v.get("score"), v.get("clinical_annotation_count"),
              links.link_csv(v.get("associated_drugs"), lambda s: links.drug_url(name=s)))
             for v in pgv],
            ROW_CAP, noun="variants"))

    # PharmGKB per-publication variant annotations — the raw evidence layer
    # (one finding per row, each a plain-English sentence + PMID), keyed by drug.
    pva = b.get("pharmgkb_var_annotation") or []
    if pva:
        L.append("\n### PharmGKB variant annotations {#pharmgkb-var-annotations}\n")
        L.append("Per-publication pharmacogenomic findings for this gene "
                 "(significant associations).\n")
        L.append(capped_table(
            ["Variant", "Drug(s)", "Category", "Finding", "PMID"],
            [(links.maybe_link(v.get("variant") or "", _variant_link(v.get("variant"))),
              links.link_csv((v.get("drugs") or "").replace(";", ", "),
                             lambda s: links.drug_url(name=s)),
              v.get("category"), v.get("sentence"),
              links.maybe_link(f"PMID:{v['pmid']}",
                               f"https://pubmed.ncbi.nlm.nih.gov/{v['pmid']}/")
              if v.get("pmid") else "")
             for v in pva],
            ROW_CAP, noun="published findings"))

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
        L.append("CPIC / DPWG genotype-guided dosing for this gene (drug × pharmacogene).\n")
        # Same column structure as the drug page's mirror of this table; only the
        # partner column differs (there it's Gene(s), here it's Drug).
        L.append(capped_table(
            ["Guideline", "Source", "Drug", "Dosing?", "Recommendation?"],
            [((g.get("name") or "")[:70],
              g.get("source"),
              links.link_csv(g.get("chemical_names"), lambda s: links.drug_url(name=s)),
              "yes" if g.get("has_dosing_info") else "",
              "yes" if g.get("has_recommendation") else "")
             for g in pgg_sorted],
            ROW_CAP, noun="dosing guidelines"))
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
    L.append("Via MONDO — disease-level, not drug-specific.\n")
    L.append(capped_table(["NCT", "Phase", "Status", "Title"],
                          # phase_label (audit #8): biobtree emits 'NaN' for trials with
                          # no interventional phase; render it as 'Not specified', not
                          # a leaked 'nan'/'NAN'.
                          [(t["id"], phase_label(t.get("phase")), t.get("status"),
                            (t.get("title") or "").strip()) for t in ct],
                          ROW_CAP, total=b.get("disease_trial_count"), noun="trials"))
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
        L.append("\n### FANTOM5 promoters {#fantom5-promoters}\n")
        L.append(capped_table(["Promoter ID", "TPM avg", "Samples expressed"],
                              [(p["id"], p.get("tpm_average"), p.get("samples_expressed")) for p in fp],
                              10, noun="alternative TSS"))
    # Tissue name as plain text — the UBERON/CL id is shown in its own column;
    # we don't link every id (consistency: links are reserved for primary
    # targets, not decorative on every row).
    # Bgee expression score (0-100, higher = more highly expressed; derived from
    # the per-condition rank). We drop the raw Rank column — it's inversely
    # redundant with the score and the pair reads as a contradiction (high score
    # next to a mid-looking rank) when it isn't.
    L.append("\n### Top tissues by expression {#tissue-expression}\n")
    # Bgee tissues are score-ranked with a long tail (median gene ~243 entities);
    # show the top 100 by expression score.
    L.append(capped_table(["Tissue", "Anatomy ID", "Expression score", "Quality"],
                          [(t.get("tissue") or "", t.get("anatomy_id") or "",
                            t.get("score"), t.get("quality")) for t in b.get("top_tissues", [])],
                          100, total=b.get("tissue_count"),
                          noun="tissues by Bgee expression score (0-100, higher = more expressed)"))
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
                             "yes" if e.get("is_marker") else "no",
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
        # intOGen's role is inferred from the statistical *pattern* of somatic
        # mutations across tumour cohorts — a computational call that can diverge
        # from a gene's established/curated role (e.g. a known tumour suppressor
        # may carry an "activating" label in a given cohort). Read it alongside
        # the curated evidence above, not in place of it.
        L.append("*intOGen's role is a computational inference from cohort "
                 "mutation patterns and may diverge from the gene's curated/"
                 "literature role; read it alongside the evidence above.*")
    return "\n\n".join(L)


def r_diseases(b):
    L = ["## Disease associations", ""]
    L.append(f"**OMIM:** gene `{', '.join(b.get('gene_omim', []))}` | "
             f"disease phenotypes: {', '.join(b.get('disease_omim', [])[:40])}")
    L.append("\n### GenCC curated gene-disease {#gencc}\n")
    # Collapse to one row per disease (audit #13: GenCC has many submissions per
    # gene–disease pair — FANCB ×4); keep the strongest classification.
    from atlas.render_common import gencc_rank
    _gc, _gc_n = {}, {}
    for g in b.get("gencc", []):
        d = g.get("disease") or ""
        if not d:
            continue
        _gc_n[d] = _gc_n.get(d, 0) + 1     # submissions per disease (surfaced as Records)
        if d not in _gc or gencc_rank(g.get("classification")) > gencc_rank(_gc[d].get("classification")):
            _gc[d] = g
    _gc_rows = sorted(_gc.values(),
                      key=lambda g: -gencc_rank(g.get("classification")))
    # "Records" mirrors the disease-page GenCC table (submission count after the
    # one-row-per-disease dedup), so both sides are transparent about the collapse.
    L.append(capped_table(["Disease", "Classification", "Inheritance", "Records"],
                          [(links.maybe_link(g.get("disease"), links.disease_url(name=g.get("disease"))),
                            g.get("classification"), g.get("inheritance"),
                            str(_gc_n.get(g.get("disease") or "", 0)) if _gc_n.get(g.get("disease") or "", 0) > 1 else "")
                           for g in _gc_rows],
                          ROW_CAP, noun="curated gene-disease records"))

    # ClinGen Gene-Disease Validity — expert-panel curated relationship strength.
    # Distinct from GenCC: ClinGen is the canonical authority for gene-disease
    # validity classifications used by clinical variant labs (ACMG/AMP guidelines).
    cgv = b.get("clingen_validity") or []
    if cgv:
        L.append("\n### ClinGen Gene-Disease Validity {#clingen}\n")
        L.append("Expert-panel classifications — Definitive > Strong > Moderate > "
                 "Limited > Disputed > Refuted.\n")
        L.append(capped_table(["Disease", "Classification", "Inheritance"],
                              [(links.maybe_link(c.get("disease"), links.disease_url(name=c.get("disease"))),
                                c.get("classification"), c.get("moi"))
                               for c in cgv],
                              ROW_CAP, noun="validity classifications"))

    # PanelApp (Genomics England) — diagnostic gene panels the gene is on, with
    # its rating per panel. Green = diagnostic-grade; the clinical "which test
    # would include this gene" signal that ClinGen/GenCC (per-disease) don't give.
    pa = b.get("panelapp") or []
    if pa:
        green = b.get("panelapp_green") or 0
        L.append("\n### Genomics England panels (PanelApp) {#panelapp}\n")
        L.append(f"Clinical diagnostic gene panels including this gene "
                 f"({len(pa)} panels, {green} green/diagnostic-grade). "
                 "Confidence: green = diagnostic, amber = moderate.\n")
        L.append(capped_table(["Panel", "Rating"],
                              [(p.get("panel"), p.get("confidence")) for p in pa],
                              ROW_CAP, noun="panels"))
    L.append(f"\n**Mondo ({len(b.get('mondo', []))}):** "
             + ", ".join(f"{links.maybe_link(m.get('name'), links.disease_url(mondo_id=m['id'], name=m.get('name')))} ({m['id']})"
                         for m in b.get("mondo", [])[:ROW_CAP]))
    L.append(more_line(len(b.get("mondo", [])), ROW_CAP))
    L.append(f"\n**Orphanet ({len(b.get('orphanet', []))}):** "
             + ", ".join(f"{o.get('name') or ''} ({o['id']})" for o in b.get("orphanet", [])[:ROW_CAP]))
    L.append(more_line(len(b.get("orphanet", [])), ROW_CAP))
    # No relevance/frequency ordering on the gene→HPO route, so don't claim
    # "(top)" — these are the first 30 by HPO id. (Frequency-ranked clinical
    # features live on the disease pages.)
    L.append("\n### HPO phenotypes {#hpo}\n")
    L.append(capped_table(["HPO", "Term"], [(h["id"], h.get("name")) for h in b.get("hpo", [])],
                          ROW_CAP, total=b.get("hpo_total"), noun="phenotypes (HPO-id order)"))
    L.append("\n### GWAS associations {#gwas-assoc}\n")
    L.append(capped_table(["Study", "Trait", "p-value"],
                          [(g["id"], g.get("trait"), pval(g.get("p_value"))) for g in b.get("gwas", [])],
                          30, total=b.get("gwas_total"), noun="associations"))

    # EFO canonical trait names — normalized vocabulary for the free-text
    # GWAS-catalog traits above. The same disease can appear under multiple
    # synonymous GWAS-catalog spellings; EFO collapses them. Useful for
    # AI agents doing cross-resource joins on trait IDs.
    efo = b.get("efo_traits") or []
    if efo:
        L.append("\n### EFO canonical traits (from GWAS) {#efo}\n")
        L.append(capped_table(["EFO ID", "Trait name"],
                              [(e["id"], e.get("name")) for e in efo],
                              ROW_CAP, noun="canonical traits"))

    # MeSH disease descriptors — NLM's controlled disease vocabulary.
    # Tree numbers (e.g. C04.700.600) classify into MeSH categories
    # (C04=Neoplasms, C16=Congenital, C18=Nutritional/Metabolic, etc.) —
    # useful for grouping diseases at a coarser level.
    # Drop nameless descriptors — biobtree occasionally returns a MeSH id with
    # no name (e.g. `C538339`); a row with a bare id misleads LLMs.
    mesh = [m for m in (b.get("mesh_descriptors") or []) if (m.get("name") or "").strip()]
    if mesh:
        L.append("\n### MeSH disease descriptors {#mesh}\n")
        L.append(capped_table(["Descriptor", "Name", "Tree numbers"],
                              [(m["id"],
                                m["name"] + (" *(supp.)*" if m.get("is_supplementary") else ""),
                                "; ".join(m.get("tree_numbers") or []))
                               for m in mesh],
                              ROW_CAP, noun="disease descriptors"))
    return "\n".join(L)


# --- Human Protein Atlas (§13) — derived renderers reading the shared bundle,
# slotted by pipeline.render_all into the Protein / Gene-structure / Disease zones.

def _pv(s):
    try:
        return float(s)
    except (TypeError, ValueError):
        return 9.9


def r_hpa_protein(bundle):
    """HPA protein facts → Protein zone: subcellular location (IHC-imaging, shown
    alongside the UniProt CC narrative), protein classes, protein-evidence tier."""
    h = (bundle.get("13") or {}).get("hpa") or {}
    if not h:
        return ""
    L = ["## Human Protein Atlas {#hpa-protein}", ""]
    cls = h.get("protein_classes") or []
    if cls:                                          # noise already filtered (s13) — show all
        L.append(f"**Protein class:** {', '.join(cls)}")
    if h.get("protein_evidence"):
        L.append(f"\n**Protein evidence:** {h['protein_evidence']}")
    loc, seen = [], set()                            # main + additional, one merged list
    for x in (h.get("subcellular_main") or []) + (h.get("subcellular_additional") or []):
        if x and x not in seen:
            seen.add(x)
            loc.append(x)
    if loc:
        L.append(f"\n**Subcellular location (HPA, imaging):** {'; '.join(loc)}")
    if h.get("secretome_location"):
        L.append(f"\n**Secretome:** {h['secretome_location']}")
    # Antibody-staining reliability — qualifies how trustworthy the HPA protein
    # data above is (IH = immunohistochemistry/tissue, IF = immunofluorescence/
    # localization; enhanced > supported > approved > uncertain).
    rel = [f"{lab} {h[k]}" for k, lab in (("reliability_ih", "IH"),
                                          ("reliability_if", "IF")) if h.get(k)]
    if rel:
        L.append(f"\n**Antibody reliability (HPA):** {' · '.join(rel)}")
    return "\n".join(L)


def r_hpa_expression(bundle):
    """HPA per-tissue/cell expression (nTPM, + IHC protein level when present).
    Renders as a SUBSECTION of "Expression profiles" (an H4 after demotion,
    alongside the Bgee/FANTOM5/SCXA expression blocks), not its own top-level
    section — HPA RNA expression is expression data, so it belongs there."""
    b = bundle.get("13") or {}
    exp = b.get("hpa_expression") or []
    spec = (b.get("hpa") or {}).get("rna_tissue_specificity")
    profiled = b.get("hpa_expression_profiled") or 0
    if not exp:
        # Nothing detected above the HPA threshold — state it concisely rather than
        # render dozens of nTPM=0 rows (the ZNF735 case). Elide entirely if HPA
        # profiled nothing at all.
        if not profiled:
            return ""
        return ("### HPA expression {#hpa-expression}\n\n"
                f"RNA tissue specificity (HPA): **{spec or 'Not detected'}** — no "
                f"expression detected above threshold across {profiled:,} profiled "
                f"tissues/cells.")
    L = ["### HPA expression {#hpa-expression}", ""]
    if spec:
        L.append(f"RNA tissue specificity: **{spec}**\n")
    if any(r.get("protein_level") for r in exp):
        L.append(capped_table(["Tissue / cell", "Axis", "nTPM", "Protein (IHC)"],
                              [(r.get("entity"), r.get("axis"), r.get("ntpm"), r.get("protein_level"))
                               for r in exp],
                              ROW_CAP, total=b.get("hpa_expression_total"), noun="HPA entities by nTPM"))
    else:
        L.append(capped_table(["Tissue / cell", "Axis", "nTPM"],
                              [(r.get("entity"), r.get("axis"), r.get("ntpm")) for r in exp],
                              ROW_CAP, total=b.get("hpa_expression_total"), noun="HPA entities by nTPM"))
    return "\n".join(L)


def r_hpa_cancer(bundle):
    """HPA cancer prognostics + RNA cancer specificity → Disease zone, alongside
    intOGen drivers and CIViC."""
    b = bundle.get("13") or {}
    h = b.get("hpa") or {}
    path = b.get("hpa_pathology") or []
    spec = h.get("rna_cancer_specificity")
    if not path and not spec:
        return ""
    L = ["## HPA cancer prognostics {#hpa-cancer}", ""]
    if spec:
        L.append(f"**RNA cancer specificity:** {spec}\n")
    if path:
        L.append("Human Protein Atlas survival analysis — favorable = higher "
                 "expression predicts better survival.\n")
        L.append(capped_table(["Cancer", "Prognostic", "p-value"],
                              [(p.get("cancer"), p.get("prognostic_type"), p.get("p_value"))
                               for p in sorted(path, key=lambda p: _pv(p.get("p_value")))],
                              ROW_CAP, noun="prognostic cancer types"))
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
