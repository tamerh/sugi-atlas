#!/usr/bin/env python3
"""Deterministic markdown renderer for drug-section bundles — NO model.
Mirrors gene/render.py + disease/render.py: one r_* fn per section + a RENDER
dict. Every fact comes verbatim from the bundle."""
from atlas.render_common import table, fnum, is_ontology_id, display_name, phase_label, more_line
from atlas.civic import therapy_label, LEGEND as CIVIC_LEGEND
from atlas.page import links

# Display cap for curated, naturally bounded tables (see the gene renderer's note).
# The bioactivity table is collector-capped (top 100 distinct by potency); CIViC /
# indications / clinical annotations render in full (already collector-bounded).
ROW_CAP = 60


def _i(n):
    return f"{n:,}" if isinstance(n, int) else (n if n is not None else "")


# ChEBI roles mix pharmacological roles with chemical / biological / application
# ones ("environmental contaminant", "xenobiotic", "human metabolite", "solvent").
# Split them so the medical roles read as a drug-class signal and the rest don't
# masquerade as one (atorvastatin's only roles were "contaminant, xenobiotic").
_NON_PHARMA_ROLE = ("contaminant", "xenobiotic", "pollutant", "metabolite",
                    "solvent", "fertilizer", "fertiliser", "food", "dye",
                    "reagent", "nmr ", "fuel", "detergent", "buffer",
                    "antifoaming", "propellant", "refrigerant", "marine",
                    "algal", "plant ", "indicator", "flavour", "flavor",
                    "fragrance", "cosmetic")


def _split_chebi_roles(roles):
    """(pharmacological_roles, other_roles) — other = chemical/environmental/
    metabolite roles matched by _NON_PHARMA_ROLE keywords."""
    pharma, other = [], []
    for r in roles or []:
        (other if any(k in (r or "").lower() for k in _NON_PHARMA_ROLE)
         else pharma).append(r)
    return pharma, other


def _render_chebi_roles(roles):
    """Two lines: medical roles first, other (non-pharmacological) roles second.
    Either is omitted when empty."""
    pharma, other = _split_chebi_roles(roles)
    out = []
    if pharma:
        out.append(f"\n**Pharmacological roles (ChEBI):** {', '.join(pharma)}.")
    if other:
        out.append(f"\n**Other ChEBI roles (chemical / environmental):** "
                   f"{', '.join(other)}.")
    return out


def r_drug_ids(b):
    L = ["## Drug identity and classification", ""]
    rows = [
        ("ChEMBL ID", b.get("chembl_id")),
        ("Name", display_name(b.get("canonical_name"))),  # audit #12: de-SHOUT
        ("Type", b.get("molecule_type")),
        ("Max phase", b.get("max_phase")),
        ("FDA approved", None if b.get("is_fda_approved") is None
         else ("yes" if b.get("is_fda_approved") else "no")),
        ("PubChem CID", b.get("pubchem_cid")),
        ("ChEBI", b.get("chebi_id")),
        ("ATC", ", ".join(b.get("atc_codes") or [])),
        ("Molecular formula", b.get("molecular_formula")),
        ("Molecular weight", b.get("molecular_weight")),
        ("InChIKey", b.get("inchi_key")),
    ]
    L.append(table(["Field", "Value"], [(k, v) for k, v in rows if v not in (None, "")]))
    if b.get("smiles"):
        L.append(f"\n**SMILES:** `{b['smiles']}`")
    if b.get("iupac_name"):
        L.append(f"\n**IUPAC name:** {b['iupac_name']}")
    if b.get("chebi_definition"):
        L.append(f"\n**ChEBI definition:** {b['chebi_definition']}")
    # ChEBI functional roles (folded in from the former standalone Pharmacology
    # section; ATC already lives in the table above). Split pharmacological from
    # chemical/environmental roles so the latter don't read as a drug class.
    L += _render_chebi_roles(b.get("chebi_roles") or [])
    an = b.get("alt_names") or []
    if an:
        L.append(f"\n**Also known as:** {', '.join(an[:12])}")
    # salt-form navigation
    if b.get("parent_chembl"):
        L.append(f"\n*Salt/anhydrous form of parent* `{b['parent_chembl']}`.")
    if b.get("child_chembls"):
        L.append(f"\n*Parent form; salt/anhydrous children:* "
                 + ", ".join(f"`{c}`" for c in b["child_chembls"]))
    # Patent footprint — SureChEMBL compound mentions (folded in from the
    # former standalone Patent literature section). We report distinct patent
    # *families* (the honest dedup metric — one invention across many
    # jurisdictions) alongside the raw mention count, which is typically
    # dominated by one promiscuous/reference structure (quantified here).
    # (assignee + CPC/IPC landscape is gated on a representative sample — needs
    # biobtree #27 date-sort/facets.)
    pt = b.get("patent_total") or 0
    if pt:
        bd = b.get("patent_compound_breakdown") or []
        n_rec = len(bd) or len(b.get("patent_compound_ids") or [])
        fam = b.get("patent_family_total") or 0
        dom = ""
        if bd and bd[0].get("patent_count"):
            top = bd[0]["patent_count"]
            pct = round(100 * top / pt) if pt else 0
            if len(bd) > 1 and pct >= 60:
                dom = (f" One matched structure accounts for {_i(top)} ({pct}%) "
                       f"of the total.")
        # Lead with the family count (the honest dedup metric the caveat itself
        # endorses); demote the raw mention figure to the parenthetical.
        head = (f"**{_i(fam)} distinct patent families** "
                f"({_i(pt)} SureChEMBL compound mentions)" if fam
                else f"**{_i(pt)} SureChEMBL compound mentions**")
        L.append(f"\n**Patent coverage:** {head}, from {n_rec} matched compound "
                 f"structure(s).{dom} Mentions count patents naming the compound "
                 f"(not distinct inventions), so promiscuous / reference molecules "
                 f"inflate the mention figure — families are the dedup metric.")
    return "\n".join(L)


def r_targets(b):
    L = ["## Targets", ""]
    pt = b.get("primary_targets") or []
    if pt:
        src = b.get("primary_source") or "gtopdb"
        label = "GtoPdb curated mechanism" if src == "gtopdb" else "ChEMBL bioactivity"
        L.append(f"**Primary targets ({label}):** the *Cancer dependency* column is "
                 f"the DepMap CRISPR fitness signal (% of screened cell lines "
                 f"dependent on the target).\n")
        def _dep(t):
            pct = t.get("dep_pct")
            if pct in (None, ""):
                return ""
            tags = []
            if t.get("dep_selective"):
                tags.append("strongly selective")
            if t.get("dep_common_essential"):
                tags.append("common-essential")
            return f"{fnum(pct)}%" + (f" ({', '.join(tags)})" if tags else "")
        def _action(a):
            # GtoPdb returns the literal string 'None'/'Unknown' for ligands with
            # no defined functional action — render blank, not a fake null cell.
            a = (a or "").strip()
            return "" if a.lower() in ("none", "unknown") else a
        L.append(table(["Gene", "Target", "Action", "pAffinity", "Cancer dependency", "UniProt"],
                       [(links.maybe_link(t.get("gene_symbol") or "", links.gene_url(symbol=t.get("gene_symbol"), hgnc_id=t.get("hgnc_id"))),
                         t.get("target_name") or "",
                         _action(t.get("action")), fnum(t.get("affinity")) if t.get("affinity") not in (None, "") else "",
                         _dep(t), t.get("uniprot") or "") for t in pt]))
    bc = b.get("bioactivity_target_count") or 0
    if bc:
        # Count only — the sample of names dumped unsorted, off-target screening
        # hits (aspirin → estrogen receptor) with duplicates; the same targets
        # appear, deduped + with potencies, in the ChEMBL bioactivities table.
        L.append(f"\n**Broader ChEMBL bioactivity targets: {_i(bc)}** "
                 "(assay-derived screening hits — not curated targets; per-target "
                 "potencies are in the ChEMBL bioactivities table below).")
    if not pt and not bc:
        L.append("*No target linkage available.*")
    return "\n".join(L)


def r_bioactivity(b):
    L = ["## Bioactivity", ""]
    ca = b.get("activities") or []
    if not ca:
        L.append("*No ChEMBL bioactivity rows at pChembl ≥ 5 "
                 "(expected for biologics / antibodies).*")
        return "\n".join(L)
    L.append(f"**ChEMBL activities: {_i(b.get('potent_count'))} potent at "
             f"pChembl ≥ 5 of {_i(b.get('activity_total'))} total. Top {_i(len(ca))} "
             f"distinct by potency (10 = 0.1 nM, 6 = 1 µM):**\n")
    # Three target columns — the assayed PROTEIN (name + UniProt, linked) and its
    # encoding GENE (symbol, linked to the Atlas page) kept distinct, not conflated.
    # value+unit merged. Activity ID retained as machine-readable provenance (it
    # resolves at the ChEMBL API to the assay + source paper).
    def _val(r):
        return f"{fnum(r.get('value'))} {r.get('unit') or ''}".strip()
    def _protein(r):
        # non-human targets carry an organism (sheep/cattle COX-1 orthologs) —
        # label it so the blank Gene cell reads as "ortholog", not a missing map.
        nm, org = r.get("protein_name") or "", r.get("organism")
        return f"{nm} ({org})" if (nm and org) else nm
    L.append(table(["Protein", "UniProt", "Gene", "pChembl", "Type", "Value", "Activity ID"],
                   [(_protein(r),
                     links.maybe_link(r.get("uniprot") or "",
                                      links.uniprot_url(r.get("uniprot"))),
                     links.maybe_link(r.get("target_symbol") or "",
                                      links.gene_url(symbol=r["target_symbol"]) if r.get("target_symbol") else None),
                     fnum(r.get("pchembl")), r.get("type"), _val(r), r.get("id")) for r in ca]))
    return "\n".join(L)


def r_indications(b):
    inds = b.get("indications") or []
    # Keep only rows with a real disease name (drop raw EFO/MeSH/MP ids), then
    # dedup by (name, MONDO) — the source repeats the same disease across xrefs.
    named, seen = [], set()
    for i in inds:
        nm = (i.get("name") or "").strip()
        if not nm or is_ontology_id(nm):
            continue
        key = (nm.lower(), i.get("mondo_id") or "")
        if key in seen:
            continue
        seen.add(key)
        named.append(i)
    # Tier on the per-indication `approved` flag (computed in s04 via
    # atlas.indication): phase 4, or an anticancer drug at phase 3 vs a cancer
    # (imatinib→CML). Heavily-trialed drugs (aspirin's phase-2/3 cancers) stay
    # investigational — they must NOT read as approved.
    # Investigational tier floors at phase 2 to match the disease-side indication
    # index (batch._build_indication_index drops phase ≤1 as exploratory noise —
    # 90% fail). Keeping the same floor on both sides means a drug↔disease trial
    # link shown here also appears on that disease's page, and vice-versa.
    approved, trials = [], []
    for i in named:
        if i.get("approved"):
            approved.append(i)
        elif 2 <= (i.get("max_phase") or 0) <= 3:
            trials.append(i)

    def _tbl(items, col0):
        return table([col0, "Phase", "MONDO", "EFO"],
                     [(links.maybe_link(i.get("name"),
                                        links.disease_url(mondo_id=i.get("mondo_id"), name=i.get("name"))),
                       i.get("max_phase"), i.get("mondo_id") or "", i.get("efo_id") or "")
                      for i in sorted(items, key=lambda i: -(i.get("max_phase") or 0))])

    L = ["## Indications", ""]
    if approved:
        L += [f"**{_i(len(approved))} approved indication"
              f"{'s' if len(approved) != 1 else ''}.** FDA phase 4, plus an anticancer "
              "drug's labelled cancer uses (which ChEMBL often logs at phase 3).", "",
              _tbl(approved, "Indication")]
    if trials:
        L += ["",
              f"**{_i(len(trials))} disease{'s' if len(trials) != 1 else ''} in clinical "
              f"trials** (phase 2–3, investigational — *not* approved indications). Highest "
              "ChEMBL trial phase per disease; a non-cancer approved use is occasionally "
              "logged at phase 3 here.", "",
              _tbl(trials, "Disease (in trials)")]
    if not named:
        n = b.get("indication_count") or 0
        L.append(f"**{_i(n)} indication record{'s' if n != 1 else ''}** "
                 f"carr{'y' if n != 1 else 'ies'} no mapped disease name "
                 "(EFO/MeSH-only); none shown.")
    else:
        dropped = len(inds) - len(named)
        if dropped > 0:
            L.append(f"\n*{dropped} further indication record"
                     f"{'s' if dropped != 1 else ''} had no mapped disease name "
                     "(EFO/MeSH-only) or were duplicates, and are omitted.*")
    return "\n".join(L)


def r_related_molecules(b):
    L = ["## Related molecules", ""]
    rm = b.get("related_molecules") or []
    if not rm:
        L.append("*No competitor molecules sharing a primary target "
                 "(ChEMBL phase ≥2 or PubChem drug-class).*")
        return "\n".join(L)
    # Brief method/datasets note — this same-mechanism view is distinctive to
    # Sugi Atlas, so spell out how it's built.
    L.append("*Molecules sharing ≥1 of this drug's curated primary targets, "
             "merged from two biobtree sources and ranked by shared-target "
             "count, then clinical phase: **ChEMBL** clinical-stage candidates "
             "(development phase ≥2) and **PubChem** drug-class bioactivity "
             "(approved / known drugs acting on the target). Deduplicated by drug "
             "name; the drug's own salt forms are excluded. Note: for a drug with "
             "few primary targets a shared-target match can reflect off-target / "
             "promiscuous binding rather than the same therapeutic mechanism — "
             "the phase ordering surfaces bona-fide therapeutics first.*")
    L.append(f"\n**{_i(b.get('competitor_count'))} molecules share ≥1 primary "
             f"target. Top {len(rm)} by shared-target count:**\n")

    def _status(r):
        ph = r.get("phase") or 0
        if ph:
            return f"Phase {ph}" + (" (approved)" if ph >= 4 or r.get("fda") else "")
        return "Approved" if r.get("fda") else "—"

    L.append(table(["Molecule", "Source", "Status", "Shared targets"],
                   [(links.maybe_link(r.get("name") or "—", links.drug_url(name=r.get("name"))),
                     " + ".join(r.get("sources") or []),
                     _status(r),
                     ", ".join(r.get("shared_targets") or [])) for r in rm]))
    return "\n".join(L)


def r_target_pathways(b):
    L = ["## Target pathways", ""]
    tg = b.get("target_genes") or []
    tp = b.get("top_pathways") or []
    go = b.get("top_go_bp") or []
    if not tg and not tp and not go:
        L.append("*No target-pathway data for this drug "
                 "(no mapped target genes).*")
        return "\n".join(L)
    tg_md = ", ".join(links.maybe_link(g, links.gene_url(symbol=g)) for g in tg)
    L.append(f"**Aggregated over {len(tg)} target gene(s): {tg_md}.**")
    if tp:
        L += ["", "### Top Reactome pathways {#target-reactome}", "",
              f"{_i(b.get('pathway_count'))} total, by targets touching each:", "",
              table(["Pathway", "Targets", "Genes"],
                    [(p.get("name") or p.get("id") or "",
                      p.get("gene_count"),
                      ", ".join(links.maybe_link(g, links.gene_url(symbol=g)) for g in (p.get("genes") or [])))
                     for p in tp])]
    if go:
        L += ["", "### Dominant GO biological processes {#target-go}", "",
              table(["GO term", "Targets"],
                    [(g.get("name") or g.get("id") or "",
                      g.get("target_count")) for g in go])]
    return "\n".join(L)


def r_pharmacogenomics(b):
    L = ["## Pharmacogenomics", ""]
    g = b.get("guidelines") or []
    pa = b.get("pharmgkb_chemical_id")
    ca, va = b.get("clinical_annotation_count"), b.get("variant_annotation_count")
    if not g and not pa:
        L.append("*No PharmGKB pharmacogenomic data curated for this drug.*")
        return "\n".join(L)
    if g:
        # Sort CPIC > DPWG > CPNDS > others, then by name — same ordering the gene
        # page uses for its mirror of this table.
        order = {"CPIC": 0, "DPWG": 1, "CPNDS": 2}
        g = sorted(g, key=lambda r: (order.get(r.get("source"), 99), r.get("name") or ""))
        L.append(f"**PharmGKB dosing guidelines ({b.get('guideline_count')}) — CPIC / "
                 f"DPWG genotype-guided dosing for this drug (drug × pharmacogene):**\n")
        L.append(table(["Guideline", "Source", "Gene(s)", "Dosing?", "Recommendation?"],
                       [((r.get("name") or r.get("id") or "")[:70],
                         r.get("source"),
                         links.link_csv(r.get("genes"), lambda s: links.gene_url(symbol=s)),
                         "yes" if r.get("has_dosing") else "",
                         "yes" if r.get("has_recommendation") else "") for r in g[:ROW_CAP]]))
        L.append(more_line(len(g), ROW_CAP))
    # Per-variant CLINICAL annotations for this drug — surfaced directly now, not
    # just counted: variant × gene × association type × phenotype × evidence level
    # (resolved via the drug's PGx genes; see s09). We don't tease the deeper
    # per-publication variant-annotation count: those records aren't accessible in
    # biobtree yet, and our aim is to aggregate data, not point users elsewhere.
    cl = b.get("clinical_annotations") or []

    def _variant_link(v):
        return f"https://www.ncbi.nlm.nih.gov/snp/{v}" if v and v.startswith("rs") else None

    if cl:
        L += ["", f"**Pharmacogenomic clinical annotations ({_i(len(cl))})** — variant × "
              "gene × association × phenotype, ranked by PharmGKB evidence level "
              "(1 strongest → 4):", "",
              table(["Gene", "Variant", "Association", "Phenotype", "Level"],
                    [(links.maybe_link(r.get("gene") or "",
                                       links.gene_url(symbol=r["gene"]) if r.get("gene") else None),
                      links.maybe_link(r.get("variant") or "", _variant_link(r.get("variant"))),
                      r.get("type"), r.get("phenotypes"), r.get("level")) for r in cl])]
    elif pa and (ca or va):
        counts = (f"{_i(ca) or 0} clinical and {_i(va) or 0} variant annotation(s) "
                  f"for this drug (see [PharmGKB](https://www.pharmgkb.org/chemical/{pa}))")
        L.append(f"\nPharmGKB also curates {counts}." if g else
                 f"*No CPIC/DPWG dosing guideline, but PharmGKB curates {counts}.*")
    elif not g and pa:
        L.append("*No CPIC/DPWG dosing guideline or drug-level clinical/variant "
                 "annotations in PharmGKB for this molecule.*")
    return "\n".join(L)


def r_clinical_evidence(b):
    L = ["## Clinical evidence (CIViC)", ""]
    ce = b.get("civic_evidence") or []
    if not ce:
        L.append("*No CIViC predictive evidence "
                 "(expected for non-precision-medicine drugs).*")
        return "\n".join(L)
    etc = b.get("civic_evidence_type_counts") or {}
    extra = ", ".join(f"{n} {k.lower()}" for k, n in etc.items() if k != "Predictive")
    L.append(f"**Variant × therapy × indication "
             f"({_i(b.get('civic_association_total'))} predictive associations from "
             f"{_i(b.get('civic_predictive_total'))} curated evidence items"
             + (f"; also {extra}" if extra else "") + "):** " + CIVIC_LEGEND + "\n")
    # Same column order as the gene page (Variant | Therapy | Indication | Effect);
    # the therapy is linked too (it's often a partner drug, not this one).
    L.append(table(["Variant", "Therapy", "Indication", "Effect", "Level", "CIViC"],
                   [(r["profile"],
                     links.maybe_link(therapy_label(r["therapy"]), links.drug_url(name=therapy_label(r["therapy"]))),
                     links.maybe_link(r["disease"], links.disease_url(name=r["disease"])),
                     r["significance"],
                     f"CIViC {r['level']}" if r.get("level") else "",
                     f"EID{r['evidence_id']}"
                     + (f" +{r['n']-1}" if r.get("n", 1) > 1 else ""))
                    for r in ce]))
    L.append(more_line(b.get("civic_association_total", 0), len(ce), "by evidence level"))
    if any(r.get("n", 1) > 1 for r in ce):
        L.append("\n*CIViC column: a representative evidence item (EID); "
                 "+N = additional CIViC items supporting the same association.*")
    return "\n".join(L)


def r_clinical_trials(b):
    L = ["## Clinical trials", "",
         f"**Total trials: {_i(b.get('trial_count'))}.**"]
    pc = b.get("phase_counts") or {}
    if pc:
        # The phase/status distribution + top list are over a sampled subset when
        # the drug has more trials than we page through (high-trial drugs).
        sampled = b.get("sampled_trials") or 0
        total = b.get("trial_count") or 0
        note = (f" (phase/status distribution below is over {_i(sampled)} sampled "
                f"trials of the {_i(total)} total)" if total and sampled and total > sampled
                else "")
        L += ["", "### Phase distribution {#trial-phases}", "",
              *([note.strip() + ".", ""] if note else []),
              table(["Phase", "Trials"],
                    [(k, _i(v)) for k, v in sorted(pc.items(), key=lambda kv: -kv[1])])]
    tt = b.get("top_trials") or []
    if tt:
        L += ["", "### Top trials by phase / activity {#top-trials}", "",
              table(["NCT", "Phase", "Status", "Title"],
                    [((t.get("id") or ""),
                      phase_label(t.get("phase")), t.get("status"), (t.get("title") or "").strip())
                     for t in tt])]
    return "\n".join(L)


def r_pharmacology(b):
    body = []
    roles = b.get("chebi_roles") or []
    if roles:
        body.append(f"**ChEBI roles:** {', '.join(roles)}.")
    atc = b.get("atc_codes") or []
    if atc:
        codes = ", ".join(atc)
        body.append(f"**ATC classification:** {codes} "
                    f"(WHO ATC level names are licensing-restricted).")
    if not body:
        body.append("*No ChEBI role or ATC classification available.*")
    return "## Pharmacology\n\n" + "\n\n".join(body)


def r_patent_literature(b):
    L = ["## Patent literature", ""]
    total = b.get("patent_total") or 0
    if total:
        L.append(f"**Patent mentions (SureChEMBL): {_i(total)}** across "
                 f"{len(b.get('patent_compound_ids') or [])} patent_compound record(s). "
                 f"Counts attach to the compound, so promiscuous molecules score high.")
    else:
        L.append("*No SureChEMBL patent_compound records "
                 "(expected for biologics).*")
    return "\n".join(L)


def r_salt_forms(b):
    parent, childs = b.get("parent_chembl"), b.get("child_chembls") or []
    if not parent and not childs:
        return ("## Salt forms & parent\n\n*Single canonical ChEMBL molecule "
                "(no salt/anhydrous forms linked).*")
    L = ["## Salt forms & parent", ""]
    if parent:
        L.append(f"This molecule (`{b.get('chembl_id')}`) is a salt/anhydrous form "
                 f"of parent `{parent}`.")
    if childs:
        L.append(f"Parent form; salt/anhydrous children: "
                 + ", ".join(f"`{c}`" for c in childs) + ".")
    return "\n".join(L)


# §6 (ChEBI roles) and §11 (patents) are folded into §1; §12 (salt forms)
# duplicates §1's salt/parent lines — all three are dropped as standalone
# sections to avoid 1-line orphans (see render_all). Their collectors still run
# so the data is available to fold in.
RENDER = {
    "1": r_drug_ids,
    "2": r_targets,
    "3": r_bioactivity,
    "4": r_indications,
    "5": r_clinical_trials,
    "7": r_related_molecules,
    "8": r_target_pathways,
    "9": r_pharmacogenomics,
    "10": r_clinical_evidence,
}


def render_all(bundles):
    """Drug page body in the FROZEN canonical H2 order (docs/PAGE_CONTRACT.md):
    Identifiers → Targets → Indications & clinical → Pharmacology → Related
    molecules. (Summary wrapped by assemble_page; Related appended after.)"""
    from atlas.render_common import demote, emit_canonical, with_heading_id
    # Fold §6 ChEBI roles + §11 patent metrics into the §1 bundle so the
    # "Identifiers" section carries them (§12 salt forms already covered by §1's
    # own parent/child lines).
    b1 = dict(bundles.get("1") or {})
    b6 = bundles.get("6") or {}
    b11 = bundles.get("11") or {}
    b1["chebi_roles"] = b6.get("chebi_roles")
    b1["patent_total"] = b11.get("patent_total")
    b1["patent_family_total"] = b11.get("patent_family_total")
    b1["patent_compound_ids"] = b11.get("patent_compound_ids")
    b1["patent_compound_breakdown"] = b11.get("patent_compound_breakdown")
    merged = dict(bundles)
    merged["1"] = b1

    def S(s, anchor):
        return with_heading_id(demote(RENDER[s](merged[s])), anchor) if s in merged else ""

    def join(*parts):
        return "\n\n".join(p for p in parts if p and p.strip())

    spec = [
        ("Identifiers", "identifiers", S("1", "drug-ids"), None),
        ("Targets", "targets",
         join(S("2", "primary-targets"), S("3", "bioactivity"), S("8", "target-pathways")),
         "No curated protein targets or measured bioactivity."),
        ("Indications & clinical", "indications",
         join(S("4", "indication-list"), S("5", "clinical-trials"), S("10", "civic")),
         "No labelled indications, trials, or CIViC evidence."),
        ("Pharmacology", "pharmacology", S("9", "pharmacogenomics"),
         "No pharmacogenomic data."),
        ("Related molecules", "related-molecules", S("7", "related-mol"),
         "No competitor molecules sharing a primary target."),
    ]
    return emit_canonical(spec)
