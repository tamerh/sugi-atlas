#!/usr/bin/env python3
"""Deterministic markdown renderer for drug-section bundles — NO model.
Mirrors gene/render.py + disease/render.py: one r_* fn per section + a RENDER
dict. Every fact comes verbatim from the bundle."""
from atlas.render_common import table
from atlas.civic import therapy_label


def _i(n):
    return f"{n:,}" if isinstance(n, int) else (n if n is not None else "")


def r_drug_ids(b):
    L = ["## Drug identity & classification", ""]
    rows = [
        ("ChEMBL ID", b.get("chembl_id")),
        ("Name", b.get("canonical_name")),
        ("Type", b.get("molecule_type")),
        ("Max phase", b.get("max_phase")),
        ("FDA approved", b.get("is_fda_approved")),
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
    # Drug class — ChEBI functional roles (folded in from the former
    # standalone Pharmacology section; ATC already lives in the table above).
    roles = b.get("chebi_roles") or []
    if roles:
        L.append(f"\n**Drug class (ChEBI roles):** {', '.join(roles)}.")
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
    # former standalone Patent literature section). The total is a count of
    # patents *mentioning* the matched structures, not distinct inventions —
    # and is typically dominated by one promiscuous/reference structure, which
    # we quantify here. Richer landscape (assignee, CPC/IPC class, distinct
    # families, jurisdiction/timeline) is gated on biobtree #25/#26/#27.
    pt = b.get("patent_total") or 0
    if pt:
        bd = b.get("patent_compound_breakdown") or []
        n_rec = len(bd) or len(b.get("patent_compound_ids") or [])
        dom = ""
        if bd and bd[0].get("patent_count"):
            top = bd[0]["patent_count"]
            pct = round(100 * top / pt) if pt else 0
            if len(bd) > 1 and pct >= 60:
                dom = (f" One matched structure accounts for {_i(top)} ({pct}%) "
                       f"of the total.")
        L.append(f"\n**Patent mentions (SureChEMBL):** {_i(pt)} across {n_rec} "
                 f"matched compound structure(s).{dom} Counts are patent mentions "
                 f"of the compound (not distinct inventions), so promiscuous / "
                 f"reference molecules inflate the total.")
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
            return f"{pct}%" + (f" ({', '.join(tags)})" if tags else "")
        L.append(table(["Gene", "Target", "Action", "pAffinity", "Cancer dependency", "UniProt"],
                       [(t.get("gene_symbol") or "", t.get("target_name") or "",
                         t.get("action") or "", t.get("affinity") or "",
                         _dep(t), t.get("uniprot") or "") for t in pt]))
    bc = b.get("bioactivity_target_count") or 0
    if bc:
        names = ", ".join(t.get("name") or t.get("chembl_target_id")
                          for t in (b.get("bioactivity_targets") or [])[:10])
        L.append(f"\n**Broader ChEMBL bioactivity targets: {bc}** "
                 f"(assay-derived). Sample: {names}.")
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
             f"pChembl ≥ 5 of {_i(b.get('activity_total'))} total. Top 30 by "
             f"potency (10 = 0.1 nM, 6 = 1 µM):**\n")
    L.append(table(["pChembl", "Type", "Value", "Unit", "Activity ID"],
                   [(r.get("pchembl"), r.get("type"), r.get("value"),
                     r.get("unit"), r.get("id")) for r in ca]))
    return "\n".join(L)


def r_indications(b):
    L = ["## Indications", "",
         f"**{_i(b.get('indication_count'))} indications "
         f"({_i(b.get('approved_count'))} at max phase 4 / approved).**"]
    inds = b.get("indications") or []
    if inds:
        L += ["", table(["Indication", "Max phase", "MONDO", "EFO"],
                        [(i.get("name") or i.get("efo_id") or i.get("mesh_id") or "",
                          i.get("max_phase"), i.get("mondo_id") or "",
                          i.get("efo_id") or "") for i in inds[:40]])]
    return "\n".join(L)


def r_related_molecules(b):
    L = ["## Related molecules", ""]
    rm = b.get("related_molecules") or []
    if not rm:
        L.append("*No phase-≥2 competitor molecules sharing a primary target.*")
        return "\n".join(L)
    L.append(f"**{_i(b.get('competitor_count'))} phase-≥2 molecules share ≥1 "
             f"primary target. Top {len(rm)} by shared-target count:**\n")
    L.append(table(["Molecule", "ChEMBL", "Max phase", "Shared targets"],
                   [(r.get("name") or r.get("id"), r.get("id"), r.get("phase"),
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
    L.append(f"**Aggregated over {len(tg)} target gene(s): {', '.join(tg)}.**")
    if tp:
        L += ["", f"**Top Reactome pathways ({_i(b.get('pathway_count'))} total), "
              f"by targets touching each:**", "",
              table(["Pathway", "Targets", "Genes"],
                    [(p.get("name") or p.get("id") or "",
                      p.get("gene_count"), ", ".join(p.get("genes") or []))
                     for p in tp])]
    if go:
        L += ["", "**Dominant GO biological processes:**", "",
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
        L.append(f"**PharmGKB dosing guidelines ({b.get('guideline_count')}) — CPIC / "
                 f"DPWG genotype-guided dosing for this drug (drug × pharmacogene):**\n")
        L.append(table(["Guideline", "Source", "Gene(s)", "Dosing", "Recommendation"],
                       [((r.get("name") or r.get("id") or "")[:70],
                         r.get("source"), r.get("genes"),
                         "yes" if r.get("has_dosing") else "",
                         "yes" if r.get("has_recommendation") else "") for r in g]))
    elif pa:
        L.append("*No CPIC/DPWG dosing guideline, but PharmGKB curates "
                 "pharmacogenomic annotations for this drug:*")
    # Coverage line — gene-keyed clinical/variant annotations live on the gene
    # pages; surface the counts + a link to the PharmGKB chemical record.
    if pa and (ca or va):
        L.append(f"\nPharmGKB also curates "
                 f"{_i(ca) or 0} clinical and {_i(va) or 0} variant annotation(s) "
                 f"for this drug (gene-keyed; see "
                 f"[PharmGKB](https://www.pharmgkb.org/chemical/{pa})).")
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
    L.append(f"**Variant × indication × effect "
             f"({_i(b.get('civic_association_total'))} predictive associations from "
             f"{_i(b.get('civic_predictive_total'))} curated evidence items"
             + (f"; also {extra}" if extra else "") + "):**\n")
    L.append(table(["Variant", "Indication", "Effect", "Therapy", "Level", "CIViC"],
                   [(r["profile"], r["disease"], r["significance"],
                     therapy_label(r["therapy"]),
                     f"CIViC {r['level']}" if r.get("level") else "",
                     f"EID{r['evidence_id']}"
                     + (f" +{r['n']-1}" if r.get("n", 1) > 1 else ""))
                    for r in ce]))
    more = (b.get("civic_association_total") or 0) - len(ce)
    if more > 0:
        L.append(f"\n*+{more} more predictive associations (showing top {len(ce)} by level).*")
    return "\n".join(L)


def r_clinical_trials(b):
    L = ["## Clinical trials", "",
         f"**Total trials: {_i(b.get('trial_count'))}.**"]
    pc = b.get("phase_counts") or {}
    if pc:
        L += ["", "**Phase distribution:**", "",
              table(["Phase", "Trials"],
                    [(k, _i(v)) for k, v in sorted(pc.items(), key=lambda kv: -kv[1])])]
    tt = b.get("top_trials") or []
    if tt:
        L += ["", "**Top trials by phase / activity:**", "",
              table(["NCT", "Phase", "Status", "Title"],
                    [((t.get("id") or ""),
                      t.get("phase"), t.get("status"), (t.get("title") or "")[:65])
                     for t in tt])]
    return "\n".join(L)


def r_pharmacology(b):
    body = []
    roles = b.get("chebi_roles") or []
    if roles:
        body.append(f"**Drug class (ChEBI roles):** {', '.join(roles)}.")
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
    # Fold §6 ChEBI roles + §11 patent metrics into the §1 bundle so the
    # "Drug identity & classification" section carries them (§12 salt forms is
    # already covered by §1's own parent/child lines, so it's simply dropped).
    b1 = dict(bundles.get("1") or {})
    b6 = bundles.get("6") or {}
    b11 = bundles.get("11") or {}
    b1["chebi_roles"] = b6.get("chebi_roles")
    b1["patent_total"] = b11.get("patent_total")
    b1["patent_compound_ids"] = b11.get("patent_compound_ids")
    b1["patent_compound_breakdown"] = b11.get("patent_compound_breakdown")
    merged = dict(bundles)
    merged["1"] = b1
    return "\n\n".join(RENDER[sid](merged[sid]) for sid in RENDER if sid in merged)
