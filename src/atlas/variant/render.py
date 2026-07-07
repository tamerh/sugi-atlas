"""Deterministic markdown renderer for a variant record — NO model. Mirrors the
gene/drug/disease renderers: every fact verbatim from the collected record.
"""
from atlas.page import links
from atlas.render_common import table

# ClinVar review status → gold-star tier (the standard 0-4 confidence scale).
_STARS = {
    "practice guideline": 4,
    "reviewed by expert panel": 3,
    "criteria provided, multiple submitters, no conflicts": 2,
    "criteria provided, single submitter": 1,
    "criteria provided, conflicting classifications": 1,
    "no assertion criteria provided": 0,
    "no classification provided": 0,
}


def _stars(review_status):
    n = _STARS.get((review_status or "").strip().lower(), 0)
    return ("★" * n + "☆" * (4 - n)) + f" ({n}/4)"


def _label(v):
    """Human page label: 'ACTA1 p.Gln248Lys (c.742C>A)'."""
    g = v.get("gene_symbol") or ""
    p, c = v.get("hgvs_p"), v.get("hgvs_c")
    core = f"{g} {p}" if p else f"{g} {c}" if c else g
    return core + (f" ({c})" if p and c else "")


def declarative(v):
    """Answer-first lead — the self-contained factual block an assistant lifts to
    answer 'is GENE pChange pathogenic?'. Includes the in-silico + population
    differentiators when present."""
    label = _label(v)
    cls = v.get("classification") or "classified"
    gene = v.get("gene_symbol")
    rs = v.get("review_status") or ""
    lead = f"**{label}** is classified **{cls}** in {gene} (ClinVar, {_stars(rs)}"
    cons = v.get("consensus")
    if cons:
        lead += f", {cons['n']} submitter" + ("s" if cons["n"] != 1 else "")
    if v.get("last_evaluated"):
        lead += f"; last evaluated {v['last_evaluated']}"
    lead += ")."
    tail = []
    am = v.get("alphamissense")
    if am and am.get("class"):
        tail.append(f"AlphaMissense: {am['class'].replace('_', ' ')} ({am['score']})")
    g = v.get("gnomad")
    if g:
        tail.append(g["band"])
    cond = (v.get("conditions") or [{}])[0].get("name")
    if cond:
        tail.append(f"associated with {cond}")
    return lead + (" " + "; ".join(s[0].upper() + s[1:] for s in tail) + "." if tail else "")


def declarative_plain(v):
    """The lead sentence, markdown-stripped — for the frontmatter/meta
    description (SEO + AI-citation snippet)."""
    import re
    s = declarative(v)
    s = re.sub(r"\[([^\]]+)\]\([^)]*\)", r"\1", s)   # link → label
    s = re.sub(r"[*`_]", "", s)
    return re.sub(r"\s+", " ", s).strip()


def render_body(v):
    L = ["## Summary", "", declarative(v), ""]
    # At a glance
    L.append("**At a glance:** "
             + " · ".join(filter(None, [
                 v.get("classification"),
                 f"review {_stars(v.get('review_status'))}",
                 (f"rsID {v['rsid']}" if v.get("rsid") else None),
                 (f"{v['submitter_count']} submitter"
                  + ("s" if v.get("submitter_count") != 1 else "")
                  if v.get("submitter_count") else None),
                 (f"{len(v['conditions'])} linked condition"
                  + ("s" if len(v['conditions']) != 1 else "")
                  if v.get("conditions") else None),
                 (f"AlphaMissense {v['alphamissense']['class'].replace('_', ' ')}"
                  if v.get("alphamissense") and v["alphamissense"].get("class") else None),
                 (v["gnomad"]["band"] if v.get("gnomad") else None),
             ])))

    # Identity
    L += ["", "## Identity {#identity}", "",
          table(["Field", "Value"], [
              ("Gene", links.maybe_link(v.get("gene_symbol"),
                                        links.gene_url(symbol=v.get("gene_symbol"), hgnc_id=v.get("hgnc_id")))),
              ("Protein change (HGVS p.)", v.get("hgvs_p")),
              ("Coding change (HGVS c.)", v.get("hgvs_c")),
              ("dbSNP", (f"[{v['rsid']}](https://www.ncbi.nlm.nih.gov/snp/{v['rsid']}/)"
                         if v.get("rsid") else None)),
              ("Variant type", v.get("variant_type")),
              ("Location", (f"chr{v['chromosome']}:{v['start']}-{v['stop']} ({v['assembly']})"
                            if v.get("chromosome") else None)),
              ("ClinVar", f"[VCV{v['variation_id']}](https://www.ncbi.nlm.nih.gov/clinvar/variation/{v['variation_id']}/)"),
          ])]
    exprs = v.get("hgvs_expressions") or []
    if exprs:
        L.append("\n**All HGVS expressions:** " + ", ".join(f"`{e}`" for e in exprs))

    # Computational & population evidence — the cross-source concordance readout
    conc = v.get("concordance") or {}
    if conc.get("verdict"):
        L += ["", "## Computational & population evidence {#evidence}", "",
              f"**Concordance:** {conc['verdict']}.", ""]
        L += [f"- {ln}" for ln in conc.get("lines", [])]
        L.append("\n*Independent of the ClinVar clinical classification: AlphaMissense "
                 "(Cheng et al. 2023, in-silico missense pathogenicity) and gnomAD v4 "
                 "population frequency. A prediction, not a clinical determination.*")

    # Clinical significance + consensus + per-submitter table
    cons = v.get("consensus")
    L += ["", "## Clinical significance {#significance}", "",
          f"**{v.get('classification')}** — review status {_stars(v.get('review_status'))} "
          f"*({v.get('review_status')})*"
          + (f", last evaluated {v['last_evaluated']}." if v.get("last_evaluated") else ".")]
    if cons:
        L.append(f"\n**Submitter consensus:** {cons['verdict']}.")
    vcep = v.get("vcep") or []
    if vcep:
        c = vcep[0]
        L.append(f"\n**ClinGen expert panel:** {c.get('assertion')} — {c.get('panel')} "
                 f"({c.get('disease')}). *Expert-panel curated (ACMG); the highest "
                 "ClinVar review tier.*")
    subs = v.get("submissions") or []
    if subs:
        L += ["", "### Submitter classifications {#submitters}", "",
              table(["Submitter", "Classification", "Review", "Method", "Date"],
                    [(s.get("submitter"), s.get("classification"),
                      s.get("review_status"), s.get("method"), s.get("date")) for s in subs])]

    # Conditions
    conds = v.get("conditions") or []
    if conds:
        L += ["", "## Associated conditions {#conditions}", "",
              ", ".join(links.maybe_link(c.get("name") or c.get("mondo_id"),
                                         links.disease_url(mondo_id=c.get("mondo_id"), name=c.get("name")))
                        for c in conds) + "."]
    extra = [p for p in (v.get("phenotype_list") or []) if p]
    if extra:
        L.append("\n**Reported phenotype names (ClinVar):** " + ", ".join(extra[:10]) + ".")

    # Same-residue hotspot context (deterministic, from the gene's own P/LP set)
    hs = v.get("hotspot")
    if hs and hs.get("others"):
        others = hs["others"]
        L += ["", "## Same-residue context {#hotspot}", "",
              f"Residue **{hs['position']}** carries **{len(others)} other pathogenic "
              f"ClinVar variant" + ("s" if len(others) != 1 else "")
              + "** — a recurrently-mutated position: "
              + ", ".join(f"[{o['label']}](/atlas/variant/{o['slug']}/)" for o in others[:12])
              + ". *Positional co-occurrence of independent ClinVar records, not "
              "functional proof.*"]

    L += ["", f"*Source: NCBI ClinVar (variation {v['variation_id']}), plus AlphaMissense "
          "and gnomAD as noted. Classifications reflect these databases as of the page's "
          "build date and may change. Research/reference use — not medical advice; "
          "consult the primary submitters and a clinician.*"]
    return "\n".join(L)
