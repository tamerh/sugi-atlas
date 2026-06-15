"""§14 — non-coding RNA layer. The curated function / disease / interaction /
drug associations the protein-centric sections miss — the real content layer for
otherwise-thin lncRNA / miRNA pages:

  - Rfam family + Rfam-derived GO (RNAcentral) — structured-RNA function
  - ncrna_disease (LncRNADisease v3.0 + HMDD)  — ncRNA -> disease
  - ncrna_interaction (NPInter v5)             — ncRNA -> partner (RNA/protein)
  - ncrna_drug                                  — ncRNA -> drug (resistance/target)

Routed via >>hgnc>>{ncrna_disease,ncrna_interaction,ncrna_drug} and
>>hgnc>>rnacentral (+ >>rnacentral>>go). These are curated, SYMBOL-keyed
associations (not positional genomic overlap), so they SURVIVE the non-coding
scrub (pipeline._scrub_noncoding) — unlike the positional variant/disease blocks.
Empty for most protein-coding genes; the renderers elide."""
from atlas.biobtree import map_all, xref_counts
from atlas.gene.sections.base import Section

# Page caps on the fetch (×100 rows). Interaction sets are huge (MALAT1 ~2,968),
# so bound the fetch; the true totals come from the hgnc xref counts (1 call,
# already paid), and the render shows a capped sample of the fetched rows.
_DISEASE_CAP = 3
_INTERACTION_CAP = 2
_DRUG_CAP = 2


def _dedup(rows, keyfn):
    seen, out = set(), []
    for r in rows:
        k = keyfn(r)
        if not k or k in seen:
            continue
        seen.add(k)
        out.append(r)
    return out


def _locus_group(hgnc_entry):
    return (((hgnc_entry or {}).get("Attributes") or {}).get("Hgnc") or {}).get("locus_group") or ""


def collect(a):
    bundle = {"section": "14_ncrna", "symbol": a.symbol}
    hid = a.hgnc_id
    # The ncRNA layer is only "about this gene" when the gene IS a non-coding RNA.
    # A protein-coding gene only appears in these datasets as a partner/target of
    # some ncRNA — orientation-inverted and out of place on its page — so skip it
    # entirely (also skips the bulk of the corpus, a real collection-cost win).
    if not hid or _locus_group(a.hgnc_entry) == "protein-coding gene":
        return bundle
    totals = xref_counts(a.hgnc_entry or {})   # exact ncrna_* counts — already fetched
    self_sym = (a.symbol or "").strip().upper()

    # Rfam family + Rfam-derived GO (RNAcentral). Distinct Rfam across the gene's
    # URS records; lncRNAs typically carry none (the projection drops the empty
    # rfam columns — read defensively). GO is the Rfam-derived functional layer.
    rfam = []
    for t in map_all(hid, ">>hgnc>>rnacentral"):
        rid = (t.get("rfam_id") or "").strip()
        if rid:
            rfam.append({"rfam_id": rid, "rfam_description": t.get("rfam_description"),
                         "rna_type": t.get("rna_type")})
    bundle["rfam"] = _dedup(rfam, lambda r: r["rfam_id"])
    bundle["go"] = _dedup(
        [{"id": t.get("id"), "type": t.get("type"), "name": t.get("name")}
         for t in map_all(hid, ">>hgnc>>rnacentral>>go") if t.get("id")],
        lambda r: r["id"])

    # ncRNA -> disease (LncRNADisease v3.0 + HMDD). Map projection carries the
    # display fields directly (no per-record entry fetch — the #45/#47 lesson).
    dis = _dedup(map_all(hid, ">>hgnc>>ncrna_disease", cap=_DISEASE_CAP),
                 lambda r: (r.get("disease_name"), r.get("causality")))
    bundle["diseases"] = [{"disease_name": r.get("disease_name"),
                           "causality": r.get("causality"),
                           "validated_method": r.get("validated_method"),
                           "category": r.get("ncrna_category")} for r in dis]
    bundle["disease_total"] = totals.get("ncrna_disease", len(dis))

    # ncRNA -> partner (NPInter): RNA-RNA / RNA-protein interactions. A record may
    # list this gene on either side, so display the OTHER side as the partner
    # (with its type) — never this gene interacting with itself.
    def _orient(r):
        if (r.get("ncrna_name") or "").strip().upper() == self_sym:
            return r.get("partner_name"), r.get("partner_type")
        return r.get("ncrna_name"), r.get("ncrna_type")   # this gene is the partner
    inter_rows = []
    for r in map_all(hid, ">>hgnc>>ncrna_interaction", cap=_INTERACTION_CAP):
        pname, ptype = _orient(r)
        inter_rows.append({"partner_name": pname, "partner_type": ptype,
                           "level": r.get("level"), "datasource": r.get("datasource")})
    inter = _dedup(inter_rows, lambda r: (r.get("partner_name"), r.get("partner_type")))
    bundle["interactions"] = inter
    bundle["interaction_total"] = totals.get("ncrna_interaction", len(inter))

    # ncRNA -> drug (resistance / target).
    drg = _dedup(map_all(hid, ">>hgnc>>ncrna_drug", cap=_DRUG_CAP),
                 lambda r: (r.get("drug_name"), r.get("relation"), r.get("condition")))
    bundle["drugs"] = [{"drug_name": r.get("drug_name"), "relation": r.get("relation"),
                        "effect": r.get("effect"), "condition": r.get("condition")}
                       for r in drg]
    bundle["drug_total"] = totals.get("ncrna_drug", len(drg))
    return bundle


SECTION = Section(
    id="14", name="ncrna",
    description=("Non-coding RNA layer: Rfam family + GO (RNAcentral), ncRNA->disease "
                 "(LncRNADisease/HMDD), ncRNA->partner (NPInter), ncRNA->drug"),
    needs=("symbol", "hgnc_id", "hgnc_entry"),
    produces=("rfam", "go", "diseases", "disease_total", "interactions",
              "interaction_total", "drugs", "drug_total"),
    datasets=("rnacentral", "go", "ncrna_disease", "ncrna_interaction", "ncrna_drug"),
    chains=(">>hgnc>>rnacentral", ">>hgnc>>rnacentral>>go", ">>hgnc>>ncrna_disease",
            ">>hgnc>>ncrna_interaction", ">>hgnc>>ncrna_drug"),
    collect_fn=collect,
)
