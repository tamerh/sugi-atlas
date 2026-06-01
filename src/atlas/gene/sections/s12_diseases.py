"""§12 — diseases: Mendelian (OMIM/GenCC/MONDO/Orphanet), phenotypes (HPO),
complex disease (GWAS). Multi-hop disease routes are the showcase of this section."""
from atlas.biobtree import map_all, xref_counts
from atlas.gene.sections.base import Section

CHAINS = (
    ">>hgnc>>mim",                                # gene MIM
    ">>hgnc>>clinvar>>mondo>>mim",                # disease-phenotype MIMs (3-hop)
    ">>hgnc>>gencc",                              # curated gene-disease
    ">>hgnc>>clinvar>>mondo", ">>hgnc>>gencc>>mondo",
    ">>hgnc>>clinvar>>orphanet", ">>hgnc>>clinvar>>mondo>>orphanet",
    ">>hgnc>>hpo", ">>hgnc>>clinvar>>mondo>>hpo",
    ">>hgnc>>gwas",
    ">>hgnc>>gwas>>gwas_study", ">>hgnc>>clinvar>>mondo>>gwas_study",
    ">>hgnc>>gwas>>efo",                          # EFO-canonical GWAS trait names
    ">>hgnc>>clinvar>>mondo>>mesh", ">>hgnc>>gencc>>mondo>>mesh",   # MeSH disease categories
    ">>hgnc>>intogen",                            # cancer-driver classification (LoF/Act)
    ">>hgnc>>civic",                              # CIViC gene-level curated narrative
    ">>hgnc>>clingen_gene_validity",              # ClinGen expert-panel gene-disease validity
)
DATASETS = ("mim", "gencc", "mondo", "orphanet", "hpo", "gwas", "gwas_study", "efo",
            "mesh", "intogen", "civic", "clinvar", "hgnc", "clingen_gene_validity")

def collect(a):
    bundle = {"section": "12_diseases", "symbol": a.symbol}
    xc = xref_counts(a.hgnc_entry)

    # Cancer-significance signals — one row per gene at most. Genes outside
    # CIViC / intOGen coverage (e.g. structural proteins like TTN) return
    # nothing; the renderer suppresses the block in that case.
    intogen_rows = map_all(a.hgnc_id, ">>hgnc>>intogen")
    if intogen_rows:
        r = intogen_rows[0]
        bundle["intogen"] = {
            "symbol": r.get("symbol"),
            "role": r.get("role"),   # 'Act' = oncogene-like, 'LoF' = tumor-suppressor-like
            "cancer_types": r.get("cancer_types"),   # comma-separated abbreviations
        }
    else:
        bundle["intogen"] = None
    civic_rows = map_all(a.hgnc_id, ">>hgnc>>civic")
    if civic_rows:
        r = civic_rows[0]
        bundle["civic"] = {
            "id": r.get("id"),
            "name": r.get("name"),
            "description": r.get("description"),   # curated paragraph
        }
    else:
        bundle["civic"] = None

    bundle["gene_omim"] = [f"MIM:{t['id']}" for t in map_all(a.hgnc_id, ">>hgnc>>mim")]
    # Disease-phenotype OMIM IDs come through MONDO (gene>>mim is only the gene).
    bundle["disease_omim"] = [f"MIM:{t['id']}" for t in
                              map_all(a.hgnc_id, ">>hgnc>>clinvar>>mondo>>mim")]
    bundle["gencc"] = [{"disease": t.get("disease_title"),
                        "classification": t.get("classification_title"),
                        "inheritance": t.get("moi_title")}
                       for t in map_all(a.hgnc_id, ">>hgnc>>gencc")]

    # MONDO: union curated (gencc) + clinical (clinvar) routes.
    mondo = {}
    for ch in (">>hgnc>>clinvar>>mondo", ">>hgnc>>gencc>>mondo"):
        for t in map_all(a.hgnc_id, ch):
            mondo[t["id"]] = t.get("name") or mondo.get(t["id"])
    bundle["mondo"] = [{"id": k, "name": v} for k, v in mondo.items()]

    # Orphanet: union direct (clinvar) + via-mondo (carries names).
    orph = {}
    for ch in (">>hgnc>>clinvar>>orphanet", ">>hgnc>>clinvar>>mondo>>orphanet"):
        for t in map_all(a.hgnc_id, ch):
            orph[t["id"]] = t.get("name") or orph.get(t["id"])
    bundle["orphanet"] = [{"id": f"Orphanet:{k}", "name": v} for k, v in orph.items()]

    # HPO: gene-level (hgnc>>hpo) + disease-level (via mondo).
    hpo = {}
    for ch in (">>hgnc>>hpo", ">>hgnc>>clinvar>>mondo>>hpo"):
        for t in map_all(a.hgnc_id, ch):
            hpo[t["id"]] = t.get("name") or hpo.get(t["id"])
    bundle["hpo"] = [{"id": k, "name": v} for k, v in hpo.items()]
    bundle["hpo_total"] = xc.get("hpo", len(hpo))

    gwas = map_all(a.hgnc_id, ">>hgnc>>gwas")
    bundle["gwas"] = [{"id": t["id"], "trait": t.get("disease_trait"),
                       "p_value": t.get("p_value")} for t in gwas]
    bundle["gwas_total"] = xc.get("gwas", len(gwas))

    # EFO trait names — canonical vocabulary for the GWAS traits above.
    # Same disease can appear under multiple GWAS Catalog free-text names
    # (e.g. "Type 2 diabetes" / "Type 2 diabetes mellitus" / "Diabetes,
    # type 2"); EFO normalizes them. AI agents prefer EFO's IDs over
    # GWAS-catalog free text for cross-resource joins.
    efo_seen = {}
    for t in map_all(a.hgnc_id, ">>hgnc>>gwas>>efo"):
        eid = t.get("id")
        if eid and eid not in efo_seen:
            efo_seen[eid] = t.get("name") or eid
    bundle["efo_traits"] = [{"id": eid, "name": name}
                            for eid, name in efo_seen.items()]

    # gwas_study ids: gene-mapped (gwas>>gwas_study) + disease-associated
    # (mondo>>gwas_study). The disease route catches genes whose GWAS aren't
    # gene-mapped directly (e.g. AR: hgnc>>gwas empty) but reach studies via disease.
    gws = {}
    for ch in (">>hgnc>>gwas>>gwas_study", ">>hgnc>>clinvar>>mondo>>gwas_study"):
        for t in map_all(a.hgnc_id, ch):
            gws[t["id"]] = 1
    bundle["gwas_studies"] = list(gws)

    # MeSH disease descriptors — union of clinvar and gencc disease routes.
    # `descriptor_class=1` is a Main MeSH heading (D-code); supplementary
    # concepts (C-code, `is_supplementary=true`) are a smaller fallback set.
    # Tree numbers (e.g. C04.700.600) classify the descriptor into the MeSH
    # category hierarchy — useful for rollup to top-level disease groupings.
    mesh = {}
    for ch in (">>hgnc>>clinvar>>mondo>>mesh", ">>hgnc>>gencc>>mondo>>mesh"):
        for t in map_all(a.hgnc_id, ch):
            mesh[t["id"]] = t  # later route wins on dupes (gencc carries class info too)
    def _msort(t):
        # Main descriptors first (descriptor_class=1), then alphabetical by name.
        return (0 if t.get("descriptor_class") == "1" else 1,
                t.get("descriptor_name") or "")
    rows = sorted(mesh.values(), key=_msort)
    bundle["mesh_descriptors"] = [{
        "id": r["id"],
        "name": r.get("descriptor_name"),
        "tree_numbers": [tn for tn in (r.get("tree_numbers") or "").split(";") if tn],
        "is_supplementary": r.get("is_supplementary") == "true",
    } for r in rows]

    # ClinGen gene-disease validity — expert-panel curated relationship strength.
    # Multiple rows per gene (one per disease). `classification` is the headline:
    # Definitive > Strong > Moderate > Limited > Disputed > Refuted > No Known.
    bundle["clingen_validity"] = [{
        "id": r.get("id"),
        "disease": r.get("disease_label"),
        "classification": r.get("classification"),
        "moi": r.get("moi"),  # AD / AR / XL / etc.
    } for r in map_all(a.hgnc_id, ">>hgnc>>clingen_gene_validity")]
    return bundle

SECTION = Section(
    id="12", name="diseases",
    description=("Disease associations — OMIM (gene + disease phenotype), GenCC, MONDO, "
                 "Orphanet, HPO phenotypes, GWAS (gene-mapped + disease-associated), "
                 "MeSH descriptors (NLM disease vocabulary with tree-number category paths)"),
    needs=("hgnc_id", "hgnc_entry"),
    produces=("gene_omim", "disease_omim", "gencc", "mondo", "orphanet", "hpo",
              "gwas", "gwas_studies", "efo_traits", "mesh_descriptors", "intogen", "civic",
              "clingen_validity"),
    datasets=DATASETS, chains=CHAINS, collect_fn=collect,
)
