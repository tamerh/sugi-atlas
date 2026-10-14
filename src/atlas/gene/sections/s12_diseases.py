"""§12 — diseases: Mendelian (OMIM/GenCC/MONDO/Orphanet), phenotypes (HPO),
complex disease (GWAS). Multi-hop disease routes are the showcase of this section."""
from atlas.biobtree import map_all, xref_counts, bbmap
from atlas.gene.sections.base import Section

# Disease links via ClinVar only from this gene's own (likely) pathogenic records —
# not VUS/benign/conflicting ones, and not copy-number records whose multi-megabase
# span merely covers the gene (C1orf159 → OX40 deficiency via a 5.6-Mb deletion;
# AKAP8 → long QT from VUS). 2,252 gene pages listed such diseases.
_CV_PATHOGENIC = ('[germline_classification.contains("athogenic") && '
                  '!germline_classification.contains("onflicting") && '
                  '!type.contains("copy number")]')
# A disease reached through ClinVar needs ≥2 of this gene's own (likely)
# pathogenic, non-region records: one record's condition list can carry unrelated
# conditions (CAPN3 → type 2 diabetes, migraine), and large deletions typed
# "Deletion" (AK1 → HHT) slip past the copy-number filter.
_CV_MIN_RECORDS = 2
_CV_RECORD_PAGES = 10             # ≤ ~1,000 P/LP records sampled per gene


def _clinvar_mondo(hgnc_id):
    """{MONDO id: name} supported by ≥_CV_MIN_RECORDS of the gene's own
    pathogenic, non-region ClinVar records (per-record clinvar>>mondo, batched)."""
    from atlas.gene.sections.s06_variants import _is_region_cnv
    recs = [r for r in map_all(hgnc_id, ">>hgnc>>clinvar" + _CV_PATHOGENIC, cap=_CV_RECORD_PAGES)
            if r.get("id") and not _is_region_cnv(r.get("name"))]
    counts, names = {}, {}
    ids = [r["id"] for r in recs]
    for i in range(0, len(ids), 100):
        page = None
        for _ in range(20):
            resp = bbmap(",".join(ids[i:i + 100]), ">>clinvar>>mondo", page)
            for m in resp.get("mappings") or []:
                for t in set(m.get("targets") or []):
                    mid, _, nm = t.partition("|")
                    counts[mid] = counts.get(mid, 0) + 1
                    names.setdefault(mid, nm)
            pg = resp.get("pagination") or {}
            page = pg.get("next_token")
            if not pg.get("has_next") or not page:
                break
    return {m: names[m] for m, n in counts.items() if n >= _CV_MIN_RECORDS}


def _via_mondo(mondo_ids, tail):
    """Targets of `>>mondo>>…` for a set of MONDO ids (batched)."""
    ids, out = sorted(mondo_ids), []
    for i in range(0, len(ids), 50):
        out += map_all(",".join(ids[i:i + 50]), ">>mondo" + tail)
    return out


CHAINS = (
    ">>hgnc>>mim",                                # gene MIM
    ">>hgnc>>clinvar" + _CV_PATHOGENIC + ">>mondo>>mim",                # disease-phenotype MIMs (3-hop)
    ">>hgnc>>gencc",                              # curated gene-disease
    ">>hgnc>>clinvar" + _CV_PATHOGENIC + ">>mondo", ">>hgnc>>gencc>>mondo",
    ">>hgnc>>clinvar" + _CV_PATHOGENIC + ">>orphanet", ">>hgnc>>clinvar" + _CV_PATHOGENIC + ">>mondo>>orphanet",
    ">>hgnc>>hpo", ">>hgnc>>clinvar" + _CV_PATHOGENIC + ">>mondo>>hpo",
    ">>hgnc>>gwas",
    ">>hgnc>>gwas>>gwas_study", ">>hgnc>>clinvar" + _CV_PATHOGENIC + ">>mondo>>gwas_study",
    ">>hgnc>>gwas>>efo",                          # EFO-canonical GWAS trait names
    ">>hgnc>>clinvar" + _CV_PATHOGENIC + ">>mondo>>mesh", ">>hgnc>>gencc>>mondo>>mesh",   # MeSH disease categories
    ">>hgnc>>intogen",                            # cancer-driver classification (LoF/Act)
    ">>hgnc>>civic",                              # CIViC gene-level curated narrative
    ">>hgnc>>clingen_gene_validity",              # ClinGen expert-panel gene-disease validity
    ">>hgnc>>panelapp_gene",                      # Genomics England clinical gene panels
)
DATASETS = ("mim", "gencc", "mondo", "orphanet", "hpo", "gwas", "gwas_study", "efo",
            "mesh", "intogen", "civic", "clinvar", "hgnc", "clingen_gene_validity",
            "panelapp_gene")

# PanelApp diagnostic-grade rating order (green = high-confidence diagnostic,
# amber = moderate, red = low — usually filtered out upstream).
_PANEL_RANK = {"green": 0, "amber": 1, "red": 2}

def _pval(p):
    try:
        v = float(p)
        return v if v == v else 1.0
    except (TypeError, ValueError):
        return 1.0


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
    cv_mondo = _clinvar_mondo(a.hgnc_id)
    bundle["disease_omim"] = [f"MIM:{t['id']}" for t in _via_mondo(cv_mondo, ">>mim")]
    bundle["gencc"] = [{"disease": t.get("disease_title"),
                        "classification": t.get("classification_title"),
                        "inheritance": t.get("moi_title")}
                       for t in map_all(a.hgnc_id, ">>hgnc>>gencc")]

    # MONDO: union curated (gencc) + clinical (clinvar) routes.
    mondo = dict(cv_mondo)
    for t in map_all(a.hgnc_id, ">>hgnc>>gencc>>mondo"):
        mondo[t["id"]] = t.get("name") or mondo.get(t["id"])
    bundle["mondo"] = [{"id": k, "name": v} for k, v in mondo.items()]

    # Orphanet: union direct (clinvar) + via-mondo (carries names).
    orph = {}
    for t in _via_mondo(cv_mondo, ">>orphanet"):
        orph[t["id"]] = t.get("name") or orph.get(t["id"])
    bundle["orphanet"] = [{"id": f"Orphanet:{k}", "name": v} for k, v in orph.items()]

    # HPO: gene-level (hgnc>>hpo) + disease-level (via mondo).
    hpo = {}
    for t in map_all(a.hgnc_id, ">>hgnc>>hpo") + _via_mondo(cv_mondo, ">>hpo"):
        hpo[t["id"]] = t.get("name") or hpo.get(t["id"])
    bundle["hpo"] = [{"id": k, "name": v} for k, v in hpo.items()]
    bundle["hpo_total"] = xc.get("hpo", len(hpo))

    gwas = map_all(a.hgnc_id, ">>hgnc>>gwas")
    # Strongest association first — the render shows the top 30, and ID order hid
    # the lead hits (97% of gene GWAS tables were unsorted).
    gwas.sort(key=lambda t: _pval(t.get("p_value")))
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
    for t in map_all(a.hgnc_id, ">>hgnc>>gwas>>gwas_study") + _via_mondo(cv_mondo, ">>gwas_study"):
        gws[t["id"]] = 1
    bundle["gwas_studies"] = list(gws)

    # MeSH disease descriptors — union of clinvar and gencc disease routes.
    # `descriptor_class=1` is a Main MeSH heading (D-code); supplementary
    # concepts (C-code, `is_supplementary=true`) are a smaller fallback set.
    # Tree numbers (e.g. C04.700.600) classify the descriptor into the MeSH
    # category hierarchy — useful for rollup to top-level disease groupings.
    mesh = {}
    for t in _via_mondo(cv_mondo, ">>mesh") + map_all(a.hgnc_id, ">>hgnc>>gencc>>mondo>>mesh"):
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
    } for r in map_all(a.hgnc_id, ">>hgnc>>clingen_gene_validity")
        if not (r.get("disease_label") or "").startswith("obsolete ")]  # skip obsolete disease names

    # PanelApp (Genomics England) — clinical diagnostic gene panels this gene sits
    # on, with the gene's rating on each (green = diagnostic-grade). One row per
    # (panel, gene); ranked diagnostic-grade first, then panel name.
    panels = [{"panel": (r.get("panel_name") or "").strip(),
               "confidence": (r.get("confidence") or "").strip().lower()}
              for r in map_all(a.hgnc_id, ">>hgnc>>panelapp_gene")
              if (r.get("panel_name") or "").strip()]
    panels.sort(key=lambda p: (_PANEL_RANK.get(p["confidence"], 3), p["panel"]))
    bundle["panelapp"] = panels
    bundle["panelapp_green"] = sum(1 for p in panels if p["confidence"] == "green")
    return bundle

SECTION = Section(
    id="12", name="diseases",
    description=("Disease associations — OMIM (gene + disease phenotype), GenCC, MONDO, "
                 "Orphanet, HPO phenotypes, GWAS (gene-mapped + disease-associated), "
                 "MeSH descriptors (NLM disease vocabulary with tree-number category paths)"),
    needs=("hgnc_id", "hgnc_entry"),
    produces=("gene_omim", "disease_omim", "gencc", "mondo", "orphanet", "hpo",
              "gwas", "gwas_studies", "efo_traits", "mesh_descriptors", "intogen", "civic",
              "clingen_validity", "panelapp", "panelapp_green"),
    datasets=DATASETS, chains=CHAINS, collect_fn=collect,
)
