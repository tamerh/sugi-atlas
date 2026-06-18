"""§5 — genes → proteins: per-cohort-gene id mapping (HGNC, Ensembl, UniProt,
protein name). REUSE wrapper — fans gene §1 (gene_ids) + §3 (protein_ids) over
the cohort and aggregates.

No new chains/datasets — everything comes from gene §1 + §3 bundles via the
cohort fan-out helper."""
from atlas.section import Section
from atlas.biobtree import map_all
from atlas.disease.cohort import fan
from atlas.gene.sections import s01_gene_ids, s03_protein_ids

CHAINS   = (">>hgnc>>ensembl>>uniprot",         # reused via gene collectors
            ">>mondo>>doid>>alliance_disease")  # Alliance/DO gene→disease associations
DATASETS = ("hgnc", "ensembl", "uniprot", "doid", "alliance_disease")

_EVIDENCE_KEYS = ("gwas", "gencc", "clinvar", "civic_evidence", "definitional")


def _classify(ev: dict) -> str:
    """Partition a gene's evidence flags into a single bucket name."""
    g = bool(ev.get("gwas"))
    gc = bool(ev.get("gencc"))
    cv = bool(ev.get("clinvar"))
    ci = bool(ev.get("civic_evidence"))
    # `definitional` genes are a last-resort fallback that ONLY exists when no
    # structured route resolved (see anchors._definitional_genes) — so such a
    # gene never carries another flag, and gets its own bucket.
    if ev.get("definitional") and not (g or gc or cv or ci):
        return "definitional"
    n = sum((g, gc, cv, ci))
    if n >= 3:
        return "multi_evidence"
    if g and gc and n == 2:
        return "gwas_and_gencc"
    if g and cv and n == 2:
        return "gwas_and_clinvar"
    if g and not (gc or cv or ci):
        return "gwas_only"
    if ci and not (g or gc or cv):
        return "civic_only"
    # Two-evidence combos not listed above (e.g. gencc+clinvar, gwas+civic)
    # collapse into multi_evidence; single-source non-GWAS/non-CIViC into
    # multi_evidence too so the partition stays exhaustive.
    if n >= 2:
        return "multi_evidence"
    return "multi_evidence"


def collect(a):
    g1_bundles = fan(s01_gene_ids.SECTION.collect_fn, a.cohort)
    g3_bundles = fan(s03_protein_ids.SECTION.collect_fn, a.cohort)
    g3_by = {b.get("symbol"): b for b in g3_bundles}

    summary = {"gwas_only": 0, "gwas_and_gencc": 0, "gwas_and_clinvar": 0,
               "civic_only": 0, "multi_evidence": 0, "definitional": 0}
    genes = []
    seen_canonical = set()
    from atlas.page.uniprot_cc import first_sentence
    function_lines = []   # cohort_function_summary: per-gene 1-sentence function
    molecular_basis = []  # per-gene UniProt CC: function/family/location/disease
    for b1 in g1_bundles:
        sym = b1.get("symbol")
        b3 = g3_by.get(sym, {})
        ev_raw = a.cohort_evidence.get(b1.get("hgnc_id")) or {}
        ev = {k: bool(ev_raw.get(k)) for k in _EVIDENCE_KEYS}
        canonical = b3.get("canonical_uniprot")
        if canonical:
            seen_canonical.add(canonical)
        # protein_name + function now flow from gene §3 (CC block, post the
        # 2026-05-31 biobtree refresh).
        protein_name = b3.get("protein_name")
        function_text = (b3.get("cc") or {}).get("function") or ""
        genes.append({
            "symbol": sym,
            "hgnc_id": b1.get("hgnc_id"),
            "ensembl_id": b1.get("ensembl_id"),
            "hgnc_name": b1.get("name"),
            "canonical_uniprot": canonical,
            "all_uniprots": b3.get("uniprot_all", []),
            "protein_name": protein_name,
            "function": function_text,
            "evidence": ev,
        })
        if function_text:
            # First sentence per gene — the cohort summary tabulates these
            # so a reader sees "what each disease gene actually does" at a
            # glance, not just HGNC ids.
            function_lines.append({
                "symbol": sym,
                "protein_name": protein_name or sym,
                "function_lead": first_sentence(function_text, max_chars=240),
            })
        # Molecular-basis fields — the rest of the UniProt CC block §3 already
        # fetched (we used to keep only `function`). Powers the disease's
        # "Molecular basis" block for its causal gene(s). `disease` is UniProt's
        # OWN curated gene→disorder statement (reporting, not inference).
        cc = b3.get("cc") or {}
        mb = {
            "symbol": sym, "protein_name": protein_name or sym, "uniprot": canonical,
            "family": first_sentence(cc.get("similarity") or "", max_chars=140),
            "function": first_sentence(function_text, max_chars=300),
            "subcellular": first_sentence(cc.get("subcellular_location") or "", max_chars=120),
            "cofactor": first_sentence(cc.get("cofactor") or "", max_chars=120),
            "disease": (cc.get("disease") or "").strip(),
        }
        if mb["function"] or mb["disease"] or mb["family"]:
            molecular_basis.append(mb)
        summary[_classify(ev)] += 1

    # Alliance / Disease Ontology gene→disease associations (via mondo→doid→
    # alliance_disease). A distinct curated source from the cohort above, with the
    # is_implicated_in (causal-ish) vs is_marker_for (biomarker) distinction nothing
    # else here carries. Human only (the rare model-organism rows are noise);
    # implicated takes precedence over marker for a gene in both.
    implicated, marker = set(), set()
    if a.mondo_id:
        for r in map_all(a.mondo_id, ">>mondo>>doid>>alliance_disease"):
            sym = r.get("gene_symbol")
            if not sym or (r.get("species") or "") != "Homo sapiens":
                continue
            (implicated if r.get("association_type") == "is_implicated_in" else marker).add(sym)
    marker -= implicated

    return {
        "section": "05_genes_proteins",
        "mondo_id": a.mondo_id,
        "genes": genes,
        "gene_count": len(genes),
        "protein_count": len(seen_canonical),
        "evidence_summary": summary,
        "cohort_function_summary": function_lines,
        "molecular_basis": molecular_basis,
        "alliance_implicated": sorted(implicated),
        "alliance_marker": sorted(marker),
    }


SECTION = Section(
    id="5", name="genes_proteins",
    description=("Per-cohort-gene identifier mapping (HGNC → Ensembl → UniProt + "
                 "protein name + evidence tier + Mendelian overlap flag). "
                 "REUSE wrapper over gene §1/§3."),
    needs=("cohort", "cohort_evidence"),
    produces=("genes", "protein_count", "gene_count", "evidence_summary",
              "cohort_function_summary", "molecular_basis",
              "alliance_implicated", "alliance_marker"),
    datasets=DATASETS, chains=CHAINS, collect_fn=collect,
)
