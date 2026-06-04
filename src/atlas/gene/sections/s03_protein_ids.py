"""§3 — protein_ids: UniProt + RefSeq proteins, InterPro/Pfam domains,
antibodies, UniProt sequence features, BRENDA EC + (NEW 2026-05-31)
UniProt CC narratives + named isoforms.

CC ('comments') is the curated free-text description set: function,
subunit, subcellular_location, tissue_specificity, disease, ptm, etc.
Closes the audit's #1 content gap. Per-canonical-uniprot, single entry call."""
from collections import Counter
from atlas.biobtree import map_all
from atlas.gene.sections.base import Section
from atlas.page.uniprot_cc import fetch_cc, strip_evidence_codes

CHAINS = (
    ">>ensembl>>uniprot",
    '>>ensembl>>refseq[type=="protein"]',
    ">>uniprot>>interpro",
    ">>uniprot>>pfam",
    ">>uniprot>>antibody",
    ">>uniprot>>ufeature",
    ">>uniprot>>brenda",
    ">>uniprot>>brenda>>brenda_kinetics",
    ">>uniprot>>rhea",                 # catalyzed reactions (substrate → product)
    ">>hgnc>>entrez",                  # anchor-resolved: NCBI gene id + summary
    ">>hgnc>>clingen_dosage",          # anchor-resolved: haplo/triplo scores
    ">>hgnc>>depmap",                  # anchor-resolved: CRISPR fitness
    ">>entrez>>generif",               # per-gene PMID-anchored claims
)
DATASETS = ("uniprot", "ensembl", "refseq", "interpro", "pfam", "antibody",
            "ufeature", "brenda", "rhea", "entrez", "clingen_dosage", "depmap", "generif")

def collect(a):
    bundle = {
        "section": "03_protein_ids", "symbol": a.symbol,
        "hgnc_id": a.hgnc_id,
        "reviewed_uniprot": list(a.reviewed_uniprots),
        "canonical_uniprot": a.canonical_uniprot,
        "ensembl_id": a.ensembl_id,
    }

    allu = map_all(a.ensembl_id, ">>ensembl>>uniprot") if a.ensembl_id else []
    bundle["uniprot_all"] = [t["id"] for t in allu]
    bundle["uniprot_count"] = len(allu)

    if a.ensembl_id:
        nps = map_all(a.ensembl_id, '>>ensembl>>refseq[type=="protein"]')
        bundle["refseq_protein"] = [{"id": t["id"], "mane": t.get("is_mane_select") == "true"}
                                    for t in nps]
        bundle["refseq_protein_count"] = len(nps)

    # domains/families + antibody + UniProt sequence features, unioned across
    # all reviewed products. (ufeature was previously leaking ortholog features
    # — fixed upstream, BIOBTREE_ISSUES.md #11 RESOLVED — so no post-filter
    # needed anymore.)
    interpro, pfam, antibody = {}, set(), 0
    ufeatures = []
    brenda_ec = []
    brenda_kinetics = []
    rhea = {}                          # normalized equation -> {id, equation}
    for u in a.reviewed_uniprots:
        for t in map_all(u, ">>uniprot>>interpro"):
            interpro[t["id"]] = {"id": t["id"], "name": t.get("short_name"),
                                 "type": t.get("type")}
        pfam.update(t["id"] for t in map_all(u, ">>uniprot>>pfam"))
        antibody += len(map_all(u, ">>uniprot>>antibody"))
        for t in map_all(u, ">>uniprot>>ufeature", cap=100):
            ufeatures.append({"uniprot": u, "id": t["id"], "type": t.get("type"),
                              "description": t.get("description"),
                              "begin": t.get("location_begin"),
                              "end": t.get("location_end")})
        # BRENDA enzyme classification — EC number + name + summary stats.
        # Non-enzyme proteins (TFs, inhibitors) return nothing; the bundle
        # list stays empty and the render block elides.
        for t in map_all(u, ">>uniprot>>brenda"):
            brenda_ec.append({"uniprot": u, "ec": t.get("id"),
                              "name": t.get("recommended_name"),
                              "organism_count": t.get("organism_count"),
                              "substrate_count": t.get("substrate_count"),
                              "inhibitor_count": t.get("inhibitor_count"),
                              "km_count": t.get("km_count"),
                              "kcat_count": t.get("kcat_count")})
        # Per-substrate kinetics (Km/kcat). TIER 2: the projection carries no
        # per-row organism, so min–max spans an aggregate across organisms /
        # conditions — show the range, not a single point value. Keep only rows
        # with a real measurement; rank best-characterized (most Km values) first.
        def _i(v):
            try:
                return int(v)
            except (TypeError, ValueError):
                return 0
        for t in map_all(u, ">>uniprot>>brenda>>brenda_kinetics", cap=150):
            kmn, kcatn = _i(t.get("km_count")), _i(t.get("kcat_count"))
            if not (kmn or kcatn):
                continue
            sub = (t.get("id") or "").split("|")[-1]
            # Drop malformed substrate names (reaction-direction fragments like
            # "-0.0032 {NADPH}") — keep clean chemical names only.
            if not sub or sub.startswith("-") or "{" in sub:
                continue
            brenda_kinetics.append({
                "substrate": sub, "km_count": kmn, "kcat_count": kcatn,
                "min_km": t.get("min_km"), "max_km": t.get("max_km")})
        # Rhea — the catalyzed reaction(s), substrate → product. biobtree now
        # projects the human-readable `equation` (BIOBTREE_ISSUES #29). A Rhea
        # group is one reaction in several forms: an undirected master (UN/BI,
        # "A = B") plus LR/RL directional variants ("A => B" / "B => A"). Collapse
        # the arrows AND the two reversed orderings to a single reaction, and
        # prefer the undirected master's equation for display.
        for t in map_all(u, ">>uniprot>>rhea"):
            eq = (t.get("equation") or "").strip()
            if not eq:
                continue
            norm = " ".join(eq.replace("<=>", "=").replace("=>", "=")
                              .replace("<=", "=").split())
            sides = [s.strip() for s in norm.split("=")]
            key = " = ".join(sorted(sides)) if len(sides) == 2 else norm
            pref = (t.get("direction") or "").upper() in ("UN", "BI", "")
            prev = rhea.get(key)
            if prev is None or (pref and not prev["pref"]):
                rhea[key] = {"id": t.get("id"), "equation": norm, "pref": pref}
    bundle["rhea_reactions"] = [{"id": r["id"], "equation": r["equation"]}
                                for r in list(rhea.values())[:12]]
    brenda_kinetics.sort(key=lambda r: -r["km_count"])
    bundle["brenda_kinetics"] = brenda_kinetics[:15]
    bundle["brenda_kinetics_total"] = len(brenda_kinetics)
    bundle["interpro"] = list(interpro.values())
    bundle["pfam"] = sorted(pfam)
    bundle["antibody_count"] = antibody
    bundle["ufeature_counts"] = dict(Counter(f["type"] for f in ufeatures))
    bundle["ufeatures"] = ufeatures
    bundle["brenda_ec"] = brenda_ec

    # UniProt CC narratives + named isoforms (one entry call on canonical
    # accession). Returns {} for non-protein-coding genes (no canonical
    # uniprot) or unreviewed accessions — bundle keys stay empty + renderer
    # elides cleanly.
    cc_blob = fetch_cc(a.canonical_uniprot) if a.canonical_uniprot else {}
    raw_comments = cc_blob.get("comments") or {}
    # Strip evidence codes once at the bundle layer so every downstream
    # consumer (render, declarative-lead, JSON-LD) sees clean text.
    bundle["cc"] = {k: strip_evidence_codes(v) for k, v in raw_comments.items()
                    if isinstance(v, str) and v.strip()}
    bundle["isoforms"] = cc_blob.get("isoforms") or []
    bundle["protein_name"] = cc_blob.get("name")  # primary UniProt name
    bundle["alternative_names"] = cc_blob.get("alternative_names") or []

    # NCBI Entrez summary — independent narrative complementary to UniProt CC.
    # Resolved at anchor time so zero extra cost here; just pass through.
    bundle["ncbi_summary"] = a.ncbi_summary or ""
    bundle["entrez_id"] = a.entrez_id

    # ClinGen dosage sensitivity (1 row/gene). Empty {} for un-curated genes
    # (~80% of human protein-coding genome — ClinGen prioritizes clinical relevance).
    bundle["clingen_dosage"] = a.clingen_dosage or {}

    # DepMap CRISPR fitness summary — drives target-quality reasoning in §10
    # (we also expose here so the gene page's protein-level block is self-contained).
    bundle["depmap"] = a.depmap or {}

    # GeneRIFs — NCBI per-gene PMID-anchored claims. ~hundreds per popular gene;
    # cap at top-40 (insertion order from biobtree ~ chronological). Each entry
    # carries gene_id + text only; full PMID list comes from entry() but we don't
    # call it here — the id format `geneid_pmid_idx` already embeds the PMID.
    # Drop HuGE Navigator stock phrases ("Observational study of gene-disease
    # association. (HuGE Navigator)" and kin) — templated, near-identical, and on
    # well-studied genes they crowd out the mechanistic RIFs. Filter BEFORE the
    # cap so 40 substantive findings survive (TP53 had 27/40 as this boilerplate).
    generifs = []
    if a.entrez_id:
        for r in map_all(a.entrez_id, ">>entrez>>generif", cap=1):
            text = r.get("text") or ""
            if "HuGE Navigator" in text:
                continue
            pmid = (r.get("id") or "").split("_")[1] if "_" in (r.get("id") or "") else None
            generifs.append({"id": r.get("id"), "pmid": pmid, "text": text})
            if len(generifs) >= 40:
                break
    bundle["generifs"] = generifs
    bundle["generif_count_in_page"] = len(generifs)

    return bundle

SECTION = Section(
    id="3", name="protein_ids",
    description="UniProt accessions (canonical+all), RefSeq proteins, InterPro/Pfam domains, antibodies, UniProt features",
    needs=("hgnc_id", "ensembl_id", "reviewed_uniprots", "canonical_uniprot"),
    produces=("reviewed_uniprot", "uniprot_all", "refseq_protein", "interpro",
              "pfam", "antibody_count", "ufeatures", "ufeature_counts",
              "brenda_ec", "brenda_kinetics", "brenda_kinetics_total", "rhea_reactions",
              "cc", "isoforms", "protein_name",
              "alternative_names", "ncbi_summary", "entrez_id",
              "clingen_dosage", "depmap", "generifs"),
    datasets=DATASETS, chains=CHAINS, collect_fn=collect,
    # refseq_protein follows the same REVIEWED-only fluctuation as
    # refseq_mrna (BIOBTREE_ISSUES.md #11 — see §2 shrinkable note).
    # ufeatures shrinks when UniProt re-curates feature annotations
    # (e.g. demotes "Probable" → "By similarity" and drops them).
    shrinkable=("refseq_protein", "ufeatures"),
)
