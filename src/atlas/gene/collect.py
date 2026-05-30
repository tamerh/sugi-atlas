#!/usr/bin/env python3
"""Deterministic biobtree collector — prototype (gene section 01_gene_ids).

Given a gene symbol, runs the canonical chain plan against the LOCAL biobtree
API (no model), returning a structured bundle + the api_calls log. Proves that
the data a section needs can be gathered deterministically in code.
"""
import json, sys, re, urllib.parse, urllib.request

API = "http://127.0.0.1:8000/api"
CALLS = []  # reproducibility log, same shape as the page frontmatter

def _get(path, params):
    qs = urllib.parse.urlencode(params)
    url = f"{API}/{path}?{qs}"
    with urllib.request.urlopen(url, timeout=15) as r:
        body = json.load(r)
    CALLS.append({"path": path, "params": params})
    return body

def search(term, source=None):
    p = {"i": term}
    if source: p["s"] = source
    return _get("search", p)

def entry(identifier, source):
    return _get("entry", {"i": identifier, "s": source})

def bbmap(ids, chain, page=None):
    params = {"i": ids, "m": chain}
    if page:
        params["p"] = page  # pagination cursor param is `p` (NOT `page`; FastAPI
                            # silently drops unknown params -> looked like a stuck cursor)
    return _get("map", params)

def rows(resp):
    """search/data rows -> list of dicts keyed by schema columns."""
    cols = resp.get("schema", "").split("|")
    return [dict(zip(cols, r.split("|"))) for r in resp.get("data", [])]

def map_targets(resp):
    """map -> flat list of target dicts keyed by schema columns."""
    cols = resp.get("schema", "").split("|")
    out = []
    for m in (resp.get("mappings") or []):
        for t in (m.get("targets") or []):
            # some fields embed an escaped pipe (\|), e.g. bgee_evidence ids
            # "ENSG..\|UBERON:.." — split on real separators only.
            parts = [p.replace("\x00", "|") for p in t.replace("\\|", "\x00").split("|")]
            out.append(dict(zip(cols, parts)))
    return out

def map_all(ids, chain, cap=60):
    """All target dicts across pages, deduped. Pagination uses the `p=` cursor
    param (see bbmap) and works correctly. Dedupe + stop-on-no-new make it
    loop-safe; `cap` bounds total pages (cap*100 rows) for very large sets."""
    out, page, n, seen = [], None, 0, set()
    while True:
        resp = bbmap(ids, chain, page)
        new = 0
        for t in map_targets(resp):
            key = t.get("id") or tuple(t.values())
            if key in seen:
                continue
            seen.add(key); out.append(t); new += 1
        pg = resp.get("pagination", {}) or {}
        nxt = pg.get("next_token")
        if not pg.get("has_next") or n >= cap or not nxt or new == 0:
            break
        page, n = nxt, n + 1
    return out

def resolve_hgnc(symbol):
    """Resolve a gene symbol -> (hgnc_id, hgnc_entry) robustly.

    Filter the search to the hgnc dataset (the unfiltered top-50 can omit the
    hgnc row for high-xref genes), then disambiguate ambiguous symbols (e.g.
    'AR' matches amphiregulin AND androgen receptor) by the exact approved
    symbol from each candidate's hgnc entry.
    """
    cand = [r["id"] for r in rows(search(symbol, source="hgnc"))
            if re.match(r"HGNC:\d+$", r.get("id", ""))]
    if not cand:  # fallback: unfiltered
        cand = [r["id"] for r in rows(search(symbol)) if r.get("dataset") == "hgnc"]
    if not cand:
        sys.exit(f"no HGNC row for {symbol}")
    if len(cand) == 1:
        return cand[0], entry(cand[0], "hgnc")
    for cid in cand[:8]:
        he = entry(cid, "hgnc")
        syms = he.get("Attributes", {}).get("Hgnc", {}).get("symbols", [])
        if symbol.upper() in [s.upper() for s in syms]:
            return cid, he
    return cand[0], entry(cand[0], "hgnc")  # last resort

def collect_gene_ids(symbol):
    bundle = {"section": "01_gene_ids", "symbol": symbol}

    # 1-2. resolve HGNC id + entry (core identifiers + xref count table = gold)
    hgnc_id, he = resolve_hgnc(symbol)
    bundle["hgnc_id"] = hgnc_id
    attrs = he.get("Attributes", {}).get("Hgnc", {})
    bundle["name"] = (attrs.get("names") or [None])[0]
    bundle["hgnc"] = {
        "symbol": (attrs.get("symbols") or [None])[0],
        "name": (attrs.get("names") or [None])[0],
        "location": attrs.get("location"),
        "locus_type": attrs.get("locus_type"),
        "status": attrs.get("status"),
        "aliases": attrs.get("aliases", []),
    }
    bundle["xref_counts"] = {r.split("|")[0]: int(r.split("|")[1])
                             for r in he.get("xrefs", {}).get("data", [])}

    # 3. cross-refs to other id systems
    ens = map_all(hgnc_id, ">>hgnc>>ensembl")
    bundle["ensembl_id"] = ens[0]["id"] if ens else None
    bundle["mim"] = [t["id"] for t in map_all(hgnc_id, ">>hgnc>>mim")]
    bundle["entrez"] = [t["id"] for t in map_all(hgnc_id, ">>hgnc>>entrez")]

    # 4. ensembl entry (biotype/genome) if we have an id
    if bundle["ensembl_id"]:
        ee = entry(bundle["ensembl_id"], "ensembl")
        ea = ee.get("Attributes", {}).get("Ensembl", {})
        bundle["ensembl"] = {"biotype": ea.get("biotype"), "genome": ea.get("genome")}

    return bundle


def collect_transcripts(symbol):
    """Section 02_transcripts: Ensembl transcripts (+biotype, count), RefSeq
    (mark MANE Select), CCDS, and exons of the canonical/MANE transcript."""
    bundle = {"section": "02_transcripts", "symbol": symbol}
    hgnc_id, _ = resolve_hgnc(symbol)
    ens = map_all(hgnc_id, ">>hgnc>>ensembl")
    ensembl_id = ens[0]["id"] if ens else None
    bundle["ensembl_id"] = ensembl_id
    if not ensembl_id:
        return bundle

    tr = map_all(ensembl_id, ">>ensembl>>transcript")
    bundle["ensembl_transcripts"] = [{"id": t["id"], "biotype": t.get("biotype")} for t in tr]
    bundle["ensembl_transcript_count"] = len(tr)

    # curated mRNA accessions (filtered to skip predicted XM_/XR_ noise). NOTE:
    # full mRNA list (paginated). MANE Select is also fetched via its own
    # filter below so it is never missed.
    rs = map_all(ensembl_id, '>>ensembl>>refseq[type=="mRNA"]')
    bundle["refseq_mrna"] = [t["id"] for t in rs]
    bundle["refseq_mrna_count"] = len(rs)

    bundle["ccds"] = [t["id"] for t in map_all(ensembl_id, ">>ensembl>>ccds")]

    # canonical transcript = MANE-Select mRNA -> its Ensembl ENST -> exons.
    # Fetched via a dedicated is_mane_select filter (NOT by scanning the full
    # refseq list, which is page-capped and would miss it for big genes).
    mane = map_all(ensembl_id, ">>ensembl>>refseq[is_mane_select==true]")
    mane_nm = next((t["id"] for t in mane if t.get("type") == "mRNA"), None)
    canonical = None
    if mane_nm:
        ct = map_all(mane_nm, ">>refseq>>transcript")
        canonical = ct[0]["id"] if ct else None
    if not canonical and tr:
        canonical = tr[0]["id"]  # fallback: first transcript
    bundle["mane_select_refseq"] = mane_nm
    bundle["canonical_transcript"] = canonical
    if canonical:
        ex = map_all(canonical, ">>transcript>>exon")
        bundle["canonical_exons"] = [{"id": e["id"], "start": e.get("start"),
                                      "end": e.get("end")} for e in ex]
        bundle["canonical_exon_count"] = len(ex)
    return bundle


def collect_protein_ids(symbol):
    """Section 03_protein_ids: all UniProt (mark canonical reviewed), RefSeq NP_
    proteins, InterPro + Pfam domains/families, antibody availability."""
    bundle = {"section": "03_protein_ids", "symbol": symbol}
    hgnc_id, _ = resolve_hgnc(symbol)
    bundle["hgnc_id"] = hgnc_id

    # reviewed (Swiss-Prot) products = HGNC's uniprot xrefs. Usually 1, but
    # dual-product genes have >1 (e.g. CDKN2A -> P42771 p16INK4a + Q8N726 p14ARF),
    # so domains must be gathered across ALL of them, not just the first.
    reviewed = [t["id"] for t in map_all(hgnc_id, ">>hgnc>>uniprot")]
    bundle["reviewed_uniprot"] = reviewed
    bundle["canonical_uniprot"] = reviewed[0] if reviewed else None

    ens = map_all(hgnc_id, ">>hgnc>>ensembl")
    ensembl_id = ens[0]["id"] if ens else None
    bundle["ensembl_id"] = ensembl_id

    # all UniProt accessions (reviewed + unreviewed) via Ensembl
    allu = map_all(ensembl_id, ">>ensembl>>uniprot") if ensembl_id else []
    bundle["uniprot_all"] = [t["id"] for t in allu]
    bundle["uniprot_count"] = len(allu)

    # RefSeq NP_ proteins
    if ensembl_id:
        nps = map_all(ensembl_id, '>>ensembl>>refseq[type=="protein"]')
        bundle["refseq_protein"] = [{"id": t["id"], "mane": t.get("is_mane_select") == "true"}
                                    for t in nps]
        bundle["refseq_protein_count"] = len(nps)

    # domains/families + antibody + UniProt sequence features, unioned across
    # all reviewed products. ufeature is keyed by source accession and biobtree
    # returns ortholog features too — filter id startswith "{u}_" per species.
    interpro, pfam, antibody = {}, set(), 0
    ufeatures = []
    for u in reviewed:
        for t in map_all(u, ">>uniprot>>interpro"):
            interpro[t["id"]] = {"id": t["id"], "name": t.get("short_name"),
                                 "type": t.get("type")}
        pfam.update(t["id"] for t in map_all(u, ">>uniprot>>pfam"))
        antibody += len(map_all(u, ">>uniprot>>antibody"))
        for t in map_all(u, ">>uniprot>>ufeature", cap=100):
            if not t["id"].startswith(u + "_"):
                continue  # ortholog feature
            ufeatures.append({"uniprot": u, "id": t["id"], "type": t.get("type"),
                              "description": t.get("description"),
                              "begin": t.get("location_begin"),
                              "end": t.get("location_end")})
    bundle["interpro"] = list(interpro.values())
    bundle["pfam"] = sorted(pfam)
    bundle["antibody_count"] = antibody
    from collections import Counter
    bundle["ufeature_counts"] = dict(Counter(f["type"] for f in ufeatures))
    bundle["ufeatures"] = ufeatures
    return bundle


def collect_structure(symbol):
    """Section 04_structure: experimental PDB structures (method + resolution)
    and AlphaFold predicted model (pLDDT), for the reviewed protein product(s).
    (PDB list fully paginated)."""
    bundle = {"section": "04_structure", "symbol": symbol}
    hgnc_id, _ = resolve_hgnc(symbol)
    bundle["hgnc_id"] = hgnc_id
    reviewed = [t["id"] for t in map_all(hgnc_id, ">>hgnc>>uniprot")]
    bundle["reviewed_uniprot"] = reviewed

    pdb, af = {}, []
    for u in reviewed:
        for t in map_all(u, ">>uniprot>>pdb"):
            pdb[t["id"]] = {"id": t["id"], "method": t.get("method"),
                            "resolution": t.get("resolution")}
        # Every UniProt protein has an AlphaFold model (id AF-<acc>-F1). biobtree's
        # alphafold map is EMPTY for very large/fragmented proteins (see
        # BIOBTREE_ISSUES.md), so construct the id and attach pLDDT when present.
        # F2..Fn for huge proteins (>2700 aa) are not enumerated.
        m = map_all(u, ">>uniprot>>alphafold")
        af.append({"id": f"AF-{u}-F1", "uniprot": u,
                   "plddt": m[0].get("global_metric") if m else None,
                   "fraction_plddt_very_high": m[0].get("fraction_plddt_very_high") if m else None})
    bundle["pdb"] = list(pdb.values())
    bundle["pdb_count"] = len(pdb)
    bundle["alphafold"] = af
    return bundle


def canonical_transcript(ensembl_id):
    """MANE-Select Ensembl transcript for a gene (fallback: first transcript)."""
    if not ensembl_id:
        return None
    mane = map_all(ensembl_id, ">>ensembl>>refseq[is_mane_select==true]")
    mane_nm = next((t["id"] for t in mane if t.get("type") == "mRNA"), None)
    if mane_nm:
        ct = map_all(mane_nm, ">>refseq>>transcript")
        if ct:
            return ct[0]["id"]
    tr = map_all(ensembl_id, ">>ensembl>>transcript")
    return tr[0]["id"] if tr else None


def xref_counts(entry_resp):
    return {r.split("|")[0]: int(r.split("|")[1])
            for r in entry_resp.get("xrefs", {}).get("data", [])}


def collect_variants(symbol):
    """Section 06_variants: ClinVar (total + per-class breakdown + top-30
    pathogenic), SpliceAI (total + top-30), AlphaMissense (total + top-30
    likely-pathogenic). Totals from xref tables (exact); per-class breakdown now
    fully paginated (real counts). 'condition' needs clinvar>>mondo (not collected)."""
    bundle = {"section": "06_variants", "symbol": symbol}
    hgnc_id, he = resolve_hgnc(symbol)
    bundle["hgnc_id"] = hgnc_id
    xc = xref_counts(he)
    bundle["clinvar_total"] = xc.get("clinvar", 0)
    bundle["spliceai_total"] = xc.get("spliceai", 0)

    classes = ["Pathogenic", "Likely pathogenic", "Uncertain significance",
               "Likely benign", "Benign"]
    breakdown, top_path = {}, []
    for cls in classes:
        rs = map_all(hgnc_id, f'>>hgnc>>clinvar[germline_classification=="{cls}"]')
        breakdown[cls] = len(rs)  # floor (capped ~100)
        if cls in ("Pathogenic", "Likely pathogenic"):
            for t in rs:
                if len(top_path) >= 30:
                    break
                top_path.append({"id": t["id"], "hgvs": t.get("name"),
                                 "classification": t.get("germline_classification")})
    bundle["clinvar_breakdown"] = breakdown
    bundle["top_pathogenic"] = top_path

    sp = sorted(map_all(hgnc_id, ">>hgnc>>spliceai"),
                key=lambda t: float(t.get("score") or 0), reverse=True)
    bundle["top_spliceai"] = [{"id": t["id"], "effect": t.get("effect"),
                               "score": t.get("score")} for t in sp[:30]]

    ens = map_all(hgnc_id, ">>hgnc>>ensembl")
    ct = canonical_transcript(ens[0]["id"] if ens else None)
    bundle["canonical_transcript"] = ct
    if ct:
        bundle["alphamissense_total"] = xref_counts(entry(ct, "transcript")).get("alphamissense", 0)
        am = sorted(map_all(ct, '>>transcript>>alphamissense[am_class=="likely_pathogenic"]'),
                    key=lambda t: float(t.get("am_pathogenicity") or 0), reverse=True)
        bundle["top_alphamissense"] = [{"id": t["id"], "variant": t.get("protein_variant"),
                                        "am_pathogenicity": t.get("am_pathogenicity")} for t in am[:30]]
    # dbSNP rsIDs via ENTREZ (direct hgnc>>dbsnp is unbacked — dbsnp xrefs entrez,
    # not hgnc; see BIOBTREE_ISSUES.md). Large set -> sampled, not exhaustive.
    dbs = map_all(hgnc_id, ">>hgnc>>entrez>>dbsnp", cap=2)
    bundle["dbsnp_sample"] = [{"id": t["id"], "pos": f"{t.get('chromosome')}:{t.get('position')}",
                               "change": f"{t.get('ref_allele')}>{t.get('alt_allele')}"} for t in dbs[:30]]
    bundle["dbsnp_sampled"] = len(dbs)
    return bundle


def collect_orthologs(symbol):
    """Section 05_orthologs: orthologous genes in model organisms (Ensembl
    Compara), id + symbol + organism."""
    bundle = {"section": "05_orthologs", "symbol": symbol}
    hgnc_id, _ = resolve_hgnc(symbol)
    ens = map_all(hgnc_id, ">>hgnc>>ensembl")
    ensembl_id = ens[0]["id"] if ens else None
    bundle["ensembl_id"] = ensembl_id
    orths = map_all(ensembl_id, ">>ensembl>>ortholog") if ensembl_id else []
    bundle["orthologs"] = [{"id": t["id"], "symbol": t.get("name"),
                            "organism": t.get("genome")} for t in orths]
    bundle["ortholog_count"] = len(orths)
    paras = map_all(ensembl_id, ">>ensembl>>paralog") if ensembl_id else []
    bundle["paralogs"] = [{"id": t["id"], "symbol": t.get("name")} for t in paras]
    bundle["paralog_count"] = len(paras)
    return bundle


def _resolve(symbol):
    """Common anchors: (hgnc_id, hgnc_entry, ensembl_id, canonical_uniprot)."""
    hgnc_id, he = resolve_hgnc(symbol)
    ens = map_all(hgnc_id, ">>hgnc>>ensembl")
    cu = map_all(hgnc_id, ">>hgnc>>uniprot")
    return hgnc_id, he, (ens[0]["id"] if ens else None), (cu[0]["id"] if cu else None)


def collect_pathways(symbol):
    """Section 07: Reactome pathways, MSigDB gene sets, GO terms (BP/MF/CC)."""
    bundle = {"section": "07_pathways", "symbol": symbol}
    hgnc_id, he, ensembl_id, uni = _resolve(symbol)
    xc = xref_counts(he)
    # all reviewed products (dual-product genes like CDKN2A → p16 + p14ARF carry
    # distinct Reactome/GO; the single canonical uniprot misses the other product).
    unis = [t["id"] for t in map_all(hgnc_id, ">>hgnc>>uniprot")]
    # Reactome: union every reviewed uniprot + the gene-level ensembl route, dedupe.
    rx = {}
    for u in unis:
        for t in map_all(u, ">>uniprot>>reactome"):
            rx[t["id"]] = t.get("name") or rx.get(t["id"])
    for t in (map_all(ensembl_id, ">>ensembl>>reactome") if ensembl_id else []):
        rx[t["id"]] = t.get("name") or rx.get(t["id"])
    bundle["reactome"] = [{"id": k, "name": v} for k, v in rx.items()]
    bundle["reactome_count"] = len(bundle["reactome"])
    msig = map_all(hgnc_id, ">>hgnc>>msigdb")
    bundle["msigdb"] = [{"id": t["id"], "name": t.get("standard_name"),
                         "collection": t.get("collection")} for t in msig]
    bundle["msigdb_total"] = xc.get("msigdb", len(msig))
    # GO: union every reviewed uniprot + the gene-level ensembl route. UniProt-GOA
    # vs Ensembl annotation diverge ~20%, and per-product GO differs (CDKN2A
    # p14ARF mitophagy terms absent from p16/ensembl), dedupe by GO id.
    go_map = {}
    for u in unis:
        for t in map_all(u, ">>uniprot>>go"):
            go_map[t["id"]] = t
    for t in (map_all(ensembl_id, ">>ensembl>>go") if ensembl_id else []):
        go_map[t["id"]] = t
    grouped = {"biological_process": [], "molecular_function": [], "cellular_component": []}
    for t in go_map.values():
        grouped.setdefault(t.get("type"), []).append({"id": t["id"], "name": t.get("name")})
    bundle["go"] = grouped
    bundle["go_counts"] = {k: len(v) for k, v in grouped.items()}
    return bundle


def collect_interactions(symbol):
    """Section 08: PPIs (STRING/IntAct/BioGRID partners + counts) and similarity
    (ESM2/Diamond). NOTE: terminal >>uniprot gives partner IDs only — per-edge
    confidence scores need the intermediate interaction record (refinement)."""
    bundle = {"section": "08_interactions", "symbol": symbol}
    _, he, _, uni = _resolve(symbol)
    def _f(x):
        try: return float(x)
        except (TypeError, ValueError): return 0.0
    # interaction RECORDS (NOT >>...>>uniprot, which collapses to bare partner ids)
    # carry the scores/evidence — see "general fields at map, detail before the
    # terminal hop". STRING score is confidence x1000.
    st = map_all(uni, ">>uniprot>>string_interaction") if uni else []
    st.sort(key=lambda t: _f(t.get("score")), reverse=True)
    bundle["string"] = [{"partner": t.get("uniprot_b"), "score": t.get("score")} for t in st[:30]]
    bundle["string_count"] = len(st)
    ia = map_all(uni, ">>uniprot>>intact") if uni else []
    ia.sort(key=lambda t: _f(t.get("confidence_score")), reverse=True)
    bundle["intact"] = [{"a": t.get("protein_a_gene"), "b": t.get("protein_b_gene"),
                         "type": t.get("interaction_type"),
                         "score": t.get("confidence_score")} for t in ia[:30]]
    bundle["intact_count"] = len(ia)
    bg = map_all(uni, ">>uniprot>>biogrid_interaction") if uni else []
    bundle["biogrid"] = [{"partner": t.get("interactor_b_symbol"),
                          "method": t.get("experimental_system")} for t in bg[:30]]
    bundle["biogrid_count"] = len(bg)
    def partners(ds):
        return [t["id"] for t in (map_all(uni, f">>uniprot>>{ds}>>uniprot") if uni else [])]
    bundle["esm2_similar"] = partners("esm2_similarity")
    bundle["diamond_similar"] = partners("diamond_similarity")
    # SIGNOR signaling relationships (directed, with effect + mechanism)
    sig = map_all(uni, ">>uniprot>>signor") if uni else []
    bundle["signor"] = [{"a": t.get("entity_a"), "b": t.get("entity_b"),
                         "effect": t.get("effect"), "mechanism": t.get("mechanism")} for t in sig]
    bundle["signor_count"] = len(sig)
    return bundle


def collect_tf_regulation(symbol):
    """Section 09: CollecTRI downstream targets / upstream regulators + JASPAR
    motifs. is_tf inferred from having downstream targets or motifs."""
    bundle = {"section": "09_tf_regulation", "symbol": symbol}
    hgnc_id, he, _, uni = _resolve(symbol)
    # direction-filtered (tf_gene/target_gene) and fully paginated, so even
    # high-degree TFs (TP53: 1207 CollecTRI records) get complete targets.
    down = map_all(hgnc_id, f'>>hgnc>>collectri[tf_gene=="{symbol}"]')
    up = map_all(hgnc_id, f'>>hgnc>>collectri[target_gene=="{symbol}"]')
    bundle["downstream_targets"] = [{"target": r.get("target_gene"),
                                     "regulation": r.get("regulation")} for r in down]
    bundle["downstream_count"] = len(down)
    bundle["upstream_regulators"] = [{"regulator": r.get("tf_gene"),
                                      "regulation": r.get("regulation")} for r in up]
    bundle["jaspar_motifs"] = [{"id": t["id"], "name": t.get("name"), "class": t.get("class"),
                                "family": t.get("family")}
                               for t in (map_all(uni, ">>uniprot>>jaspar") if uni else [])]
    bundle["is_transcription_factor"] = bool(down or bundle["jaspar_motifs"])
    return bundle


def collect_drugs(symbol):
    """Section 10: ChEMBL targets + targeting molecules, PharmGKB. is_target from
    having ChEMBL targets. Clinical-trial enumeration left as a refinement."""
    bundle = {"section": "10_drugs", "symbol": symbol}
    hgnc_id, he, _, uni = _resolve(symbol)
    targets = map_all(uni, ">>uniprot>>chembl_target") if uni else []
    bundle["chembl_targets"] = [{"id": t["id"], "title": t.get("title"), "type": t.get("type")}
                                for t in targets]
    # Targeting DRUGS only (development phase >=1) — the raw chembl_molecule edge
    # returns thousands of screening compounds ordered by ID (useless); filtering
    # to phased drugs + ranking by highestDevelopmentPhase surfaces the real ones.
    def _phase(d): return int(d) if (d or "").isdigit() else 0
    drugs = {}
    for t in targets:
        for m in map_all(t["id"], '>>chembl_target>>chembl_molecule[highestDevelopmentPhase>=1]', cap=10):
            drugs[m["id"]] = {"id": m["id"], "name": m.get("name"), "type": m.get("type"),
                              "phase": m.get("highestDevelopmentPhase")}
    bundle["molecules"] = sorted(drugs.values(), key=lambda d: _phase(d["phase"]), reverse=True)
    bundle["molecule_count"] = len(drugs)  # drugs in development/approved (not screening hits)
    bundle["pharmgkb"] = [{"id": t["id"], "vip": t.get("is_vip"),
                           "cpic_guideline": t.get("has_cpic_guideline")}
                          for t in map_all(hgnc_id, ">>hgnc>>pharmgkb_gene")]
    # binding affinities (BindingDB can be huge -> sample first pages)
    bd = map_all(uni, ">>uniprot>>bindingdb", cap=2) if uni else []
    bundle["bindingdb_sample"] = [{"ligand": t.get("ligand_name"), "ki": t.get("ki"),
                                   "ic50": t.get("ic50")} for t in bd[:30]]
    bundle["bindingdb_sampled"] = len(bd)
    # Clinical trials via the DISEASE route (gene -> MONDO -> clinical_trials) —
    # biobtree's intended pattern (clinical_trials xrefs MONDO/ChEMBL, no gene->
    # drug->trials test exists). The chembl_molecule>>clinical_trials path is
    # avoided: ChEMBL target->molecule is bioactivity-based, so off-target approved
    # drugs (Levodopa@phase4 vs EGFR) and their unrelated trials swamp it. This is
    # trials for the gene's associated diseases (not drug-specific, but clean).
    trials = {}
    for chain in (">>hgnc>>gencc>>mondo>>clinical_trials",
                  ">>hgnc>>clinvar>>mondo>>clinical_trials"):
        for t in map_all(hgnc_id, chain, cap=2):
            trials[t["id"]] = {"id": t["id"], "title": t.get("brief_title"),
                               "phase": t.get("phase"), "status": t.get("overall_status")}
    def _ph(t):
        p = (t.get("phase") or "").upper().replace("PHASE", "").strip()
        return int(p) if p.isdigit() else 0
    bundle["disease_trials"] = sorted(trials.values(), key=_ph, reverse=True)
    bundle["disease_trial_count"] = len(trials)
    bundle["is_drug_target"] = bool(targets)
    return bundle


def collect_expression(symbol):
    """Section 11: Bgee gene-level expression summary + single-cell datasets
    (scxa). NOTE: per-tissue scores need bgee_evidence (refinement)."""
    bundle = {"section": "11_expression", "symbol": symbol}
    _, he, ensembl_id, uni = _resolve(symbol)
    bgee = map_all(ensembl_id, ">>ensembl>>bgee") if ensembl_id else []
    if bgee:
        b = bgee[0]
        bundle["bgee"] = {"breadth": b.get("expression_breadth"),
                          "present_calls": b.get("total_present_calls"),
                          "max_expression_score": b.get("max_expression_score")}
    # per-tissue expression (Bgee evidence) — top by expression score
    def _num(x):
        try: return float(x)
        except (TypeError, ValueError): return 0.0
    tissues = map_all(ensembl_id, ">>ensembl>>bgee>>bgee_evidence") if ensembl_id else []
    tissues.sort(key=lambda t: _num(t.get("expression_score")), reverse=True)
    bundle["top_tissues"] = [{"tissue": t.get("anatomical_entity_name"),
                              "score": t.get("expression_score"), "rank": t.get("expression_rank"),
                              "quality": t.get("call_quality")} for t in tissues[:30]]
    bundle["tissue_count"] = len(tissues)
    # FANTOM5 CAGE gene-level expression (TPM + breadth) — complements Bgee.
    f5 = map_all(ensembl_id, ">>ensembl>>fantom5_gene") if ensembl_id else []
    if f5:
        x = f5[0]
        bundle["fantom5"] = {"tpm_average": x.get("tpm_average"), "tpm_max": x.get("tpm_max"),
                             "samples_expressed": x.get("samples_expressed"),
                             "breadth": x.get("expression_breadth")}
    # FANTOM5 CAGE promoters (alternative TSS usage) — sorted by activity.
    def _num(x):
        try: return float(x)
        except (TypeError, ValueError): return 0.0
    proms = map_all(ensembl_id, ">>ensembl>>fantom5_promoter") if ensembl_id else []
    proms.sort(key=lambda t: _num(t.get("tpm_average")), reverse=True)
    bundle["fantom5_promoters"] = [{"id": t["id"], "tpm_average": t.get("tpm_average"),
                                    "samples_expressed": t.get("samples_expressed")}
                                   for t in proms]
    bundle["single_cell_datasets"] = [{"id": t["id"], "description": t.get("description"),
                                       "cells": t.get("number_of_cells")}
                                      for t in (map_all(ensembl_id, ">>ensembl>>scxa") if ensembl_id else [])]
    return bundle


def collect_diseases(symbol):
    """Section 12: Mendelian (OMIM/GenCC/Mondo/Orphanet), phenotypes (HPO),
    complex-disease (GWAS)."""
    bundle = {"section": "12_diseases", "symbol": symbol}
    hgnc_id, he, _, _ = _resolve(symbol)
    xc = xref_counts(he)
    bundle["gene_omim"] = [f"MIM:{t['id']}" for t in map_all(hgnc_id, ">>hgnc>>mim")]
    # disease OMIM phenotype ids (gene>>mim is only the gene; disease MIMs come
    # through MONDO): >>hgnc>>clinvar>>mondo>>mim
    bundle["disease_omim"] = [f"MIM:{t['id']}" for t in map_all(hgnc_id, ">>hgnc>>clinvar>>mondo>>mim")]
    bundle["gencc"] = [{"disease": t.get("disease_title"), "classification": t.get("classification_title"),
                        "inheritance": t.get("moi_title")} for t in map_all(hgnc_id, ">>hgnc>>gencc")]
    # MONDO: union clinvar + gencc routes (curated + clinical)
    mondo = {}
    for ch in (">>hgnc>>clinvar>>mondo", ">>hgnc>>gencc>>mondo"):
        for t in map_all(hgnc_id, ch):
            mondo[t["id"]] = t.get("name") or mondo.get(t["id"])
    bundle["mondo"] = [{"id": k, "name": v} for k, v in mondo.items()]
    # Orphanet: union clinvar route (ids) + mondo route (carries names)
    orph = {}
    for ch in (">>hgnc>>clinvar>>orphanet", ">>hgnc>>clinvar>>mondo>>orphanet"):
        for t in map_all(hgnc_id, ch):
            orph[t["id"]] = t.get("name") or orph.get(t["id"])
    bundle["orphanet"] = [{"id": f"Orphanet:{k}", "name": v} for k, v in orph.items()]
    # HPO: gene-level (hgnc>>hpo) + disease-level phenotypes (via mondo)
    hpo = {}
    for ch in (">>hgnc>>hpo", ">>hgnc>>clinvar>>mondo>>hpo"):
        for t in map_all(hgnc_id, ch):
            hpo[t["id"]] = t.get("name") or hpo.get(t["id"])
    bundle["hpo"] = [{"id": k, "name": v} for k, v in hpo.items()]
    bundle["hpo_total"] = xc.get("hpo", len(hpo))
    gwas = map_all(hgnc_id, ">>hgnc>>gwas")
    bundle["gwas"] = [{"id": t["id"], "trait": t.get("disease_trait"), "p_value": t.get("p_value")}
                      for t in gwas]
    bundle["gwas_total"] = xc.get("gwas", len(gwas))
    # gwas_study ids: gene-mapped (gwas>>gwas_study) + disease-associated
    # (mondo>>gwas_study). The latter catches genes whose GWAS aren't gene-mapped
    # directly (e.g. AR: >>hgnc>>gwas empty) but are reachable via their diseases.
    gws = {}
    for ch in (">>hgnc>>gwas>>gwas_study", ">>hgnc>>clinvar>>mondo>>gwas_study"):
        for t in map_all(hgnc_id, ch):
            gws[t["id"]] = 1
    bundle["gwas_studies"] = list(gws)
    return bundle


SECTIONS = {
    "1": collect_gene_ids,
    "2": collect_transcripts,
    "3": collect_protein_ids,
    "4": collect_structure,
    "5": collect_orthologs,
    "6": collect_variants,
    "7": collect_pathways,
    "8": collect_interactions,
    "9": collect_tf_regulation,
    "10": collect_drugs,
    "11": collect_expression,
    "12": collect_diseases,
}

if __name__ == "__main__":
    sym = sys.argv[1] if len(sys.argv) > 1 else "TP53"
    sec = sys.argv[2] if len(sys.argv) > 2 else "1"
    b = SECTIONS[sec](sym)
    print(json.dumps(b, indent=2))
    print(f"\n--- {len(CALLS)} api calls ---", file=sys.stderr)
    for c in CALLS:
        print(f"  {c['path']}({c['params']})", file=sys.stderr)
