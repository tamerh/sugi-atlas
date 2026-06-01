"""Build the deterministic one-sentence lead for an Atlas gene page.

The first thing AI agents extract as the gene's definition. Pure function of
the bundle (no LLM, no network) — every fact in the sentence is sourced from
biobtree via the §1 / §3 section bundles, so there's no hallucination surface.

Sits above the LLM exec summary so it dominates the first-extracted text.
Format is plain Markdown (no blockquote) so AI clients see clean prose:

    **TP53** (tumor protein p53, HGNC:11998) is a protein-coding gene on
    chromosome 17p13.1, encoding the reviewed UniProt protein P04637.

The function is intentionally tolerant — missing optional fields are dropped
silently rather than crashing or emitting "None". The output is always a
single, well-formed English sentence ending in a period.
"""

# Ensembl biotype -> human-readable gene class. The fallback for unmapped
# biotypes is to scan the HGNC locus_type string for keywords, then default
# to "gene".
BIOTYPE_HUMAN = {
    "protein_coding": "protein-coding gene",
    "lncRNA": "long non-coding RNA gene",
    "miRNA": "microRNA gene",
    "snoRNA": "small nucleolar RNA gene",
    "snRNA": "small nuclear RNA gene",
    "rRNA": "ribosomal RNA gene",
    "scaRNA": "small Cajal body RNA gene",
    "misc_RNA": "non-coding RNA gene",
    "Mt_tRNA": "mitochondrial tRNA gene",
    "Mt_rRNA": "mitochondrial rRNA gene",
    "pseudogene": "pseudogene",
    "processed_pseudogene": "processed pseudogene",
    "unprocessed_pseudogene": "unprocessed pseudogene",
    "transcribed_unprocessed_pseudogene": "unprocessed pseudogene",
    "IG_V_gene": "immunoglobulin V gene",
    "IG_C_gene": "immunoglobulin C gene",
    "IG_J_gene": "immunoglobulin J gene",
    "IG_D_gene": "immunoglobulin D gene",
    "TR_V_gene": "T-cell receptor V gene",
    "TR_J_gene": "T-cell receptor J gene",
    "TR_C_gene": "T-cell receptor C gene",
    "TR_D_gene": "T-cell receptor D gene",
}

def _gene_class(b1):
    """Map biotype (preferred) / locus_type (fallback) to human-readable class."""
    biotype = ((b1.get("ensembl") or {}).get("biotype") or "").strip()
    if biotype in BIOTYPE_HUMAN:
        return BIOTYPE_HUMAN[biotype]
    locus = ((b1.get("hgnc") or {}).get("locus_type") or "").lower()
    if "protein product" in locus:
        return "protein-coding gene"
    if "long non-coding" in locus or "lncrna" in locus:
        return "long non-coding RNA gene"
    if "microrna" in locus or "mirna" in locus:
        return "microRNA gene"
    if "pseudogene" in locus:
        return "pseudogene"
    return "gene"

def declarative_sentence(bundle):
    """Compose the lead sentence from a full bundle dict {section_id: bundle_dict}.

    Sources used:
      bundle["1"] - HGNC symbol + name + location + locus_type, ensembl biotype
      bundle["3"] - canonical_uniprot (optional, only for protein-coding)

    Returns a single sentence string (no trailing newline)."""
    b1 = bundle.get("1") or {}
    b3 = bundle.get("3") or {}

    symbol  = b1.get("symbol") or "?"
    hgnc    = b1.get("hgnc") or {}
    name    = hgnc.get("name") or ""
    hgnc_id = b1.get("hgnc_id")
    location = hgnc.get("location")
    klass   = _gene_class(b1)

    parens = []
    if name and name != symbol:
        parens.append(name)
    if hgnc_id:
        parens.append(hgnc_id)
    head = f"**{symbol}** ({', '.join(parens)})" if parens else f"**{symbol}**"

    where = f" on chromosome {location}" if location else ""
    sentence = f"{head} is a {klass}{where}"

    # Append protein clause for protein-coding genes only.
    canon = b3.get("canonical_uniprot")
    if canon and klass == "protein-coding gene":
        protein_name = b3.get("protein_name")
        if protein_name:
            sentence += f", encoding **{protein_name}** ({canon})"
        else:
            sentence += f", encoding the reviewed UniProt protein {canon}"

    sentence += "."

    # Second sentence: lead with UniProt's curated FUNCTION when available.
    # This was the audit's #1 content gap — biobtree's 2026-05-31 refresh
    # finally exposed the CC block on uniprot entries. We surface the lead
    # sentence here so AI agents extract a real function description, not
    # just identifier facts.
    from atlas.page.uniprot_cc import first_sentence
    function = (b3.get("cc") or {}).get("function")
    if function:
        sentence += " " + first_sentence(function)
        if not sentence.endswith("."):
            sentence += "."

    # Answer-first precision-oncology verdict — promote the top CIViC predictive
    # association (§10) into the lead. The most clinically-loaded fact, otherwise
    # buried ~1000 lines down in the §10 table. Deterministic, CC0 (CIViC).
    from atlas.civic import predictive_verdict
    b10 = bundle.get("10") or {}
    verdict = predictive_verdict(b10.get("civic_evidence") or [])
    if verdict:
        total = b10.get("civic_association_total") or 0
        more = (f"; {total - 1} further curated variant–drug associations are listed below"
                if total > 1 else "")
        sentence += f" In precision oncology, {verdict}{more}."

    # Functional-genomics verdict — DepMap cancer-dependency + ClinGen dosage
    # (§3), surfaced only when notable (a dependency, or sufficient-evidence
    # haploinsufficiency).
    dep = _dependency_clause(b3)
    if dep:
        sentence += " " + dep

    return sentence


def _dependency_clause(b3):
    """Concise DepMap + ClinGen-dosage verdict; '' unless notable."""
    dm = b3.get("depmap") or {}
    cd = b3.get("clingen_dosage") or {}
    bits = []
    pct = dm.get("pct_dependent")
    try:
        p = float(pct)
    except (TypeError, ValueError):
        p = None
    if p is not None:
        if dm.get("common_essential") == "true":
            bits.append(f"a common-essential gene (DepMap: required in {pct}% of cancer cell lines)")
        elif dm.get("strongly_selective") == "true" or p >= 10:
            bits.append(f"a selective cancer dependency (DepMap: {pct}% of cell lines)")
    if cd.get("haplo_score") == "3":
        bits.append("haploinsufficient (ClinGen: sufficient evidence)")
    return ("It is " + " and ".join(bits) + ".") if bits else ""
