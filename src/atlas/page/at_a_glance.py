"""Deterministic "At a glance" flag-sheet for Atlas gene pages.

A scannable bullet list of the gene's binary / verdict facts — the one-line
notes that were previously buried mid-page (drug-target flag, TF flag, DepMap
dependency, ClinGen dosage, MANE Select, CIViC/intOGen verdicts). Pure dict
lookups over the collector bundle: no LLM, no network, no hallucination
surface.

Sits directly under the declarative lead sentence. The lead is the prose
"definition + single headline verdict"; this block is the at-a-glance fact
sheet you scan, distinct on purpose (no sentence synthesis, just labelled
facts). Bullets with no backing data are dropped; if nothing qualifies the
whole block elides (returns "").
"""
import re

from atlas.page import evidence

# ClinGen dosage score scale (same mapping the §-detail block uses).
_DOSAGE_SCALE = {
    "3": "sufficient evidence", "2": "emerging evidence",
    "1": "little evidence", "0": "no evidence",
    "30": "autosomal recessive", "40": "dosage sensitivity unlikely",
}

# Gene-disease validity classification strength (best first). Shared by
# ClinGen Gene-Disease Validity and GenCC.
_GD_RANK = {"definitive": 0, "strong": 1, "moderate": 2, "supportive": 3,
            "limited": 4, "disputed": 5, "refuted": 6}

# Minimum CollecTRI downstream targets before the TF flag is headline-worthy.
# Below this the flag is too weak to promote (it still appears in §9 Regulation).
_TF_MIN_TARGETS = 10


def _truthy(v):
    return v in (True, "true", "True", "1", 1, "yes", "Yes")


def _symbol(bundle):
    for sid in ("3", "2", "8", "6", "12", "4", "5"):
        s = (bundle.get(sid) or {}).get("symbol")
        if s:
            return s
    return ""


def _sense_gene_bullet(symbol):
    """For an antisense lncRNA (SYMBOL-AS1/-AS2…), orient the reader to its sense
    partner and quote the SENSE gene's biology — strictly attributed as the sense
    gene's, never imported into this transcript's evidence. Manifest-gated (only
    when the sense gene is a built page). '' otherwise."""
    m = re.match(r"^(.+?)-AS\d*$", symbol or "")
    if not m:
        return ""
    sense = m.group(1)
    from atlas.page import links
    url = links.gene_url(symbol=sense)
    if not url:
        return ""
    c = evidence.components_for("gene", sense)
    bits = []
    if c.get("gwas_count"):
        bits.append(f"{c['gwas_count']:,} GWAS associations")
    if c.get("drug_count"):
        bits.append(f"{c['drug_count']:,} ChEMBL molecules")
    if c.get("variant_count"):
        bits.append(f"{c['variant_count']:,} ClinVar variants")
    carries = (f" — the sense gene carries {', '.join(bits[:2])} "
               "(the sense gene's biology, not this transcript's)") if bits else ""
    return f"**Sense gene:** antisense to {links.maybe_link(sense, url)}{carries}"


def _gene_disease(b12) -> str:
    """Top curated gene–disease relationship (ClinGen validity + GenCC),
    ranked by classification strength, with a count of the rest. APOE's real
    headline (Alzheimer disease 2, Definitive) lives here, not in the flag set."""
    cands = []
    for c in (b12.get("clingen_validity") or []):
        cands.append((c.get("disease"), c.get("classification"), "ClinGen"))
    for g in (b12.get("gencc") or []):
        cands.append((g.get("disease"), g.get("classification"), "GenCC"))
    cands = [c for c in cands if c[0] and c[1]]
    if not cands:
        return ""
    cands.sort(key=lambda c: _GD_RANK.get((c[1] or "").strip().lower(), 9))
    # Dedup by disease (this is one gene's page) — GenCC/ClinGen repeat the same
    # gene-disease relationship across submissions; count distinct relationships,
    # not raw rows. Sort-first keeps the strongest classification per disease.
    seen, uniq = set(), []
    for c in cands:
        d = (c[0] or "").strip().lower()
        if d in seen:
            continue
        seen.add(d)
        uniq.append(c)
    cands = uniq
    disease, classification, source = cands[0]
    n_more = len(cands) - 1
    more = (f" — +{n_more} more curated relationship{'s' if n_more != 1 else ''}"
            if n_more > 0 else "")
    return f"**Gene–disease (curated):** {disease} ({classification}, {source}){more}"


def at_a_glance(bundle) -> str:
    """Compose the `## At a glance` markdown block from a full bundle dict.

    Sources: §2 (MANE), §3 (DepMap + ClinGen dosage), §9 (TF flag),
    §10 (drug-target flag + CIViC predictive total), §12 (intOGen driver)."""
    b2 = bundle.get("2") or {}
    b3 = bundle.get("3") or {}
    b6 = bundle.get("6") or {}
    b9 = bundle.get("9") or {}
    b10 = bundle.get("10") or {}
    b12 = bundle.get("12") or {}

    bullets = []

    # Affirmative non-coding line — the positional variant/disease/trial blocks
    # are scrubbed for non-coding genes (see pipeline._scrub_noncoding), so say
    # so explicitly rather than leaving a suspiciously empty page.
    noncoding = bundle.get("_noncoding")
    if noncoding:
        # Name the overlapping protein-coding gene(s) when known (from the
        # positional variant/splice records, §6) — so "positional, from an
        # overlapping gene" is concrete + navigable, not vague.
        from atlas.page import links
        overlap = (bundle.get("6") or {}).get("overlap_genes") or []
        linked = [links.maybe_link(g, links.gene_url(symbol=g)) for g in overlap[:4]]
        ov = (f" The variant/disease data at this locus belongs to the overlapping "
              f"protein-coding gene{'s' if len(linked) != 1 else ''} "
              f"{', '.join(linked)}, not this transcript.") if linked else ""
        bullets.append(f"**Gene type:** non-coding ({noncoding}) — no protein "
                       f"product, so no protein-based drug-target data (Atlas drug "
                       f"coverage is protein-target-based; RNA-targeting therapeutics "
                       f"aside). Variant/disease associations are omitted (they would "
                       f"be positional).{ov}")
        sense = _sense_gene_bullet(_symbol(bundle))
        if sense:
            bullets.append(sense)

    # Curated gene–disease relationship — often the real headline (esp. for
    # Mendelian / complex-disease genes the cancer/drug flags miss).
    gd = _gene_disease(b12)
    if gd:
        bullets.append(gd)

    # GWAS associations — common-disease genetics, complements the Mendelian
    # gene–disease line above.
    gwas = b12.get("gwas_total") or 0
    if gwas:
        bullets.append(f"**GWAS associations:** {gwas:,}"
                       + evidence.rank_clause("gene", "gwas_count", gwas))

    # Clinical variant burden (ClinVar) — total + the pathogenic / likely-
    # pathogenic floor. High clinical signal, populated for most disease genes.
    cv_total = b6.get("clinvar_total") or 0
    if cv_total:
        bd = b6.get("clinvar_breakdown") or {}
        parts = []
        if bd.get("Pathogenic"):
            parts.append(f"{bd['Pathogenic']} pathogenic")
        if bd.get("Likely pathogenic"):
            parts.append(f"{bd['Likely pathogenic']} likely-pathogenic")
        detail = f" — {', '.join(parts)}" if parts else ""
        bullets.append(f"**Clinical variants (ClinVar):** {cv_total:,} total{detail}"
                       + evidence.rank_clause("gene", "variant_count", cv_total))

    # Phenotypes (HPO) — clinical feature count across associated conditions.
    hpo = b12.get("hpo_total") or 0
    if hpo:
        bullets.append(f"**Phenotypes (HPO):** {hpo:,}")

    # Druggable target — drug-discovery headline flag.
    if _truthy(b10.get("is_drug_target")):
        mc = b10.get("molecule_count") or 0
        extra = f" — {mc:,} molecules with ChEMBL bioactivity" if mc else ""
        bullets.append(f"**Druggable target:** yes{extra}"
                       + evidence.rank_clause("gene", "drug_count", mc))

    # Precision-oncology evidence (CIViC) — curated variant–drug associations.
    civ_total = b10.get("civic_association_total") or 0
    if civ_total:
        bullets.append(
            f"**Precision-oncology evidence (CIViC):** {civ_total} curated "
            f"variant–drug association{'s' if civ_total != 1 else ''}")

    # Cancer driver classification (intOGen).
    intogen = b12.get("intogen") or {}
    if intogen:
        role = {
            "Act": "activating (oncogene-like)",
            "LoF": "loss-of-function (tumor-suppressor-like)",
            "ambiguous": "ambiguous (mixed evidence)",
        }.get(intogen.get("role") or "", intogen.get("role") or "")
        cancers = [c for c in (intogen.get("cancer_types") or "").split(",") if c]
        scope = f" across {len(cancers)} cancer types" if cancers else ""
        if role:
            bullets.append(f"**Cancer driver (intOGen):** {role}{scope}")

    # Cancer dependency (DepMap CRISPR fitness) — only when notable. A low
    # pct_dependent (e.g. 0.2%) means "not a dependency" — noise as a headline,
    # so we suppress it unless common-essential, strongly-selective, or ≥10%.
    dm = b3.get("depmap") or {}
    pct = dm.get("pct_dependent")
    try:
        notable = float(pct) >= 10
    except (TypeError, ValueError):
        notable = False
    notable = (notable or dm.get("common_essential") == "true"
               or dm.get("strongly_selective") == "true")
    if pct not in (None, "") and notable:
        tags = []
        if dm.get("common_essential") == "true":
            tags.append("common-essential")
        if dm.get("strongly_selective") == "true":
            tags.append("strongly selective")
        tag = f" ({', '.join(tags)})" if tags else ""
        bullets.append(
            f"**Cancer dependency (DepMap):** dependent in {pct}% of "
            f"screened cell lines{tag}")

    # Dosage sensitivity (ClinGen).
    cd = b3.get("clingen_dosage") or {}
    if cd:
        haplo = cd.get("haplo_score", "")
        triplo = cd.get("triplo_score", "")
        bullets.append(
            f"**Dosage sensitivity (ClinGen):** haploinsufficiency "
            f"{_DOSAGE_SCALE.get(haplo, 'unscored')}, triplosensitivity "
            f"{_DOSAGE_SCALE.get(triplo, 'unscored')}")

    # Transcription factor flag (CollecTRI) — promoted only when it has a
    # meaningful regulon. A 1–2 target flag (e.g. TTN, APOE) is a source
    # false-positive and stays in §9 Regulation, not the headline set.
    if _truthy(b9.get("is_transcription_factor")):
        n = b9.get("downstream_count") or 0
        if n >= _TF_MIN_TARGETS:
            bullets.append(
                f"**Transcription factor:** yes — {n:,} downstream targets (CollecTRI)")

    # MANE Select transcript — the reference transcript for clinical reporting.
    mane = b2.get("mane_select_refseq")
    if mane and mane != "None":
        bullets.append(f"**MANE Select transcript:** `{mane}`")

    # Interaction hub (STRING) — surfaced ONLY when corpus-relatively notable: a
    # raw partner count means little without the distribution, so the rank IS the
    # reason for the bullet. Elides for ordinary genes.
    inter = (bundle.get("8") or {}).get("string_count") or 0
    hub = evidence.rank_clause("gene", "interaction_count", inter)
    if hub:
        bullets.append(f"**Interaction hub:** {inter:,} STRING partners{hub}")

    # Notable callout — deterministic anomaly: heavily sequenced (many ClinVar
    # variants) yet no curated precision-oncology actionability. A "studied but
    # not yet actionable" observation no single source states.
    if cv_total >= 50 and civ_total == 0:
        bullets.append(f"**Notable:** {cv_total:,} clinical variants but no curated "
                       f"precision-oncology (CIViC) evidence yet")

    if not bullets:
        return ""
    # Bold label, not a "## " header — this is an intro block, not a peer
    # data section.
    return "**At a glance**\n\n" + "\n".join(f"- {b}" for b in bullets)
