"""Deterministic "At a glance" flag-sheet for Atlas disease pages.

Disease analog of atlas.page.at_a_glance (genes). A scannable bullet digest of
the disease's headline numbers — prevalence, cohort size, GWAS / ClinVar /
phenotype counts, trial activity, precision-medicine evidence. Pure dict
lookups over the collector bundle: no LLM, no network, no hallucination
surface.

Sits directly under the declarative lead as a bold "**At a glance**" intro
block (not a "## " section). Bullets with no backing data drop; the whole
block elides if nothing qualifies.
"""


def _format_int(n):
    try:
        return f"{int(n):,}"
    except (TypeError, ValueError):
        return None


def _prevalence(b1) -> str:
    """Best Orphanet prevalence row — widest, validated, point-prevalence first
    (same ranking as the JSON-LD epidemiology field). '' for diseases with no
    Orphanet epidemiology (most cancers / common conditions)."""
    prevs = b1.get("prevalences") or []
    if not prevs:
        return ""
    def _rank(p):
        geo = p.get("geographic") or ""
        return (0 if p.get("validation_status") == "Validated" else 1,
                0 if p.get("prevalence_type") == "Point prevalence" else 1,
                0 if geo == "Worldwide" else (1 if geo in ("Europe", "Americas") else 2))
    best = sorted(prevs, key=_rank)[0]
    parts = [best.get("prevalence_class") or ""]
    if best.get("geographic"):
        parts.append(f"({best['geographic']})")
    if best.get("validation_status") == "Validated":
        parts.append("[Orphanet-validated]")
    label = " ".join(p for p in parts if p).strip()
    return f"**Prevalence:** {label}" if label else ""


def at_a_glance(bundle) -> str:
    """Compose the `**At a glance**` block from a full disease bundle dict.

    Sources: §1 (class + prevalence + HPO + xref_counts), §2 (GWAS),
    §3 (ClinVar), §5 (cohort genes), §13 (trials + CIViC)."""
    b1  = bundle.get("1")  or {}
    b2  = bundle.get("2")  or {}
    b5  = bundle.get("5")  or {}
    b13 = bundle.get("13") or {}
    xc  = b1.get("xref_counts") or {}

    bullets = []

    # Classification — only the informative "cancer" flag. (Orphanet's
    # disorder_type is generic filler — "Disease"/"Category" — so it's skipped.)
    if b1.get("is_cancer"):
        bullets.append("**Classification:** Cancer")

    # Prevalence (Orphanet) — the clinical headline for rare disease.
    prev = _prevalence(b1)
    if prev:
        bullets.append(prev)

    # Cohort genes (§5) — size of the associated-gene set this page aggregates.
    gc = b5.get("gene_count") or 0
    if gc:
        bullets.append(f"**Cohort genes:** {_format_int(gc)}")

    # GWAS associations (§2) — common-variant genetics.
    ga = b2.get("assoc_total") or 0
    if ga:
        bullets.append(f"**GWAS associations:** {_format_int(ga)}")

    # ClinVar variants — accurate xref total from the Mondo entry (§1
    # xref_counts), NOT §3's clinvar_total which is a paginated fetch floor
    # (caps at ~600). Same accurate count the declarative lead uses.
    cv = xc.get("clinvar") or 0
    if cv:
        bullets.append(f"**ClinVar variants:** {_format_int(cv)}")

    # HPO phenotypes (§1) — clinical feature count.
    hp = b1.get("phenotype_count") or 0
    if hp:
        bullets.append(f"**Phenotypes (HPO):** {_format_int(hp)}")

    # Clinical trials — the VALIDATED §13 count (title-matched), NOT the raw
    # xref count (biobtree's mondo→clinical_trials edge is contaminated). 0 is a
    # valid answer for a rare disease, so don't fall back to the raw count.
    tr = b13.get("trial_count") or 0
    if tr:
        bullets.append(f"**Clinical trials:** {_format_int(tr)}")

    # Precision-medicine evidence (CIViC, §13) — subtype–drug associations.
    civ = b13.get("civic_association_total") or 0
    if civ:
        bullets.append(
            f"**Precision-medicine evidence (CIViC):** {_format_int(civ)} "
            f"subtype–drug association{'s' if civ != 1 else ''}")

    if not bullets:
        return ""
    return "**At a glance**\n\n" + "\n".join(f"- {b}" for b in bullets)
