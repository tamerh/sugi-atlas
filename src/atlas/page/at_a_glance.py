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

# ClinGen dosage score scale (same mapping the §-detail block uses).
_DOSAGE_SCALE = {
    "3": "sufficient evidence", "2": "emerging evidence",
    "1": "little evidence", "0": "no evidence",
    "30": "autosomal recessive", "40": "dosage sensitivity unlikely",
}


def _truthy(v):
    return v in (True, "true", "True", "1", 1, "yes", "Yes")


def at_a_glance(bundle) -> str:
    """Compose the `## At a glance` markdown block from a full bundle dict.

    Sources: §2 (MANE), §3 (DepMap + ClinGen dosage), §9 (TF flag),
    §10 (drug-target flag + CIViC predictive total), §12 (intOGen driver)."""
    b2 = bundle.get("2") or {}
    b3 = bundle.get("3") or {}
    b9 = bundle.get("9") or {}
    b10 = bundle.get("10") or {}
    b12 = bundle.get("12") or {}

    bullets = []

    # Druggable target — drug-discovery headline flag.
    if _truthy(b10.get("is_drug_target")):
        mc = b10.get("molecule_count") or 0
        extra = f" — {mc:,} molecules with ChEMBL bioactivity" if mc else ""
        bullets.append(f"**Druggable target:** yes{extra}")

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

    # Cancer dependency (DepMap CRISPR fitness).
    dm = b3.get("depmap") or {}
    pct = dm.get("pct_dependent")
    if pct not in (None, ""):
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

    # Transcription factor flag (CollecTRI).
    if _truthy(b9.get("is_transcription_factor")):
        n = b9.get("downstream_count") or 0
        extra = (f" — {n:,} downstream target{'s' if n != 1 else ''} (CollecTRI)"
                 if n else "")
        bullets.append(f"**Transcription factor:** yes{extra}")

    # MANE Select transcript — the reference transcript for clinical reporting.
    mane = b2.get("mane_select_refseq")
    if mane and mane != "None":
        bullets.append(f"**MANE Select transcript:** `{mane}`")

    if not bullets:
        return ""
    return "## At a glance\n\n" + "\n".join(f"- {b}" for b in bullets)
