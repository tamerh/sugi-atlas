"""Deterministic "At a glance" flag-sheet for Atlas drug pages.

Drug analog of atlas.page.at_a_glance (genes). A scannable bullet digest of the
drug's headline facts — development status, modality, ATC class, target /
indication / trial counts, precision-oncology evidence, chemistry. Pure dict
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


def _status(b1) -> str:
    """Development status from max_phase / FDA flag."""
    if b1.get("is_fda_approved") or b1.get("max_phase") == 4:
        return "Approved (max clinical phase 4)"
    mp = b1.get("max_phase")
    try:
        mp = int(float(mp))
    except (TypeError, ValueError):
        return ""
    return f"Max clinical phase {mp} (not approved)" if mp else ""


def at_a_glance(bundle) -> str:
    """Compose the `**At a glance**` block from a full drug bundle dict.

    Sources: §1 (status + modality + ATC + chemistry), §2 (targets),
    §4 (indications), §5 (trials), §10 (CIViC)."""
    b1  = bundle.get("1")  or {}
    b2  = bundle.get("2")  or {}
    b4  = bundle.get("4")  or {}
    b5  = bundle.get("5")  or {}
    b10 = bundle.get("10") or {}

    bullets = []

    # Development status.
    status = _status(b1)
    if status:
        bullets.append(f"**Status:** {status}")

    # Modality (molecule type).
    mtype = (b1.get("molecule_type") or "").strip()
    if mtype:
        bullets.append(f"**Modality:** {mtype}")

    # ATC therapeutic class.
    atc = b1.get("atc_codes") or []
    if atc:
        more = f" (+{len(atc) - 1} more)" if len(atc) > 1 else ""
        bullets.append(f"**ATC class:** {atc[0]}{more}")

    # Targets (§2 primary targets) — count + top symbols.
    prim = b2.get("primary_targets") or []
    genes = [t.get("gene_symbol") for t in prim if t.get("gene_symbol")]
    if genes:
        shown = ", ".join(genes[:3])
        extra = f" ({shown}…)" if len(genes) > 3 else f" ({shown})"
        bullets.append(f"**Targets:** {len(genes)}{extra}")

    # Indications (§4).
    ind = b4.get("indication_count") or 0
    if ind:
        bullets.append(f"**Indications:** {_format_int(ind)} condition"
                       f"{'s' if ind != 1 else ''}")

    # Clinical trials (§5).
    tr = b5.get("trial_count") or 0
    if tr:
        bullets.append(f"**Clinical trials:** {_format_int(tr)}")

    # Precision-oncology evidence (CIViC, §10) — variant–indication associations.
    civ = b10.get("civic_association_total") or 0
    if civ:
        bullets.append(
            f"**Precision-oncology evidence (CIViC):** {_format_int(civ)} "
            f"variant–indication association{'s' if civ != 1 else ''}")

    # Chemistry (§1) — molecular weight + formula (small molecules).
    mw = b1.get("molecular_weight")
    formula = b1.get("molecular_formula")
    if mw or formula:
        bits = []
        if mw:
            bits.append(f"{mw} Da")
        if formula:
            bits.append(formula)
        bullets.append(f"**Chemistry:** {' · '.join(bits)}")

    if not bullets:
        return ""
    return "**At a glance**\n\n" + "\n".join(f"- {b}" for b in bullets)
