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
from atlas.page import evidence


def _format_int(n):
    try:
        return f"{int(n):,}"
    except (TypeError, ValueError):
        return None


def _status(b1, b5=None) -> str:
    """Development status from max_phase / FDA flag / Phase-4 trial signal."""
    from atlas.indication import has_phase4_trial
    if b1.get("is_fda_approved") or b1.get("max_phase") == 4:
        return "Approved (max clinical phase 4)"
    mp = b1.get("max_phase")
    try:
        mp = int(float(mp))
    except (TypeError, ValueError):
        return ""
    # ChEMBL under-phases some approved drugs (e.g. non-oncology oligonucleotides
    # like inclisiran); a registered Phase-4 trial only exists post-approval.
    if mp >= 3 and has_phase4_trial(b5):
        return f"Approved (registered Phase 4 trials; ChEMBL max clinical phase {mp})"
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
    status = _status(b1, b5)
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
        # rank on the broader ChEMBL bioactivity-target breadth (the evidence
        # component), which is what the corpus distribution is built from.
        tgt = b2.get("bioactivity_target_count") or len(genes)
        bullets.append(f"**Targets:** {len(genes)}{extra}"
                       + evidence.rank_clause("drug", "target_count", tgt))

    # Indications (§4).
    ind = b4.get("indication_count") or 0
    if ind:
        bullets.append(f"**Indications:** {_format_int(ind)} condition"
                       f"{'s' if ind != 1 else ''}"
                       + evidence.rank_clause("drug", "indication_count", ind))

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

    # Cross-entity reach over the curated (GtoPdb) targets — how broadly the
    # drug's molecular targets touch disease genetics vs where it's actually used,
    # and how crowded each target is. Reverse-index reads; elide off-corpus.
    from atlas.page import links
    land = links.drug_target_landscape(bundle)
    reach = land["reach"]
    if reach:
        ind_urls = set()
        for i in (b4.get("indications") or []):
            u = links.disease_url(mondo_id=i.get("mondo_id"), name=i.get("name"))
            if u:
                ind_urls.add(u)
        nt = land["n_targets"]
        tail = ""
        if ind_urls:
            overlap = len(reach & ind_urls)
            tail = (f"; indicated or in trials for {len(ind_urls)} "
                    f"({overlap} on both)")
        bullets.append(
            f"**Target biology:** {nt} curated target{'s' if nt != 1 else ''} "
            f"implicated as cohort gene{'s' if nt != 1 else ''} in {len(reach):,} "
            f"disease{'s' if len(reach) != 1 else ''} corpus-wide{tail}")
        comp = sorted(((s, nd - 1) for s, _d, nd in land["per_target"] if nd > 1),
                      key=lambda x: -x[1])
        if comp:
            lead_s, lead_n = comp[0]
            extra = "; ".join(f"{s}: {n}" for s, n in comp[1:4])
            more = f" ({extra} others)" if extra else ""
            bullets.append(
                f"**Competitive landscape:** target {lead_s} is also targeted by "
                f"{lead_n} other curated drug{'s' if lead_n != 1 else ''} in the "
                f"corpus{more}")

    # Notable callout — substantial trial activity without approval (a drug in
    # active development, surfaced as a deterministic observation).
    from atlas.indication import molecule_approved
    if tr >= 10 and not molecule_approved(b1, b5):
        bullets.append(f"**Notable:** {_format_int(tr)} clinical trials, not yet approved")

    if not bullets:
        return ""
    return "**At a glance**\n\n" + "\n".join(f"- {b}" for b in bullets)
