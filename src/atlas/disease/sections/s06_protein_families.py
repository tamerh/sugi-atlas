"""§6 — protein-family classification: cohort proteins classified by druggable
families (Kinase / GPCR / Ion channel / Nuclear receptor / Protease / TF / ...)
via InterPro. REUSE wrapper over gene §3 + custom classifier.

Druggable vs difficult split drives §16 (druggability pyramid)."""
from collections import Counter

from atlas.section import Section
from atlas.disease.cohort import fan
from atlas.gene.sections import s03_protein_ids
from atlas.ora import enrich, background

CHAINS   = (">>uniprot>>interpro",)  # reused via gene §3
DATASETS = ("uniprot", "interpro", "pfam")

# Families that are pharmacologically tractable (small molecule or biologic).
_DRUGGABLE = {
    "Kinase", "GPCR", "Ion channel", "Nuclear receptor", "Protease",
    "Phosphatase", "Transporter", "Enzyme (other)", "Complement",
    "Antibody/Immunoglobulin",
}
_DIFFICULT = {"Transcription factor", "Scaffold/PPI"}
_UNKNOWN = {"Other/Unknown"}


def _classify(interpro_names: list[str], has_ec: bool) -> str:
    """Apply the keyword cascade — first hit wins. See module docstring for
    rationale; the order matters (kinase before protease before TF).

    InterPro short_names use compressed CamelCase tokens joined by `_`/`/`/`-`
    (e.g. `Prot_kinase_dom`, `K_chnl_dom`, `Znf_C2H2_type`, `Nucl_hrmn_rcpt`).
    We normalise separators to spaces so substring checks land on tokens; we
    also include common biobtree abbreviations (chnl/hrmn/rcpt/znf) so the
    keyword list covers what actually arrives from biobtree."""
    blob = " ".join(n for n in interpro_names if n).lower()
    blob = blob.replace("_", " ").replace("/", " ").replace("-", " ")
    # pad so naked-token checks (e.g. " hd ", " tf ") get word boundaries on
    # the edges of the string too.
    blob = f" {blob} "

    def has(*kws):
        return any(k in blob for k in kws)

    if has("kinase"):
        return "Kinase"
    if has("gpcr", "g protein coupled", "rhodopsin", "7tm "):
        return "GPCR"
    if has("ion channel", "ion transport", "ionic channel",
           " ion trans ", " k chnl ", " na chnl ", " ca chnl ", " cl chnl ",
           " trpc channel", "trp dom", "trpm ", "trpv ", "trpa ",
           " chnl dom"):
        return "Ion channel"
    if has("nuclear receptor", "nucl hrmn rcpt", "nuclear hrmn rcpt",
           "znf hrmn rcpt", "steroid receptor"):
        return "Nuclear receptor"
    if has("protease", "peptidase", "trypsin dom", "merops"):
        return "Protease"
    if has("phosphatase"):
        return "Phosphatase"
    if has("transporter", "solute carrier", "abc transporter",
           " abc tran ", " mfs ", " mfs dom", " mct ", " slc ",
           " mfs trans"):
        return "Transporter"
    # TF: explicit phrase + InterPro shorthands (Znf_*, HD, homeobox, bHLH,
    # tumour_suppressor for p53). The `_TF_` shorthand shows up inside compound
    # names like `p53 like tf dna bd` after normalisation.
    if has("transcription factor", " tf dna ", " tf bd ", " tf ",
           "zinc finger", " znf ", "homeobox", "homeodomain", " hd ",
           "bhlh", "tumour suppressor", "tumor suppressor",
           "p53 ", "runt "):
        return "Transcription factor"
    if has("complement", "factor h", " c3 ", " c5 ", " c9 ",
           "sushi scr ccp", "anaphylatoxin", "macpf"):
        return "Complement"
    if has("immunoglobulin", "antibody", " ig like ", " ig dom",
           " ig sub "):
        return "Antibody/Immunoglobulin"
    if has("scaffold", " ww dom", " pdz", " sh3 ", " sh2 ",
           " ph domain", " ph dom", " bar dom", " ankyrin",
           " arm rpt", " wd40"):
        return "Scaffold/PPI"
    if has_ec:
        return "Enzyme (other)"
    return "Other/Unknown"


def collect(a):
    g3_bundles = fan(s03_protein_ids.SECTION.collect_fn, a.cohort)

    assignments = []
    for b in g3_bundles:
        interpro = b.get("interpro") or []
        names = [(d.get("name") or "") for d in interpro]
        pfam = list(b.get("pfam") or [])
        brenda = b.get("brenda_ec") or []
        has_ec = bool(brenda)
        ec = brenda[0].get("ec") if brenda else None

        family = _classify(names, has_ec)
        druggable = family in _DRUGGABLE

        assignments.append({
            "symbol": b.get("symbol"),
            "hgnc_id": b.get("hgnc_id"),
            "canonical_uniprot": b.get("canonical_uniprot"),
            "interpro_names": names,
            "pfam_ids": pfam,
            "assigned_family": family,
            "druggable": druggable,
            "ec": ec,
        })

    family_counts = Counter(x["assigned_family"] for x in assignments)
    druggable_count = sum(1 for x in assignments if x["druggable"])
    unknown_count = sum(1 for x in assignments
                        if x["assigned_family"] in _UNKNOWN)
    difficult_count = len(assignments) - druggable_count - unknown_count
    total = len(assignments)
    druggable_fraction = round(druggable_count / total, 2) if total else 0.0

    # Over-representation vs a genome-wide family background (same _classify over
    # all protein-coding genes): surfaces that the cohort is enriched for e.g.
    # Kinases, and demotes the uninformative catch-all Other/Unknown bucket that
    # otherwise leads a raw-count table. Counts are kept; fold + FDR added.
    fam_universe_n, fam_sizes = background("family")
    family_enrichment = {it["id"]: {"fold": it["fold"], "fdr": it["fdr"], "K": it["K"]}
                         for it in enrich(
                             [{"id": fam, "k": cnt, "K": fam_sizes.get(fam, 0)}
                              for fam, cnt in family_counts.items()],
                             cohort_n=total, universe_n=fam_universe_n)}

    return {
        "section": "06_protein_families",
        "mondo_id": a.mondo_id,
        "family_assignments": assignments,
        "family_counts": family_counts,
        "family_enrichment": family_enrichment,
        "druggable_count": druggable_count,
        "difficult_count": difficult_count,
        "unknown_count": unknown_count,
        "druggable_fraction": druggable_fraction,
    }


SECTION = Section(
    id="6", name="protein_families",
    description=("Cohort proteins classified by druggable family (Kinase, GPCR, "
                 "Ion channel, Nuclear receptor, Protease, Phosphatase, Enzyme, "
                 "TF, Scaffold) via InterPro. Druggable vs difficult split."),
    needs=("cohort",),
    produces=("family_assignments", "family_counts", "family_enrichment",
              "druggable_count", "difficult_count"),
    datasets=DATASETS, chains=CHAINS, collect_fn=collect,
)
