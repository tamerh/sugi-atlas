"""§7 — expression context: per-cohort-gene tissue + single-cell expression
(Bgee, FANTOM5, SCXA). REUSE wrapper over gene §11.

Fans gene §11 over the cohort, then summarizes:
 - per-gene breadth/presence + top tissues,
 - breadth distribution (narrow/moderate/broad/unknown),
 - cohort-wide tissue counts (how many genes mention each tissue),
 - count of "low confidence" genes (no bgee + no fantom5 + no scxa).
"""
from collections import Counter

from atlas.section import Section
from atlas.disease.cohort import fan
from atlas.gene.sections import s11_expression

CHAINS   = (">>hgnc>>ensembl>>bgee", ">>hgnc>>ensembl>>scxa")  # via gene §11
DATASETS = ("hgnc", "ensembl", "bgee", "scxa")

_TOP_TISSUES_PER_GENE = 3
_COHORT_TOP_TISSUES = 20


def _breadth_bucket(b):
    """Bin a bgee breadth integer into the four reporting buckets."""
    if b is None:
        return "unknown"
    try:
        n = int(b)
    except (TypeError, ValueError):
        return "unknown"
    if n <= 5:
        return "narrow (1-5 tissues)"
    if n <= 20:
        return "moderate (6-20)"
    return "broad (>20)"


def collect(a):
    g11_bundles = fan(s11_expression.SECTION.collect_fn, a.cohort)

    per_gene = []
    breadth_distribution = {
        "narrow (1-5 tissues)": 0,
        "moderate (6-20)": 0,
        "broad (>20)": 0,
        "unknown": 0,
    }
    cohort_tissue_counts: Counter = Counter()
    no_expression_count = 0

    for b in g11_bundles:
        bgee = b.get("bgee") or {}
        fantom5 = b.get("fantom5") or {}
        scxa = b.get("single_cell_datasets") or []
        top_t_raw = b.get("top_tissues") or []

        # bgee exposes both a categorical label ("narrow"/"broad"/"ubiquitous")
        # and a per-tissue present_calls integer. The integer is what the
        # downstream "1-5 / 6-20 / >20" buckets need, so we surface it as
        # bgee_breadth; the categorical label is left to gene §11 consumers.
        bgee_breadth = None
        if bgee:
            try:
                bgee_breadth = int(bgee.get("present_calls"))
            except (TypeError, ValueError):
                bgee_breadth = None
        fantom5_breadth = fantom5.get("breadth") if fantom5 else None
        scxa_present = bool(scxa)

        # Top tissues by score — gene §11 already sorted bgee_evidence by
        # expression_score desc. Fall back to empty list if none.
        top_tissues = [t.get("tissue") for t in top_t_raw[:_TOP_TISSUES_PER_GENE]
                       if t.get("tissue")]

        per_gene.append({
            "symbol": b.get("symbol"),
            "hgnc_id": b.get("hgnc_id"),
            "bgee_breadth": bgee_breadth,
            "fantom5_breadth": fantom5_breadth,
            "scxa_present": scxa_present,
            "top_tissues": top_tissues,
        })

        breadth_distribution[_breadth_bucket(bgee_breadth)] += 1

        # Per-gene set so a gene contributes at most once per tissue name.
        # sorted() makes Counter insertion order deterministic → most_common()
        # ties break stably across runs (set() iteration order otherwise varies
        # by process hash-seed, churning rendered row order).
        for tissue in sorted(set(top_tissues)):
            cohort_tissue_counts[tissue] += 1

        if not bgee and not fantom5 and not scxa_present:
            no_expression_count += 1

    top_cohort_tissues = Counter(dict(cohort_tissue_counts.most_common(_COHORT_TOP_TISSUES)))

    return {
        "section": "07_expression_context",
        "mondo_id": a.mondo_id,
        "per_gene_expression": per_gene,
        "breadth_distribution": breadth_distribution,
        "cohort_tissue_counts": top_cohort_tissues,
        "no_expression_count": no_expression_count,
    }


SECTION = Section(
    id="7", name="expression_context",
    description=("Per-cohort-gene Bgee tissue expression + FANTOM5 + SCXA "
                 "single-cell. Tissue-specificity ranking; flag genes with "
                 "no expression evidence as low confidence."),
    needs=("cohort",),
    produces=("per_gene_expression", "breadth_distribution",
              "cohort_tissue_counts", "no_expression_count"),
    datasets=DATASETS, chains=CHAINS, collect_fn=collect,
)
