"""§13 — Human Protein Atlas. Gene/protein-level summary (subcellular location,
protein classes, RNA tissue/cancer specificity, protein-evidence tier), per-
tissue/cell expression (nTPM + IHC protein level), and cancer pathology
prognostics. Routed via >>hgnc>>hpa and its child datasets. Empty for non-coding
genes (HPA is protein-focused) — the renderers elide.

This is a SHARED-bundle section: it's collected once, and pipeline.render_all
threads its bundle into derived renderers that slot HPA facts into the Protein
(subcellular/classes), Gene-structure (expression), and Disease (cancer
prognostics) zones — so HPA isn't fetched three times."""
import re

from atlas.biobtree import map_all, entry
from atlas.gene.sections.base import Section

CHAINS = (">>hgnc>>hpa", ">>hgnc>>hpa>>hpa_expression", ">>hgnc>>hpa>>hpa_pathology")
DATASETS = ("hpa", "hpa_expression", "hpa_pathology")

# HPA dumps ~50 protein_classes/gene, most of them noise: per-method membrane/
# secreted topology predictions (MEMSAT3/Phobius/SignalP/… — the SAME call from
# a dozen tools), COSMIC mutation-tally subclasses, and protein-evidence labels
# (we already surface the evidence tier as its own field). Drop those; keep the
# real functional/annotation classes (Kinases, Transporters, Cancer-related, …).
_CLASS_NOISE = re.compile(
    r"predicted by|MEMSAT|SPOCTOPUS|Phobius|SCAMPI|THUMBUP|TMHMM|DeepSig|SignalP"
    r"|^COSMIC |Evidence at protein level|neXtProt|Protein evidence \("
    # also drop HPA's verbose cancer-/disease-AREA groupings ("Cancers of the
    # digestive system", "Head and neck cancers", "Immune system diseases") — the
    # gene's actual cancer/disease associations live in the Disease zone. Keep the
    # functional + gene-flag classes (Kinases, FDA approved drug targets, …).
    r"|^Cancers$|^Cancers of|cancers$|diseases$", re.I)


def _num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return 0.0


def collect(a):
    bundle = {"section": "13_hpa", "symbol": a.symbol}
    hits = map_all(a.symbol, ">>hgnc>>hpa")
    if not hits:
        return bundle                          # no HPA record (non-coding / unmapped)
    hid = hits[0].get("id")
    try:                                       # full attrs live on the entry node
        at = (entry(hid, "hpa").get("Attributes") or {}).get("Hpa") or {}
    except Exception:
        at = hits[0]
    bundle["hpa"] = {
        "protein_classes": [c for c in (at.get("protein_classes") or [])
                            if c and not _CLASS_NOISE.search(c)],
        "protein_evidence": at.get("protein_evidence"),
        "subcellular_main": at.get("subcellular_main") or [],
        "subcellular_additional": at.get("subcellular_additional") or [],
        "rna_tissue_specificity": at.get("rna_tissue_specificity"),
        "rna_tissue_distribution": at.get("rna_tissue_distribution"),
        "rna_cancer_specificity": at.get("rna_cancer_specificity"),
        "secretome_location": at.get("secretome_location"),
        "top_tissues": at.get("top_tissues") or [],   # ["Placenta|61.8", ...]
    }
    # Per-tissue/cell expression — nTPM, highest first; IHC protein_level when present.
    exp = [{"entity": r.get("entity_name"), "axis": r.get("axis"),
            "ntpm": r.get("ntpm"), "protein_level": r.get("protein_level"),
            "reliability": r.get("reliability")}
           for r in map_all(a.symbol, ">>hgnc>>hpa>>hpa_expression")]
    exp.sort(key=lambda r: -_num(r.get("ntpm")))
    bundle["hpa_expression"] = exp
    bundle["hpa_expression_total"] = len(exp)
    # Cancer pathology — keep the PROGNOSTIC cancers (favorable/unfavorable + p);
    # the non-prognostic rows are just "expression detected", no signal.
    bundle["hpa_pathology"] = [
        {"cancer": r.get("cancer"), "prognostic_type": r.get("prognostic_type"),
         "p_value": r.get("p_value")}
        for r in map_all(a.symbol, ">>hgnc>>hpa>>hpa_pathology")
        if r.get("is_prognostic") == "true"]
    return bundle


SECTION = Section(
    id="13", name="hpa",
    description=("Human Protein Atlas: subcellular location, protein classes, RNA "
                 "tissue/cancer specificity, per-tissue expression (nTPM), and cancer "
                 "prognostic markers, via >>hgnc>>hpa (+ hpa_expression, hpa_pathology)"),
    needs=("symbol",),
    produces=("hpa", "hpa_expression", "hpa_expression_total", "hpa_pathology"),
    datasets=DATASETS,
    chains=CHAINS,
    collect_fn=collect,
)
