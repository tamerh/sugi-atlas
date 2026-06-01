"""§3 — variant details: top GWAS variants → dbSNP attributes + tier
classification (Tier 1 coding / Tier 2 splice-UTR / Tier 3 regulatory /
Tier 4 intronic-intergenic)."""
import re, html
from atlas.section import Section
from atlas.biobtree import map_all, entry

CHAINS   = (">>mondo>>gwas>>dbsnp", ">>mondo>>clinvar")
DATASETS = ("mondo", "gwas", "dbsnp", "clinvar")

# ClinVar germline-classification ranking — clinically actionable first.
# Order matters: more-specific phrases first (substring match, first hit wins),
# so "conflicting … of pathogenicity" doesn't get caught by the "pathogenic" key.
_CV_RANK = [("conflicting", 4),
            ("pathogenic/likely pathogenic", 0), ("likely pathogenic", 1),
            ("pathogenic", 0), ("risk", 2), ("drug response", 2),
            ("uncertain", 5), ("likely benign", 6), ("benign", 6),
            ("not provided", 7)]


def _cv_rank(classification):
    c = (classification or "").lower()
    for key, r in _CV_RANK:
        if key in c:
            return r
    return 5

TIER1 = {"missense_variant", "missense", "frameshift_variant", "frameshift",
         "stop_gained", "nonsense", "stop_lost", "start_lost",
         "inframe_insertion", "inframe_deletion"}
TIER2_EXACT = {"5_prime_UTR_variant", "3_prime_UTR_variant"}
TIER3_EXACT = {"regulatory_region_variant", "TF_binding_site_variant"}

T1 = "Tier 1: coding"
T2 = "Tier 2: splice/UTR"
T3 = "Tier 3: regulatory"
T4 = "Tier 4: intronic/intergenic"

def _classify(consequence):
    if not consequence:
        return T4
    c = consequence.strip()
    if c in TIER1:
        return T1
    if "splice" in c or c in TIER2_EXACT:
        return T2
    if c in TIER3_EXACT or "regulatory" in c:
        return T3
    return T4

_MAF_RX = re.compile(r"([0-9]*\.?[0-9]+)")

def _parse_maf(raf):
    """Parse risk_allele_frequency like '0.70 (HapMap CEU)' -> 0.70.
    Returns None if not parseable. MAF = min(raf, 1-raf)."""
    if not raf or raf in ("NR", "NA", "."):
        return None
    m = _MAF_RX.search(str(raf))
    if not m:
        return None
    try:
        v = float(m.group(1))
    except ValueError:
        return None
    if v < 0 or v > 1:
        return None
    return min(v, 1.0 - v)

def _maf_bucket(maf):
    if maf is None:
        return "unknown"
    if maf >= 0.05:
        return "common (>=0.05)"
    if maf >= 0.01:
        return "low_freq (0.01-0.05)"
    return "rare (<0.01)"

def collect(a):
    # 1) All GWAS associations for the disease (already cheap, paginated).
    gwas_rows = map_all(a.mondo_id, ">>mondo>>gwas", cap=10)

    # 2) Fetch per-GWAS entry to get p_value, rsid, consequence (context),
    #    mapped_gene, risk_allele_frequency. The map-row only has chr/p_value
    #    string + mapped_gene, no rsid or context — entry() is the source.
    #    Bound to 50 calls: sort the rows by p_value string ascending (smaller
    #    is more significant). The string form '4.000000e-08' sorts
    #    lexicographically the same way as float for typical GWAS p-values
    #    (all 'NeXX' with same width), but parse to float to be safe.
    def _pv(r):
        try:
            return float(r.get("p_value") or "1")
        except ValueError:
            return 1.0
    gwas_rows.sort(key=_pv)

    # GWAS rows aren't unique by rsid (one rsid can recur across studies). We
    # take the smallest-p_value row per rsid, then keep TOP 50 unique rsids.
    # Sort first ensures the iteration below picks the best p-value per rsid.
    # We still cap entry() calls below via the dbsnp_budget.
    top = gwas_rows  # entry() iteration is bounded by dbsnp_budget + break

    top_variants = []
    tier_counts = {T1: 0, T2: 0, T3: 0, T4: 0}
    maf_distribution = {"common (>=0.05)": 0, "low_freq (0.01-0.05)": 0,
                         "rare (<0.01)": 0, "unknown": 0}
    consequence_distribution = {}

    seen_rs = set()
    dbsnp_budget = 50
    UNIQUE_CAP = 50
    for g in top:
        if len(top_variants) >= UNIQUE_CAP:
            break
        gid = g.get("id")
        if not gid:
            continue
        try:
            ge = entry(gid, "gwas")
        except Exception:
            continue
        gattrs = (ge.get("Attributes") or {}).get("Gwas") or {}
        rsid = gattrs.get("snp_id") or ""
        if not rsid or rsid in seen_rs:
            continue
        seen_rs.add(rsid)

        consequence = gattrs.get("context") or ""
        gene_symbol = gattrs.get("mapped_gene") or ""
        pvalue = gattrs.get("p_value")
        chrom = str(gattrs.get("chr_id") or "")
        pos = gattrs.get("chr_pos")
        raf = gattrs.get("risk_allele_frequency")
        maf = _parse_maf(raf)

        # dbSNP entry: fills alleles + corroborates chrom/pos. Bounded.
        alleles = ""
        if dbsnp_budget > 0:
            try:
                de = entry(rsid, "dbsnp")
                dattrs = (de.get("Attributes") or {}).get("Dbsnp") or {}
                ref = dattrs.get("ref_allele") or ""
                alt = dattrs.get("alt_allele") or ""
                if ref or alt:
                    alleles = f"{ref}>{alt}" if ref and alt else (ref or alt)
                if not chrom:
                    chrom = str(dattrs.get("chromosome") or "")
                if not pos:
                    pos = dattrs.get("position")
                # dbsnp has is_common as a coarse MAF hint (>=0.01 typically),
                # but only adopt when GWAS lacked frequency info.
                if maf is None and dattrs.get("is_common"):
                    maf = 0.05  # sentinel mapping to "common" bucket
                dbsnp_budget -= 1
            except Exception:
                pass

        tier = _classify(consequence)
        tier_counts[tier] += 1
        maf_distribution[_maf_bucket(maf)] += 1
        ckey = consequence or "unknown"
        consequence_distribution[ckey] = consequence_distribution.get(ckey, 0) + 1

        top_variants.append({
            "rsid": rsid,
            "chrom": chrom,
            "pos": pos,
            "alleles": alleles,
            "maf": maf,
            "consequence": consequence,
            "gene_symbol": gene_symbol,
            "pvalue": pvalue,
            "tier": tier,
        })

    # ── ClinVar germline variants — the rare/coding genetic evidence (the
    # Mendelian counterpart to the GWAS tiering above). Many diseases have no
    # GWAS but a rich ClinVar set (e.g. ataxia-telangiectasia → 100+ ATM
    # variants). Classification distribution + top variants, clinically
    # actionable (Pathogenic/Likely pathogenic) first.
    cv_class_counts = {}
    cv_variants, seen_cv = [], set()
    for r in map_all(a.mondo_id, ">>mondo>>clinvar", cap=5):
        cid = r.get("id")
        if not cid or cid in seen_cv:
            continue
        seen_cv.add(cid)
        cls = (r.get("germline_classification") or "not provided").strip()
        cv_class_counts[cls] = cv_class_counts.get(cls, 0) + 1
        cv_variants.append({
            "id": cid,
            "hgvs": html.unescape(r.get("name") or ""),   # biobtree HTML-escapes '>'
            "gene": r.get("gene_symbol"),
            "classification": cls,
            "review_status": r.get("review_status"),
        })
    cv_variants.sort(key=lambda v: (_cv_rank(v["classification"]), v["gene"] or ""))

    return {
        "section": "03_variant_details",
        "mondo_id": a.mondo_id,
        "top_variants": top_variants,
        "tier_counts": tier_counts,
        "maf_distribution": maf_distribution,
        "consequence_distribution": consequence_distribution,
        "clinvar_total": len(seen_cv),
        "clinvar_class_counts": cv_class_counts,
        "clinvar_variants": cv_variants[:30],
    }

SECTION = Section(
    id="3", name="variant_details",
    description=("Top GWAS variants with dbSNP details (rsID / chr / pos / "
                 "alleles / MAF / functional consequence) + genetic-evidence "
                 "tier classification."),
    needs=("mondo_id",),
    produces=("top_variants", "tier_counts", "maf_distribution",
              "consequence_distribution", "clinvar_total",
              "clinvar_class_counts", "clinvar_variants"),
    datasets=DATASETS, chains=CHAINS, collect_fn=collect,
)
