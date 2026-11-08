"""§3 — variant details: top GWAS variants → dbSNP attributes + tier
classification (Tier 1 coding / Tier 2 splice-UTR / Tier 3 regulatory /
Tier 4 intronic-intergenic)."""
import re, html
from collections import Counter
from urllib.parse import unquote
from atlas.section import Section
from atlas.biobtree import map_all, entry
from atlas.render_common import clinvar_stars, clinvar_id_num

CHAINS   = (">>mondo>>gwas>>dbsnp", ">>mondo>>clinvar", ">>mondo>>dbvar")
DATASETS = ("mondo", "gwas", "dbsnp", "clinvar", "dbvar")

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

# ── dbVar structural variants ───────────────────────────────────────────────
# dbVar is per-SUBMISSION, not per-locus: Duchenne's 1,066 records are 685
# overlapping deletions of the same region. Rendering the raw rows gives 15
# near-identical multi-Mb deletions, so records are clustered into loci by
# (chromosome, svtype, >=50% reciprocal overlap). "Records" then carries the
# submission count behind each locus. Same discipline as the ClinVar block:
# clinically actionable first, honest counts, no per-row source links.
_DBVAR_RANK = [("conflicting", 3), ("pathogenic/likely pathogenic", 0),
               ("likely pathogenic", 1), ("pathogenic", 0),
               ("uncertain", 4), ("likely benign", 5), ("benign", 5),
               ("not provided", 6)]
_DBVAR_OVERLAP = 0.5
DBVAR_CAP = 10          # loci rendered per page


def _dbvar_rank(sig):
    s = (sig or "").lower()
    for key, r in _DBVAR_RANK:
        if key in s:
            return r
    return 4


def _dbvar_int(v):
    try:
        return int(float(v))
    except (TypeError, ValueError):
        return 0


def dbvar_size(n):
    if not n or n <= 0:
        return "—"
    if n >= 1_000_000:
        return f"{n / 1_000_000:.1f} Mb"
    if n >= 1_000:
        return f"{n / 1_000:.0f} kb"
    return f"{n} bp"


def _dbvar_cluster(records, frac=_DBVAR_OVERLAP):
    """Group overlapping dbVar records into loci. Same chromosome + same svtype
    + >=frac reciprocal overlap merges; the locus span is the union. Each locus
    keeps its submission count, its majority significance (or an honest 'mixed'
    label), and the best record's ClinVar accession."""
    loci = []
    for r in sorted(records, key=lambda x: (x["chrom"], x["start"])):
        span = max(1, r["end"] - r["start"])
        for c in loci:
            if c["chrom"] != r["chrom"] or c["svtype"] != r["svtype"]:
                continue
            ov = min(c["end"], r["end"]) - max(c["start"], r["start"])
            if ov > 0 and ov / min(span, max(1, c["end"] - c["start"])) >= frac:
                c["start"] = min(c["start"], r["start"])
                c["end"] = max(c["end"], r["end"])
                c["members"].append(r)
                break
        else:
            loci.append({"chrom": r["chrom"], "start": r["start"], "end": r["end"],
                         "svtype": r["svtype"], "members": [r]})
    for c in loci:
        best = sorted(c["members"], key=lambda m: (_dbvar_rank(m["sig"]),
                                                   -(m["end"] - m["start"])))
        ranks = Counter(_dbvar_rank(m["sig"]) for m in c["members"])
        cnt = ranks[min(ranks)]
        c["n"] = len(c["members"])
        c["sig"] = (best[0]["sig"] if cnt * 2 >= c["n"]
                    else f"mixed ({cnt}/{c['n']} pathogenic/likely pathogenic)")
        c["clinvar"] = best[0]["clinvar"]
        # Union span; insertions have start == end, so fall back to svlen.
        c["size"] = c["end"] - c["start"] or abs(_dbvar_int(c["members"][0].get("svlen")))
        del c["members"]                      # keep the cached bundle small
    loci.sort(key=lambda c: (-c["n"], c["chrom"], c["start"]))
    return loci


def _collect_dbvar(a):
    """dbVar structural variants for this disease. The authoritative total comes
    from the Mondo entry's xref counts (no call when zero); the paginated fetch
    is disclosed as a floor when it falls short."""
    total = int((a.xref_counts or {}).get("dbvar") or 0)
    empty = {"dbvar_total": 0, "dbvar_retrieved": 0, "dbvar_plp": 0,
             "dbvar_loci": [], "dbvar_svtype_counts": {}}
    if total <= 0:
        return empty
    cap = min(40, total // 100 + 2)           # ~100 rows/page + headroom
    rows = map_all(a.mondo_id, ">>mondo>>dbvar", cap=cap)
    seen, recs = set(), []
    for r in rows:
        vid = r.get("id")
        if not vid or vid in seen:
            continue
        seen.add(vid)
        recs.append({
            "id": vid,
            "chrom": str(r.get("chromosome") or ""),
            "start": _dbvar_int(r.get("start")),
            "end": _dbvar_int(r.get("end")),
            "svlen": r.get("svlen"),
            "svtype": (r.get("svtype") or "?").strip(),
            "sig": unquote(r.get("clinical_significance") or "not provided").strip(),
            "origin": r.get("origin"),
            "clinvar": (r.get("clinvar_accessions") or "").split(",")[0] or None,
        })
    if not recs:
        return empty
    return {
        "dbvar_total": total,
        "dbvar_retrieved": len(recs),
        "dbvar_plp": sum(1 for r in recs if _dbvar_rank(r["sig"]) <= 1),
        "dbvar_loci": _dbvar_cluster(recs),
        "dbvar_svtype_counts": dict(Counter(r["svtype"] for r in recs).most_common()),
    }


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
    # (Likely) pathogenic, non-region-CNV ClinVar records per gene for THIS
    # disease — the evidence weight behind the ClinVar-only causal tier (a single
    # stray submission put SGCB on sialidosis; a deletion's bystander, UQCRB, on
    # Leri pleonosteosis).
    from atlas.gene.sections.s06_variants import _is_region_cnv, _mane_transcript
    plp = {}
    for v in cv_variants:
        c = (v["classification"] or "").lower()
        if ("pathogenic" in c and "conflicting" not in c and "benign" not in c
                and v.get("gene") and not _is_region_cnv(v["hgvs"])):
            plp[v["gene"]] = plp.get(v["gene"], 0) + 1
    cv_variants.sort(key=lambda v: (_cv_rank(v["classification"]),
                                    -clinvar_stars(v["review_status"]), clinvar_id_num(v["id"])))
    # Attach the MANE Select consequence to the rows we actually render (after
    # the sort — enriching the pre-sort slice put them on the wrong variants).
    cv_top = cv_variants[:30]
    for v in cv_top:
        try:
            at = (entry(v["id"], "clinvar").get("Attributes") or {}).get("Clinvar") or {}
        except Exception:
            continue
        rec = _mane_transcript(at)
        if rec:
            v["consequence"] = rec.get("consequence")

    dbvar = _collect_dbvar(a)

    return {
        "section": "03_variant_details",
        "mondo_id": a.mondo_id,
        "top_variants": top_variants,
        "tier_counts": tier_counts,
        "maf_distribution": maf_distribution,
        "consequence_distribution": consequence_distribution,
        "clinvar_total": len(seen_cv),
        "clinvar_class_counts": cv_class_counts,
        "clinvar_variants": cv_top,
        "clinvar_plp_by_gene": plp,
        **dbvar,
    }

SECTION = Section(
    id="3", name="variant_details",
    description=("Top GWAS variants with dbSNP details (rsID / chr / pos / "
                 "alleles / MAF / functional consequence) + genetic-evidence "
                 "tier classification; ClinVar germline variants; dbVar "
                 "structural variants grouped into loci."),
    needs=("mondo_id", "xref_counts"),
    produces=("top_variants", "tier_counts", "maf_distribution",
              "consequence_distribution", "clinvar_total",
              "clinvar_class_counts", "clinvar_variants", "clinvar_plp_by_gene",
              "dbvar_total", "dbvar_retrieved", "dbvar_plp", "dbvar_loci",
              "dbvar_svtype_counts"),
    datasets=DATASETS, chains=CHAINS, collect_fn=collect,
)
