"""§3 — variant details: top GWAS variants → dbSNP attributes + tier
classification (Tier 1 coding / Tier 2 splice-UTR / Tier 3 regulatory /
Tier 4 intronic-intergenic).

NEW collector. Skeleton."""
from atlas.section import Section

CHAINS   = (">>mondo>>gwas>>dbsnp", ">>mondo>>gwas_study>>gwas>>dbsnp")
DATASETS = ("mondo", "gwas", "gwas_study", "dbsnp")

def collect(a):
    bundle = {"section": "03_variant_details", "mondo_id": a.mondo_id,
              "_todo": "Implement: TOP 50 GWAS variants -> dbsnp -> tier"}
    return bundle

SECTION = Section(
    id="3", name="variant_details",
    description=("Top GWAS variants with dbSNP details (rsID / chr / pos / "
                 "alleles / MAF / functional consequence) + genetic-evidence "
                 "tier classification."),
    needs=("mondo_id",),
    produces=("top_variants", "tier_counts", "maf_distribution",
              "consequence_distribution"),
    datasets=DATASETS, chains=CHAINS, collect_fn=collect,
)
