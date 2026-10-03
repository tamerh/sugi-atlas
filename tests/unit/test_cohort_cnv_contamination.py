"""Disease cohort: a gene whose ONLY ClinVar support is a large multi-gene region
CNV is a bystander, not an associated gene — it must not seed the cohort (audit: a
272 kb chr1 deletion put LDLRAP1 on the SELENON / rigid-spine-MD disease, polluting
its pathways/druggability with LDL-clearance)."""
from atlas.disease import anchors


def test_is_large_cnv():
    assert anchors._is_large_cnv("NC_000001.10:g.(?_25870190)_(26142209_?)del")   # 272 kb, fuzzy
    assert anchors._is_large_cnv("NC_000001.11:g.1000_900000dup")                 # 899 kb
    assert not anchors._is_large_cnv("NM_206926.2(SELENON):c.1A>G")               # transcript SNV
    assert not anchors._is_large_cnv("NC_000001.11:g.100_110del")                 # 10 bp
    assert not anchors._is_large_cnv("NC_000001.11:g.100_199del")                 # 99 bp < 200 kb
    assert not anchors._is_large_cnv("")


def test_cnv_only_gene_is_dropped_localized_kept(monkeypatch):
    def fake_map_all(ident, chain, **k):
        if chain == ">>mondo>>clinvar":
            return [
                {"gene_symbol": "SELENON", "name": "NM_206926.2(SELENON):c.1A>G"},      # localized
                {"gene_symbol": "SELENON", "name": "NC_000001.10:g.(?_26120000)_(26130000_?)del"},  # 10kb, localized-ish (<200kb)
                {"gene_symbol": "LDLRAP1", "name": "NC_000001.10:g.(?_25870190)_(26142209_?)del"},   # 272kb CNV-only bystander
            ]
        if chain == ">>hgnc":
            return {"LDLRAP1": [{"id": "HGNC:18640"}], "SELENON": [{"id": "HGNC:15999"}]}.get(ident, [])
        return []
    monkeypatch.setattr(anchors, "map_all", fake_map_all)
    drop = anchors._clinvar_cnv_only_genes("MONDO:0011271")
    assert drop == {"HGNC:18640"}        # LDLRAP1 (only a large CNV) dropped; SELENON kept (has localized)


def test_gene_with_any_localized_record_not_dropped(monkeypatch):
    # A gene appearing in BOTH a large CNV and a transcript record stays.
    def fake_map_all(ident, chain, **k):
        if chain == ">>mondo>>clinvar":
            return [
                {"gene_symbol": "BIGGENE", "name": "NC_000002.12:g.(?_1)_(9000000_?)del"},  # 9 Mb CNV
                {"gene_symbol": "BIGGENE", "name": "NM_1.2(BIGGENE):c.500G>T"},             # localized → keep
            ]
        return [{"id": "HGNC:1"}] if chain == ">>hgnc" else []
    monkeypatch.setattr(anchors, "map_all", fake_map_all)
    assert anchors._clinvar_cnv_only_genes("MONDO:x") == set()
