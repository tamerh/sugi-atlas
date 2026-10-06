"""§2 — targets. Primary = GtoPdb curated mechanism targets (gene + action +
pAffinity; covers antibodies). Secondary = the raw ChEMBL bioactivity target
set (count + sample). Each primary target is annotated with its DepMap
cancer-dependency signal (>>hgnc>>depmap) — the target-tractability view:
"is this drug's target a cancer dependency?". Gene symbols become /atlas/gene/
links once links.py exists (post-build pass)."""
from atlas.biobtree import map_all, entry
from atlas.section import Section


def _depmap(hgnc_id):
    """DepMap per-gene CRISPR fitness aggregate; {} on miss (non-cancer-screened
    genes / no hgnc)."""
    if not hgnc_id:
        return {}
    r = map_all(hgnc_id, ">>hgnc>>depmap", cap=1)
    return r[0] if r else {}


def _hgnc_symbol(hgnc_id):
    """HGNC id -> primary symbol (the chembl_mechanism>>hgnc projection carries no
    symbol). None on miss."""
    try:
        h = (entry(hgnc_id, "hgnc").get("Attributes") or {}).get("Hgnc") or {}
        syms = h.get("symbols") or []
        return syms[0] if syms else None
    except Exception:
        return None


def _mechanisms(chembl_id, child_ids=()):
    """ChEMBL curated mechanism-of-action (drug_mechanism table) → (moa rows, target
    genes, target organisms).

    - Queries the parent AND its salt forms: ChEMBL often keys the mechanism on the
      salt (imatinib's lives on CHEMBL1642 imatinib mesylate), so a parent-only
      lookup left ~800 drug pages with no curated mechanism.
    - Genes: the >>chembl_mechanism>>hgnc edge only resolves nucleic-acid targets
      (inclisiran → PCSK9 mRNA), so protein targets are followed via
      chembl_target>>uniprot>>hgnc — SINGLE PROTEIN targets only (a family/complex,
      e.g. diazepam's GABA-A receptor group, stays a named target and is not
      expanded into every member gene). Adalimumab → TNF; imatinib → ABL1/PDGFRB/KIT.
    - A target with no human gene is a pathogen/non-human target: its organism is
      returned, for any target type (isoniazid → Mycobacterium tuberculosis,
      telaprevir's NS3/NS4A complex → hepatitis C virus)."""
    if not chembl_id:
        return [], [], []
    ids = [chembl_id] + [c for c in (child_ids or ()) if c and c != chembl_id][:5]
    moa, seen_moa = [], set()
    for mid in ids:
        for r in map_all(mid, ">>chembl_molecule>>chembl_mechanism", cap=2):
            key = (r.get("mechanism_of_action"), r.get("action_type"), r.get("target_name"))
            if key in seen_moa:
                continue
            seen_moa.add(key)
            moa.append({"mechanism_of_action": r.get("mechanism_of_action"),
                        "action_type": r.get("action_type"),
                        "target_name": r.get("target_name"),
                        "target_type": r.get("target_type")})
    genes, seen, organisms = [], set(), []

    def _add_gene(hid):
        if hid and hid not in seen:
            seen.add(hid)
            genes.append({"hgnc_id": hid, "gene_symbol": _hgnc_symbol(hid)})

    seen_tgt = set()
    for mid in ids:
        for r in map_all(mid, ">>chembl_molecule>>chembl_mechanism>>hgnc", cap=2):
            _add_gene(r.get("id"))                       # nucleic-acid targets
        for t in map_all(mid, ">>chembl_molecule>>chembl_mechanism>>chembl_target", cap=2):
            tid = t.get("id")
            if not tid or tid in seen_tgt:
                continue
            seen_tgt.add(tid)
            if (t.get("type") or "").upper() == "SINGLE PROTEIN":
                hg = [x.get("id") for x in map_all(tid, ">>chembl_target>>uniprot>>hgnc", cap=1)
                      if (x.get("id") or "").startswith("HGNC:")]
                if hg:
                    for h in hg[:1]:
                        _add_gene(h)
                    continue
            # Organism for any non-human target, whatever its type — a viral PROTEIN
            # COMPLEX (telaprevir → HCV NS3/NS4A) is as much a pathogen target as a
            # single protein. Human families/complexes resolve to Homo sapiens and are
            # skipped. (biobtree v2.13.0 now emits these targets' metadata, #61.)
            for tx in map_all(tid, ">>chembl_target>>taxonomy", cap=1):
                org = (tx.get("name") or "").strip()
                if org and org != "Homo sapiens" and org not in organisms:
                    organisms.append(org)
    return moa, genes, organisms


def collect(a):
    primary = []
    for t in a.targets:
        dm = _depmap(t.hgnc_id)
        primary.append({
            "gene_symbol": t.gene_symbol, "hgnc_id": t.hgnc_id, "uniprot": t.uniprot,
            "target_name": t.target_name, "target_type": t.target_type,
            "source": t.source, "action": t.action, "affinity": t.affinity,
            "dep_pct": dm.get("pct_dependent"),
            "dep_selective": dm.get("strongly_selective") == "true",
            "dep_common_essential": dm.get("common_essential") == "true",
        })
    bioactivity = [{"chembl_target_id": t.get("chembl_target_id"),
                    "name": t.get("name"), "type": t.get("type")}
                   for t in a.bioactivity_targets]
    mechanisms, mechanism_genes, mechanism_organisms = _mechanisms(
        a.chembl_id, getattr(a, "child_chembls", None) or ())
    from atlas.predict import resolve_schembl
    return {
        "section": "02_targets",
        "smiles": a.smiles,   # Sugi Predict cross-link (search fallback)
        "predict_schembl": resolve_schembl(a.smiles),   # direct compound page when in the atlas
        "primary_targets": primary,
        "primary_source": (a.targets[0].source if a.targets else None),
        "bioactivity_target_count": len(bioactivity),
        "bioactivity_targets": bioactivity[:30],
        "mechanisms": mechanisms,
        "mechanism_genes": mechanism_genes,
        "mechanism_organisms": mechanism_organisms,
    }


SECTION = Section(
    id="2", name="targets",
    description=("Primary mechanism targets (GtoPdb-curated: gene + action + "
                 "pAffinity; covers antibodies) annotated with DepMap cancer-"
                 "dependency + secondary ChEMBL bioactivity target set"),
    needs=("targets", "bioactivity_targets"),
    produces=("smiles", "predict_schembl", "primary_targets", "primary_source",
              "bioactivity_target_count", "bioactivity_targets", "mechanisms", "mechanism_genes",
              "mechanism_organisms"),
    datasets=("gtopdb_ligand", "gtopdb_interaction", "gtopdb", "uniprot",
              "hgnc", "chembl_target", "chembl_mechanism", "depmap", "taxonomy"),
    chains=(">>gtopdb_ligand>>gtopdb_interaction>>gtopdb>>uniprot>>hgnc",
            ">>chembl_molecule>>chembl_target", ">>hgnc>>depmap",
            ">>chembl_molecule>>chembl_mechanism",
            ">>chembl_molecule>>chembl_mechanism>>hgnc",
            ">>chembl_molecule>>chembl_mechanism>>chembl_target",
            ">>chembl_target>>uniprot>>hgnc", ">>chembl_target>>taxonomy"),
    collect_fn=collect,
)
