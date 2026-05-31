#!/usr/bin/env python3
"""Drug anchors — resolve a drug name (or ChEMBL id) to the ID set, chemistry
descriptors, target list, and indication list the 12 drug sections need, ONCE.

Anchored at a ChEMBL molecule id. More involved than the gene anchor because:
- chemistry descriptors (SMILES / InChIKey / formula / MW) are NOT on the
  chembl_molecule entry — they live on the pubchem + chebi dataset entries
  (both open-licensed), reached via >>chembl_molecule>>pubchem / >>chebi.
- the pharmacological mode of action is a *target list* (chembl_target →
  uniprot → hgnc; the direct chembl_target>>hgnc edge returns 0, so we hop
  through uniprot), each needing a gene symbol for the §2 cross-link.
- drug-class semantics come from ChEBI `roles` (open) rather than the
  licensing-restricted WHO ATC name table; the raw ATC code is link-out only.
- salt forms link via parent/childs (a child entry carries `parent`, the
  parent carries `childs`).

All probed live 2026-05-31 — see SPEC_drug_entity.md §15 + the data layer
confirmed during spec iteration.
"""
import re, sys
from dataclasses import dataclass, field
from typing import Optional, Tuple, Dict

from atlas.biobtree import search, entry, rows, map_all, bbmap, map_targets


@dataclass(frozen=True)
class TargetAnchor:
    chembl_target_id: Optional[str] # CHEMBL1862 (None for gtopdb-only targets)
    target_type: str                # 'SINGLE PROTEIN' | 'PROTEIN COMPLEX' | ...
    target_name: str                # target label
    uniprot: Optional[str]          # P00519
    gene_symbol: Optional[str]      # ABL1
    hgnc_id: Optional[str]          # HGNC:76
    source: str = "gtopdb"          # 'gtopdb' (curated mechanism) | 'chembl' (bioactivity)
    action: Optional[str] = None    # GtoPdb: Inhibition / Agonist / Antibody ...
    affinity: Optional[str] = None  # GtoPdb pAffinity value (where measured)


@dataclass(frozen=True)
class IndicationRecord:
    efo_id: Optional[str]
    mesh_id: Optional[str]
    mondo_id: Optional[str]         # cross-walked from efo (fallback mesh)
    name: Optional[str]             # mondo canonical name where resolved
    max_phase: int                  # phase for THIS indication
    slug: Optional[str]             # cross-link key to Atlas disease page


@dataclass(frozen=True)
class DrugAnchors:
    name: str                       # caller-given input
    chembl_id: str                  # CHEMBL941
    chembl_entry: dict              # full entry for traceability
    canonical_name: str             # 'IMATINIB' from entry.name
    molecule_type: str              # 'Small molecule' | 'Antibody' | ...
    max_phase: int                  # 0..4
    atc_codes: Tuple[str, ...]      # raw codes only (e.g. ('L01EA01',))
    alt_names: Tuple[str, ...]      # filtered brand/generic/INN names
    parent_chembl: Optional[str]    # set if this is a salt/anhydrous form
    child_chembls: Tuple[str, ...]  # set if this is the parent form
    targets: Tuple[TargetAnchor, ...]            # primary: GtoPdb curated (or chembl fallback)
    bioactivity_targets: Tuple[dict, ...]        # secondary: raw chembl_target bioactivity set
    indications: Tuple[IndicationRecord, ...]
    xref_counts: Dict[str, int]
    # Chemistry descriptors — from pubchem + chebi entries (small-mol only).
    pubchem_cid: Optional[str]
    chebi_id: Optional[str]
    inchi_key: Optional[str]
    smiles: Optional[str]
    iupac_name: Optional[str]
    molecular_formula: Optional[str]
    molecular_weight: Optional[str]
    chebi_definition: Optional[str]
    chebi_roles: Tuple[str, ...]    # decoded role names (drug-class semantics)
    is_fda_approved: Optional[bool]


# Chemistry IUPAC strings pollute chembl_molecule.altNames; keep human-facing
# brand/generic/INN/dev-code names, drop the long systematic-chemistry strings.
_CHEM_NAME_RX = re.compile(r"[0-9]\-|piperazin|pyrimidin|benzamid|\bN\-\[|methyl\-", re.I)


def _alt_name_ok(n: str) -> bool:
    if not n or len(n) > 40:
        return False
    return not _CHEM_NAME_RX.search(n)


_HTML_TAG_RX = re.compile(r"<[^>]+>")


def _strip_html(s):
    """ChEBI definitions/role-names embed markup (<em>, <small><sup>, ...).
    Strip tags so renders/JSON-LD carry clean text."""
    if not s:
        return s
    return _HTML_TAG_RX.sub("", s).strip()


def _mol_attrs(en: dict) -> dict:
    ch = (en.get("Attributes") or {}).get("Chembl") or {}
    return ch.get("molecule") if isinstance(ch.get("molecule"), dict) else ch


def _phase(v) -> int:
    try:
        return int(v or 0)
    except (TypeError, ValueError):
        return 0


def resolve_chembl(name_or_id: str) -> Tuple[str, dict]:
    """name_or_id -> (chembl_id, entry). Direct fetch for CHEMBL ids; otherwise
    search the chembl_molecule dataset and take the highest-xref-count hit
    (the canonical reviewed molecule beats screening duplicates)."""
    if re.match(r"^CHEMBL\d+$", name_or_id):
        return name_or_id, entry(name_or_id, "chembl_molecule")
    cand = [r for r in rows(search(name_or_id, source="chembl_molecule"))
            if (r.get("id") or "").startswith("CHEMBL")]
    if not cand:
        raise LookupError(f"no chembl_molecule row for {name_or_id!r}")
    cand.sort(key=lambda r: int(r.get("xref_count") or 0), reverse=True)
    cid = cand[0]["id"]
    return cid, entry(cid, "chembl_molecule")


def _batch_map(ids, chain) -> Dict[str, list]:
    """bbmap over comma-joined ids → {input_id: [target_id, ...]}, paginated.
    Keeps the per-input keying that map_all/map_targets flatten away."""
    out: Dict[str, list] = {}
    page = None
    for _ in range(10):
        resp = bbmap(",".join(ids), chain, page)
        for m in (resp.get("mappings") or []):
            inp = (m.get("input") or "").strip()
            for t in (m.get("targets") or []):
                out.setdefault(inp, []).append(t.split("|", 1)[0])
        pg = resp.get("pagination", {}) or {}
        if not pg.get("has_next") or not pg.get("next_token"):
            break
        page = pg.get("next_token")
    return out


def _symbols_for(hgnc_ids) -> Dict[str, str]:
    """hgnc id → primary symbol (one entry per unique id; bounded by target count)."""
    out: Dict[str, str] = {}
    for h in sorted(set(h for h in hgnc_ids if h)):
        try:
            syms = ((entry(h, "hgnc").get("Attributes") or {}).get("Hgnc") or {}).get("symbols") or []
            if syms:
                out[h] = syms[0]
        except Exception:
            continue
    return out


def _gtopdb_ligand_id(canonical_name: str) -> Optional[str]:
    """Resolve a drug's GtoPdb ligand id by name. The chembl_molecule>>gtopdb_ligand
    forward edge is unwired (BIOBTREE_ISSUES #18 — xref count exists, no
    traversal), so we name-search instead. Prefer an exact case-insensitive
    name match, else the highest-xref hit."""
    res = rows(search(canonical_name, source="gtopdb_ligand"))
    if not res:
        return None
    exact = [r for r in res if (r.get("name") or "").lower() == canonical_name.lower()]
    pick = (exact or sorted(res, key=lambda r: int(r.get("xref_count") or 0), reverse=True))[0]
    return pick.get("id")


def _gtopdb_ligand_for(chembl_id: str, canonical_name: str) -> Optional[str]:
    """Resolve the GtoPdb ligand id for a drug. Prefer the authoritative
    ID-join `>>chembl_molecule>>gtopdb_ligand` (works for small molecules that
    carry a gtopdb xref — Imatinib→5687, Olaparib→7519; unambiguous, no
    name collisions). Fall back to name search for drugs without that xref —
    notably **antibodies**, whose GtoPdb ligands carry no ChEMBL id, so there's
    no shared key (BIOBTREE_ISSUES #18(a), a documented biologics-coverage
    gap, not a bug)."""
    rows_ = map_all(chembl_id, ">>chembl_molecule>>gtopdb_ligand", cap=1)
    if rows_ and rows_[0].get("id"):
        return rows_[0]["id"]
    return _gtopdb_ligand_id(canonical_name)


def _gtopdb_targets(chembl_id: str, canonical_name: str) -> Tuple[TargetAnchor, ...]:
    """Curated mechanism targets via GtoPdb — works for BOTH small molecules
    and antibodies (ChEMBL bioactivity edges miss antibodies entirely).

      gtopdb_ligand →gtopdb_interaction→ gtopdb(target) →uniprot→ hgnc

    The interaction row carries target_name + action + affinity; the gtopdb
    target node cross-walks to the human UniProt (the interaction maps to
    human + rodent orthologs — we keep the one that resolves to an HGNC id)."""
    lig = _gtopdb_ligand_for(chembl_id, canonical_name)
    if not lig:
        return ()
    inter = map_all(lig, ">>gtopdb_ligand>>gtopdb_interaction")
    # group by gtopdb target id (the "{target}_{ligand}" interaction id).
    # GUARD: until the gtopdb re-update lands (BIOBTREE_ISSUES #18(b), fixed
    # upstream — bare-numeric synonym tokens like "7519" from "AT-7519" no
    # longer collide with the ligand-ID namespace), the interaction edge can
    # leak another ligand's rows by id-substring. Keep only rows whose trailing
    # ligand id matches. Harmless (no-op) once the re-update ships.
    tinfo: Dict[str, dict] = {}
    for r in inter:
        parts = (r.get("id") or "").split("_")
        if len(parts) < 2 or parts[-1] != lig:
            continue
        tid = parts[0]
        if tid and tid not in tinfo:
            tinfo[tid] = {"name": r.get("target_name") or "",
                          "type": r.get("type") or "",
                          "action": r.get("action") or None,
                          "affinity": r.get("affinity") or None}
    if not tinfo:
        return ()
    gt2uni = _batch_map(list(tinfo), ">>gtopdb>>uniprot")
    all_uni = sorted({u for us in gt2uni.values() for u in us})
    uni2hgnc = _batch_map(all_uni, ">>uniprot>>hgnc") if all_uni else {}
    sym = _symbols_for(h[0] for h in uni2hgnc.values() if h)
    out = []
    for tid, info in tinfo.items():
        uni = hg = None
        for u in gt2uni.get(tid, []):     # pick the human ortholog (resolves to hgnc)
            if uni2hgnc.get(u):
                uni, hg = u, uni2hgnc[u][0]
                break
        out.append(TargetAnchor(
            chembl_target_id=None, target_type=info["type"], target_name=info["name"],
            uniprot=uni, hgnc_id=hg, gene_symbol=sym.get(hg) if hg else None,
            source="gtopdb", action=info["action"], affinity=info["affinity"]))
    return tuple(out)


def _resolve_targets(canonical_name: str, chembl_id: str):
    """(primary_targets, bioactivity_targets).

    Primary = GtoPdb curated mechanism targets (action + affinity; covers
    antibodies). Fallback to gene-resolved ChEMBL targets only when GtoPdb has
    no ligand. Secondary = the raw chembl_target bioactivity set (id/name/type,
    NOT gene-resolved — cheap; feeds §2 'broader targets' + §3 grouping)."""
    primary = _gtopdb_targets(chembl_id, canonical_name)

    bio = map_all(chembl_id, ">>chembl_molecule>>chembl_target")
    bioactivity = tuple({"chembl_target_id": t.get("id"), "name": t.get("title"),
                         "type": t.get("type")} for t in bio if t.get("id"))

    # Fallback: no GtoPdb curated targets → gene-resolve the top ChEMBL targets.
    # chembl_target's xref is to UniProt (components.acc), so the gene is reached
    # via chembl_target>>uniprot>>hgnc — biobtree's documented pattern, not a gap.
    if not primary and bioactivity:
        top = [b["chembl_target_id"] for b in bioactivity[:25]]
        t2uni = _batch_map(top, ">>chembl_target>>uniprot")
        uni2hgnc = _batch_map(sorted({u for us in t2uni.values() for u in us}),
                              ">>uniprot>>hgnc")
        sym = _symbols_for(h[0] for h in uni2hgnc.values() if h)
        resolved = []
        for b in bioactivity[:25]:
            uni = (t2uni.get(b["chembl_target_id"]) or [None])[0]
            hg = (uni2hgnc.get(uni) or [None])[0] if uni else None
            resolved.append(TargetAnchor(
                chembl_target_id=b["chembl_target_id"], target_type=b["type"] or "",
                target_name=b["name"] or "", uniprot=uni, hgnc_id=hg,
                gene_symbol=sym.get(hg) if hg else None, source="chembl"))
        primary = tuple(resolved)
    return primary, bioactivity


def _resolve_indications(mol: dict) -> Tuple[IndicationRecord, ...]:
    """entry.indications[] = [{highestDevelopmentPhase, efo, mesh}]. Dedupe by
    efo (many rows repeat), cross-walk efo→mondo (fallback mesh→mondo), keep
    the max phase seen per indication."""
    from atlas.disease.slug import slugify as disease_slug

    raw = mol.get("indications") or []
    by_efo: Dict[str, dict] = {}
    for ind in raw:
        efo = ind.get("efo")
        mesh = ind.get("mesh")
        key = efo or mesh
        if not key:
            continue
        ph = _phase(ind.get("highestDevelopmentPhase"))
        cur = by_efo.get(key)
        if cur is None or ph > cur["phase"]:
            by_efo[key] = {"efo": efo, "mesh": mesh, "phase": ph}

    out = []
    for rec in by_efo.values():
        mondo_id = name = None
        # cross-walk to mondo via efo, then mesh
        for src_id, src in ((rec["efo"], "efo"), (rec["mesh"], "mesh")):
            if not src_id:
                continue
            mres = map_all(src_id, f">>{src}>>mondo", cap=1)
            if mres:
                mondo_id = mres[0].get("id")
                name = mres[0].get("name")
                break
        out.append(IndicationRecord(
            efo_id=rec["efo"], mesh_id=rec["mesh"], mondo_id=mondo_id,
            name=name, max_phase=rec["phase"],
            slug=disease_slug(name) if name else None,
        ))
    # highest-phase indications first
    return tuple(sorted(out, key=lambda r: -r.max_phase))


def _resolve_chemistry(chembl_id: str):
    """pubchem + chebi descriptors. Returns a dict (empty for biologics with
    no small-mol CID)."""
    chem = {"pubchem_cid": None, "chebi_id": None, "inchi_key": None,
            "smiles": None, "iupac_name": None, "molecular_formula": None,
            "molecular_weight": None, "chebi_definition": None,
            "chebi_roles": (), "is_fda_approved": None}

    pc = map_all(chembl_id, ">>chembl_molecule>>pubchem", cap=1)
    if pc:
        chem["pubchem_cid"] = pc[0].get("id")
        fda = (pc[0].get("is_fda_approved") or "").lower()
        chem["is_fda_approved"] = True if fda == "true" else (False if fda == "false" else None)
        try:
            pa = (entry(chem["pubchem_cid"], "pubchem").get("Attributes") or {}).get("Pubchem") or {}
            chem.update(inchi_key=pa.get("inchi_key"), smiles=pa.get("smiles"),
                        iupac_name=pa.get("iupac_name"),
                        molecular_formula=pa.get("molecular_formula"),
                        molecular_weight=pa.get("molecular_weight"))
        except Exception:
            pass

    cb = map_all(chembl_id, ">>chembl_molecule>>chebi", cap=1)
    if cb:
        chem["chebi_id"] = cb[0].get("id")
        try:
            ca = (entry(chem["chebi_id"], "chebi").get("Attributes") or {}).get("Chebi") or {}
            chem["chebi_definition"] = _strip_html(ca.get("definition"))
            chem.setdefault("inchi_key") or chem.update(inchi_key=chem["inchi_key"] or ca.get("inchi_key"))
            if not chem["smiles"]:
                chem["smiles"] = ca.get("smiles")
            if not chem["molecular_formula"]:
                chem["molecular_formula"] = ca.get("formula")
            # decode role ids -> names (open-licensed drug-class semantics)
            roles = []
            for rid in (ca.get("roles") or []):
                try:
                    rn = _strip_html(((entry(rid, "chebi").get("Attributes") or {}).get("Chebi") or {}).get("name"))
                    if rn:
                        roles.append(rn)
                except Exception:
                    continue
            chem["chebi_roles"] = tuple(roles)
        except Exception:
            pass
    return chem


def resolve(name_or_id: str) -> DrugAnchors:
    """Drug name (or ChEMBL id) → DrugAnchors. ~15-25 biobtree calls."""
    chembl_id, en = resolve_chembl(name_or_id)
    mol = _mol_attrs(en)
    xc = {r.split("|")[0]: int(r.split("|")[1])
          for r in (en.get("xrefs", {}).get("data") or [])}
    canonical = mol.get("name") or chembl_id

    chem = _resolve_chemistry(chembl_id)
    targets, bioactivity_targets = _resolve_targets(canonical, chembl_id)

    return DrugAnchors(
        name=name_or_id,
        chembl_id=chembl_id,
        chembl_entry=en,
        canonical_name=canonical,
        molecule_type=mol.get("type") or "",
        max_phase=_phase(mol.get("highestDevelopmentPhase")),
        atc_codes=tuple(mol.get("atcClassification") or ()),
        alt_names=tuple(n for n in (mol.get("altNames") or []) if _alt_name_ok(n)),
        parent_chembl=mol.get("parent"),
        child_chembls=tuple(mol.get("childs") or ()),
        targets=targets,
        bioactivity_targets=bioactivity_targets,
        indications=_resolve_indications(mol),
        xref_counts=xc,
        pubchem_cid=chem["pubchem_cid"],
        chebi_id=chem["chebi_id"],
        inchi_key=chem["inchi_key"],
        smiles=chem["smiles"],
        iupac_name=chem["iupac_name"],
        molecular_formula=chem["molecular_formula"],
        molecular_weight=chem["molecular_weight"],
        chebi_definition=chem["chebi_definition"],
        chebi_roles=chem["chebi_roles"],
        is_fda_approved=chem["is_fda_approved"],
    )


if __name__ == "__main__":
    a = resolve(sys.argv[1] if len(sys.argv) > 1 else "Imatinib")
    print(f"{a.canonical_name} ({a.chembl_id}) type={a.molecule_type} phase={a.max_phase}")
    print(f"  ATC={a.atc_codes}  roles={a.chebi_roles}")
    print(f"  chem: CID={a.pubchem_cid} ChEBI={a.chebi_id} formula={a.molecular_formula} "
          f"MW={a.molecular_weight} inchikey={a.inchi_key}")
    print(f"  fda_approved={a.is_fda_approved}  alt_names={a.alt_names[:6]}")
    print(f"  parent={a.parent_chembl} childs={a.child_chembls}")
    print(f"  primary targets={len(a.targets)} (src={a.targets[0].source if a.targets else '-'}):")
    for t in a.targets[:6]:
        print(f"      {t.gene_symbol or t.target_name[:30]:30} {t.action or ''} aff={t.affinity or '-'} ({t.uniprot or '?'})")
    print(f"  bioactivity targets (chembl, raw)={len(a.bioactivity_targets)}")
    print(f"  indications={len(a.indications)}  (sample: "
          + ", ".join(f"{i.name or i.efo_id}(p{i.max_phase})" for i in a.indications[:5]) + ")")
