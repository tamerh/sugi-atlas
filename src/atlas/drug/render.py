#!/usr/bin/env python3
"""Deterministic markdown renderer for drug-section bundles — NO model.
Mirrors gene/render.py + disease/render.py: one r_* fn per section + a RENDER
dict. Every fact comes verbatim from the bundle."""
from atlas.render_common import table


def _i(n):
    return f"{n:,}" if isinstance(n, int) else (n if n is not None else "")


def r_drug_ids(b):
    L = ["## Drug identifiers", ""]
    rows = [
        ("ChEMBL ID", b.get("chembl_id")),
        ("Name", b.get("canonical_name")),
        ("Type", b.get("molecule_type")),
        ("Max phase", b.get("max_phase")),
        ("FDA approved", b.get("is_fda_approved")),
        ("PubChem CID", b.get("pubchem_cid")),
        ("ChEBI", b.get("chebi_id")),
        ("ATC", ", ".join(b.get("atc_codes") or [])),
        ("Molecular formula", b.get("molecular_formula")),
        ("Molecular weight", b.get("molecular_weight")),
        ("InChIKey", b.get("inchi_key")),
    ]
    L.append(table(["Field", "Value"], [(k, v) for k, v in rows if v not in (None, "")]))
    if b.get("smiles"):
        L.append(f"\n**SMILES:** `{b['smiles']}`")
    if b.get("iupac_name"):
        L.append(f"\n**IUPAC name:** {b['iupac_name']}")
    if b.get("chebi_definition"):
        L.append(f"\n**ChEBI definition:** {b['chebi_definition']}")
    an = b.get("alt_names") or []
    if an:
        L.append(f"\n**Also known as:** {', '.join(an[:12])}")
    # salt-form navigation
    if b.get("parent_chembl"):
        L.append(f"\n*Salt/anhydrous form of parent* `{b['parent_chembl']}`.")
    if b.get("child_chembls"):
        L.append(f"\n*Parent form; salt/anhydrous children:* "
                 + ", ".join(f"`{c}`" for c in b["child_chembls"]))
    return "\n".join(L)


def r_bioactivity(b):
    L = ["## Bioactivity", ""]
    ca = b.get("activities") or []
    if not ca:
        L.append("*No ChEMBL bioactivity rows at pChembl ≥ 5 "
                 "(expected for biologics / antibodies).*")
        return "\n".join(L)
    L.append(f"**ChEMBL activities: {_i(b.get('potent_count'))} potent at "
             f"pChembl ≥ 5 of {_i(b.get('activity_total'))} total. Top 30 by "
             f"potency (10 = 0.1 nM, 6 = 1 µM):**\n")
    L.append(table(["pChembl", "Type", "Value", "Unit", "Activity ID"],
                   [(r.get("pchembl"), r.get("type"), r.get("value"),
                     r.get("unit"), r.get("id")) for r in ca]))
    return "\n".join(L)


def r_clinical_trials(b):
    L = ["## Clinical trials", "",
         f"**Total trials: {_i(b.get('trial_count'))}.**"]
    pc = b.get("phase_counts") or {}
    if pc:
        L += ["", "**Phase distribution:**", "",
              table(["Phase", "Trials"],
                    [(k, _i(v)) for k, v in sorted(pc.items(), key=lambda kv: -kv[1])])]
    tt = b.get("top_trials") or []
    if tt:
        L += ["", "**Top trials by phase / activity:**", "",
              table(["NCT", "Phase", "Status", "Title"],
                    [(f"[{t['id']}](https://clinicaltrials.gov/study/{t['id']})" if t.get("id") else "",
                      t.get("phase"), t.get("status"), (t.get("title") or "")[:65])
                     for t in tt])]
    return "\n".join(L)


def r_pharmacology(b):
    body = []
    roles = b.get("chebi_roles") or []
    if roles:
        body.append(f"**Drug class (ChEBI roles):** {', '.join(roles)}.")
    atc = b.get("atc_codes") or []
    if atc:
        links = ", ".join(
            f"[{c}](https://www.whocc.no/atc_ddd_index/?code={c})" for c in atc)
        body.append(f"**ATC classification:** {links} "
                    f"(WHO ATC level names are licensing-restricted — code links out).")
    if not body:
        body.append("*No ChEBI role or ATC classification available.*")
    return "## Pharmacology\n\n" + "\n\n".join(body)


def r_patent_literature(b):
    L = ["## Patent literature", ""]
    total = b.get("patent_total") or 0
    if total:
        L.append(f"**Patent mentions (SureChEMBL): {_i(total)}** across "
                 f"{len(b.get('patent_compound_ids') or [])} patent_compound record(s). "
                 f"Counts attach to the compound, so promiscuous molecules score high.")
    else:
        L.append("*No SureChEMBL patent_compound records "
                 "(expected for biologics).*")
    return "\n".join(L)


RENDER = {
    "1": r_drug_ids,
    "3": r_bioactivity,
    "5": r_clinical_trials,
    "6": r_pharmacology,
    "11": r_patent_literature,
}


def render_all(bundles):
    return "\n\n".join(RENDER[sid](bundles[sid]) for sid in RENDER if sid in bundles)
