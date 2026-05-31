"""§1 — drug identifiers + chemistry. Pure anchor read (no extra biobtree
calls): ChEMBL / PubChem CID / ChEBI / ATC ids, molecule type + max phase,
salt-form parent/child chain, filtered alt-names, and the chemistry block
(SMILES / InChIKey / IUPAC / formula / MW) + ChEBI one-line definition. The
chemistry block is populated only for small molecules (antibodies elide it)."""
from atlas.section import Section


def collect(a):
    return {
        "section": "01_drug_ids",
        "chembl_id": a.chembl_id,
        "canonical_name": a.canonical_name,
        "molecule_type": a.molecule_type,
        "max_phase": a.max_phase,
        "is_fda_approved": a.is_fda_approved,
        "pubchem_cid": a.pubchem_cid,
        "chebi_id": a.chebi_id,
        "atc_codes": list(a.atc_codes),
        "alt_names": list(a.alt_names),
        "parent_chembl": a.parent_chembl,
        "child_chembls": list(a.child_chembls),
        "inchi_key": a.inchi_key,
        "smiles": a.smiles,
        "iupac_name": a.iupac_name,
        "molecular_formula": a.molecular_formula,
        "molecular_weight": a.molecular_weight,
        "chebi_definition": a.chebi_definition,
    }


SECTION = Section(
    id="1", name="drug_ids",
    description=("Drug identifiers (ChEMBL / PubChem CID / ChEBI / ATC), molecule "
                 "type + max phase, salt-form chain, alt-names, and the chemistry "
                 "block (SMILES / InChIKey / IUPAC / formula / MW) from pubchem+chebi"),
    needs=("chembl_id", "canonical_name", "molecule_type", "max_phase",
           "pubchem_cid", "chebi_id", "atc_codes", "smiles", "inchi_key"),
    produces=("chembl_id", "canonical_name", "molecule_type", "max_phase",
              "pubchem_cid", "chebi_id", "atc_codes", "alt_names",
              "parent_chembl", "child_chembls", "inchi_key", "smiles",
              "molecular_formula", "molecular_weight", "chebi_definition"),
    datasets=("chembl_molecule", "pubchem", "chebi"),
    chains=(">>chembl_molecule>>pubchem", ">>chembl_molecule>>chebi"),
    collect_fn=collect,
)
