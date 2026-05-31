"""§12 — salt forms & parent. Render-only navigation over the anchor's
parent/child ChEMBL linkage. ChEMBL treats salt/anhydrous forms as distinct
molecule IDs (CHEMBL941 IMATINIB parent ↔ CHEMBL1642 IMATINIB MESYLATE child);
this section makes the relationship explicit + navigable."""
from atlas.section import Section


def collect(a):
    return {
        "section": "12_salt_forms",
        "chembl_id": a.chembl_id,
        "canonical_name": a.canonical_name,
        "parent_chembl": a.parent_chembl,
        "child_chembls": list(a.child_chembls),
    }


SECTION = Section(
    id="12", name="salt_forms",
    description="Parent / salt-form (child) ChEMBL linkage for navigation",
    needs=("parent_chembl", "child_chembls"),
    produces=("parent_chembl", "child_chembls"),
    datasets=("chembl_molecule",),
    chains=(),
    collect_fn=collect,
)
