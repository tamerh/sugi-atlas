"""§15 — MeSH pharmacology. NLM's curated definition of the drug (the MeSH scope
note), its MeSH Pharmacological Action classes ("Anticoagulants", "Tyrosine
Kinase Inhibitors") and its broader MeSH heading(s) (the chemical class it sits
under: "4-Hydroxycoumarins", "Piperazines").

Resolution is an EXACT (case-insensitive, whitespace-collapsed) match of a drug
name to the MeSH record name — never fuzzy, never via entry terms/synonyms. No
id chain is usable: biobtree's pubchem→mesh edge is itself a keyword join and
links the wrong record for common drugs (warfarin → "Coumarins", a class), and
CTD's MeSH-keyed InChIKeys are sparse. Corpus diagnostics on exact-name hits:
where an id-grade signal exists (SCR registry number = PubChem UNII, CTD InChIKey
= drug InChIKey, mesh→pubchem contains the drug's CID) it agrees, and the rare
disagreements are salt/stereo variants of the same drug — so exact name is the
gate, plus a category guard (below).

Names tried, in order: the drug's own canonical name (and the caller's input
name), then — only if those miss — the parent compound's name (salt pages), then
the salt forms' names (parent pages: imatinib's only MeSH record is the
descriptor "Imatinib Mesylate").

Category guard: the record must be a chemical — a Descriptor with a D-tree
(Chemicals and Drugs) number, or a Supplementary Concept whose mapped heading is
one. That keeps same-name disease / protocol SCRs out.

Supplementary Concept "notes" are often MeSH indexing chatter rather than a
definition ("structure in first source", "RN given refers to parent cpd",
"was minor descriptor") — those clauses are stripped by clean_note()."""
import re

from atlas.biobtree import search, entry, rows, map_all
from atlas.drug.anchors import _mol_attrs
from atlas.section import Section

_PARENT_CHAIN = ">>mesh>>meshparent"     # descriptor → its broader tree heading(s)
_CHILD_CAP = 8                            # salt forms tried on a parent page
_BROADER_CAP = 6
_MESH_ID = re.compile(r"^[DC]\d{6,9}$")
_CHEMBL_ID = re.compile(r"^CHEMBL\d+$", re.I)

# One `;`-separated clause of a MeSH note that is indexing metadata, not content.
_JUNK_CLAUSE = re.compile(
    r"""^(?:
        rn\b                                            # RN given refers to …
      | (?:amino\s+acid\s+)?(?:sequence|structure)      # structure (given) in first source /
        (?:\s*&\s*rn)?                                  #   Structure & RN given in …
        (?:\s+(?:given|shown))?(?:\s+in\b.*)?\.?$       #   in Negwer, 5th ed …
      | index\s+medicus\b                               # INDEX MEDICUS search ETHANOLAMINES (72)
      | request\s+from\b                                # request from searcher 3/76
      | was\s+(?:minor|major|mh|see|heading|indexed|ep|y)\b  # was minor descriptor / was EP to …
      | .*\b(?:minor|major)\s+descriptor\b              # (Chloride was) minor descriptor (75-80)
      | .*\bwas\s+see\b                                 # Y 6047 was see CLOTIAZEPAM 1978-94
      | file\s+maintained\b                             # File maintained to MORPHOLINES (66-86)
      | note\s*:                                        # Note: tradenames that start with …
      | .*\(\d{2}-\d{2}\)\.?$                           # any clause ending in an MeSH year span
      | on[\s-]?line\s*&                                # on-line & Index Medicus search …
      | see\s+(?:also|under)\b
      | use\s+.+\bto\s+search\b
      | do\s+not\s+confuse\b
      | inchikey\b                                      # InChIKey: LIRY…
      | synonyms?\b                                     # Synonyms 3024 CERM … refers to di-HCl
      | \d+\(\d+\):\d+                                  # citation tail 386(11):1013-1025
      | .*\b(?:is|are)\s+(?:the|an?)\b.*isomers?        # "cerivastatin is the ((E)-(+))-isomer",
        (?:,?\s*respectively)?\.?$                      #   "X and Y are the (R)- and (S)-isomers"
      | .*\bmay\s+also\s+refers?\s+to\b                 # "Polcortolon may also refers to …"
      | .*\bCA\s+Vol\b                                  # "N1 is from CA Vol 90 Form Index"
      | .*\bN1\s+(?:is|in)\b                            # "N1 in Chemline is same as synonym 8"
    )""", re.I | re.X)
_CITATION = re.compile(r"\b(?:19|20)\d{2}\s+(?:jan|feb|mar|apr|may|jun|jul|aug|sep|"
                       r"oct|nov|dec)\b", re.I)
_NOT_INCLUDED = re.compile(r"\bnot included\b", re.I)
# The same "… in first source" chatter embedded mid-clause after a comma:
# "A recombinant form of TFPI, RN & amino acid sequence in first source".
_INLINE_SOURCE = re.compile(r",\s*(?:rn\s*&\s*)?(?:amino\s+acid\s+)?(?:sequence|structure)"
                            r"(?:\s*&\s*rn)?\s+(?:given\s+)?in\s+first\s+source\b", re.I)


def _norm(s):
    return " ".join((s or "").split()).casefold()


def clean_note(note, supplementary=True):
    """MeSH scope note / SCR note → display text, or '' when nothing but
    indexing chatter remains. Only Supplementary Concept notes are filtered: their
    junk `;`-clauses are dropped, and when none are the text is unchanged (bar a
    closing period and an initial capital on an all-lowercase first word — SCR
    notes start "a carbocyclic nucleoside …"). A Descriptor's scope note is a
    curated definition and is kept verbatim (a clause filter could eat a real
    defining sentence such as "… is the d-isomer of …")."""
    note = " ".join((note or "").split())
    if not note:
        return ""
    if supplementary:
        keep = [_INLINE_SOURCE.sub("", c).strip() for c in note.split(";")]
        keep = [c for c in keep if c and not _JUNK_CLAUSE.match(c)
                and not _CITATION.search(c) and not _NOT_INCLUDED.search(c)]
        note = "; ".join(keep)
    text = note.strip().rstrip(";,").strip()
    if len(text) < 10 or len(text.split()) < 2:
        return ""
    first = text.split(" ", 1)[0].rstrip(",:")
    if first.isalpha() and first.islower():
        text = text[0].upper() + text[1:]
    if text[-1] not in ".!?":
        text += "."
    return text


def _mesh_attrs(mid):
    try:
        return (entry(mid, "mesh").get("Attributes") or {}).get("Mesh") or {}
    except Exception:
        return {}


def _is_chemical_descriptor(m):
    return any(str(t).startswith("D") for t in (m.get("tree_numbers") or []))


def _split_ids(v):
    if isinstance(v, (list, tuple)):
        v = ";".join(str(x) for x in v)
    return [x for x in re.split(r"[;,\s]+", str(v or "")) if _MESH_ID.match(x)]


def _resolve_record(mid, m):
    """A matched MeSH record → the bundle fields, or None when it fails the
    chemical-category guard."""
    supp = str(m.get("is_supplementary")).lower() == "true"
    if supp:
        mapped = [(i, _mesh_attrs(i)) for i in _split_ids(m.get("heading_mapped_to"))]
        if not any(_is_chemical_descriptor(a) for _i, a in mapped):
            return None                          # disease / protocol / organism SCR
        broader = [a.get("descriptor_name") for _i, a in mapped if a.get("descriptor_name")]
    else:
        if not _is_chemical_descriptor(m):
            return None
        try:
            broader = [r.get("descriptor_name") for r in map_all(mid, _PARENT_CHAIN)
                       if r.get("descriptor_name")]
        except Exception:
            broader = []
    # Descriptor PAs come as names; SCR PAs as descriptor ids → resolve.
    actions = []
    for pa in m.get("pharmacological_actions") or []:
        pa = str(pa).strip()
        nm = (_mesh_attrs(pa).get("descriptor_name") if _MESH_ID.match(pa) else pa)
        if nm and nm not in actions:
            actions.append(nm)
    return {
        "mesh_id": mid,
        "mesh_name": m.get("descriptor_name"),
        "mesh_record_type": "supplementary" if supp else "descriptor",
        "definition": clean_note(m.get("scope_note"), supplementary=supp),
        "pharmacological_actions": actions,
        "broader_headings": list(dict.fromkeys(broader))[:_BROADER_CAP],
    }


def _match_name(name):
    """The single chemical MeSH record whose name equals `name` exactly, or None.
    Descriptor beats Supplementary Concept; two same-name records of the best
    kind are ambiguous → None (never guess)."""
    want = _norm(name)
    if not want:
        return None
    try:
        hits = rows(search(name, "mesh"))
    except Exception:
        return None
    ids = sorted({r["id"] for r in hits
                  if r.get("id") and _norm(r.get("name")) == want})
    found = []
    for mid in ids:
        m = _mesh_attrs(mid)
        if _norm(m.get("descriptor_name")) != want:   # re-check on the record itself
            continue
        rec = _resolve_record(mid, m)
        if rec:
            found.append(rec)
    if not found:
        return None
    desc = [r for r in found if r["mesh_record_type"] == "descriptor"]
    best = desc or found
    return best[0] if len(best) == 1 else None


def _chembl_name(chembl_id):
    try:
        return (_mol_attrs(entry(chembl_id, "chembl_molecule")).get("name") or "").strip()
    except Exception:
        return ""


def _candidates(a):
    """(via, chembl_id, name) in priority order; parent/salt names fetched
    lazily so the common own-name hit costs no extra lookups."""
    seen = set()
    own = [a.canonical_name]
    if a.name and not _CHEMBL_ID.match(a.name):
        own.append(a.name)
    for n in own:
        if n and _norm(n) not in seen:
            seen.add(_norm(n))
            yield "name", a.chembl_id, n
    if a.parent_chembl:
        n = _chembl_name(a.parent_chembl)
        if n and _norm(n) not in seen:
            seen.add(_norm(n))
            yield "parent", a.parent_chembl, n
    for cid in sorted(a.child_chembls or ())[:_CHILD_CAP]:
        n = _chembl_name(cid)
        if n and _norm(n) not in seen:
            seen.add(_norm(n))
            yield "salt", cid, n


def collect(a):
    out = {"section": "15_mesh_pharmacology", "mesh_id": None}
    for via, cid, name in _candidates(a):
        rec = _match_name(name)
        if rec:
            out.update(rec, matched_via=via, matched_chembl=cid, matched_name=name)
            break
    return out


SECTION = Section(
    id="15", name="mesh_pharmacology",
    description=("MeSH definition (scope note), Pharmacological Action classes and "
                 "broader MeSH heading(s), from the MeSH record whose name EXACTLY "
                 "matches the drug (or its parent / salt form); nothing on no match"),
    needs=("canonical_name", "parent_chembl", "child_chembls"),
    produces=("mesh_id", "mesh_name", "mesh_record_type", "definition",
              "pharmacological_actions", "broader_headings", "matched_via",
              "matched_chembl", "matched_name"),
    datasets=("mesh", "chembl_molecule"),
    chains=(_PARENT_CHAIN,),
    collect_fn=collect,
)
