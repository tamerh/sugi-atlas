#!/usr/bin/env python3
"""Pathway anchor resolution — a Reactome pathway id (or name) → the shared facts
every section reuses: name, member genes (with symbols), parent/child hierarchy,
and the GO xref. Resolved ONCE per pathway (mirrors drug/disease anchors).

Reactome chains (verified live):
  R-HSA-… >>reactome>>ensembl        -> member genes (id|name|biotype|genome)
  R-HSA-… >>reactome>>reactomeparent -> parent pathway (id|name|tax_id|is_disease_pathway)
  R-HSA-… >>reactome>>reactomechild  -> subpathways
  R-HSA-… >>reactome>>go             -> GO term xref
The reactome entry carries only {name, tax_id} — there is no Reactome-supplied
prose definition, so the page is structural (members + hierarchy), not narrative.
"""
from dataclasses import dataclass, field

from atlas.biobtree import entry, map_all, search, rows


@dataclass
class PathwayAnchors:
    reactome_id: str
    name: str
    members: list = field(default_factory=list)     # [{ensg, symbol, biotype}]
    parent: dict = field(default_factory=dict)       # {id, name}
    children: list = field(default_factory=list)     # [{id, name}]
    go_id: str = None


def _resolve_name(name):
    """Pathway name → R-HSA id (first reactome search hit). Returns None on miss."""
    for r in rows(search(name, source="reactome")):
        if (r.get("id") or "").upper().startswith("R-HSA-"):
            return r["id"]
    return None


def resolve(id_or_name):
    rid = (id_or_name if (id_or_name or "").upper().startswith("R-HSA-")
           else _resolve_name(id_or_name))
    if not rid:
        return PathwayAnchors(reactome_id="", name=id_or_name or "")
    at = ((entry(rid, "reactome").get("Attributes") or {}).get("Reactome") or {})
    name = at.get("name") or rid

    members = [{"ensg": t["id"], "symbol": t.get("name"), "biotype": t.get("biotype")}
               for t in map_all(rid, ">>reactome>>ensembl") if t.get("id")]
    # human, protein-coding first then by symbol — deterministic
    members.sort(key=lambda m: (m.get("biotype") != "protein_coding",
                                (m.get("symbol") or m["ensg"]).upper()))

    par = map_all(rid, ">>reactome>>reactomeparent")
    parent = {"id": par[0]["id"], "name": par[0].get("name")} if par else {}
    children = [{"id": t["id"], "name": t.get("name")}
                for t in map_all(rid, ">>reactome>>reactomechild") if t.get("id")]
    go = map_all(rid, ">>reactome>>go")
    go_id = go[0]["id"] if go else None
    return PathwayAnchors(reactome_id=rid, name=name, members=members,
                          parent=parent, children=children, go_id=go_id)
