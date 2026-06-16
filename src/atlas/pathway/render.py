#!/usr/bin/env python3
"""Pathway page body — deterministic markdown from the collect bundle. Structural
(members + hierarchy), no narrative: Reactome supplies no prose definition.

Member genes link to their Atlas gene pages (manifest-gated); parent/children
link to their Atlas pathway pages. The reverse edge — "pathways containing this
gene" on each gene page — is produced by the mesh (links.related_targets), not
here.
"""
from atlas.page import links
from atlas.render_common import table, capped_table

ROW_CAP = 100


def _reactome_url(rid):
    return f"https://reactome.org/content/detail/{rid}" if rid else None


def declarative_sentence(b):
    """One-line lead: name, id, member count, parent, Reactome link."""
    name = b.get("name") or b.get("reactome_id") or "?"
    rid = b.get("reactome_id")
    head = f"**{name}**" + (f" (`{rid}`)" if rid else "")
    n = b.get("member_count") or 0
    s = f"{head} is a Reactome pathway with {n:,} member gene{'s' if n != 1 else ''} (human)"
    par = b.get("parent") or {}
    if par.get("name"):
        s += f"; a subpathway of {links.maybe_link(par['name'], links.pathway_url(reactome_id=par.get('id'), name=par.get('name')))}"
    nc = b.get("child_count") or 0
    if nc:
        s += f", itself comprising {nc} subpathway{'s' if nc != 1 else ''}"
    s += "."
    return s


def r_hierarchy(b):
    par = b.get("parent") or {}
    children = b.get("children") or []
    if not par.get("name") and not children:
        return ""
    L = ["## Pathway hierarchy {#hierarchy}", ""]
    if par.get("name"):
        L.append(f"**Parent pathway:** "
                 + links.maybe_link(par["name"], links.pathway_url(reactome_id=par.get("id"), name=par.get("name"))))
    if children:
        links_ = ", ".join(links.maybe_link(c.get("name"),
                                            links.pathway_url(reactome_id=c.get("id"), name=c.get("name")))
                           for c in children)
        L.append(f"\n**Subpathways ({len(children)}):** {links_}")
    return "\n".join(L)


def r_members(b):
    members = b.get("members") or []
    if not members:
        return ""
    L = [f"## Member genes ({b.get('member_count', len(members))}) {{#member-genes}}", "",
         "Genes whose products participate in this pathway (Reactome → Ensembl).", ""]
    L.append(capped_table(["Gene", "Ensembl", "Biotype"],
                          [(links.maybe_link(m.get("symbol") or m["ensg"], links.gene_url(symbol=m.get("symbol"))),
                            m["ensg"], (m.get("biotype") or "").replace("_", " "))
                           for m in members],
                          ROW_CAP, total=b.get("member_count"), noun="member genes"))
    return "\n".join(L)


def summary_block(b):
    """The Summary-section content: declarative lead + GO + Reactome source. Used
    by assemble_page (pathway branch) as the `## Summary` body, mirroring the
    gene/disease/drug declarative lead."""
    parts = [declarative_sentence(b)]
    go = b.get("go_id")
    if go:
        # Id in a code span, link the word (never the bare ontology id as a label).
        parts.append(f"**Gene Ontology:** `{go}` "
                     f"([QuickGO](https://www.ebi.ac.uk/QuickGO/term/{go}))")
    rid = b.get("reactome_id")
    if rid:
        parts.append(f"*Source: [Reactome {rid}]({_reactome_url(rid)}).*")
    return "\n\n".join(parts)


def body(b):
    """The data sections (hierarchy + members) — the Summary is added separately
    by assemble_page so the page matches the gene/disease/drug shape."""
    out = []
    for r in (r_hierarchy(b), r_members(b)):
        if r:
            out.append(r)
    return "\n\n".join(out)


def render_all(b):
    """Full standalone body (summary + sections) — used by the standalone
    collect→render path and tests. The batch path uses summary_block()/body()."""
    return summary_block(b) + ("\n\n" + body(b) if body(b) else "")
