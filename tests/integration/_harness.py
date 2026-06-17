"""Corpus-level integration checks over a BUILT dist (the dense set).

The unit suite guards the logic; this suite guards the actual rendered pages —
a release gate to run after a dense build, before committing to the full corpus.
Every check here corresponds to an invariant we rely on (the frozen page
contract, frontmatter schema, data-quality guards, mesh-link integrity,
JSON-LD validity), so a regression in any of them fails loudly against real
output.

Point it at a dist with ATLAS_INTEGRATION_DIST (default ./dist at the repo root,
the local gitignored build dir `atlas.sh` writes to). The whole suite skips
cleanly if no dist is present, so a plain `pytest` on a machine without a build
still runs the unit tests and skips these.

    pytest -m integration                      # corpus checks (needs a dist)
    pytest -m "not integration"                # unit only
    ATLAS_INTEGRATION_DIST=/tmp/x pytest -m integration
"""
import glob
import json
import os
import re

import pytest
import yaml

_REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
DIST = os.environ.get("ATLAS_INTEGRATION_DIST", os.path.join(_REPO, "dist"))
ATLAS = os.path.join(DIST, "atlas")

# The FROZEN page contract (docs/PAGE_CONTRACT.md) — anchor ids in order, per
# entity. Hardcoded here independently of the generator, so a code change that
# drifts the H2 set/order fails this gate.
H2_IDS = {
    "gene":    ["summary", "identifiers", "gene-structure", "protein",
                "function", "disease", "drugs", "related"],
    "disease": ["summary", "clinical", "identifiers", "family", "genetics",
                "genes", "function", "drugs", "trials", "related"],
    "drug":    ["summary", "identifiers", "targets", "indications",
                "pharmacology", "related-molecules", "related"],
    # Pathway is a lighter type: hierarchy (no parent for top-level terms) and
    # related (members may be off-corpus) legitimately elide, so its H2s are a
    # SUBSEQUENCE of this frozen order, not an exact match (see test_contract).
    "pathway": ["summary", "hierarchy", "member-genes", "related"],
}
ID_RE = {
    "gene":    re.compile(r"^[A-Za-z0-9][A-Za-z0-9._@\-]*$"),  # HGNC symbol (incl. @ for gene clusters, e.g. PCDHA@)
    "disease": re.compile(r"^MONDO:\d+$"),
    "drug":    re.compile(r"^CHEMBL\d+$"),
    "pathway": re.compile(r"^R-HSA-\d+$"),
}

# Frozen section-H3 id allow-list per entity (the ids render_all assigns). A
# page's H3 ids must be a SUBSET of these — a new/renamed id that drifts from
# the contract fails the gate.
H3_IDS = {
    "gene": {"gene-ids", "transcripts", "expression", "regulation",
             "functional-genomics", "generif", "orthologs", "protein-ids",
             "structure", "residue-map", "pathways", "gene-sets", "interactions", "cancer",
             "variants", "disease-assoc", "drug-data",
             "ncrna-function", "ncrna-disease", "ncrna-interactions", "ncrna-drugs",
             "hpa-protein", "hpa-cancer"},   # hpa-expression is now an H4 under §expression
    "pathway": set(),       # pathway sections are all H2 (no H3s)
    "disease": {"clinical-description", "epidemiology", "symptoms", "disease-ids", "gwas", "variant-tiers", "mendelian",
                "cohort-genes", "protein-families", "expression", "interactions",
                "structural", "pathways", "indicated", "mechanism-alignment",
                "molecular-basis", "alliance-genes", "drug-targets", "bioactivity",
                "pharmacogenomics", "tractability", "druggability", "undrugged",
                "clinical-trials"},
    "drug": {"drug-ids", "primary-targets", "bioactivity", "target-pathways",
             "indication-list", "clinical-trials", "civic", "pharmacogenomics",
             "related-mol"},
}

# Frozen H4 id allow-list per entity. H4s are the table-block titles inside a
# section (e.g. gene §drug-data). A page's H4 ids must be a SUBSET of these.
H4_IDS = {
    "gene": {  # §transcripts / §protein-ids / §structure
             "canonical-exons", "isoforms", "domains", "uniprot-cc", "brenda",
             "enzyme-ec", "rhea", "pdb", "alphafold", "sabdab",
             # §variants
             "clinvar", "top-pathogenic", "clingen-variants", "spliceai", "alphamissense", "dbsnp",
             # §pathways / §interactions / §regulation
             "reactome", "reactome-categories", "go-categories",
             "string", "intact", "signor", "corum", "cellphonedb", "interactome-enrichment",
             "collectri", "jaspar", "mirdb",
             # §expression / §disease-assoc / §ncrna
             "fantom5-promoters", "tissue-expression", "scxa", "hpa-expression",
             "cross-species-homologs", "ncrna-go",
             "gencc", "clingen", "hpo", "gwas-assoc", "efo", "mesh",
             # §drug-data
             "chembl-molecules", "civic", "civic-variants", "pharmgkb-clinical", "pharmgkb-variants",
             "pharmgkb-var-annotations",
             "pharmgkb-guidelines", "gtopdb", "bindingdb", "chembl-bioactivity",
             "pubchem-bioassay", "ctd", "chembl-assays", "cellosaurus", "gene-trials"},
    "disease": {  # clinical / genetics
              "prevalence", "hpo-features", "gwas-associations", "gwas-studies",
              "clinvar-variants", "tier-distribution", "maf-distribution",
              "consequences", "top-variants",
              # genes & proteins
              "dual-evidence", "somatic-drivers", "gencc-validity", "orphanet-linkage",
              "omim-shared", "evidence-partition", "cohort-genes-full", "cohort-function",
              "family-distribution", "family-assignment", "breadth-distribution",
              "cohort-tissues", "per-gene-tissue", "hub-genes", "intra-cohort-edges",
              "cohort-pdb", "cohort-alphafold",
              # function / therapeutics / trials
              "cohort-pathways", "go-enrichment", "approved-drug-genes", "top-targets", "cohort-drugs",
              "cohort-bioactivity", "cohort-enzymes", "screening-signal", "cohort-pgx",
              "trial-phases", "top-trials", "trial-drugs", "trial-intervention-drugs", "civic"},
    # drug pages are flat (one table per H3 section); only the two multi-table
    # sections (Target pathways, Clinical trials) carry H4s.
    "drug": {"mechanism", "target-reactome", "target-go", "trial-phases", "top-trials"},
    "pathway": set(),       # pathway sections carry no H4 table-blocks
}

_H2_LINE = re.compile(r"^## (.+?)(?:\s*\{#([a-z0-9-]+)\})?\s*$")
_H3_LINE = re.compile(r"^### (.+?)(?:\s*\{#([a-z0-9-]+)\})?\s*$")
_H4_LINE = re.compile(r"^#### (.+?)(?:\s*\{#([a-z0-9-]+)\})?\s*$")
_ANY_HEADING = re.compile(r"^(#{2,6}) (.+?)(?:\s*\{#([a-z0-9-]+)\})?\s*$")
_TABLE_ROW = re.compile(r"^\s*\|.*\|\s*$")
_CELL_SPLIT = re.compile(r"(?<!\\)\|")     # unescaped pipe = real column divider
_INTERNAL_LINK = re.compile(r"\]\((/atlas/(gene|disease|drug|pathway)/([^/)#]+)/[^)]*)\)")


def _maybe_skip():
    if not os.path.isdir(ATLAS) or not glob.glob(os.path.join(ATLAS, "*", "*", "index.md")):
        pytest.skip(f"no built dist at {ATLAS} (set ATLAS_INTEGRATION_DIST)",
                    allow_module_level=True)


class Page:
    def __init__(self, entity, slug, path):
        self.entity, self.slug, self.path = entity, slug, path
        self.raw = open(path).read()
        parts = self.raw.split("---\n", 2)
        if len(parts) >= 3:
            try:
                self.fm = yaml.safe_load(parts[1]) or {}
            except yaml.YAMLError:
                self.fm = {"__yaml_error__": True}
            self.body = parts[2]
        else:
            self.fm, self.body = {}, self.raw
        jp = os.path.join(os.path.dirname(path), "entity.jsonld")
        self.jsonld_raw = open(jp).read() if os.path.exists(jp) else None

    @property
    def h2(self):
        """[(label, id_or_None)] for every '## ' heading, in source order."""
        return [(m.group(1).strip(), m.group(2))
                for line in self.body.splitlines()
                if (m := _H2_LINE.match(line))]

    @property
    def h3(self):
        return [(m.group(1).strip(), m.group(2))
                for line in self.body.splitlines()
                if (m := _H3_LINE.match(line))]

    @property
    def h4(self):
        return [(m.group(1).strip(), m.group(2))
                for line in self.body.splitlines()
                if (m := _H4_LINE.match(line))]

    @property
    def url(self):
        return f"https://sugi.bio/atlas/{self.entity}/{self.slug}/"

    def h2_blocks(self):
        """{h2_id: block_text} — content from each '## …{#id}' to the next '## '."""
        blocks, cur_id, buf = {}, None, []
        for line in self.body.splitlines():
            m = _H2_LINE.match(line)
            if m:
                if cur_id is not None:
                    blocks[cur_id] = "\n".join(buf)
                cur_id, buf = m.group(2), []
            else:
                buf.append(line)
        if cur_id is not None:
            blocks[cur_id] = "\n".join(buf)
        return blocks

    def tables(self):
        """List of tables; each a list of rows, each row a list of cells. Splits
        on UNESCAPED '|' only (a cell may legitimately contain '\\|')."""
        out, cur = [], []
        for line in self.body.splitlines():
            if _TABLE_ROW.match(line):
                cells = _CELL_SPLIT.split(line.strip().strip("|"))
                cur.append(cells)
            elif cur:
                out.append(cur)
                cur = []
        if cur:
            out.append(cur)
        return out

    def jsonld(self):
        return json.loads(self.jsonld_raw) if self.jsonld_raw else None

    def inline_jsonld(self):
        tag = '<script type="application/ld+json">'
        if tag not in self.body:
            return None
        return json.loads(self.body.split(tag, 1)[1].split("</script>", 1)[0].strip())


def _load():
    pages = []
    for path in glob.glob(os.path.join(ATLAS, "*", "*", "index.md")):
        entity = os.path.basename(os.path.dirname(os.path.dirname(path)))
        slug = os.path.basename(os.path.dirname(path))
        pages.append(Page(entity, slug, path))
    return pages


_URL_PARTS = re.compile(r"/atlas/(gene|disease|drug|pathway)/([^/)#]+)/")


def page_exists(url):
    """True if an /atlas/<entity>/<slug>/ url has a page dir in the dist."""
    m = _URL_PARTS.search(url or "")
    return bool(m) and os.path.isdir(os.path.join(ATLAS, m.group(1), m.group(2)))


def norm(s):
    """Lowercase + whitespace-collapse, for label/title comparison."""
    return re.sub(r"\s+", " ", (s or "").strip().lower())


def report(violations, limit=25):
    """A readable assertion message: count + the first `limit` offenders."""
    n = len(violations)
    head = "\n".join(f"  - {v}" for v in violations[:limit])
    more = f"\n  … +{n - limit} more" if n > limit else ""
    return f"{n} violation(s):\n{head}{more}"
