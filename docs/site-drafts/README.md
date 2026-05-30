# Site-file drafts (pre-launch)

These files will eventually be deployed to the **`biobtree-content`** Hugo
site (the place Atlas's generated pages get served from). They live here as
drafts so we can iterate on them alongside the codebase, then copy them into
`biobtree-content/static/` manually right before launch.

## Files

- `llms.txt` — site orientation for LLM agents (Anthropic / Perplexity read this).
- `robots.txt` — bot allowlist, with the AI fetchers explicitly named so they
  don't fall back to the wildcard rule.

## Why not deploy now

We're still building substantial content (drugs, diseases, more biobtree
datasets, per-fact provenance anchors, Hugo theme work for JSON-LD discovery
hints). Sending the AI bots toward the site too early means they crawl an
incomplete corpus and cache that view. Launch the site files alongside the
finished V1 corpus.

## Pre-launch workflow

1. Re-read both files when V1 content is complete (gene corpus done, exec
   summary models locked in, Hugo theme updated).
2. Confirm the file is consistent with `docs/research/NEXT.md` (URLs, scope,
   citation guidance).
3. Copy into `biobtree-content/static/llms.txt` and
   `biobtree-content/static/robots.txt`.
4. Verify with: `curl -I https://sugi.bio/llms.txt` and
   `curl -I https://sugi.bio/robots.txt`.
5. After deploy, watch `chatgpt_report.py` for new AI-fetch patterns.

## Related work still pending in biobtree-content

Tracked under "Pre-launch / cross-repo work" in
[`docs/research/NEXT.md`](../research/NEXT.md):

- Hugo theme adds `<link rel="alternate" type="application/ld+json" href="entity.jsonld">`
  to the page head (machine clients then discover the sidecar without
  scraping the body).
- Hugo `config.toml`: enable per-page `<lastmod>` in the sitemap (drives
  the Perplexity freshness signal).
- Hugo serves `Last-Modified` HTTP header from page `generated_at` frontmatter.
- nginx logrotate review (already extended to ~3 years).
