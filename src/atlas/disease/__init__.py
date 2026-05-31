"""Atlas disease-entity pipeline (collector + render + workflow).

## Slug stability — caller-supplied override

By default `slugify(canonical_name)` derives the URL-safe slug from Mondo's
canonical_name field. If Mondo recurates ("endometrial cancer" → "endometrial
carcinoma") Atlas would silently change the dist subdir and the published
URL — bad for stable links + bad for body_gate snapshot keys (a renamed slug
looks like a fresh `first_run`).

Mitigation: callers can pin the slug. Both the standalone pipeline and the
Enju workflow accept it:

  # Pipeline (standalone)
  from atlas.pipeline import run_disease
  # Caller passes the disease record manually — slug derives from
  # canonical_name today; future enhancement = add a `slug` argument that
  # overrides. (Open follow-up; not yet wired through run_disease().)

  # Enju workflow — the `diseases` parameter is a list<record> with both
  # `name` and `slug` fields. Pass the desired slug explicitly:
  diseases:
    - {name: "Endometrial Cancer", slug: "endometrial-cancer"}    # pin
    - {name: "Bipolar Disorder"}                                  # derive

  The task scripts (bin/run_disease_backlog.py, src/atlas/disease/tasks/*)
  honor an explicit `slug` field when present; fall back to slugify() when
  absent.

When you change a slug after a snapshot exists, the new slug looks like a
first_run because the snapshot is keyed on the old slug name. Either rename
the snapshot file (snapshots/disease/<old>.json → <new>.json) or accept the
first_run verdict + refresh.
"""
