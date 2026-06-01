"""Atomic file writes for published artifacts.

A crash, kill, or parallel clobber mid-write would otherwise leave a truncated
page.md / entity.jsonld / manifest.json in the published dist that an AI
crawler may fetch. Write to a temp file in the same directory, then
os.replace() — atomic on POSIX (rename within a filesystem) — so a reader
ever sees only the complete old or complete new file.
"""
import json
import os


def write_text(path: str, text: str) -> None:
    d = os.path.dirname(path) or "."
    os.makedirs(d, exist_ok=True)
    tmp = f"{path}.tmp.{os.getpid()}"
    with open(tmp, "w") as f:
        f.write(text)
    os.replace(tmp, path)


def write_json(path: str, obj, indent=0, sort_keys=False) -> None:
    write_text(path, json.dumps(obj, indent=indent, sort_keys=sort_keys))
