#!/usr/bin/env python3
"""Atlas CLI — placeholder entry point.

Today: thin shims, mostly redirecting to per-module __main__ blocks. As the
pipeline solidifies this will grow into a proper subcommand dispatcher used by
the Enju workflows.

  python -m atlas.cli gene collect TP53
  python -m atlas.cli gene render  TP53 all
"""
import sys
from atlas.gene import collect, render

USAGE = "usage: python -m atlas.cli gene {collect,render} <symbol> [args...]"

def main(argv=None):
    a = argv or sys.argv[1:]
    if len(a) < 3 or a[0] != "gene":
        print(USAGE); sys.exit(2)
    _, op, *rest = a
    sys.argv = [op] + rest  # delegate
    if op == "collect": collect.main() if hasattr(collect, "main") else _no_main("collect")
    elif op == "render": render.main() if hasattr(render, "main") else _no_main("render")
    else: print(USAGE); sys.exit(2)

def _no_main(name):
    print(f"atlas.gene.{name} has no main() yet — run with `python -m atlas.gene.{name}` for now")
    sys.exit(2)

if __name__ == "__main__":
    main()
