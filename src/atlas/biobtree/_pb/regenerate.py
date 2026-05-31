#!/usr/bin/env python3
"""Regenerate Python protobuf stubs from biobtree's proto definitions.

Pulls from /data/biobtree/src/pbuf/{app,attr}.proto and emits stubs into
this directory:
    attr_pb2.py, app_pb2.py, app_pb2_grpc.py

Run from project root: `python -m atlas.biobtree._pb.regenerate`
Or directly:           `python src/atlas/biobtree/_pb/regenerate.py`

When biobtree adds a new dataset attribute or RPC, re-run this. The
generated files are committed (so an end-user clone works without
regenerating); this script exists so we can refresh deterministically.

Adapted from bioyoda/sugi-agent/modules/integrations/biobtree_pb/regenerate_protobuf.py.
Two import rewrites needed because protoc emits `import attr_pb2` (top-level)
while we import the modules as a package; `from . import attr_pb2` is the
correct form for relative-package import.
"""
import subprocess
import sys
from pathlib import Path

PROTO_SRC = Path("/data/biobtree/src/pbuf")
OUTPUT_DIR = Path(__file__).parent.resolve()


def regenerate():
    attr_proto = PROTO_SRC / "attr.proto"
    app_proto = PROTO_SRC / "app.proto"
    if not attr_proto.exists() or not app_proto.exists():
        print(f"ERROR: proto files not found under {PROTO_SRC}", file=sys.stderr)
        sys.exit(1)

    print(f"Proto source: {PROTO_SRC}")
    print(f"Output dir:   {OUTPUT_DIR}")

    cmd = [
        sys.executable, "-m", "grpc_tools.protoc",
        f"-I{PROTO_SRC}",
        f"--python_out={OUTPUT_DIR}",
        f"--grpc_python_out={OUTPUT_DIR}",
        str(attr_proto), str(app_proto),
    ]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        print(f"protoc failed:\n{r.stderr}", file=sys.stderr)
        sys.exit(1)

    # Rewrite imports to package-relative form.
    for filepath, old, new in [
        (OUTPUT_DIR / "app_pb2.py",       "import attr_pb2", "from . import attr_pb2"),
        (OUTPUT_DIR / "app_pb2_grpc.py",  "import app_pb2",  "from . import app_pb2"),
    ]:
        text = filepath.read_text()
        if old in text and new not in text:
            filepath.write_text(text.replace(old, new))
            print(f"  rewrote imports in {filepath.name}")

    print("done.")


if __name__ == "__main__":
    regenerate()
