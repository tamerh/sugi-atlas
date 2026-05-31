"""Pluggable transports for biobtree's API.

Each transport implements three primitives — `search`, `entry`, `bbmap` —
returning identical Python-dict shapes (same as biobtree's REST JSON
schema). The dispatcher in `atlas.biobtree.client` picks one based on the
`ATLAS_BIOBTREE_TRANSPORT` env var.

Available:
  urllib       — stdlib urllib.request, no connection reuse. Default.
                 The reproducibility baseline.
  urllib_pool  — urllib3 PoolManager with HTTP keep-alive. ~15-20% faster
                 than urllib by avoiding per-call TCP setup. Same wire
                 format (JSON over HTTP/1.1), same response shape.
  grpc         — (reserved) gRPC over HTTP/2 with binary protobuf.
                 Requires regenerated stubs at atlas.biobtree._pb;
                 needs a per-dataset field-rename map to produce
                 REST-compatible dict shapes. Not yet implemented.

A/B testing:
  ATLAS_BIOBTREE_TRANSPORT=urllib       python -m atlas.disease.corpus run ...
  ATLAS_BIOBTREE_TRANSPORT=urllib_pool  python -m atlas.disease.corpus run ...

Validation: tests/integration/test_transport_parity.py compares
representative search/entry/bbmap outputs across active transports.
"""
