"""gRPC transport — talks to biobtree on :7776 with HTTP/2 keep-alive.

Wire layer only: calls Search/Mapping/Entry RPCs, returns the proto-native
dict via MessageToDict(preserving_proto_field_name=True). The dispatcher's
rows()/map_targets() in client.py detect the shape and adapt — see
_grpc_adapter.py for the REST-shape projection (uses biobtree's own
compact_fields config as the schema source of truth).

Why MessageToDict instead of working with proto messages directly:
- collectors are written against dicts (REST shape), not protobuf objects
- one boundary for the shape adapter, not a sprinkling of `.HasField`/.GetX()
  through every collector
- the cost (proto -> dict serialization) is small vs the network call

gRPC channel options match bioyoda's working setup: keep-alive pings + no
ping cap so the channel survives the per-disease build (~11k calls over
several minutes).
"""
import grpc
from google.protobuf.json_format import MessageToDict

from atlas.biobtree._pb import app_pb2, app_pb2_grpc
from atlas.biobtree._transports import _grpc_adapter as _adapter

# Shared CALLS log (dispatcher imports this).
CALLS = []
API = "127.0.0.1:7776"

_OPTIONS = [
    ("grpc.keepalive_time_ms", 30_000),
    ("grpc.keepalive_timeout_ms", 10_000),
    ("grpc.keepalive_permit_without_calls", 1),
    ("grpc.http2.max_pings_without_data", 0),
    # 32 MiB — biobtree mapping responses can grow large (string_interaction etc.)
    ("grpc.max_receive_message_length", 32 * 1024 * 1024),
]

# One channel + stub, module-level. Created lazily on first call so importing
# this module doesn't connect.
_channel = None
_stub = None


def _get_stub():
    global _channel, _stub
    if _stub is None:
        _channel = grpc.insecure_channel(API, options=_OPTIONS)
        _stub = app_pb2_grpc.BiobtreeServiceStub(_channel)
    return _stub


def _to_dict(msg) -> dict:
    return MessageToDict(msg, preserving_proto_field_name=True)


def search(term: str, source: str = None) -> dict:
    req = app_pb2.SearchRequest(terms=[term])
    if source:
        req.dataset = source
    CALLS.append({"path": "search", "params": {"i": term, **({"s": source} if source else {})}})
    return _adapter.search_to_rest(_to_dict(_get_stub().Search(req)))


def entry(identifier: str, source: str) -> dict:
    req = app_pb2.EntryRequest(identifier=identifier, dataset=source)
    CALLS.append({"path": "entry", "params": {"i": identifier, "s": source}})
    return _adapter.entry_to_rest(_to_dict(_get_stub().Entry(req)))


def bbmap(ids: str, chain: str, page: str = None) -> dict:
    # REST sends a comma-joined string for `i`; the proto is `repeated string terms`.
    terms = [t for t in ids.split(",") if t]
    req = app_pb2.MappingRequest(terms=terms, query=chain)
    if page:
        req.page = page
    CALLS.append({"path": "map", "params": {"i": ids, "m": chain, **({"p": page} if page else {})}})
    return _adapter.mapping_to_rest(_to_dict(_get_stub().Mapping(req)))
