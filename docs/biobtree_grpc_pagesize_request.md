# Honor `mode="lite"` (or add explicit `page_size`) in gRPC Mapping/Search

## Summary

The gRPC `Mapping` and `Search` handlers ignore the request's `mode` field and
always route through the full-mode code path with `maxMappingResult = 30` /
`resultPageSize = 10`. The REST lite handlers use the much larger
`maxMappingResultLite = 150` / `resultPageSizeLite = 50`. The proto field
already exists; only the server routing is missing.

Fixing this is a small change (a few-line branch in `service/grpc.go`)
that closes the gap between gRPC and REST throughput.

## Where this lives in the code

`src/service/grpc.go:99` — `Mapping`:

```go
func (g *biobtreegrpc) Mapping(ctx context.Context, in *pbuf.MappingRequest) (*pbuf.MappingResponse, error) {
    ...
    // Note: gRPC lite mode currently returns full results
    // The REST API lite mode uses the new pipe-delimited format
    // TODO: Update protobuf definitions for new lite format if gRPC lite is needed

    // Full mode (default)
    res, err := g.service.MapFilter(in.Terms, in.Query, in.Page)
    ...
}
```

`src/service/mapfilter.go:25` — REST lite path uses the larger limit:

```go
fullResult, err := s.MapFilterWithLimit(ids, mapFilterQuery, page, s.maxMappingResultLite)
```

`src/service/service.go:127,144`:

```go
s.maxMappingResult = 30          // gRPC currently uses this
s.maxMappingResultLite = 150     // 5x larger; REST lite uses this
```

The proto already declares the mode field (`src/pbuf/app.proto`):

```proto
message MappingRequest {
    repeated string terms = 1;
    string dataset = 2;
    string query = 3;
    string page = 4;
    string mode = 5;  // "full" (default) or "lite"  ← parsed but ignored by gRPC
}
```

## Requested change

Inside `grpc.Mapping`, branch on `in.Mode`:

```go
limit := g.service.MaxMappingResult()
if in.Mode == "lite" {
    limit = g.service.MaxMappingResultLite()
}
res, err := g.service.MapFilterWithLimit(in.Terms, in.Query, in.Page, limit)
```

Same shape change for `grpc.Search` (use `resultPageSizeLite` when
`in.Mode == "lite"`).

**Important:** we are **not** asking for the pipe-encoded `ResultLite` /
`MapFilterResultLite` payload over gRPC. The full proto attribute messages
are exactly what we want — they carry strictly more information than the
compact format. Only the page-size limit needs to follow the `mode` field.

## Optional alternative

Add an explicit `page_size` (uint32) to `MappingRequest` / `SearchRequest`
so clients can pick any limit without overloading `mode`. Either approach
unblocks the use case below; we have no preference.

## Why this matters (downstream context)

We built a gRPC client that reconstructs the REST `{schema, data}` shape
locally in Python — driven entirely by biobtree's own
`conf/source{1,2}.dataset.json` (compact_fields + dataset id-to-name maps).
The adapter is ~200 LOC and survives biobtree dataset additions without
edits.

On an end-to-end Breast Cancer atlas page (~12k biobtree calls) we measure:

| transport       | wall clock | per-call cost | total calls |
| --------------- | ---------: | ------------: | ----------: |
| REST urllib     |      ~42 s |         3.3ms |      12,577 |
| REST urllib3 KA |      ~38 s |         3.0ms |      12,577 |
| gRPC (current)  |      ~38 s |         2.9ms |      13,195 |

gRPC IS measurably faster per call (binary protocol, HTTP/2 keep-alive),
but biobtree returns ~5% more pages over gRPC due to the smaller per-page
limit — which cancels the per-call savings.

With the page-size fix above, the dominant `map` path would round-trip
~5× less and the end-to-end build should drop to ~10-15s — a real 3-4×
speedup for any client doing fan-out work over biobtree.

## Acceptance / verification

A simple parity probe: same query over REST `/api/map?mode=lite` and gRPC
`Mapping(mode="lite")` returns the same number of targets per page (within
biobtree's normal post-filtering tolerance).

---

*Filed by the sugi-atlas project. Atlas adapter and parity tests live at
`src/atlas/biobtree/_transports/` and `tests/integration/test_transport_parity.py`
in the sugi-atlas repo.*
