"""pytest configuration for Atlas.

Unit tests live under tests/unit/ and exercise PURE functions only — no
biobtree, no network, no filesystem (beyond tmp_path). The integration-shaped
checks live in src/atlas/validation/ (coverage, stress, body_gate snapshot
diffs) — those need a live biobtree at 127.0.0.1:8000 and run via the
respective `python -m atlas.validation.<name>` entry points, not pytest.

Run unit tests:
    pip install -e ".[dev]"
    pytest                       # runs tests/unit/
    pytest tests/unit/test_body_gate.py
"""
