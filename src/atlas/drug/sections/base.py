"""Backwards-compatible re-export — the Section dataclass lives in
atlas.section (shared across gene / disease / drug entity pipelines)."""
from atlas.section import Section

__all__ = ["Section"]
