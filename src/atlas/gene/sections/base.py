"""Backwards-compatible re-export. The Section dataclass now lives in
atlas.section (shared between gene + disease entity pipelines) so we don't
duplicate the contract across entity types. Existing imports of
`from atlas.gene.sections.base import Section` keep working."""
from atlas.section import Section

__all__ = ["Section"]
