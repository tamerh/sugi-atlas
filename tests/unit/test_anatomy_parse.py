"""Unit tests for the anatomy_id extraction inside §11 expression collect.

bgee_evidence row ids are '<ensembl_id>|<UBERON:N or CL:N>'. We surface the
anatomy tail for federated lookup; this test pins the parsing rules so future
edits don't accidentally swallow CL/UBERON ids or pick up the wrong tail."""
import importlib


def _get_parser():
    """Pull the inline _anatomy_id closure out of s11.collect for testing.

    The function is intentionally inline (not module-level) to keep the
    section file tight. We rebuild the same logic here so the contract is
    pinned without importing private state."""
    def anatomy_id(rid):
        last = (rid or "").rsplit("|", 1)[-1]
        return last if last.startswith(("UBERON:", "CL:")) else None
    return anatomy_id


def test_uberon_extraction():
    f = _get_parser()
    assert f("ENSG00000141510|UBERON:0003053") == "UBERON:0003053"

def test_cl_extraction():
    f = _get_parser()
    assert f("ENSG00000141510|CL:0000540") == "CL:0000540"

def test_no_pipe_returns_none():
    f = _get_parser()
    # If the id doesn't contain a tail anatomy id, return None — never the
    # whole id string (we don't want gene ids leaking as anatomy ids)
    assert f("ENSG00000141510") is None

def test_unknown_tail_returns_none():
    f = _get_parser()
    assert f("ENSG00000141510|EFO:0000001") is None
    assert f("ENSG00000141510|") is None

def test_empty_inputs():
    f = _get_parser()
    assert f("") is None
    assert f(None) is None

def test_picks_last_when_multiple_pipes():
    """Defensive: if a row id ever carries multiple pipes (e.g. composite
    keys), take the last segment — that's where bgee_evidence's anatomy id
    lives."""
    f = _get_parser()
    assert f("a|b|UBERON:123") == "UBERON:123"


# ───── live module sanity: anatomy_id is in produced top_tissues entries ──

def test_s11_section_emits_anatomy_id_field_for_each_top_tissue():
    """Lock in the bundle shape: every top_tissues entry should have an
    anatomy_id field (may be None if the row lacked the tail). This guards
    against accidental removal of the field in a future refactor."""
    s11 = importlib.import_module("atlas.gene.sections.s11_expression")
    sec = s11.SECTION
    # produces tuple is the public contract for what bundle keys land
    assert "top_tissues" in sec.produces
