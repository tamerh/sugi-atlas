"""Gene §SpliceAI: collapse the saturated positional grid into splice-region ranges
(audit: 30 adjacent Δ=1.0000 rows conveyed nothing)."""
from atlas.gene.render import _collapse_spliceai


def test_consecutive_positions_collapse_to_a_range():
    rows = [
        {"id": "17:7670604", "effect": "donor_loss", "score": "1.0000"},
        {"id": "17:7670605", "effect": "donor_loss", "score": "1.0000"},
        {"id": "17:7670606", "effect": "donor_loss", "score": "1.0000"},
        {"id": "5:100", "effect": "acceptor_gain", "score": "0.90"},   # scattered single
        {"id": "weird-id", "effect": "x", "score": "0.5"},             # unparseable
    ]
    out = _collapse_spliceai(rows)
    assert ("17:7670604–7670606", "donor_loss", "1.0000", "3") in out  # run → one range, 3 sites
    assert ("5:100", "acceptor_gain", "0.90", "") in out                     # single stays single
    assert ("weird-id", "x", "0.5", "") in out                               # loose row passes through


def test_same_positions_different_effect_do_not_merge():
    rows = [
        {"id": "1:50", "effect": "donor_gain", "score": "0.8"},
        {"id": "1:51", "effect": "acceptor_loss", "score": "0.8"},   # different effect → separate
    ]
    out = _collapse_spliceai(rows)
    labels = {r[0] for r in out}
    assert labels == {"1:50", "1:51"}            # not merged into a 1:50–51 range
