"""Sugi Atlas's biobtree client now lives in the shared `sugibiobtree` package
(single source of truth, also used by Sugi Variant). This module aliases to it,
so every `atlas.biobtree` / `atlas.biobtree.client` import — including the
transport test that reaches into client internals — resolves to the one client.
"""
import sys
from sugibiobtree import client as _client
sys.modules[__name__] = _client
