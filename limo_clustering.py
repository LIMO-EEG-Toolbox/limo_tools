"""Compatibility entry point; implementation lives in :mod:`limo.limo_clustering`."""

import sys

from _limo_compat import load_module

_module = load_module("limo_clustering")

sys.modules[__name__] = _module
