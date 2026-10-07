"""Compatibility entry point; implementation lives in :mod:`limo.limo_irls`."""

import sys

from _limo_compat import load_module

_module = load_module("limo_irls")

sys.modules[__name__] = _module
