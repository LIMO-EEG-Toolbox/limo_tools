"""Compatibility entry point; implementation lives in :mod:`limo.limo_WLS`."""

import sys

from _limo_compat import load_module

_module = load_module("limo_WLS")

sys.modules[__name__] = _module
