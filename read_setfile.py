"""Compatibility entry point; implementation lives in :mod:`limo.read_setfile`."""

import sys

from _limo_compat import load_module

_module = load_module("read_setfile")

if __name__ == "__main__":
    _module.main()
else:
    sys.modules[__name__] = _module
