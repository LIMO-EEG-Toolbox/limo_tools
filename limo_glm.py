"""Compatibility entry point; implementation lives in :mod:`limo.limo_glm`."""

import sys

from _limo_compat import load_module

_module = load_module("limo_glm")

if __name__ == "__main__":
    _module.main()
else:
    sys.modules[__name__] = _module
