"""Compatibility entry point; implementation lives in :mod:`limo.eeglab_import`."""

import sys

from _limo_compat import load_module

_module = load_module("eeglab_import")

if __name__ == "__main__":
    _module.main()
else:
    sys.modules[__name__] = _module
