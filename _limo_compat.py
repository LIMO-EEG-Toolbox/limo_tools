"""Load canonical package modules for legacy checkout scripts/imports."""

from importlib import import_module
from pathlib import Path
import sys


def load_module(name: str):
    source = str(Path(__file__).resolve().parent / "src")
    if source not in sys.path:
        sys.path.insert(0, source)
    return import_module(f"limo.{name}")
