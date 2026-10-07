import importlib
from pathlib import Path
import subprocess
import sys
import os

import pytest


@pytest.mark.parametrize("module", [
    "eeglab_import", "limo_design", "limo_glm", "limo_contrast",
    "limo_tfce", "read_setfile", "limo_WLS", "limo_irls", "limo_clustering",
])
def test_package_modules_import(module):
    assert importlib.import_module(f"limo.{module}") is not None


@pytest.mark.parametrize("module", [
    "eeglab_import", "limo_design", "limo_glm", "limo_contrast", "limo_tfce", "read_setfile",
])
def test_module_and_legacy_cli_help(module):
    repo = Path(__file__).resolve().parents[1]
    commands = {
        "eeglab_import": "limo-import", "limo_design": "limo-design",
        "limo_glm": "limo-glm", "limo_contrast": "limo-contrast",
        "limo_tfce": "limo-tfce", "read_setfile": "limo-inspect-set",
    }
    executable = Path(sys.executable).parent / (commands[module] + (".exe" if os.name == "nt" else ""))
    for args in [[sys.executable, "-m", f"limo.{module}", "--help"],
                 [sys.executable, str(repo / f"{module}.py"), "--help"],
                 [str(executable), "--help"]]:
        result = subprocess.run(args, cwd=repo, capture_output=True, text=True)
        assert result.returncode == 0, result.stderr
        assert "usage:" in result.stdout
