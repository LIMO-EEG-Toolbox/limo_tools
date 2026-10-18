# Revised GAE implementation

This development branch adds a self-contained implementation under `gae_revised/`.
It does not yet replace the original LIMO pipeline.

- [Stage 1: MATLAB post-processing](stage1/README.md): errors, scores, weights and diagnostics.
- [Stage 2: Python GAE](stage2/README.md): shared-beta, participant-held-out masked predictions.
- Validation: [MATLAB/Octave](stage1/VALIDATION.md) and [Python](stage2/VALIDATION.md).

## Python setup and tests

From this directory, create an isolated environment:

```powershell
python -m venv .gae_env
& '.\.gae_env\Scripts\python.exe' -m pip install -r stage2/requirements.txt
& '.\.gae_env\Scripts\python.exe' stage2/tests/test_limo_gae_predict.py
& '.\.gae_env\Scripts\python.exe' stage2/example_synthetic.py
```

Select a CUDA-enabled PyTorch build appropriate for the hardware using the
[official PyTorch installer selector](https://pytorch.org/get-started/locally/)
when GPU support is required. The function defaults to CUDA when available,
otherwise CPU; users can override with `device="cpu"` or `device="cuda:N"`.

The frozen file in `review/gae_source_snapshot/` is included solely for the
regression test against the original implementation at commit
`2e68021088066be1f8177f9b5f108f3d8b785ff5`.

## MATLAB or Octave tests

From this directory:

```matlab
addpath('stage1');
addpath(fullfile('stage1','tests'));
test_limo_gae_weights
```

Validation completed: 12 Octave test groups and 14 Python tests passed.
GPU execution and proprietary MATLAB runtime execution remain unverified.
The MATLAB/Python bridge, integration with `limo_random_robust`, and null
calibration remain subsequent steps. No automatic degrees-of-freedom change
is introduced.
