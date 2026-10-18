# GAE extraction validation — 18 October 2026

**14/14 automated Python tests passed** on CPU (exit code 0).
The actual GCN training and prediction implementation was executed, not mocked.

Tested environment: Windows; Python 3.12.14; NumPy 2.5.4;
PyTorch 2.14.1+cpu; PyTorch Geometric 2.8.1. Dependencies are isolated in
the project-local `.gae_env` environment.

| Test | Result |
|---|---|
| Wiki beta axes and fitted-response dimensions | Passed with distinct axis-marker values |
| Output shape, fold coverage, metadata and input preservation | Passed |
| Singleton axes and explicit one-participant in-sample mode | Passed |
| Disjoint graph batching | Batched outputs match separate graphs; changing one graph does not change another |
| Target-channel masking and channel chunks | Targets hidden in their own graphs; output channel order preserved |
| Empty training mask protection | At least one masked node even at very low mask probability |
| Held-out participant leakage | Changing all beta maps of held-out participants leaves that fold's trained parameters exactly unchanged |
| Reproducibility and caller state | Repeated predictions exactly equal on tested CPU; caller Torch RNG/determinism restored |
| Participant order invariance | Exact equality when IDs and explicit fold assignments accompany the permutation |
| Learning smoke test | Average late loss below average initial loss in every fold on a predictable synthetic dataset |
| Regression against original code | Exact single-model parameter and derived-deviation agreement under matched initialization/masks/full batches |
| Optional modes and input failures | Unmasked/minibatch path finite; malformed arrays/graphs/folds/settings rejected |

Two additional tests verify automatic CUDA preference/CPU fallback and explicit
CPU/CUDA/index overrides using simulated hardware availability, plus an actual
CPU training run with the device argument omitted. GPU kernel execution has
not been tested; the local validation environment uses CPU-only PyTorch.
Requested and selected device settings are now recorded separately.

Some related assertions share a test method, giving 14 test methods overall.

The standalone synthetic example also ran with input/output shape
`[2, 12, 4, 6]`, finite predictions, and final loss below initial loss for
each fold. Its saved output is `synthetic_result.json`. The example is a
short training smoke run, not a reconstruction-quality benchmark or evidence
of artifact-detection validity. A dependency emits a TorchScript deprecation
warning during import; it did not affect execution or test results.

## Source and dimensional scope

The model/training code was extracted from the reviewed `gae_weighting.py`
snapshot at commit `2e68021088066be1f8177f9b5f108f3d8b785ff5`.
The regression test loads only the original class and training/scoring
definitions, avoiding the original module's environment/CLI side effects.

The linked LIMO wiki documents first-level beta arrays as channel × frame ×
parameter, and fitted responses as channel × frame × trial. It does not
specify fixed sample counts. Tests therefore use explicit axis markers and
small synthetic arrays; no actual participant recordings were downloaded.

## Deliberately pending

- MATLAB-to-Python conversion and full `limo_gae_weights` integration.
- Real-data model performance, GPU execution and large-array memory behavior.
- `limo_random_robust` integration, treatment of learned weights during
  resampling, and repeated-simulation null calibration.

No GitHub changes, scientific tuning, or degrees-of-freedom changes were made.
