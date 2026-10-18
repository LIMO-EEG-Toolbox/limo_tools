# Shared-beta GAE: focused extraction

This is the existing newer GAE, refactored to return predictions for
`stage1/limo_gae_weights.m`. It does not calculate weights, run inference,
read MATLAB files, or modify LIMO. A single public entry point is provided:

```python
from limo_gae_predict import limo_gae_predict

result = limo_gae_predict(
    betas, adjacency,
    folds=4, masked=True, epochs=60, seed=0,  # GPU if available, otherwise CPU
)
reconstructed = result["reconstructed"]
prediction_info = result["prediction_info"]
```

GPU selection is automatic by default (`device="auto"`). Override it with
`device="cpu"` to force CPU, `device="cuda"` to require the default CUDA GPU,
or `device="cuda:1"` to select a particular GPU. An explicit unavailable GPU
raises an error; only automatic selection falls back to CPU. Both requested
and selected devices are recorded in the returned settings.

CUDA use requires a CUDA-enabled PyTorch installation and a compatible NVIDIA
driver/GPU. The isolated environment used for these tests contains CPU-only
PyTorch, so its automatic selection uses CPU. The synthetic test example
intentionally forces CPU for reproducible reference results.

## Dimensions from the LIMO wiki

Source: [information and dimensions about files on drive](https://github.com/LIMO-EEG-Toolbox/limo_tools/wiki/information-and-dimensions-about-files-on-drive), read 18 October 2026.

| Data | Dimensions |
|---|---|
| One participant's beta array | channel × frame × parameter |
| Stacked participant beta arrays | channel × frame × parameter × participant |
| GAE input and output | parameter × frame × channel × participant |
| One participant's design X | trial × parameter |
| One participant's Yr / Yhat | channel × frame × trial |
| Adjacency | channel × channel |

Explicit conversion after stacking aligned beta arrays:

```python
# Each item must retain channel x frame x parameter, including singleton axes.
stacked = np.stack(subject_betas, axis=3)
betas = stacked.transpose(2, 1, 0, 3)
# Convert predictions back to the stacked on-disk axis convention if needed:
stacked_predictions = reconstructed.transpose(2, 1, 0, 3)
```

MATLAB equivalent: `B = permute(cat(4, subjectBetas{:}), [3 2 1 4]);`.
Do not squeeze singleton dimensions. Python requires an explicitly 4D array;
the future MATLAB bridge must preserve trailing singleton axes during conversion.
Channel identities, frame grids, parameter meanings and participant ordering
must be aligned by the caller before stacking. The GAE validates sizes, not
scientific alignment. The wiki supplies axis conventions, not fixed numerical
sizes; tests use small synthetic 2 × 7 × 4 × 6 arrays and singleton examples.
Bootstrap replicates and repeated-measure axes must not be treated as independent
participants. Time-frequency reshaping is outside this stage.

## Preserved behavior and targeted changes

Preserved from `gae_weighting.py` at commit
`2e68021088066be1f8177f9b5f108f3d8b785ff5`: one model shared across all betas
within each participant fold; GCNConv encoder; linear/ReLU decoder; SmoothL1
loss; Adam; channel masking; default 4 folds, 60 epochs, 32 hidden units,
one GCN layer, dropout 0.3, mask probability 0.3, learning rate 0.001 and
weight decay 0.0001. `batch_size=None` retains full-batch optimization.

Changes:

- Return predictions in the same axis order as input instead of reducing to deviations.
- Guarantee a nonempty training mask in each batch (source could compute an empty loss).
- Reject nonfinite inputs, float32 overflow, incompatible adjacency and invalid folds/options.
- Record explicit folds, participant identities, per-fold seeds and loss histories.
- Seed each fold separately and restore the caller's PyTorch random/determinism state.
- Optionally process masked channel predictions in chunks to limit memory use.
- Remove the source's global `KMP_DUPLICATE_LIB_OK` environment override.

These changes do not promise byte-for-byte agreement with a complete old
multi-fold run: the new per-fold seeding differs from the old single global
random stream. A regression test does verify identical single-model training
and derived deviations with equal initialization, masks, and full-batch inputs.

## Options and outputs

- `prediction`: `"out-of-fold"` (default), or `"in-sample"` for explicit diagnostics.
- `folds`: 4 by default, 2 through participant count. Stride assignment matches
  the source; the seed does not shuffle folds.
- `fold_ids`: optional nonnegative integer label per participant, at least two
  unique labels; overrides `folds`. Folds are never split within participants.
- `participant_ids`: unique nonempty strings. Supply these with explicit folds
  for order-invariant reruns; training participants are sorted by these IDs.
- `masked`: true by default. Each returned target channel is predicted with
  that channel hidden. Such a map combines distinct masked predictions; it is
  not one simultaneous unmasked reconstruction.
- `epochs`, `hidden_dim`, `depth`, `dropout`, `mask_probability`, `learning_rate`,
  `weight_decay`: explicit training settings with source defaults above.
- `batch_size`: None (full batch) or positive number of beta maps per step.
- `channel_batch_size`: 16 (number of separate target-channel masks per inference batch).
- `seed`: 0; each fold gets seed + its ordinal, recorded in metadata.
- `device`: `"auto"` (CUDA if available), `"cpu"`, or a CUDA device string.
- `deterministic`: true; unsupported deterministic device operations raise
  rather than silently switching to nondeterminism. CPU testing is the baseline.
- `return_models`: false; if true, return CPU models keyed by fold label.

`reconstructed` is float32 with exactly the input shape. `prediction_info`
contains layout, original participant ordering, zero-based train/test indices,
fold assignments, per-fold training-map counts and seeds, loss histories,
resolved settings and dependency versions. It is JSON-serializable. Identical
CPU runs in the tested environment are reproducible; cross-device/version bitwise
equivalence is not promised. The seeded context modifies global PyTorch state
temporarily, so do not run simultaneous training calls in threads sharing one process.

The binary symmetric adjacency may contain diagonal values; they are removed
and GCNConv adds its self-loops. No learned preprocessing or early stopping is
performed. Held-out participants never contribute training beta maps. A masked
prediction still uses the other channels of that same held-out participant.

For repeated sessions from one person, combine them appropriately before calling
or assign all related entries to the same explicit fold. The function cannot
infer shared biological identity from distinct input participant IDs.

## Run and review

From `gae_revised/`, after creating an environment and installing the requirements:

```powershell
& '.\.gae_env\Scripts\python.exe' stage2/tests/test_limo_gae_predict.py
& '.\.gae_env\Scripts\python.exe' stage2/example_synthetic.py
```

For another environment, install `requirements.txt`; select a CPU PyTorch wheel
from the official PyTorch CPU index when CUDA is not required. Exact versions
used in validation are recorded in `VALIDATION.md`.

The tests cover wiki axis mapping, shapes/singletons, disjoint graph batching,
target masking, participant holdout, unchanged held-out-fold model parameters
when test data change, seed repeatability, stable-ID ordering, finite learning,
source regression, optional unmasked/minibatch paths, and rejected inputs.

Next, after review: implement the MATLAB adapter, feed these predictions to
`limo_gae_weights`, then integrate with `limo_random_robust` and assess null
calibration. No degrees-of-freedom change or new inferential procedure is
introduced here. No real EEG data or scientific accuracy claims are involved
in the synthetic tests.
