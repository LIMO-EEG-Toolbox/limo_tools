# Stage 1 validation — 18 October 2026

## Executed successfully

- **Octave execution: 12/12 test groups passed** on GNU Octave 11.3.0
  (Windows x64), exit code 0. This executed the actual MATLAB-language
  implementation, including the missing-data and error-path assertions.
- **Independent numerical references: 5/5 checks passed.** NumPy verified
  hand-computed MAE/RMSE, design-matrix multiplication, bisquare boundary and
  midpoint values, robust-score normalization/effective n, and zero-MAD
  examples with and without a usable IQR fallback.
- **MATLAB syntax parsing: 2/2 files passed.** MISS_HIT Core 0.9.44 parsed
  `limo_gae_weights.m` and `tests/test_limo_gae_weights.m` using MATLAB 2017b
  syntax. No syntax errors were reported. Style rules were disabled.

## Execution details and remaining limitation

The MATLAB assertion suite contains 12 test groups, including errors,
weights, output dimensions, fitted-space calculations, pooled/sample-SD
calibration, missing values, permutations, scaling, validation errors, and
side-effect checks. The first Octave run passed 11 groups but the directory
side-effect assertion compared Octave-specific access-time metadata. The
test now compares names, types, sizes and modification times, excluding
`.` and `..`. No change to the numerical implementation was needed. The
full suite then passed all 12 groups.

Execution in proprietary MATLAB remains pending; Octave execution does not
guarantee MATLAB runtime compatibility. No training, statistical-null
calibration, or scientific performance claims are made by these tests.

Run the MATLAB/Octave tests using the instructions in README.md.

## Review notes

- No existing GAE code or synced reference file was modified.
- No changes were made to GitHub, model training, or `limo_random_robust`.
- The earlier report was amended to retain LIMO's bootstrap framework for
  stage 2, without prescribing an automatic degrees-of-freedom adjustment.
- Degenerate calibration returns NaN plus a status, rather than invented
  finite weights. Callers must inspect validity/status before using outputs.
