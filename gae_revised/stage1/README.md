# GAE post-processing — stage 1

`limo_gae_weights.m` is a new MATLAB implementation. It accepts original and
reconstructed arrays in **beta × frame × channel × subject** order and returns
errors, scores, weights, flags, reference statistics and provenance. It does
not train a model, read/write data, change the working directory, alter the
random stream, or run statistical inference. No additional MATLAB toolboxes
are required. MATLAB R2017b or later is the intended baseline; the numerical
suite has been executed successfully with GNU Octave 11.3.0 on Windows.

```matlab
addpath('stage1');
R = limo_gae_weights(B, Brec, 'BetaIndices', [1 2], 'Scope', 'both');

% Optional normalized copy; raw quality weights remain available:
R = limo_gae_weights(B, Brec, 'Normalization', 'mean-one');

% Fitted-response discrepancy from subject-specific trial x beta designs:
R = limo_gae_weights(B, Brec, 'ErrorDomain', 'fitted', ...
    'DesignMatrices', X, 'Scope', 'both');
```

If the existing Python GAE returns beta × subject × channel × frame,
transpose it with `permute(double(reconstructed),[1 4 3 2])` before calling.
The caller must align subject, channel, frame and beta identities. Metadata
is descriptive: there is no silent matching, reordering, or inference of axes.

## Options

| Option | Default | Alternatives / meaning |
|---|---|---|
| ErrorDomain | beta | fitted: calculate X*(B−Brec), then reduce trials and frames |
| DesignMatrices | empty | Required cell array, one finite trial × original-beta matrix per subject, in fitted mode only |
| BetaIndices, FrameIndices, ChannelIndices | all | Unique in-range integer indices applied identically to both inputs |
| Metric | mae | rmse; evaluated before aggregation across betas |
| BetaAggregation | mean | median, none; fitted mode accepts only the default mean |
| Scope | channel | subject, both; global and local weights are never multiplied together |
| Standardization | median-mad | mean-sd (sample SD), none (threshold is then in error units) |
| Reference | within-channel | pooled: calibrate each beta/aggregate over all selected channel–subject cells |
| WeightFunction | bisquare | powerlaw, none |
| Threshold | 3 | Nonnegative attenuation onset; not a significance threshold |
| TuningConstant | 4.685 | Positive bisquare taper width; explicit use with another mapping is rejected |
| Normalization | none | mean-one across valid subjects in each output row |
| MissingPolicy | error | omit: nonfinite input pairs are marked missing |
| MinValidObservations | 1 | Minimum finite frames (or trial–frame cells in fitted mode) per error |
| ReferenceMask | all subjects | Logical vector selecting reference subjects |
| MinReferenceSubjects | 3 | Minimum distinct subjects with valid reference values; at least 2 allowed |
| ScaleTolerance | 1e-12 | Relative to maximum absolute reference value |
| SubjectIDs, ChannelIDs | empty | Unique numeric/string/cellstr vectors matching full input axes |
| FrameValues | empty | Finite numeric vector matching the full frame axis |
| PredictionInfo | scalar struct | Caller-supplied model/fold/seed/input identifiers; not inferred or verified |
| ReturnReconstruction | false | Return selected reconstruction; invalid input pairs are represented as NaN |

These are development defaults, not validated artifact-detection settings.
Select the relevant betas explicitly: averaging coefficients of different
scales or meanings needs a scientific justification. MissingPolicy='omit'
also needs care: subjects can then be scored over different valid frames.
Set MinValidObservations to the required coverage rather than relying on its
technical minimum. Missing data in any selected beta invalidates the beta
aggregate; subject summaries require all selected channels to be valid.

## Numerical definition

Channel errors are standardized using median and 1.4826*MAD by default.
If this scale is degenerate, the fallback is IQR/1.349, using linear quantiles
at position 1+(n−1)*p (explicitly implemented to avoid toolbox dependencies).
An entirely constant reference and identical scored values produce z=0.
A constant reference with differing scored values, or nonconstant data with
no usable MAD/IQR, produces NaN scores and `degenerate-scale` status. Too few
reference subjects produces `insufficient-reference`. No tiny absolute
epsilon is added. Standardization='none' bypasses reference calibration.

For bisquare, e=max(z−Threshold,0), q=(1−(e/TuningConstant)^2)^2 for
e<TuningConstant, otherwise q=0. Powerlaw uses q=1/(1+e^2). Mapping 'none'
returns one for valid scores. Raw q is in [0,1]. Normalized W=q*nvalid/sum(q)
can exceed one. All-zero valid rows raise an explicit error. Invalid rows
remain NaN; they are never silently replaced by unit weights.

Subject scores are computed by averaging calibrated channel scores, then
standardizing that summary across reference subjects using the same method.
Raw subject errors are separately reported in the original error units.
With BetaAggregation='none', all processing is independent per beta.

## Outputs

- `errors.by_beta_channel`: beta × channel × subject in beta mode, empty in fitted mode.
- `errors.channel`, `scores.channel`: channel × subject, or beta × channel × subject when retaining betas.
- `errors.subject`: 1 × subject, or beta × subject; descriptive mean of channel errors.
- `weights.raw_channel`, `weights.channel`: bounded raw and optionally normalized weights, when requested.
- `scores.subject`, `weights.raw_subject`, `weights.subject`: subject outputs, when requested.
- `flags.channel/subject.attenuated`: score exceeds threshold (false with mapping 'none'); independent of normalization.
- `flags.channel/subject.zero_weight`: raw weight is zero. Invalid values have false flags; inspect `valid` too.
- `valid.channel/subject`: finite calibrated scores; `observation_count` and `selected_input` document coverage.
- `reference.channel/subject`: cells with location, scale, status, subject count and value count.
- `diagnostics`: invalid channel-score count and effective n=(sum q)^2/sum(q.^2), for requested scopes.
- `options`: resolved settings, including supplied design matrices.
- `provenance`: input shape, selected metadata, output layout, version and caller prediction information.
- `reconstructed`: optional selected reconstruction.

Effective n is a descriptive measure of weight concentration. It is **not**
an instruction to replace LIMO's degrees of freedom. Stage 2 retains
`limo_random_robust` as the intended inference framework and evaluates the
weighted bootstrap, including the role of learning weights from the cohort.

## Basic testing

From MATLAB or Octave, with this folder as the current directory:

```matlab
addpath(pwd);
addpath(fullfile(pwd,'tests'));
test_limo_gae_weights
```

The assertion suite has 12 groups: identity/singletons, hand-computed MAE
and RMSE, fitted-space multiplication and coefficient selection, beta
aggregation, known weight curves, normalization and ESS, ordering and
scale invariance, degenerate references, missing values, subject scope,
input rejection, and absence of I/O/random-state side effects.

`tests/check_numerical_references.py` independently verifies five arithmetic
fixtures with NumPy. It does not execute MATLAB and cannot establish that
the MATLAB implementation passes its tests. See VALIDATION.md for actual
execution status.
