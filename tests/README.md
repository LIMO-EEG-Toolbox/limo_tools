# Command line regression tests

Initialize a compatible EEGLAB installation and its installed dependencies in a fresh MATLAB process, then run this file explicitly:

```matlab
addpath(eeglabFolder);
eeglab nogui;
results = runtests(fullfile(limoFolder, 'tests', 'test_command_line.m'));
assertSuccess(results);
```

Validated locally with MATLAB R2025b, Statistics and Machine Learning Toolbox, and Parallel Computing Toolbox. A temporary setting disables automatic pool creation so the small TFCE bootstrap tests execute locally without launching workers. It does not disable the TFCE calculations or change the user's saved parallel preferences. The tests use a serial PSOM configuration to avoid concurrent batch workers. Full integration tests should also exercise the normal parallel configuration.

The suite creates deterministic synthetic data for 18 subjects in temporary directories. It checks caller STUDY preservation, standalone contrasts, current and legacy result filenames, filename container normalization, paired subject alignment, supplied channel vectors, subject and group ordering, and the existing trimmed mean numerical convention. Time and time frequency TFCE bootstrap outputs are compared with direct calls to the unchanged TFCE implementation.

Test fixtures replace input dialogs with functions that throw an error. An unexpected dialog is a failure, not an automatically accepted answer. These fixtures are installed only by the test suite. **Do not add `tests/fixtures` recursively to a production MATLAB path.** Tests restore their temporary path fixtures, working directory, random state, warnings, figure visibility and base STUDY binding.

The small suite does not replace integration testing with real EEG data. The complementary `sccn/eeglab_tests` LIMO workflows exercise the Wakeman and Henson dataset, OLS and WLS first level models, contrasts, all seven second level analysis groups, bootstrap calculations and TFCE.
