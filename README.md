# LInear MOdeling of MEEG data

The LInear MOdelling of MEEG data (LIMO MEEG) toolbox is a toolbox dedicated to the statistical analysis of MEEG data. 

This is the experimental Python branch. The installable distribution is
`limo-eeg`; its Python import namespace is `limo`. Python 3.10 or newer is
required. This backend provides the import, design, GLM and contrast stages
needed for a future EEGPrep STUDY integration.

## Installation

From a checkout of this branch:

```sh
python -m pip install .
```

For development and regression checks:

```sh
python -m pip install -e ".[test]"
python -m pytest -q
python -m build
```

`dist/` contains the wheel and source distribution after building. The runtime
dependencies are NumPy, SciPy and h5py. Optional extras are `[mne]` for MNE
sensor adjacency and `[plot]` for diagnostic figures. No PyPI release is
assumed by these instructions.

## First-level analysis

The source dataset must be epoched. Categorical and continuous regressors have
one row per trial, in the source dataset's trial order. They can be NumPy arrays
or paths to numeric `.txt`/`.mat` files. The design stage reorders trials and
regressors together. For a two-condition analysis:

```python
from limo.eeglab_import import export_limo_h5
from limo.limo_design import limo_design, read_hdf5_structure
from limo.limo_glm import limo_glm
from limo.limo_contrast import limo_contrast

model = export_limo_h5(
    "sub-01.set",
    cat="conditions.txt",  # trial labels 1 or 2
    defaults={
        "name": "results/sub-01",
        "analysis": "Time",
        "method": "WLS",
        "bootstrap": 0,
        "tfce": 0,
    },
)
model, results = limo_design(model)
limo_glm(model)
# Inspect the design before defining contrasts: two condition columns
# followed by the intercept in this example.
X = read_hdf5_structure(model)["LIMO"]["design"]["X"]
_, _, _, dataset = limo_contrast(model, {"C": [1, -1, 0], "F": 0})
contrast_result = read_hdf5_structure(results)[dataset]
```

`method` accepts `OLS`, `WLS` and `IRLS`. WLS computes trial weights separately
for each channel; IRLS computes weights separately for each channel and frame.
`analysis` accepts `Time`, `Frequency` and `Time-Frequency`. The latter two
read precomputed `datspec` or `dattimef` channel measures referenced by the SET
file. Complex time-frequency coefficients are converted to power once.

Outputs are `LIMO.h5` (metadata/design) and `limo_results.h5` (data, fit,
effects and contrasts). Time/frequency arrays retain channel × frame × trial
axes; time-frequency arrays retain channel × frequency × time × trial axes.
Beta arrays replace the last axis with parameters. The five values in a
first-level t-contrast's final axis are estimate, standard error, degrees of
freedom, t statistic and p-value.

The same workflow is available through installed commands. Save the defaults
above in `defaults.json`, then run:

```sh
limo-import sub-01.set --defaults defaults.json --cat conditions.txt
limo-design results/sub-01/LIMO.h5
limo-glm results/sub-01/LIMO.h5
limo-contrast results/sub-01/LIMO.h5 --contrast "1,-1,0" --test T
```

Use `--help` on any command for its options. `limo-tfce` and
`limo-inspect-set` are also installed. Module execution, for example
`python -m limo.limo_glm --help`, is supported. Existing root-script invocations
and imports continue to work from a source checkout via compatibility wrappers;
the wheel exposes the `limo` namespace and installed commands.

## Validation and current limits

The regression suite checks OLS against independent least-squares and
t-contrast calculations, WLS betas against a pinned MATLAB reference, channel
weight independence, missing-observation subsets, singleton axes, SET/FDT
sample order, and full file workflows for time, spectrum and time-frequency
analyses across OLS/WLS/IRLS. Time-frequency effect assembly is also checked
against separate frequency fits. The MATLAB WLS reference is recorded in
`tests/test_glm.py` with its source commit.

The backend remains experimental. Synthetic tests and one WLS reference do
not establish full numerical parity with MATLAB LIMO. Full real-data validation,
bootstrap/TFCE and clustering parity, native MATLAB v7.3 input coverage,
second-level inference, and clustered-ICA workflows need further work. Clustered
ICA design construction explicitly raises `NotImplementedError`.

The output files use the Python backend's HDF5 schema; they are not MATLAB
v7.3 MAT files or drop-in replacements for the MATLAB LIMO result-file layout.
An EEGPrep STUDY adapter and a dependency-aware execution/restart runner are
separate integration steps. Building the design initializes its result file;
avoid repeating that stage over a completed analysis. The GLM can reopen a
completed model without refitting it.

## Documentation
The [wiki](https://github.com/LIMO-EEG-Toolbox/limo_eeg/wiki) provides documentation on the various tools available and files created.  
We also have a full [tutorial](https://github.com/LIMO-EEG-Toolbox/limo_meeg/wiki) taking you through an analysis.

## Citation and method reporting
Published papers related to the method(s) used here are listed in the [citations.nbib file](https://github.com/LIMO-EEG-Toolbox/limo_tools/blob/master/citations.nbib). More generally, we recommended using [boilerplate texts from the wiki](https://github.com/LIMO-EEG-Toolbox/limo_tools/wiki/Reporting-results-differs-with-the-method-used).

## LIMO tutorial dataset
The tutorial uses data prepared using [EEG-BIDS](https://www.nature.com/articles/s41597-019-0104-8) avaialble here: https://openneuro.org/datasets/ds002718/versions/1.0.2.
There is also an older dataset that can be downloaded here: http://datashare.is.ed.ac.uk/handle/10283/2189. 

# LIMO versions

4.1 / 4.1.0 - Updated data workflow and BIDS compliance

4.1.1 - Updated GUI issue for contrast, 2025 figure compatibility

## Contribute

No brainer --> comment on anything you want (usage/doc/design) in free format on this [google doc](https://docs.google.com/document/d/1g6C4axnrJq5sItnXFbTQ0aR-3iJ_SHMXxfN55006eQ0/edit?usp=sharing)

Submit Python backend pull requests against the
[python branch](https://github.com/LIMO-EEG-Toolbox/limo_tools/tree/python).
MATLAB fixes use the [HotFixes branch](https://github.com/LIMO-EEG-Toolbox/limo_tools/tree/HotFixes).

Anyone is welcome to contribute ! check here [how you can get involved](https://github.com/LIMO-EEG-Toolbox/limo_eeg/blob/master/contributing.md), the [code of conduct](https://github.com/LIMO-EEG-Toolbox/limo_eeg/blob/master/code_of_conduct.md). Contributors are listed [here](https://github.com/LIMO-EEG-Toolbox/limo_eeg/blob/master/contributors.md)
