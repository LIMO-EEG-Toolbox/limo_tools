# LInear MOdeling of MEEG data

The LInear MOdelling of MEEG data (LIMO MEEG) toolbox is a Matlab toolbox dedicated to the statistical analysis of MEEG data. Once data are imported, all computations are performed within the toolbox, and can thus work for any data sets from any software (e.g. [EEGLAB](https://sccn.ucsd.edu/eeglab/), [FieldTrip](http://www.fieldtriptoolbox.org/), [BrainStorm](https://neuroimage.usc.edu/brainstorm/)). It is interfaced with [EEGLAB](https://sccn.ucsd.edu/eeglab/) (via [STUDY](https://eeglab.org/tutorials/10_Group_analysis/working_with_study_designs.html)) acting as a plug in ; it also uses topolot from EEGLAB for result visualization. In the [LIMO_FT_integration branch](https://github.com/LIMO-EEG-Toolbox/limo_tools/tree/LIMO_FT_integration), ERP data from FieldTrip (scalp and source) can also be imported and processed.

This repository (master) is the stable version of LIMO MEEG (v3). A quick overview on how this can help you achieve a fully reproducbile workflow can see in the NeuroMatch YouTube video:  

[![BIDS-EEGLAB-LIMO-Workflow](https://github.com/LIMO-EEG-Toolbox/limo_meeg/blob/master/resources/images/nm.jpg)](https://youtu.be/fx6KIOh-jk0 "BIDS-EEGLAB-LIMO-Workflow")  

## Installation

Have EEGLAB installed (because we call some functions) and LIMO in the plug-in directory.
Data of each subject must be in subject specific folders -- ideally follow the [Brain Imaging Data Structure](https://bids.neuroimaging.io/) making working with EEGLAB/FieldTrip/LIMO easier.

LIMO 3.0 has been tested with EEGLAB 2021.0. This [test script](https://github.com/sccn/eeglab-testcases/blob/master/unittesting_limo/limo_preproc_stats_hw.m) runs nightly to make sure LIMO remains stable.

## Documentation
The [wiki](https://github.com/LIMO-EEG-Toolbox/limo_eeg/wiki) provides documentation on the various tools available and files created.  
We also have a full [tutorial](https://github.com/LIMO-EEG-Toolbox/limo_meeg/wiki) taking you through an analysis.

## Citation and method reporting
Published papers related to the method(s) used here are listed in the [citations.nbib file](https://github.com/LIMO-EEG-Toolbox/limo_tools/blob/master/citations.nbib). More generally, we recommended using [boilerplate texts from the wiki](https://github.com/LIMO-EEG-Toolbox/limo_tools/wiki/Reporting-results-differs-with-the-method-used).

## LIMO tutorial dataset
The tutorial uses data prepared using [EEG-BIDS](https://www.nature.com/articles/s41597-019-0104-8) avaialble here: https://openneuro.org/datasets/ds002718/versions/1.0.2.
There is also an older dataset that can be downloaded here: http://datashare.is.ed.ac.uk/handle/10283/2189. 

## Contribute

No brainer --> comment on anything you want (usage/doc/design) in free format on this [google doc](https://docs.google.com/document/d/1g6C4axnrJq5sItnXFbTQ0aR-3iJ_SHMXxfN55006eQ0/edit?usp=sharing)

Push any changes, submit pull request against the [HotFixes branch](https://github.com/LIMO-EEG-Toolbox/limo_tools/tree/HotFixes).

Anyone is welcome to contribute ! check here [how you can get involved](https://github.com/LIMO-EEG-Toolbox/limo_eeg/blob/master/contributing.md), the [code of conduct](https://github.com/LIMO-EEG-Toolbox/limo_eeg/blob/master/code_of_conduct.md). Contributors are listed [here](https://github.com/LIMO-EEG-Toolbox/limo_eeg/blob/master/contributors.md)
