# LInear MOdeling of MEEG data

The LInear MOdelling of MEEG data (LIMO MEEG) toolbox is a Matlab toolbox dedicated to the statistical analysis of MEEG data. It has some  interfacing with EEGLAB (in particular the STUDY in the EEGLAB develop version) to act as a plug in. However, once data are imported all is performed within LIMO MEEG and the toolbox can thus work for any data sets.

This repo is the stable version of LIMO MEEG (v3) to be used with EEGLAB (https://sccn.ucsd.edu/eeglab/) but can be used with in other applications like FieldTrip (http://www.fieldtriptoolbox.org/) or BrainStorm (https://neuroimage.usc.edu/brainstorm/) for your research applications. A quick overview on how this can help you achieve a fully reproducbile workflow can see in the NeuroMatch YouTube video:  

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

Push any changes, submit pull request against the [HotFixes branch](https://github.com/LIMO-EEG-Toolbox/limo_tools/tree/HotFixes).

Anyone is welcome to contribute ! check here [how you can get involved](https://github.com/LIMO-EEG-Toolbox/limo_eeg/blob/master/contributing.md), the [code of conduct](https://github.com/LIMO-EEG-Toolbox/limo_eeg/blob/master/code_of_conduct.md). Contributors are listed [here](https://github.com/LIMO-EEG-Toolbox/limo_eeg/blob/master/contributors.md)
