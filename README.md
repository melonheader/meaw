[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.8042349.svg)](https://doi.org/10.5281/zenodo.8042349)
# meaw

Assembly of scripts to process MEA electrophysiology recordings.

## Downloading the source
First, one needs to download the source to a desired path. Example usage: 
 ```bash
 git clone https://github.com/melonheader/meaw.git
 cd meaw
 ```
To perform the default analysis, copy all the MEA readings in .csv format directly to the ./meaw.
Make sure to write div in capitals within the filenames of MEA readings.
Down the way, the instuctions are written inside of the analysis notebooks. 
1. Process.Rmd
2. prepareResults.RMD

First notebook processes raw MEA readings and assembles summary files with MEA features, including MFR (mean firing rate), Burst rate and Network spike rate.
Second reads in the processed results and prepares some default plots for data visualization.
