# Code to reproduce analysis and figures for "The 4.2 ka event is not remarkable in the context of Holocene climate variability"

This is the repository for the code used to produce the results and figures for the manuscript "The 4.2 ka event is not remarkable in the context of Holocene climate variability" by McKay, Nicholas P.; Kaufman, Darrell S.; Arcusa, Stéphanie; Kolus, Hannah; Edge, David; Erb, Michael P.; Hancock, Chris; Routson, Cody C.; Żarczyński, Maurycy; Marshall, Leah; Roberts, Georgia; Telles, Frank. This manuscript is currently in review. 

To run the key analyses and figures, start with the file "Code/workflow.R"

We tried to anonymize/standardize the paths, but they will likely take some tinkering to get to work on your local machine

## System requirements
This code should run on Mac, Windows and Linux operating systems using R v4.1.0 and above. It has been tested on multiple machines and versions of R that meet these criteria, although the most exhaustive testing has been conducted on MacOS.

## Installation guide
After installing R, install the required packages as listed in the software using `install.packages()` for CRAN packages, and `remotes::install_github()` for packages hosted on Github.
Typical install time should be less than 10 minutes.

## Usage
The file `workflow.R` shows the workflow through the script and how to reproduce the analyses and figures in the manuscript. This is a computationally heavy analysis, and expected run time on a typical desktop or laptop computer is several days. There will be progress bars.

The actR software that underlies most of this analysis can be used to perform similar analyses on other datasets. See [the documentation](https://linked.earth/actR) for details.
