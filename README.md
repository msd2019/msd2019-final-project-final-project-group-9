# Final Project for Modeling Social Data, 2019
George Austin (gia2105), Calvin Tong (cyt2113), Mia Fryer (mzf2106)

This repository has code to attempt to replicate and extend the results in Systematic Inequality and Hierarchy in Faculty Hiring Networks by Aaron Clauset, Samuel Arbesman, Daniel B. Larremore, Science Advances 12 Feb 2015, Vol. 1, No. 1

A complete report of our results is in `02_final_report.pdf`, which can be generated by cloning this repository and running `make` to execute the commands in the `Makefile` file. All data are in `data/original` and any original source code provided by the authors is in `authors_original_code/`.

To redownload the data run `make doanload_data`. The data will be downloaded to the `data/original` folder. The data can be found at the URL `http://tuvalu.santafe.edu/~aaronc/facultyhiring/replicationData_all.zip`

The `authors_original_code` directory contains a single matlab file `mvrsample.m`.
This was the only code the authors provided and preforms Markov Chain Monte
Carlo to sample the minimum violation rankings (MVRs) of the network. We did not
reproduce this portion of the paper, but the code is provided for reference.

An author created interactive visualization of the networks can be found [here](http://danlarremore.com/faculty/)


The repository is structured as follows:

1. `01_get_original_data.sh` gets the original data used by the authors and places a copy in `data/`
2. `02_final_report.Rmd` analysis both the original and new data to replicate and extend the results of the original paper, and produces the final report `02_final_report.pdf`

