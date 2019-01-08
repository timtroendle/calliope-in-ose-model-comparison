# OSE Model Comparison -- Euro Calliope

Simulation runs of Euro Calliope within the OSE Model Comparison.

This repository contains the entire scientific project, including code and report. The philosophy behind this repository is that no intermediary results are included, but all results are computed from raw data and code.

## Getting ready

You need [conda](https://conda.io/docs/index.html) to run the analysis. Using conda, you can create a conda environment from within you can run it:

    conda env create -f requirements.yml

Also, you need to manually add the following data:

* Copy Euro Calliope model files into `./data/model/` # TODO should be linked.
* Copy the load time series into `./data/load_time_series.xlsx` # TODO check in?
* Copy the 8760 steps renewables time series into `./data/res_time_series_8760h.xlsx` # TODO check in?
* Copy the NTC data into `./data/NTC.xlsx` # TODO check in?

## Run the analysis

    snakemake

This will run all analysis steps to reproduce results and eventually build the report.

You can also run certain parts only by using other `snakemake` rules; to get a list of all rules run `snakemake --list`.

To generate a PDF of the dependency graph of all steps, and if you have `dot` installed, run:

    snakemake --rulegraph | dot -Tpdf > dag.pdf

## Run the tests

    snakemake test

## Repo structure

* `report`: contains all files necessary to build the report; plots and result files are not in here but generated automatically
* `src`: contains the Python source code
* `tests`: contains the test code
* `config`: configurations used in the study
* `data`: place for raw data
* `build`: will contain all results (does not exist initially)
