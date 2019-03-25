# OSE Model Comparison -- Euro Calliope

Simulation runs of Euro Calliope within the OSE Model Comparison.

This repository contains the entire scientific project, including code and report. The philosophy behind this repository is that no intermediary results are included, but all results are computed from raw data and code.

## Getting ready

1. Clone the repo. Because `euro-calliope` is added as a git submodule, you may want to clone using `git clone --recurse-submodules <link-to-this-repo>`.

2. Create an environment to run the analysis. You need [conda](https://conda.io/docs/index.html) to run the analysis. Using conda, you can create a conda environment from within you can run it:

    `conda env create -f environment.yaml`

3. Make sure you have a Gurobi license, or install and configure another solver.

4. Provide the input data for Euro-Calliope, as defined in "Getting Ready" in  `./euro-calliope/README.md`.

5. Provide the following data manually:

* Copy the load time series into `./data/load_time_series.xlsx` # TODO check in?
* Copy the 8760 steps renewables time series into `./data/res_time_series_8760h.xlsx` # TODO check in?
* Copy the NTC data into `./data/NTC.xlsx` # TODO check in?
* Copy the generation capacities into `./data/generation_capacity.xlsx` # TODO check in?
* Copy the RES and CO2 bounds into `./data/bound_RES_and_CO2.xlsx` # TODO check in?

## Run the analysis

    snakemake --use-conda

This will run all analysis steps to reproduce results and eventually build the report.

You can also run certain parts only by using other `snakemake` rules; to get a list of all rules run `snakemake --list`.

To generate a PDF of the dependency graph of all steps, and if you have `dot` installed, run:

    snakemake --rulegraph | dot -Tpdf > dag.pdf

## Run on Euler cluster

To run on Euler, use the following command:

    snakemake --use-conda --profile config/euler [--config email=<you@provider.org>]

By providing an email address, you will be informed by mail when Snakemake finishes execution.

If you want to run on another cluster, read [snakemake's documentation on cluster execution](https://snakemake.readthedocs.io/en/stable/executable.html#cluster-execution) and take `config/euler` as a starting point.

## Run the tests

    snakemake test --use-conda

## Units and scaling

The default units within the optimisation model are `MW`, `MWh`, `EUR`, `Mt`, and `km2`, but you can scale all of these using the configuration values in `config/default.yaml`. Apart from convenience, this may be important to handle numerical issues with your solver.

## Repo structure

* `report`: contains all files necessary to build the report; plots and result files are not in here but generated automatically
* `src`: contains the Python source code
* `tests`: contains the test code
* `config`: configurations used in the study
* `data`: place for raw data
* `build`: will contain all results (does not exist initially)
