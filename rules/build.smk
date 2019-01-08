rule preprocess_load:
    message: "Preprocess load data."
    input:
        src = "src/load.py",
        raw = "data/load_time_series.xlsx"
    output: "build/model/electricity-demand.csv"
    params: year = 2015
    script: "../src/load.py"


rule preprocess_capacityfactors:
    message: "Preprocess capacityfactors of {wildcards.technology}."
    input:
        src = "src/renewables.py",
        raw = "data/res_time_series_8760h.xlsx"
    output: "build/model/capacityfactors-{technology}.csv"
    params: year = 2015
    script: "../src/renewables.py"
