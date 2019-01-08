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
    input: "data/model/capacityfactors-{technology}.csv"
    output: "build/model/capacityfactors-{technology}.csv"
    shell: "cp {input} {output}"
