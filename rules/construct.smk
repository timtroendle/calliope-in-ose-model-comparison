rule preprocess_load:
    message: "Preprocess load data."
    input:
        src = "src/construct/load.py",
        raw = "data/load_time_series.xlsx"
    output: "build/model/electricity-demand.csv"
    params: year = 2015
    script: "../src/construct/load.py"


rule preprocess_capacityfactors:
    message: "Preprocess capacityfactors of {wildcards.technology}."
    input:
        src = "src/construct/renewables.py",
        raw = "data/res_time_series_8760h.xlsx"
    output: "build/model/capacityfactors-{technology}.csv"
    params: year = 2015
    script: "../src/construct/renewables.py"


rule links:
    message: "Define links between locations based on net transfer capacities."
    input:
        src = "src/construct/links.py",
        ntc = "data/NTC.xlsx"
    output: "build/model/links.yaml"
    script: "../src/construct/links.py"


rule model:
    message: "Build entire model."
    input:
        rules.preprocess_load.output,
        expand(
            "build/model/capacityfactors-{technology}.csv",
            technology=["open-field-pv", "rooftop-pv", "wind-offshore", "wind-onshore"]
        ),
        rules.links.output,
        definition = "src/construct/model.yaml"
    output: "build/model/model.yaml"
    shell: "cp {input.definition} {output}"
