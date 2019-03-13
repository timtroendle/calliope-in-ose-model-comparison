subworkflow eurocalliope:
    workdir: "./euro-calliope"
    snakefile: "./euro-calliope/Snakefile"


localrules: copy_euro_calliope, model, pumped_hydro, raw_run_of_river_data_zipped

URL_RUNOFF_DATA_SWITZERLAND = "https://data.sccer-jasm.ch/runofriver_production/"\
                              "opsd-runofriver_production-2017-10-16.zip"


rule raw_run_of_river_data_zipped:
    message: "Download run of river generation data as zip."
    output:
        protected("data/automatic/raw-run-of-river.zip")
    shell:
        "curl -sLo {output} '{URL_RUNOFF_DATA_SWITZERLAND}'"


rule run_of_river_capacity_factors:
    message: "Create capacity factor timeseries from Swiss generation data."
    input:
        src = "src/construct/runoff.py",
        zip = rules.raw_run_of_river_data_zipped.output,
    shadow: "minimal"
    output: "build/model/run-of-river.csv"
    script: "../src/construct/runoff.py"


rule copy_euro_calliope:
    message: "Copy file ./build/model/{wildcards.definition_file}.yaml from euro-calliope."
    input: eurocalliope("build/model/{definition_file}.yaml"),
    output: "build/model/{definition_file}.yaml"
    shell: "cp {input} {output}"


rule copy_national_locations_from_euro_calliope:
    message: "Copy file locations.yaml from euro-calliope."
    input: eurocalliope("build/model/national/locations.yaml"),
    output: "build/model/locations.yaml"
    shell: "cp {input} {output}"


rule preprocess_load:
    message: "Preprocess load data."
    input:
        src = "src/construct/load.py",
        raw = "data/load_time_series.xlsx"
    output: "build/model/electricity-demand.csv"
    script: "../src/construct/load.py"


rule preprocess_capacityfactors:
    message: "Preprocess capacityfactors of {wildcards.technology}."
    input:
        src = "src/construct/renewables.py",
        raw = "data/res_time_series_8760h.xlsx"
    output: "build/model/capacityfactors-{technology}.csv"
    script: "../src/construct/renewables.py"


rule links:
    message: "Define links between locations based on net transfer capacities."
    input:
        src = "src/construct/links.py",
        ntc = "data/NTC.xlsx"
    output: "build/model/links.yaml"
    script: "../src/construct/links.py"


rule generation_capacities:
    message: "Define the bounds of the generation capacities."
    input:
        src = "src/construct/capacity.py",
        capacity = "data/generation_capacity.xlsx"
    output:
        csv = "build/input/capacity.csv",
        yaml = "build/model/capacity.yaml"
    script: "../src/construct/capacity.py"


rule pumped_hydro:
    message: "Define the charge rates of pumped hydro."
    input:
        src = "src/construct/pumped_hydro.py"
    output:
        yaml = "build/model/pumped-hydro.yaml"
    script: "../src/construct/pumped_hydro.py"


rule renewable_shares:
    message: "Ensure minimal renewable shares per country."
    input:
        src = "src/construct/renewable_shares.py",
        shares = "data/bound_RES_and_CO2.xlsx"
    output:
        csv = "build/input/renewable-shares.csv",
        yaml = "build/model/renewable-shares.yaml"
    script: "../src/construct/renewable_shares.py"


rule model:
    message: "Build entire model."
    input:
        "build/model/interest-rate.yaml",
        "build/model/renewable-techs.yaml",
        "build/model/locations.yaml",
        "build/model/link-techs.yaml",
        "build/model/storage-techs.yaml",
        rules.preprocess_load.output,
        expand(
            "build/model/capacityfactors-{technology}.csv",
            technology=["open-field-pv", "rooftop-pv", "wind-offshore", "wind-onshore"]
        ),
        rules.links.output,
        rules.generation_capacities.output.yaml,
        rules.pumped_hydro.output,
        rules.renewable_shares.output,
        rules.run_of_river_capacity_factors.output,
        legacy_techs = "src/template/legacy-tech.yaml",
        definition = "src/template/model.yaml"
    output:
        legacy_techs = "build/model/legacy-tech.yaml",
        model = "build/model/model.yaml"
    run:
        import jinja2

        shell("cp {input.legacy_techs} {output.legacy_techs}")
        with open(input.definition, "r") as template, open(output.model, "w") as result_file:
            result_file.write(jinja2.Template(template.read()).render(config=config))
