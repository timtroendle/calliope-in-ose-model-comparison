subworkflow eurocalliope:
    workdir: "./euro-calliope"
    snakefile: "./euro-calliope/Snakefile"
    configfile: "./config/default.yaml"

configfile: "./config/default.yaml"

localrules: copy_euro_calliope, copy_national_definitions_from_euro_calliope, model
wildcard_constraints:
    definition_file = "((renewable-techs)|(locations)|(hydro-capacities)|"
                      "(energy-generation-hydro-ror)|(energy-inflow-hydro-reservoir)|"
                      "(storage-techs))"

rule copy_euro_calliope:
    message: "Copy file ./build/model/{wildcards.definition_file}.yaml from euro-calliope."
    input: eurocalliope("build/model/{definition_file}.yaml"),
    output: "build/model/{definition_file}.yaml"
    shell: "cp {input} {output}"


rule copy_national_definitions_from_euro_calliope:
    message: "Copy file {input[0]} from euro-calliope."
    input: eurocalliope("build/model/national/{definition_file}.{suffix}"),
    output: "build/model/national/{definition_file}.{suffix}"
    shell: "cp {input} {output}"


rule preprocess_load:
    message: "Preprocess load data."
    input:
        src = "src/construct/load.py",
        raw = "data/load_time_series.xlsx"
    params: scaling_factor = config["scaling-factors"]["power"]
    output: "build/model/national/electricity-demand.csv"
    conda: "../envs/default.yaml"
    script: "../src/construct/load.py"


rule preprocess_capacityfactors:
    message: "Preprocess capacityfactors of {wildcards.technology}."
    input:
        src = "src/construct/renewables.py",
        raw = "data/res_time_series_8760h.xlsx"
    output: "build/model/national/capacityfactors-{technology}.csv"
    conda: "../envs/default.yaml"
    script: "../src/construct/renewables.py"


rule links:
    message: "Define links between locations based on net transfer capacities."
    input:
        src = "src/construct/links.py",
        ntc = "data/NTC.xlsx"
    params: scaling_factor = config["scaling-factors"]["power"]
    output: "build/model/national/links.yaml"
    conda: "../envs/default.yaml"
    script: "../src/construct/links.py"


rule generation_capacities:
    message: "Define the bounds of the generation capacities."
    input:
        src = "src/construct/capacity.py",
        capacity = "data/generation_capacity.xlsx"
    params: scaling_factor = config["scaling-factors"]["power"]
    output:
        csv = "build/input/capacity.csv",
        yaml = "build/model/national/capacity.yaml"
    conda: "../envs/default.yaml"
    script: "../src/construct/capacity.py"


rule renewable_shares:
    message: "Ensure minimal renewable shares per country."
    input:
        src = "src/construct/renewable_shares.py",
        shares = "data/bound_RES_and_CO2.xlsx"
    output:
        csv = "build/input/renewable-shares.csv",
        yaml = "build/model/national/renewable-shares.yaml"
    conda: "../envs/default.yaml"
    script: "../src/construct/renewable_shares.py"


rule co2_caps:
    message: "Ensure CO2 bounds in each country."
    input:
        src = "src/construct/co2_caps.py",
        caps = "data/bound_RES_and_CO2.xlsx"
    params: scaling_factor = config["scaling-factors"]["co2"]
    output:
        csv = "build/input/co2-caps.csv",
        yaml = "build/model/national/co2-caps.yaml"
    conda: "../envs/default.yaml"
    script: "../src/construct/co2_caps.py"


rule model:
    message: "Build entire model."
    input:
        "build/model/renewable-techs.yaml",
        "build/model/national/locations.yaml",
        "build/model/national/hydro-capacities.yaml",
        "build/model/national/energy-generation-hydro-ror.csv",
        "build/model/national/energy-inflow-hydro-reservoir.csv",
        "build/model/storage-techs.yaml",
        rules.preprocess_load.output,
        expand(
            "build/model/national/capacityfactors-{technology}.csv",
            technology=["open-field-pv", "rooftop-pv", "wind-offshore", "wind-onshore"]
        ),
        rules.links.output,
        rules.generation_capacities.output.yaml,
        rules.renewable_shares.output,
        rules.co2_caps.output,
        legacy_techs = "src/template/legacy-tech.yaml",
        definition = "src/template/model.yaml"
    params:
        from_date = config["from_date"],
        to_date = config["to_date"],
        resolution = config["resolution"],
        power_scaling_factor = config["scaling-factors"]["power"],
        monetary_scaling_factor = config["scaling-factors"]["monetary"],
        co2_scaling_factor = config["scaling-factors"]["co2"]
    output:
        legacy_techs = "build/model/legacy-tech.yaml",
        model = "build/model/national/model.yaml"
    run:
        import jinja2

        specific_monetary_scaling_factor = params.monetary_scaling_factor / params.power_scaling_factor
        specific_co2_scaling_factor = params.co2_scaling_factor / params.power_scaling_factor

        with open(input.legacy_techs, "r") as template, open(output.legacy_techs, "w") as result_file:
            result_file.write(jinja2.Template(template.read()).render(
                monetary_scaling_factor=specific_monetary_scaling_factor,
                co2_scaling_factor=specific_co2_scaling_factor,
            ))

        with open(input.definition, "r") as template, open(output.model, "w") as result_file:
            result_file.write(jinja2.Template(template.read()).render(
                monetary_scaling_factor=specific_monetary_scaling_factor,
                co2_scaling_factor=specific_co2_scaling_factor,
                from_date=params.from_date,
                to_date=params.to_date,
                resolution=params.resolution
            ))
