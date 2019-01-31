subworkflow eurocalliope:
    workdir: "./euro-calliope"
    snakefile: "./euro-calliope/Snakefile"


rule copy_euro_calliope:
    message: "Copy file ./model/{wildcards.definition_file}.yaml from euro-calliope."
    input: eurocalliope("model/{definition_file}.yaml"),
    output: "build/model/{definition_file}.yaml"
    group: "construction"
    shell: "cp {input} {output}"


rule preprocess_load:
    message: "Preprocess load data."
    input:
        src = "src/construct/load.py",
        raw = "data/load_time_series.xlsx"
    output: "build/model/electricity-demand.csv"
    group: "construction"
    script: "../src/construct/load.py"


rule preprocess_capacityfactors:
    message: "Preprocess capacityfactors of {wildcards.technology}."
    input:
        src = "src/construct/renewables.py",
        raw = "data/res_time_series_8760h.xlsx"
    output: "build/model/capacityfactors-{technology}.csv"
    group: "construction"
    script: "../src/construct/renewables.py"


rule location_specific_techs:
    message: "To apply renewable targets, build location specific renewable techs."
    input:
        src = "src/construct/location_specific_techs.py",
        locations = eurocalliope("model/locations.yaml")
    output: "build/model/location-specific-techs.yaml"
    group: "construction"
    script: "../src/construct/location_specific_techs.py"


rule links:
    message: "Define links between locations based on net transfer capacities."
    input:
        src = "src/construct/links.py",
        ntc = "data/NTC.xlsx"
    output: "build/model/links.yaml"
    group: "construction"
    script: "../src/construct/links.py"


rule generation_capacities:
    message: "Define the bounds of the generation capacities."
    input:
        src = "src/construct/capacity.py",
        capacity = "data/generation_capacity.xlsx"
    output:
        csv = "build/input/capacity.csv",
        yaml = "build/model/capacity.yaml"
    group: "construction"
    script: "../src/construct/capacity.py"


rule renewable_shares:
    message: "Ensure minimal renewable shares per country."
    input:
        src = "src/construct/renewable_shares.py",
        shares = "data/bound_RES_and_CO2.xlsx",
        demand = rules.preprocess_load.output[0]
    output:
        csv = "build/input/renewable-shares.csv",
        yaml = "build/model/renewable-shares.yaml"
    group: "construction"
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
        rules.location_specific_techs.output,
        rules.links.output,
        rules.generation_capacities.output.yaml,
        rules.renewable_shares.output,
        legacy_techs = "src/template/legacy_tech.yaml",
        definition = "src/template/model.yaml"
    output:
        legacy_techs = "build/model/legacy_tech.yaml",
        model = "build/model/model.yaml"
    group: "construction"
    run:
        import jinja2

        shell("cp {input.legacy_techs} {output.legacy_techs}")
        with open(input.definition, "r") as template, open(output.model, "w") as result_file:
            result_file.write(jinja2.Template(template.read()).render(config=config))
