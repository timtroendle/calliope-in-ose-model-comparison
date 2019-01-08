rule run:
    message: "Run the model."
    input:
        model = "src/simple-model.yaml",
        load = rules.preprocess_load.output,
        capacityfactors = expand(
            "build/model/capacityfactors-{technology}.csv",
            technology=["open-field-pv", "rooftop-pv", "wind-offshore", "wind-onshore"]
        )
    output: "build/results.nc"
    shell: "calliope run {input.model} --save_netcdf {output} --scenario=diw_assumptions"


rule plot:
    message: "Visualises the results."
    input:
        src = "src/analyse/vis.py",
        results = rules.run.output
    output: "build/plot.png"
    script: "../src/analyse/vis.py"


rule capacity:
    message: "Calculate the installed capacities."
    input:
        src = "src/analyse/capacity.py",
        results = rules.run.output
    output: "build/capacity.csv"
    script: "../src/analyse/capacity.py"
