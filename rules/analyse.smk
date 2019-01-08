rule run:
    message: "Run the model."
    input:
        model = rules.model.output
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
    message: "Excavate the installed capacities."
    input:
        src = "src/analyse/capacity.py",
        results = rules.run.output
    output: "build/capacity.csv"
    script: "../src/analyse/capacity.py"


rule trade:
    message: "Excavate the traded electricity."
    input:
        src = "src/analyse/trade.py",
        results = rules.run.output
    output: "build/trade.csv"
    script: "../src/analyse/trade.py"
