rule run:
    message: "Run the model."
    input:
        model = rules.model.output.model
    output: "build/output/results.nc"
    shell:
        """
        calliope run {input.model} --save_netcdf {output} --scenario=diw_assumptions,location_specific_techs
        echo -e "import calliope\nif calliope.read_netcdf('{output}').results.termination_condition != 'optimal':\n raise ValueError('non optimal')" \
        | python # see https://github.com/calliope-project/calliope/issues/182
        """


rule plot:
    message: "Visualises the results."
    input:
        src = "src/analyse/vis.py",
        results = rules.run.output
    output: "build/output/plot.png"
    script: "../src/analyse/vis.py"


rule capacity:
    message: "Excavate the installed capacities."
    input:
        src = "src/analyse/capacity.py",
        results = rules.run.output
    output:
        raw = "build/output/capacity-raw.csv",
        publish = "build/output/capacity-publish.csv"
    script: "../src/analyse/capacity.py"


rule trade:
    message: "Excavate the traded electricity."
    input:
        src = "src/analyse/trade.py",
        results = rules.run.output
    output: "build/output/trade.csv"
    script: "../src/analyse/trade.py"


rule test:
    message: "Run tests"
    input:
        "tests/test_capacity_constraints.py",
        rules.capacity.output.raw,
        rules.renewable_shares.output,
        rules.run.output
    output: "build/test-report.html"
    shell:
        "py.test ./tests/ --html={output} --self-contained-html"
