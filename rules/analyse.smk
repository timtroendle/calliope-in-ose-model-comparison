configfile: "config/default.yaml"

rule run:
    message: "Run the model for scenario {wildcards.scenario}."
    input:
        model = rules.model.output.model
    output: "build/output/{scenario}/results.nc"
    shell:
        """
        calliope run {input.model} --save_netcdf {output} --scenario={wildcards.scenario}
        echo -e "import calliope\nif calliope.read_netcdf('{output}').results.termination_condition != 'optimal':\n raise ValueError('non optimal')" \
        | python # see https://github.com/calliope-project/calliope/issues/182
        """


rule plot:
    message: "Visualises the results of scenario {wildcards.scenario}."
    input:
        src = "src/analyse/vis.py",
        results = rules.run.output
    output: "build/output/{scenario}/plot.png"
    script: "../src/analyse/vis.py"


rule capacity:
    message: "Excavate the installed capacities of scenario {wildcards.scenario}."
    input:
        src = "src/analyse/capacity.py",
        results = rules.run.output
    output:
        raw = "build/output/{scenario}/capacity-raw.csv",
        publish = "build/output/{scenario}/capacity-publish.csv"
    script: "../src/analyse/capacity.py"


rule trade:
    message: "Excavate the traded electricity of scenario {wildcards.scenario}."
    input:
        src = "src/analyse/trade.py",
        results = rules.run.output
    output: "build/output/{scenario}/trade.csv"
    script: "../src/analyse/trade.py"


rule test:
    message: "Run tests"
    input:
        "tests/test_capacity_constraints.py",
        expand("build/output/{scenario}/capacity-raw.csv", scenario=config["scenarios"]),
        expand("build/output/{scenario}/results.nc", scenario=config["scenarios"])
    output: "build/test-report.html"
    shell:
        "py.test ./tests/ --html={output} --self-contained-html"
