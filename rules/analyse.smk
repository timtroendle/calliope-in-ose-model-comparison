configfile: "config/default.yaml"
# The two rules `run_germany` and `run_europe` exist because they are vastly different
# computational requirements and in seperation they can be configured for cluster execution
# seperately. I am using the ruleorder approach here, because "negative" regular expressions
# would be needed for run_europe and they do not work in Snakemake at the moment. See:
# https://bitbucket.org/snakemake/snakemake/issues/684/wildcard_constraints-not-working-with
ruleorder: run_germany > run_europe



rule run_europe:
    message: "Run the model for scenario {wildcards.scenario}."
    input:
        model = rules.model.output.model
    output: "build/output/{scenario}/results.nc"
    shell:
        "calliope run {input.model} --save_netcdf {output} --scenario={wildcards.scenario}"


rule run_germany:
    message: "Run the model for scenario {wildcards.scenario}."
    input:
        model = rules.model.output.model
    output: "build/output/{scenario}/results.nc"
    wildcard_constraints:
        scenario = ".*(germany).*"
    shell:
        "calliope run {input.model} --save_netcdf {output} --scenario={wildcards.scenario}"


rule plot:
    message: "Visualises the results of scenario {wildcards.scenario}."
    input:
        src = "src/analyse/vis.py",
        results = "build/output/{scenario}/results.nc"
    params: scaling_factor = config["scaling-factors"]["power"]
    output: "build/output/{scenario}/plot.png"
    script: "../src/analyse/vis.py"


rule capacity:
    message: "Excavate the installed capacities of scenario {wildcards.scenario}."
    input:
        src = "src/analyse/capacity.py",
        results = "build/output/{scenario}/results.nc"
    params: scaling_factor = config["scaling-factors"]["power"]
    output:
        raw = "build/output/{scenario}/capacity-raw.csv",
        publish = "build/output/{scenario}/capacity-publish.csv"
    script: "../src/analyse/capacity.py"


rule storage_capacity:
    message: "Excavate the installed storage capacities of scenario {wildcards.scenario}."
    input:
        src = "src/analyse/storage_capacity.py",
        results = "build/output/{scenario}/results.nc"
    params: scaling_factor = config["scaling-factors"]["power"]
    output:
        raw = "build/output/{scenario}/storage-capacity-raw.csv",
        publish = "build/output/{scenario}/storage-capacity-publish.csv"
    script: "../src/analyse/storage_capacity.py"


rule trade:
    message: "Excavate the traded electricity of scenario {wildcards.scenario}."
    input:
        src = "src/analyse/trade.py",
        results = "build/output/{scenario}/results.nc"
    params: scaling_factor = config["scaling-factors"]["power"]
    output: "build/output/{scenario}/trade.csv"
    script: "../src/analyse/trade.py"


rule capacity_diff:
    message: "Show diff in installed capacities compared to baseline."
    input:
        src = "src/analyse/capacity_diff.py",
        baseline = "build/output/baseline/capacity-raw.csv",
        other = "build/output/{scenario}/capacity-raw.csv"
    output:
        "build/output/{scenario}/capacity-diff.csv"
    script: "../src/analyse/capacity_diff.py"


rule storage_capacity_diff:
    message: "Show diff in installed storage capacities compared to baseline."
    input:
        src = "src/analyse/storage_capacity_diff.py",
        baseline = "build/output/baseline/storage-capacity-raw.csv",
        other = "build/output/{scenario}/storage-capacity-raw.csv"
    output:
        "build/output/{scenario}/storage-capacity-diff.csv"
    script: "../src/analyse/storage_capacity_diff.py"


rule cost:
    message: "Excavate levelised costs from scenarios."
    input:
        src = "src/analyse/cost.py",
        results = expand("build/output/{scenario}/results.nc", scenario=config["scenarios"])
    params:
        power_scaling_factor = config["scaling-factors"]["power"],
        monetary_scaling_factor = config["scaling-factors"]["monetary"]
    output:
        "build/output/cost.csv"
    script: "../src/analyse/cost.py"


rule test:
    message: "Run tests"
    input:
        "src/analyse/test_runner.py",
        "tests/test_capacity_constraints.py",
        "tests/test_renewable_shares.py",
        "tests/test_co2_caps.py",
        expand("build/output/{scenario}/capacity-raw.csv", scenario=config["scenarios"]),
        re_shares = rules.renewable_shares.output.csv,
        results = expand("build/output/{scenario}/results.nc", scenario=config["scenarios"])
    params:
        scaling_factors = config["scaling-factors"]
    output: "build/logs/test-report.html"
    script: "../src/analyse/test_runner.py"
