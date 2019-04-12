PANDOC = "pandoc --filter pantable --filter pandoc-fignos --filter pandoc-tablenos --filter pandoc-citeproc"

include: "rules/construct.smk"
include: "rules/analyse.smk"
include: "rules/sync.smk"

configfile: "config/default.yaml"
localrules: all, report, clean


onstart:
    shell("mkdir -p build/logs")
onsuccess:
    if "email" in config.keys():
        shell("echo "" | mail -s 'ose-model-comparison succeeded' {config[email]}")
onerror:
    if "email" in config.keys():
        shell("echo "" | mail -s 'ose-model-comparison crashed' {config[email]}")


rule all:
    message: "Run entire analysis and compile report."
    input:
        "build/report.html",
        "build/output/all-results.csv",
        "build/output/all-ts-results.csv"


rule output:
    message: "Collect all results in standardised format."
    input:
        src = "src/analyse/output.py",
        scenarios = expand(
            "build/output/{scenario}/results.nc",
            scenario=[scenario for scenario in config["scenarios"]]
        )
    params: scaling_factors = config["scaling-factors"]
    output: "build/output/all-results.csv"
    script: "src/analyse/output.py"



rule ts_output:
    message: "Collect all time series results in standardised format."
    input:
        src = "src/analyse/ts_output.py",
        scenarios = expand(
            "build/output/{scenario}/results.nc",
            scenario=[scenario for scenario in config["scenarios"]]
        )
    params: scaling_factors = config["scaling-factors"]
    output: "build/output/all-ts-results.csv"
    conda: "envs/output.yaml"
    script: "src/analyse/ts_output.py"


rule report:
    message: "Compile report."
    input:
        "report/literature.bib",
        "report/main.md",
        "report/pandoc-metadata.yml",
        "build/output/baseline/plot.png",
        "build/output/baseline/capacity-publish.csv",
        "build/output/baseline/storage-capacity-publish.csv",
        "build/output/baseline/trade.csv",
        expand(
            "build/output/{scenario}/capacity-diff.csv",
            scenario=[scenario for scenario in config["scenarios"] if scenario != "baseline"]
        ),
        expand(
            "build/output/{scenario}/storage-capacity-diff.csv",
            scenario=[scenario for scenario in config["scenarios"] if scenario != "baseline"]
        ),
        "build/output/cost.csv",
        rules.test.output
    output:
        "build/report.html"
    shell:
        """
        cd ./report
        {PANDOC} --self-contained --css=report.css main.md pandoc-metadata.yml \
        --to html5 -o ../build/report.html
        """


rule clean: # removes all generated results
    shell:
        """
        rm -r ./build/*
        echo "Data downloaded to data/ has not been cleaned."
        """

