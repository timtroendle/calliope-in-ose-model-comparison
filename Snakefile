PANDOC = "pandoc --filter pantable --filter pandoc-fignos --filter pandoc-tablenos --filter pandoc-citeproc"

include: "rules/construct.smk"
include: "rules/analyse.smk"

configfile: "config/default.yaml"

rule all:
    message: "Run entire analysis and compile report."
    input: "build/report.html"


rule report:
    message: "Compile report."
    input:
        "report/literature.bib",
        "report/main.md",
        "report/pandoc-metadata.yml",
        "build/output/baseline/plot.png",
        "build/output/baseline/capacity-publish.csv",
        "build/output/capacity-diff.csv",
        "build/output/cost-diff.csv",
        "build/output/baseline/trade.csv",
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

