PANDOC = "pandoc --filter pantable --filter pandoc-fignos --filter pandoc-tablenos --filter pandoc-citeproc"


rule all:
    message: "Run entire analysis and compile report."
    input: "build/report.html"


rule run:
    message: "Run the model."
    input: "src/simple-model.yaml"
    output: "build/results.nc"
    shell: "calliope run {input} --save_netcdf {output}"


rule plot:
    message: "Visualises the results."
    input:
        src = "src/vis.py",
        results = rules.run.output
    output: "build/plot.png"
    script: "src/vis.py"


rule report:
    message: "Compile report."
    input:
        "report/literature.bib",
        "report/main.md",
        "report/pandoc-metadata.yml",
        rules.plot.output
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


rule test:
    shell:
        "py.test"
