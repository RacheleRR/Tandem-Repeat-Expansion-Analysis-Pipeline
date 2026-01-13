rule dbscan:
    input:
        manifest = rules.qc.output.tsv
    output:
        tsv = config["results"]["dbscan"]+ "/DBSCAN_expansions.tsv"
    threads: min(12, workflow.cores)
    log:
        "logs/dbscan/dbscan.log"
    shell:
        """
        mkdir -p $(dirname {output.tsv})
        bash scripts/3_DBSCAN.sh \
            -m {input.manifest} \
            -o $(dirname {output.tsv}) \
            -t {threads} \
            > {log} 2>&1
        """

