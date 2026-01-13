rule annotate:
    input:
        dbscan_dir = rules.dbscan.output.tsv
    output:
        tsv = config["results"]["annotation"]+ "/ehdn_DBSCAN_annotated.tsv"
    params:
        buildver = "hg38"
    log:
        "logs/annotate/annotate.log"
    shell:
        """
        mkdir -p $(dirname {output.tsv})
        bash scripts/4_annotate.sh \
            -i {input.dbscan_dir} \
            -o $(dirname {output.tsv}) \
            -b {params.buildver} \
            > {log} 2>&1
        """