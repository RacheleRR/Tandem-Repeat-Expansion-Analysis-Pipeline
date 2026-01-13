rule rexprt:
    input:
        annotated = rules.annotate.output.tsv
    output:
        tsv = config["results"]["rexprt"] + "/ehdn_DBSCAN_annotated_rex_RExPRTscores.txt"
    params:
        liftover = config["program"]["liftover"],
        chain = config["resources"]["liftover_chain"],
        rexprt = config["program"]["RExPRT"]
    log:
        "logs/rexprt/rexprt.log"
    shell:
        """
        mkdir -p $(dirname {output.tsv})
        bash scripts/5_run_rexptr.sh \
            -i {input.annotated} \
            -o $(dirname {output.tsv}) \
            -l {params.liftover} \
            -c {params.chain} \
            -r {params.rexprt} \
            > {log} 2>&1
        """
