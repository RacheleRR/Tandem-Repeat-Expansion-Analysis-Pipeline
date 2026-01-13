rule qc:
    input:
        ehdn_dir = rules.ehdn.output.done
    output:
        tsv  = config["results"]["qc"] + "/filtered_manifest.tsv"
    log:
        "logs/qc/qc.log"
    shell:
        """
        mkdir -p $(dirname {output.tsv})
        PYTHONPATH={workflow.basedir}/../helper:${{PYTHONPATH:-}}\
        bash scripts/2_QC.sh \
            -i {config[results][ehdn]} \
            -o $(dirname {output.tsv}) \
            > {log} 2>&1
        """
