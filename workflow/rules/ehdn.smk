rule ehdn:
    input:
        manifest = config["input"]["manifest_raw"]
    output:
        done = config["results"]["ehdn"] + "/.done"
    params:
        outdir = config["results"]["ehdn"],
        ehdn_bin = config["program"]["ehdn"],
        ref = config["resources"]["reference_genome"],
        # semantic name
        parallelize = "-p" if config["pipeline"].get("ehdn_parallelize", True) else "",
        max_parallel = config["pipeline"].get("ehdn_max_parallel", 4)
    log:
        "logs/ehdn/ehdn.log"
    shell:
        """
        mkdir -p {params.outdir}
        bash scripts/1_run_ehdn.sh \
            -m {input.manifest} \
            -r {params.ref} \
            -e {params.ehdn_bin} \
            -o {params.outdir} \
            {params.parallelize} \
            --max-parallel {params.max_parallel} \
            > {log} 2>&1
        touch {output.done}
        """
