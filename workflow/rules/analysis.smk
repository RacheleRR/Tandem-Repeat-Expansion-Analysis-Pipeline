ANALYSIS = config["results"]["analysis"]

rule prepare_data:
    input:
        ehdn = "results/annotation/ehdn_DBSCAN_annotated.tsv",
        manifest = "results/qc/filtered_manifest.tsv"
    output:
        annotated = f"{ANALYSIS}/data_prepare/ehdn_results_annotated.tsv",
        modified  = f"{ANALYSIS}/data_prepare/ehdn_modified.tsv",
        manifest  = f"{ANALYSIS}/data_prepare/manifest_ufficial.tsv"
    params:
        outdir = f"{ANALYSIS}/data_prepare",
        manifest_mode = config["manifest_mode"],
        manifest_complete = config["input"].get("manifest_complete", "")
    log:
        "logs/analysis/prepare_data.log"
    shell:
        """
        mkdir -p {params.outdir}
        Rscript scripts/6_prepare_data.r \
            {input.ehdn} \
            {params.outdir} \
            {params.manifest_mode} \
            {input.manifest} \
            "{params.manifest_complete}" \
            > {log} 2>&1
        """

rule correlation:
    input:
        annotated = rules.prepare_data.output.annotated,
        manifest  = rules.prepare_data.output.manifest
    output:
        done = f"{ANALYSIS}/correlation/{config['analysis_privacy']}/done"
    params:
        data_dir = f"{ANALYSIS}/data_prepare",
        outdir = f"{ANALYSIS}/correlation/{config['analysis_privacy']}",
        privacy = config["analysis_privacy"],
        resource_dir = config["resources"]["resource_directory"],
        extra_comparison = config["extra_comparisons"]
    log:
        f"logs/analysis/correlation_{config['analysis_privacy']}.log"
    shell:
        """
        mkdir -p {params.outdir}
        Rscript scripts/7_correlation.r \
            {params.data_dir} \
            {params.outdir} \
            {params.privacy} \
            {params.resource_dir} \
            {params.extra_comparison}\
            > {log} 2>&1
        """

rule proximity:
    input:
        annotated = rules.prepare_data.output.annotated,
        manifest  = rules.prepare_data.output.manifest
    output:
        done = f"{ANALYSIS}/proximity/{config['analysis_privacy']}/done"
    params:
        data_dir = f"{ANALYSIS}/data_prepare",
        outdir = f"{ANALYSIS}/proximity/{config['analysis_privacy']}",
        privacy = config["analysis_privacy"],
        comparisons = config["proximity_comparisons"],
        tss_window = config["tss_window"],
        sj_window = config["sj_window"],
        resource_dir = config["resources"]["resource_directory"],
        combined_count_dir = config["results"]["dbscan"]
    log:
        f"logs/analysis/proximity_{config['analysis_privacy']}.log"
    shell:
        """
        mkdir -p {params.outdir}
        Rscript scripts/8_proximity.r \
            {params.data_dir} \
            {params.outdir} \
            {params.privacy} \
            "{params.comparisons}" \
            {params.tss_window} \
            {params.sj_window} \
            {params.resource_dir} \
            {params.combined_count_dir} \
            > {log} 2>&1
        """

rule genomic_burden:
    input:
        prep_done = rules.prepare_data.output.annotated
    output:
        done = f"{ANALYSIS}/genomic_burden/{config['analysis_privacy']}_{config['analysis_purity']}/done"
    params:
        data_dir = f"{ANALYSIS}/data_prepare",
        outdir = f"{ANALYSIS}/genomic_burden/{config['analysis_privacy']}_{config['analysis_purity']}",
        privacy = config["analysis_privacy"],
        purity = config["analysis_purity"],
        excl1 = config["exclude_samples_1"],
        excl2 = config["exclude_samples_2"],
        include_cpg = config["include_cpg"]

    log:
        f"logs/analysis/genomic_burden_{config['analysis_privacy']}_{config['analysis_purity']}.log"
    shell:
        """
        mkdir -p {params.outdir}
        Rscript scripts/9_genomic_burden.r \
            {params.data_dir} \
            {params.outdir} \
            {params.privacy} \
            {params.purity} \
            "{params.excl1}" \
            "{params.excl2}" \
            {params.include_cpg} \
            > {log} 2>&1
        """

rule regression:
    input:
        burden_done = rules.genomic_burden.output.done
    output:
        done = f"{ANALYSIS}/regression/{config['analysis_privacy']}_{config['analysis_purity']}/done"
    params:
        data_dir = f"{ANALYSIS}/data_prepare",
        outdir = f"{ANALYSIS}/regression/{config['analysis_privacy']}_{config['analysis_purity']}",
        privacy = config["analysis_privacy"],
        purity = config["analysis_purity"],
        group_order = config.get("group_order", "auto"),
        include_cpg = config["include_cpg"]
    log:
        f"logs/analysis/regression_{config['analysis_privacy']}_{config['analysis_purity']}.log"
    shell:
        """
        mkdir -p {params.outdir}
        Rscript scripts/10_regression.R \
            {params.data_dir} \
            {params.outdir} \
            {params.privacy} \
            {params.purity} \
            "{params.group_order}" \
            {params.include_cpg} \
            > {log} 2>&1
        """

rule geneset_burden:
    input:
        burden_done = rules.genomic_burden.output.done
    output:
        done = f"{ANALYSIS}/geneset_burden/{config['analysis_privacy']}_{config['analysis_purity']}/done"
    params:
        data_dir = f"{ANALYSIS}/data_prepare",
        outdir = f"{ANALYSIS}/geneset_burden/{config['analysis_privacy']}_{config['analysis_purity']}",
        geneset_dir = config["resources"]["geneset_directory"],
        privacy = config["analysis_privacy"],
        purity = config["analysis_purity"],
        excl1 = config["exclude_samples_1"],
        excl2 = config["exclude_samples_2"],
        plot_option = config["plot_options"],
        custom_geneset_dir =config["custom_geneset_dir"],
        geneset_list_arg=",".join(config["geneset_list"]),
        geneset_mode= config["geneset_mode"]

    log:
        f"logs/analysis/geneset_burden_{config['analysis_privacy']}_{config['analysis_purity']}.log"
    shell:
        """
        mkdir -p {params.outdir}
        Rscript scripts/11_geneset_burden.r \
            {params.data_dir} \
            {params.outdir} \
            {params.privacy} \
            {params.purity} \
            "{params.excl1}" \
            "{params.excl2}" \
            {params.plot_option} \
            {params.geneset_dir} \
            "{params.custom_geneset_dir}"\
            "{params.geneset_list_arg}"\
            "{params.geneset_mode}"\
            > {log} 2>&1
        """
        

rule generate_report:
    input:
        corr_done = rules.correlation.output.done,
        prox_done = rules.proximity.output.done,
        geno_done = rules.genomic_burden.output.done,
        reg_done  = rules.regression.output.done,
        gene_done = rules.geneset_burden.output.done
    output:
        done = f"{ANALYSIS}/reports/{config['analysis_privacy']}_{config['analysis_purity']}/done"
    params:
        analysis_dir = f"{ANALYSIS}/data_prepare",
        outdir = f"{ANALYSIS}/reports/{config['analysis_privacy']}_{config['analysis_purity']}",
        privacy = config["analysis_privacy"],
        purity = config["analysis_purity"],
        excl1 = config["exclude_samples_1"],
        excl2 = config["exclude_samples_2"],
        geneset_dir = config["resources"]["geneset_directory"],
        custom_geneset_dir =config["custom_geneset_dir"],
        geneset_list_arg=",".join(config["geneset_list"]),
        geneset_mode= config["geneset_mode"]
    log:
        f"logs/analysis/report_{config['analysis_privacy']}_{config['analysis_purity']}.log"
    shell:
        """
        mkdir -p {params.outdir}
        Rscript scripts/12_generate_report.r \
            {params.analysis_dir} \
            {params.outdir} \
            {params.privacy} \
            {params.purity} \
            "{params.excl1}" \
            "{params.excl2}" \
            {params.geneset_dir} \
            "{params.custom_geneset_dir}"\
            "{params.geneset_list_arg}"\
            "{params.geneset_mode}"\
            > {log} 2>&1
        """

rule ora_analysis:
    input:
        modified = rules.prepare_data.output.modified,
        manifest = rules.prepare_data.output.manifest
    output:
        done = f"{ANALYSIS}/ORA/{{privacy}}_{{purity}}/done"
    params:
        privacy = config["analysis_privacy"],
        purity = config["analysis_purity"],
        gem_per_group = str(config.get("gem_per_group", True)).lower(),
        min_term_size = config.get("min_term_size", 2),
        max_term_size = config.get("max_term_size", 2500),
        data_dir = f"{ANALYSIS}/data_prepare",
        outdir = directory(f"{ANALYSIS}/ORA/{{privacy}}_{{purity}}"),
        gmt_file =config["gmt_file"]   
    log:
        "logs/analysis/ora_{privacy}_{purity}.log"
    shell:
        """
        mkdir -p {params.outdir}
        Rscript scripts/13_ORA_multiquerry.r \
            {params.data_dir} \
            {params.outdir} \
            {params.privacy} \
            {params.purity} \
            {params.gem_per_group} \
            {params.min_term_size} \
            {params.max_term_size} \
            {params.gmt_file} \
            > {log} 2>&1
        """
