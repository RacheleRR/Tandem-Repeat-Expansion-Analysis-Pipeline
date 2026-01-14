
##  Run Pipeline

This pipeline is executed using Snakemake. You may run the full workflow or specific modules depending on your needs

---

## Dry run :
Perform a dry run to inspect which rules will be executed and to check that the DEG is healthy without launching any jobs:

```bash
cd TRs-pipeline
snakemake -n
```

## Full Pipeline
Execute the complete pipeline using the specified number of cores available to you:

```bash
snakemake --cores 12
```

It is raccomended to use at least 12 cores. 

## Run Specific Modules
You may also run specific components of the pipeline by targeting intermediate or final output files.

```bash
# Just EHDN profiling + QC
snakemake --cores 48 results/qc/filtered_manifest.tsv

# Annotation only (requires DBSCAN output)
snakemake --cores 12 results/annotation/ehdn_DBSCAN_annotated.tsv

# All analysis (requires annotated data)
snakemake --cores 16 results/analysis/reports/all_mixed/done
```
