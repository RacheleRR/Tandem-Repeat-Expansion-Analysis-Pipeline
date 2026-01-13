
##  Run Pipeline
### Usage

**Dry run** (check what will be executed):
```bash
cd TRs-pipeline
snakemake -n
```

### Full Pipeline

```bash
snakemake --cores 12
```

### Run Specific Modules

```bash
# Just EHDN profiling + QC
snakemake --cores 48 results/qc/filtered_manifest.tsv

# Annotation only (requires DBSCAN output)
snakemake --cores 12 results/annotation/ehdn_DBSCAN_annotated.tsv

# All analysis (requires annotated data)
snakemake --cores 16 results/analysis/reports/all_mixed/done
```
