## Configuration
This pipeline is configured through the file:

```
config/config.yaml
```

This file controls input data, analysis behavior, group comparisons, and gene set enrichment.
Most users will only need to modify a subset of parameters described below.

## 1. Input Files
### 1.1 Required Manifest

```
input:
  manifest_raw: "/path/to/manifest.tsv"
```
The raw manifest must be a TSV file with exactly three columns:

Column name	Description
sample_id	Unique sample identifier
Status	Group label (e.g. Control, Case)
BAM_path	Absolute path to BAM file

```
Example:
sample_id	Status	BAM_path
sample_001	Control	/path/sample_001.bam
sample_002	Case	/path/sample_002.bam
```

### 1.2 Optional Manifest
```
manifest_mode: "filtered"   # or "complete"
manifest_complete: "/path/to/manifest_complete.csv"
```

- ```manifest_mode: filtered```  → only uses information from manifest_raw
- ```manifest_mode: complete```  → merges informations from manifest_complete

The manifest_complete file must contain at least:

Column	Required
sample_id	yes
Status	yes

## 2. Analysis Behavior
### 2.1 Privacy and Purity Filters
```
analysis_privacy: "all"    # "all" or "private"
analysis_purity: "mixed"   # "mixed" or "pure"
```

- ```analysis_privacy: private``` → only consider expansions unique to a single individual
- ```analysis_purity: pure``` → exclude loci with mixed group membership (useful for strict case-control analyses)


### Key Parameters to Customize 
Edit `config/config.yaml`:

```yaml
input:
  manifest_raw: "/path/to/manifest.tsv"  # Required: 3 columns (Sample_ID | Status | BAM_path)

  # Optional: additional manifest with covariates
  manifest_mode: "filtered"   # or "complete"
  manifest_complete: "/path/to/full_manifest.csv" # Required: minimum 2 columns with exactly these names (sample_id | Status )

analysis_privacy: "all"     # "all" or "private" (expansions in only one individual)
analysis_purity: "mixed"    # "mixed" or "pure" (exclude multi-group outliers)
include_cpg: TRUE # or false 

# Group order matters for regression (reference = first)
group_order: ["Control", "Case"]  # Reference group first

gmt_file: "path/to/filegmt" 


```

## CUSTOMIZATION

## Group Comparisons
**3+ groups**:
```yaml
# Manifest Status column has: SCZ, BD, Control
# Pipeline performs: SCZ vs Control, BD vs Control, SCZ vs BD

group_order: ["Control", "SCZ", "BD"]

```

example of manifest_complete.tsv
```csv
sample_id,Status,Sex,Age
sample_001,SCZ,M,45
sample_002,Control,F,50
sample_003,BD,M,38
```

## GENE SETS
### Adding New Gene Sets
 Your gene sets should be in CSV format with a header column "Gene" and one gene symbol per row and they should already be filtere by significance etc.  
 Place your custom gene set files in the `resources/genesets/custom` directory.

After adding your custom gene sets, update the `geneset_list` in `config/config.yaml` to include the names of your new gene sets (without file extensions). For example:
```yaml
geneset_list: ["brain", "schema_pval", "my_custom_geneset"]
geneset_mode: "combined"  # or "different" depending on your analysis needs
custom_geneset_dir: "data/custom_genesets/"  # Optional 

```

## Parallel EHDN Profiling
#! what if i want to do paralelization
For large cohorts (>50 samples), use parallel execution:

```bash
# Edit workflow/rules/ehdn.smk to use per-sample rules
snakemake --cores 48 --resources mem_mb=200000 all_ehdn
```
