## 🔧 Configuration
This pipeline is configured through the file:

```
config/config.yaml
```

This file controls input data, analysis behavior, group comparisons, and gene set enrichment.
Most users will only need to modify a subset of parameters described below.

For a the Complete Configuration Parameter Reference please click [here](06_Configuration_Parameter_Reference.md)
 


## 1. Minimal Required Configuration
### 1.1 Input Manifest (Required)

```
input:
  manifest_raw: "/path/to/manifest.tsv"
```

The raw manifest must be a **TSV file** with **exactly three columns**:

| Column name | Description |
|------------|-------------|
| `sample_id` | Unique sample identifier |
| `Status` | Group label (e.g. Control, Case) |
| `BAM_path` | Absolute path to BAM file |


```
Example:
sample_id	Status	BAM_path
sample_001	Control	/path/sample_001.bam
sample_002	Case	/path/sample_002.bam
```
---

### 1.2 Group Definition (Required)

```yaml
group_order: "Control,Case"
```

- The **first group is the reference group**
- Regression coefficients are interpreted relative to this group
- Group labels must exactly match the `Status` column in the manifest

---

## 2. Common Optional Customization

These parameters are **not required to be modified**, but are frequently adjusted depending on the analysis.

---

### 2.1 Parallel ExpansionHunterDenovo Execution

For large cohorts (>50 samples), it is best to parallelized ExpansionHunterDenovo profiling .

```yaml
  ehdn_parallelize: true
  ehdn_max_parallel: 8
```

- `ehdn_max_parallel`
  →  the amount of samples to run in parallele
- `ehdn_parallelize`
  → if to parallelize or not the samples 


---

### 2.1 Optional Manifest
```
manifest_mode: "filtered"   # or "complete"
manifest_complete: "/path/to/manifest_complete.csv"
```

- `manifest_mode: filtered`  
  → only information from `manifest_raw` is used

- `manifest_mode: complete`  
  → merges informations from `manifest_complete`


The `manifest_complete` file must contain **at least**:

| Column | Required |
|-------|----------|
| `sample_id` | yes |
| `Status` | yes |


---

### 2.2 Privacy and Purity Filters
```yaml
analysis_privacy: "all"     # "all" or "private"
analysis_purity: "mixed"   # "mixed" or "pure"
```

- `analysis_privacy: private`  
  → only consider expansions unique to a single individual

- `analysis_purity: pure`  
  → exclude loci with mixed group membership  
  → useful for strict case-control analyses

---

### 2.3 CpG Annotation

```yaml
include_cpg: TRUE
```

If enabled, CpG overlap information is included in downstream annotation
and regression models.

---

### 2.4 Comparisons

```yaml
proximity_comparisons: "auto"  # or specify as "Group1-Group2,Group3-Group4"
extra_comparisons: "" # e.g., "Group1-Group2" or "Case-KnownSTRs"
```


`proximity_comparisons: `
- `auto` decides & performs all valid group comparisons
- Custom comparisons can be specified as comma-separated pairs  
  (e.g. `Group1-Group2,Group3-Group4`)

  
`extra_comparisons`
- Additional comparisons specified as group pairs  
  (e.g. `SCZ-BD`)


---

### 2.5 Burden Analysis Options

```yaml
exclude_samples_1: ""
exclude_samples_2: ""
plot_options: "significant" # "significant", "all", "global_only"
```
- If you have samples that in hindeside need to be excluded from the analysis use exclude_samples_1 or 2. Excluded samples can be specified as comma-separated lists 
- What plots should be saved for the wilxocoxon/kruskal analysis
  
---

### 2.6 Gene Set Burden
If you decide that you want to  perform the anlysis with your custom geneset/s you will have to modify these parameters 

```yaml
geneset_list: ["brain", "brain_ntpm"]
geneset_mode: "basic"
custom_geneset_dir: ""
```
It is important to note that the custom Geneset: 
- Must be **CSV files** or **TSV files**
- Required column: `Gene`
- One gene symbol per row
- Files should already be filtered for relevance/significance

Custom gene sets can be placed in a user-defined directory.

```geneset_list: []``` <- describes the gene sets that will be used for the analysis for example `["brain", "brain_ntpm"]`. If geneset_list is not modified, the pipeline will use the default gene sets defined within the pipeline.

```geneset_mode: "basic"``` <- describes which genesets should be used, **basic** stands for the once already present in the pipeline, **diffrent** stands for the once that are customized and **combined** stands for when you use both the basic and the customized genesets 


---

### 2.7 Over-Representation Analysis (ORA)

```yaml
gmt_file: "path/to/file.gmt"
gem_per_group: TRUE
min_term_size: 2
max_term_size: 2500
```

- `gmt_file` must follow standard GMT format
- `gem_per_group` controls whether GEMs are placed in one general folder or  a folder will be generated per group 
- `Term size limits` control ORA filtering

---


