# Tandem Repeat Expansion Analysis Pipeline

## 🔍 Overview
This pipeline provides an end-to-end, reproducible framework for the genome-wide analysis of tandem repeat (TR) expansions from whole-genome sequencing data. It enables systematic detection, characterization, and biological interpretation of repeat expansions across multiple cohorts and phenotypic groups, with an emphasis on statistical rigor and downstream functional relevance.

## ✨ Key Capabilities
- **Genome-wide TR expansion analysis** across single or multiple cohorts
- **Cohort-aware modeling**, supporting case–control and multi-group study designs
- **Robust statistical framework** integrating burden tests, non-parametric comparisons, and regression modeling
- **Biological contextualization** of expansions via genomic annotation, regulatory proximity, and gene set enrichment
- **Flexible filtering strategies**, including private vs shared expansions and pure vs mixed group assignments
- **Integrated functional enrichment and network analysis**, with optional Cytoscape visualization
- **Fully reproducible execution**, implemented as a modular Snakemake workflow with YAML-based configuration

## 📦 Requirements
- **Cores**: Minimum 4, recommended 12+
- **Software**: Snakemake ≥7.0, R ≥4.1, Python ≥3.8, Cytoscape ≥ 3.9 (optional, for network visualization),Annovar 
- **R packages**: See `scripts/fundamental_functions.r` for auto-installation
- **Python packages**: pandas, biopython 


## 🚀 Installation
### 1.Setup
```bash
# 1. Clone repository
git clone https://github.com/yourusername/ehdn-pipeline.git
cd ehdn-pipeline

# 2. Run automated setup (downloads tools, resources, helper scripts)
bash setup/setup.sh

```
### 2. Manuali download annovar files 
register and download annovar from this website 
https://www.openbioinformatics.org/annovar/annovar_download_form.php
or this https://annovar.openbioinformatics.org/en/latest/user-guide/download/
```bash
tar xvfz annovar.latest.tar.gz 
mv annovar helper/annovar 
```

### 3. Download gmt from gprofiler 
If you desire to do the network analysis you need to go onto https://biit.cs.ut.ee/gprofiler/gost and to Data source select the data sources you would like to use for your network analysis and download the combined names gmt 

```bash
mv combined_names.gmt  resources/ 
```


## Configuration
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

## Citation

If you use this pipeline, please cite:

- **ExpansionHunterDenovo**: Dolzhenko et al. (2020) Genome Biology
- **g:Profiler**: Raudvere et al. (2019) Nucleic Acids Research
- **DBSCAN**: Ester et al. (1996) KDD-96

## License

MIT License - see LICENSE file

## Contact

For questions or issues, please open a GitHub issue or contact [your email].

## Acknowledgments

Built on methods from:
- Fazal et al. (2020) - EHDN helper scripts
- Trost et al. (2020) - BTlib utilities


## 👥 Contact

For questions or issues:
- GitHub Issues: [link]
- Email: [your email]

## 🔄 Version History

- **v1.0.0** (2026-01): Initial release
  - Complete EHDN workflow
  - Multi-group statistical analysis
  - Enrichment and network analysis
  





