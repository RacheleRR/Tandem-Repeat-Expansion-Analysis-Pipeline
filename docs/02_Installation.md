# 🚀 Installation

## Prerequisites

Before installing the pipeline, ensure the following software is available:

- Python 3.8 or higher
- Snakemake (>= 7.0)
- R studio (≥4.1)
- git-lfs (necessary to download the bigger files)
- GNU parallel (necessary to run samples in parallel in ehdn) 
- Cytoscape ≥ 3.9 (optional, for network visualization)

## 1.Setup
```bash
# 1. Clone repository
git lfs install
git clone https://github.com/yourusername/ehdn-pipeline.git
cd ehdn-pipeline

# 2. Run automated setup (downloads tools, resources, helper scripts)
bash setup/setup.sh

```
## 2. Manual download of ANNOVAR files
1. Register and download ANNOVAR from here:
   [OpenBioinformatics ANNOVAR download form](https://www.openbioinformatics.org/annovar/annovar_download_form.php)
2. Extract and move ANNOVAR files into the pipeline helper directory:

```bash
tar xvfz annovar.latest.tar.gz 
mv annovar helper/annovar

```

## 3. Optional: Download GMT files for network analysis
If you plan to run the network analysis:
1. Go to g:Profiler [https://biit.cs.ut.ee/gprofiler/gost](https://biit.cs.ut.ee/gprofiler/gost)
3. Select the desired data sources
4. Download the combined names GMT file
5. Move the GMT file into the resources/ folder:

```bash
mv combined_names.gmt  resources/

```

## 4. Optional: Custom Geneset
If you decide that you want to perform the anlysis with your custom geneset/s. 
It is important to respect some rules.

The custom Geneset:
- Must be **CSV files** or **TSV files**
- Required column: **Gene**
- One gene symbol per row
- Files should already be filtered for relevance/significance
- The Custom Genese/s need to me moved to the `resources/geneset/custom` directory

  
