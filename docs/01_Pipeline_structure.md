

## 📊 Pipeline Stages

### Stage 1: Tandem Repeat Detection (1-3)
Manifest → EHDN Profiling → QC Filtering → DBSCAN Clustering
- **Input**: BAM files and sample manifest (sample_id, status, path)
- **Process**: Genome-wide TRs detection (using anchored IRR (in-repeat read) patterns), sample-level QC, and density-based clustering
- **Output**: High-confidence TR expansion loci with outlier samples removed
- **Key QC**: ±3 SD outlier filtering; minimum two samples per locus

---

### Stage 2: Annotation & Pathogenicity Scoring (4-5)
TR Loci → ANNOVAR Annotation → RExPRT Scoring
- **Annotations**: Gene overlap, genomic region 
- **Scoring**: Pathogenicity prediction using RExPRT (hg19 liftover)

---

### Stage 3: Statistical & Genomic Analysis

#### 3.1 Correlation with Genomic Features
- Features: Fragile sites, GC content, PhastCons, phyloP
- Logistic regression across genomic bins
- **Outputs**: Feature-wise odds ratios tsv + bar plot 

#### 3.2 Proximity Analysis
- Distance-based analysis to TSSs and splice junctions
- Wilcoxon-based statistical comparisons
- **Outputs**: Density plots and summary tables

#### 3.3 Genomic Burden Analysis
- Fisher’s exact tests for regional enrichment
- Non-parametric comparisons of TR counts across groups
- **Outputs**: results tables, volcano& Manhattan & bar plots

#### 3.4 Regression Modeling
- Logistic regression for binary phenotypes
- Multinomial regression for multi-group designs
- **Predictors**: Genomic region, CpG content 
- **Outputs**: Summary tabel  & bar plot 

#### 3.5 Gene Set Burden 
- Fisher and non-parametric tests 
- Supported sets: SCHEMA, BipEx, brain-expressed genes, custom inputs
- **Outputs**:  Summary tabel and burden visualizations

---

### Stage 4: Functional Enrichment & Reporting
- Automated generation of summary tables, figures, and reports
- ORA via g:Profiler (GO, KEGG, Reactome)
- Optional network visualization with Cytoscape EnrichmentMap
