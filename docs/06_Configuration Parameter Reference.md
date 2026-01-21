## Complete Configuration Parameter Reference

This section documents **all available parameters** in `config/config.yaml`.



### Project

| Parameter | Type | Modification Required | Description |
|---------|------|----------|-------------|
| `project.name` | string | no | Project name |

---

### Input

| Parameter | Type | Modification Required | Description |
|---------|------|----------|-------------|
| `input.manifest_raw` | path | **yes** | Raw sample manifest |
| `manifest_mode` | string | no | `filtered` or `complete` |
| `manifest_complete` | path | no | Added manifest |

---

### Programs

| Parameter | Type | Modification Required | Description |
|---------|------|----------|-------------|
| `program.ehdn` | path | no | ExpansionHunterDenovo binary |
| `program.liftover` | path | no | UCSC liftOver binary |
| `program.RExPRT` | path | no | RExPRT script |

---

### Resources

| Parameter | Type | Modification Required | Description |
|---------|------|----------|-------------|
| `resources.resource_directory` | path | no | Resource root directory |
| `resources.geneset_directory` | path | no | Gene set directory |
| `resources.reference_genome` | path | no | Reference FASTA |
| `resources.liftover_chain` | path | no | LiftOver chain file |
| `resources.genome_bins` | path | no | Genome bin definitions |
| `resources.phylop` | path | no | PhyloP scores |
| `resources.fragile_sites` | path | no | Fragile site annotations |
| `resources.known_repeats` | path | no | Known STR annotations |
| `resources.refflat` | path | no | RefFlat annotation |
| `resources.introns` | path | no | Intron annotation |
| `resources.exons` | path | no | Exon annotation |
| `resources.appris_108` | path | no | APPRIS v108 |
| `resources.appris_109` | path | no | APPRIS v109 |

---

### Results

| Parameter | Type | Modification Required | Description |
|---------|------|----------|-------------|
| `results.ehdn` | path | no | EHDN output |
| `results.qc` | path | no | QC results |
| `results.dbscan` | path | no | DBSCAN results |
| `results.annotation` | path | no | Annotation results |
| `results.rexprt` | path | no | RExPRT output |
| `results.analysis` | path | no | Final analysis output |

---

### Pipeline

| Parameter | Type | Modification Required | Description |
|---------|------|----------|-------------|
| `pipeline.ehdn_parallelize` | boolean | no | Parallelize ehdn profiler   |
| `pipeline.ehdn_max_parallel` | int | no | How many samples running parallel |


---

### Analysis Parameters

| Parameter | Type | Modification Required | Description |
|---------|------|----------|-------------|
| `analysis_privacy` | string | no | Use all or private expansions |
| `analysis_purity` | string | no | Mixed or pure loci |
| `include_cpg` | boolean | no | Include CpG annotation |
| `group_order` | list | **yes** | Group order (reference first) |
| `tss_window` | int | no | size of window (distance) |
| `sj_window` | int | no | size of window (distance) |
| `proximity_comparisons` | list | no | Comparisons to be made |
| `extra_comparisons` |string | no | Custom comparison |
| `exclude_samples_1` | path | no | path to excluded samples |
| `exclude_samples_2` | path | no | path to excluded samples |
| `plot_options` |string | no |which plot is to be saved |

---
### Gene Sets 

| Parameter | Type | Modification Required | Description |
|---------|------|----------|-------------|
| `geneset_list` | list | no | Gene sets to include |
| `geneset_mode` | string | no | Gene set analysis mode |
| `custom_geneset_dir` | path | no | Custom gene set directory |

---

### ORA

| Parameter | Type | Modification Required | Description |
|---------|------|----------|-------------|
| `gmt_file` | path | conditional | GMT file for ORA |
| `gem_per_group` | boolean | no | Separate GEMs per group |
| `min_term_size` | int | no | Minimum ORA term size |
| `max_term_size` | int | no | Maximum ORA term size |

---
