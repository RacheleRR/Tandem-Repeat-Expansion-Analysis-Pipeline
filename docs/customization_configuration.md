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
