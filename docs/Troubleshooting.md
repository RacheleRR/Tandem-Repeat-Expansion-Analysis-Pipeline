

## 🔧 Troubleshooting

### No Significant Results

**Symptom**: Analysis completes but no significant findings

**Solutions**:
1. Check if you have sufficient samples per group (n ≥ 10 recommended)(`manifest_ufficial.tsv`)
2. Try relaxing `analysis_purity` to "mixed"
3. Try  `analysis_privacy: "all"` to increase power
4. Use `plot_options: "all"` to see all results regardless of significance
5. Check `ANALYSIS_SUMMARY.txt` files for details


### Cytoscape Errors

**Symptom**: `Warning: Cytoscape not running!`

**Solutions**:
1. Start Cytoscape before running ORA analysis
2. Install RCy3 package
3. Ensure EnrichmentMap app is installed in Cytoscape
4. Or disable Cytoscape integration (networks wont be generated)

### "No enrichment results"
- Check gene lists have ≥5 genes
- Verify term size filters in config (default: 2-2500)


### File Not Found Errors

**Symptom**: `Error: file.exists(x) is not TRUE`

**Solutions**:
1. Check all paths in `config.yml` are absolute or relative to project root
2. Verify resource files are downloaded
3. Check file permissions


### Patch Application Errors
- Ensure you have `sed` installed
- if the patch fails ( do this  sed -i 's/\r$//' helper/BTlib_docx_imports.patch)
