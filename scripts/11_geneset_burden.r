#!/usr/bin/env Rscript
# ============================================================================
#11_genesetburden 
# PURPOSE: Gene set burden analysis for psychiatric disease genes
# 
# INPUT:
#   - args[1]: Input directory (contains data_prepare outputs)
#   - args[2]: Output directory path
#   - args[3]: Privacy filter ("private" or "all")
#   - args[4]: Purity filter ("pure" or "mixed")
#   - args[5]: samples to Exclude list 1 file path
#   - args[6]: samples to exclude list 2 file path
#   - args[7]: Plot option ("significant", "all", "global_only")
#   - args[8]  Base geneset directory (default / built-in genesets)
#   - args[9]  Custom geneset directory ("" if unused)
#   - args[10] Geneset list to analyze ("all" or comma-separated names)
#   - args[11] Geneset mode ("basic", "combined", "different")
#
# OUTPUT:
#   - Fisher test results for each gene set
#   - Non-parametric test results (Wilcoxon/Kruskal)
#   - Volcano, Manhattan, and dot plots for fisher 
#   - Boxplots for Wilcoxon/Kruskal-Wallis 
#   - done file: Completion marker
#
# STATISTICS:
#   - Fisher's exact test per gene set
#   - Wilcoxon (2 groups) or Kruskal-Wallis (>2 groups)
#   - Dunn's post-hoc tests for significant Kruskal results
#
# Baisc GENE SETS:
#   - Brain-expressed genes (consensus)
#   - Brain-expressed genes (nTPM)
#   - SCHEMA schizophrenia genes (p-value)
#   - Bipolar disorder genes
#   - SCHEMA schizophrenia genes (q-value)
#
# GENESET MODES:
#   - basic     : built-in gene sets only
#   - combined  : built-in + user-provided gene sets
#   - different : user-provided gene sets only
#
# FILTERS:
#   - Privacy: "private" (count=1) vs "all" TRs
#   - Purity: "pure" (no mixed labels) vs "mixed" (all labels)
#   - Sample exclusions: Removes specified samples from analysis
#   - Geneset list to analyze 
#
# AUTHOR: Rachele R. Rubiu
# DATE: 12/12/2025
# VERSION: 6.0
# ============================================================================
args <- commandArgs(trailingOnly = TRUE)

cat("Parameters:\n")
print(args)
cat("\n")

input_dir  <- args[1]
output_dir <- args[2]
privacy    <- args[3]
purity     <- args[4]
excl1_path <- args[5]
excl2_path <- args[6]
plot_option <- args[7]
geneset_dir <- args[8]
custom_geneset_dir  <- args[9]   
geneset_list_arg    <- args[10] 
geneset_mode <- args[11]  

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================
# Load libraries & functions 
# ============================================
suppressPackageStartupMessages({
library(dplyr)
library(tidyr)
library(readr)
})

source("scripts/fundamental_functions.r")
source("scripts/statistics_functions.r")
source("scripts/plots_functions.r")

# ============================================
# Load data  & prepare data 
# ============================================

ehdn_modified <- read.delim(
  file.path(input_dir, "ehdn_modified.tsv")
)

manifest_ufficial <- read.delim(
  file.path(input_dir, "manifest_ufficial.tsv")
)

# Filter by privacy
if (privacy == "private") {
  tr_data <- ehdn_modified %>% filter(count == 1)
} else {
  tr_data <- ehdn_modified
}

if (purity == "pure"){
  tr_data <- tr_data %>% filter(!outlier_label == "mixed")
} else if (purity == "mixed"){
  tr_data <- tr_data 
}

read_exclusions <- function(path) {
  if (is.null(path) || path == "" || !file.exists(path)) {
    return(character(0))
  }
  read.delim(path, header = FALSE)[[1]]
}

all_exclude_samples <- unique(c(
  read_exclusions(excl1_path),
  read_exclusions(excl2_path)
))


groups_info <- extract_groups_from_manifest(manifest_ufficial)
group_names <- groups_info$groups 
group_counts <- groups_info$counts
data_name <- make_data_name(
  base = "ehdn",
  privacy = privacy,
  purity  = purity
)

cat("Using data_name:", data_name, "\n")
#===========================
#GENELISTS
#===============================
genesets_basic <- load_genesets(geneset_dir)

basic_gene_lists <- list(
  brain        = genesets_basic$genes_brain,
  brain_ntpm   = genesets_basic$genes_brain_ntpm,
  schema_pval  = genesets_basic$genes_schema_pval,
  schema_qval  = genesets_basic$genes_schema_qval,
  bipolar      = genesets_basic$genes_bipolar
)

supplementary_gene_lists <- NULL

if (geneset_mode %in% c("combined", "different")) {

  if (custom_geneset_dir == "" || !dir.exists(custom_geneset_dir)) {
    stop("Custom geneset directory must be provided for mode: ", geneset_mode)
  }

  supplementary_gene_lists <- load_custom_genesets(custom_geneset_dir)
}

gene_lists <- select_genesets(
  geneset_mode = geneset_mode,
  basic_genesets = basic_gene_lists,
  supplementary_genesets = supplementary_gene_lists
)


if (!is.null(geneset_list_arg) && geneset_list_arg != "all") {

  requested_sets <- trimws(unlist(strsplit(geneset_list_arg, ",")))

  missing <- setdiff(requested_sets, names(gene_lists))
  if (length(missing) > 0) {
    stop(
      "Requested genesets not found: ",
      paste(missing, collapse = ", "),
      "\nAvailable: ",
      paste(names(gene_lists), collapse = ", ")
    )
  }

  gene_lists <- gene_lists[requested_sets]
}

if (length(gene_lists) == 0) {
  stop("No genesets selected after filtering")
}

cat("\n=== Final genesets to analyze ===\n")
for (nm in names(gene_lists)) {
  cat(sprintf(" - %s (%d genes)\n", nm, length(gene_lists[[nm]])))
}
cat("================================\n\n")


# ============================================
# fisher analysis 
# ============================================
fisher_tabel <-create_fisher_count_tables(tr_data, data_name, gene_lists, private = FALSE, filter_sample = all_exclude_samples, group_names) 
fisher <- perform_fisher_analysis(fisher_tabel,manifest_ufficial,group_counts=group_counts) 
fisher$type <- extract_type_from_row_name(fisher$Row_Name)

plot_fisher_volcano(fisher,Row_Name= "type" ,save_path = file.path(output_dir, "fisher_volcano.png"))

# Manhattan plot  
plot_fisher_manhattan(fisher,Row_Name="type" ,save_path = file.path(output_dir, "fisher_manhattan.png"))

# Dot plot
plot_fisher_dot(fisher, Row_Name = "type", save_path = file.path(output_dir, "fisher_dot.png"))

write.table(
  fisher,
  file.path(output_dir, "fisher_results.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


# ============================================
# Wilcoxon or kruskal analyisis
# ============================================
complete_data <-create_wilcoxon_kruskal_tables(tr_data, data_name, gene_lists, private = FALSE, manifest_data = manifest_ufficial, filter_sample = all_exclude_samples,group_names = group_names) 
results_nonpara <- perform_nonparametric_analysis(complete_data,  group_col = "group")

#need to put a thing here so that if it doesnt have any signifiant p values etc it want just crush the plot function 

if (nrow(results_nonpara) > 0) {
  plot_results <- plot_wilcoxon_kruskal_flexible(
    test_results = results_nonpara,
    count_data = complete_data,
    group_col = "group",
    test_type = "wilcoxon",
    plot_option = plot_option,
    save_dir = file.path(output_dir, "plots")
)
} else {
  message("No significant non-parametric results; skipping plot.")
}


write.table(
  results_nonpara,
  file.path(output_dir, "nonparametric_results.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

done_file <- file.path(output_dir, "done")
file.create(done_file)

cat("Geneset analysis Complete!\n")
cat("All plots and results saved in:", output_dir, "\n")
