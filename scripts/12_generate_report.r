#!/usr/bin/env Rscript
#======================
# 12_Generate reports 
# PURPOSE: Generate comprehensive analysis reports and summary visualizations
# 
# INPUT:
#   - args[1]: Input directory (contains data_prepare outputs)
#   - args[2]: Output directory path
#   - args[3]: Privacy filter ("private" or "all")
#   - args[4]: Purity filter ("pure" or "mixed")
#   - args[5]: samples to Exclude list 1 file path
#   - args[6]: samples to exclude list 2 file path
#   - args[7]: Geneset directory path (default / built-in genesets)
#   - args[8]  Custom geneset directory ("" if unused)
#   - args[9]  Geneset list to analyze ("all" or comma-separated names)
#   - args[10] Geneset mode ("basic", "combined", "different")
#
# OUTPUT:
#   - Count tables: Summary of TRs for all gene sets
#   - Excel files: Gene lists organized by group and gene set
#   - Summmary tabel: List of genes associated to Disease and geneset 
#   - Summary plots: Count distributions by region, motif, contig
#   - Nucleotide composition plots
#   - done file: Completion marker
#
# FILTERS:
#   - Privacy: "private" (count=1) vs "all" TRs
#   - Purity: "pure" (no mixed labels) vs "mixed" (all labels)
#   - Sample exclusions: Removes specified samples from analysis
#
# GENESET MODES:
#   - basic     : built-in gene sets only
#   - combined  : built-in + user-provided gene sets
#   - different : user-provided gene sets only
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
geneset_dir <- args[7]
custom_geneset_dir  <- args[8]   
geneset_list_arg    <- args[9] 
geneset_mode <- args[10]  

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

#======================
# Load libraries & functions 
#======================
suppressPackageStartupMessages({
library(dplyr)
library(tidyr)
library(readr)
})
source("scripts/fundamental_functions.r")
source("scripts/statistics_functions.r")
source("scripts/plots_functions.r")


# ============================================
# Load and prepare data
# ============================================
ehdn_modified <- read.delim(
  file.path(input_dir, "ehdn_modified.tsv")
)
ehdn_results_annotated <- read.delim(
  file.path(input_dir, "ehdn_results_annotated.tsv")
)

manifest_ufficial <- read.delim(
  file.path(input_dir, "manifest_ufficial.tsv")
)

# Filter by privacy
if (privacy == "private") {
  tr_data_expanded <- ehdn_modified %>% filter(count == 1)
  tr_data <- ehdn_results_annotated %>% filter(count == 1)
} else {
  tr_data_expanded <- ehdn_modified
  tr_data <- ehdn_results_annotated

}

if (purity == "pure"){
  tr_data_expanded <- tr_data_expanded %>% filter(!outlier_label == "mixed")
  tr_data <- tr_data %>% filter(!outlier_label == "mixed")
} else if (purity == "mixed"){
    tr_data_expanded <- tr_data_expanded
  tr_data <- tr_data 
}

data_name <- make_data_name(
  base = "ehdn",
  privacy = privacy,
  purity  = purity
)

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
# general tabel 
# ============================================

count_tabel <-create_count_tables(tr_data_expanded,  data_name, gene_lists, private = FALSE,filter_sample = all_exclude_samples,group_names = group_names) 

write.table(
  count_tabel,
  file.path(output_dir, "count_table.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

# ============================================
#excel with gene lists
# ============================================
old_wd <- getwd()
setwd(output_dir)

results <- create_gene_lists_excel(
  df = tr_data_expanded,
  dataset_name = "ehdn_modified",
  gene_lists = gene_lists
)

setwd(old_wd)

# ============================================
#summary for pathogenic genes 
# ============================================
 summary_results <- create_gene_list_summary(
  df = tr_data_expanded,                    # Your dataframe
  gene_lists = gene_lists               # Your gene lists
)

# ============================================
#plots 
# ============================================
variables_to_plot <- c("region", "motif", "contig", "motif_length")
plots_raw <- list()
plots_normalized <- list()

for (var in variables_to_plot) {
  counts <- compute_counts(ehdn_results_annotated, var, manifest_ufficial)
  
  # Raw plot
  p_raw <- plot_counts(
    counts$raw,
    y_var = "n_loci",
    title = paste("Raw counts of", var),
    x_label = var
  )
  plots_raw[[var]] <- p_raw
  
  # Save raw plot
  ggsave(
    filename = file.path(output_dir, paste0("raw_counts_", var, ".png")),
    plot = p_raw,
    width = 10,
    height = 6,
    dpi = 300
  )
  
  # Normalized plot
  p_norm <- plot_counts(
    counts$normalized,
    y_var = "normalized_count",
    title = paste("Normalized counts per sample of", var),
    x_label = var
  )
  plots_normalized[[var]] <- p_norm
  
  # Save normalized plot
  ggsave(
    filename = file.path(output_dir, paste0("normalized_counts_", var, ".png")),
    plot = p_norm,
    width = 10,
    height = 6,
    dpi = 300
  )
}

# ============================================
# Nucleotide plots
# ============================================
nucleotide_counts <- compute_nucleotide_counts(ehdn_results_annotated, manifest_ufficial)

# Raw nucleotide plot
nucleotide_raw_plot <- plot_nucleotides(
  nucleotide_counts$raw, 
  y_var = "n_loci", 
  title = "Raw nucleotide counts per cohort"
)

# Save raw nucleotide plot
ggsave(
  filename = file.path(output_dir, "raw_nucleotide_counts.png"),
  plot = nucleotide_raw_plot,
  width = 10,
  height = 6,
  dpi = 300
)

# Normalized nucleotide plot
nucleotide_norm_plot <- plot_nucleotides(
  nucleotide_counts$normalized, 
  y_var = "normalized_count", 
  title = "Normalized nucleotide counts per cohort"
)

# Save normalized nucleotide plot
ggsave(
  filename = file.path(output_dir, "normalized_nucleotide_counts.png"),
  plot = nucleotide_norm_plot,
  width = 10,
  height = 6,
  dpi = 300
)

done_file <- file.path(output_dir, "done")
file.create(done_file)


cat("Report generation complete !")
cat("All plots and dataframes saved in:", output_dir, "\n")