#!/usr/bin/env Rscript
#==================
# 9_genomic_burden.R
# PURPOSE: Genomic burden analysis of tandem repeats across groups
# 
# INPUT:
#   - args[1]: Input directory (contains data_prepare outputs)
#   - args[2]: Output directory path
#   - args[3]: Privacy filter ("private" or "all")
#   - args[4]: Purity filter ("pure" or "mixed")
#   - args[5]: samples to Exclude list 1 file path
#   - args[6]: samples to exclude list 2 file path
#   - args[7]: Include CpG analysis (TRUE/FALSE)
#
# OUTPUT:
#   - Fisher test results: ORs and p-values for each genomic region
#   - Rate confrontation tables: TR rates per group
#   - Wilcoxon/Kruskal-Wallis results: Burden comparisons
#   - Volcano, Manhattan, and dot plots for Fisher 
#   - Plots for Wilcoxon/Kruskal-Wallis results
#   - done file: Completion marker
#
# STATISTICS:
#   - Fisher's exact test for each genomic region
#   - Odds ratios with confidence intervals
#   - Rate calculations (TRs per sample)
#   - Non-parametric tests (Wilcoxon for 2 groups, Kruskal-Wallis for >2)
#
# FILTERS:
#   - Privacy: "private" (count=1) vs "all" TRs
#   - Purity: "pure" (no mixed labels) vs "mixed" (all labels)
#   - Sample exclusions: Removes specified samples from analysis
#
# AUTHOR: Rachele R. Rubiu
# DATE: 12/12/2025
# VERSION: 5.0
#====================

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
include_cpg <- as.logical(args[7])

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

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
# Load data 
# ============================================
ehdn_results_annotated <- read.delim(
  file.path(input_dir, "ehdn_results_annotated.tsv")
)

manifest_ufficial <- read.delim(
  file.path(input_dir, "manifest_ufficial.tsv")
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

#=================
# prepare data 
#===================
#Filter by privacy
if (privacy == "private") {
  tr_data <- ehdn_results_annotated %>% filter(count == 1)
} else {
  tr_data <- ehdn_results_annotated
}

if (purity == "pure"){
  tr_data <- tr_data %>% filter(!outlier_label == "mixed")
} else if (purity == "mixed"){
  tr_data <- tr_data 
}

groups_info <- extract_groups_from_manifest(manifest_ufficial)
group_names <- groups_info$groups 
group_counts <- groups_info$counts


tr_data <- add_group_counts(tr_data,group_names)

# ============================================
# Genomic burden calculation
# ============================================

unique_counts <- count_unique_samples(
  tr_data, 
  manifest_ufficial
)

unique_cpg_counts <- tr_data %>%
  filter(CpG == 1) %>%
  count_unique_samples(manifest_ufficial)


fisher <- perform_fisher_analysis(unique_counts,manifest_ufficial,group_counts=group_counts) 

if (include_cpg) {
fisher_cpg <- perform_fisher_analysis(unique_cpg_counts,manifest_ufficial,group_counts =group_counts)

}

#rate confrontation 

rates <-  rate_confront(unique_counts,group_counts=group_counts)

if (include_cpg) {
rates_cpg <- rate_confront(unique_cpg_counts,group_counts=group_counts)
}

#WILCOXON TEST FOR BURDEN

full_table_complete <- create_complete_sample_table_simple(
  data = tr_data,
  manifest_df = manifest_ufficial,
  filter_samples = all_exclude_samples,
  tr_type = "all"  # Per tutti i TR
)

if (include_cpg) {
full_table_cpg_complete <- create_complete_sample_table_simple(
  data = tr_data,
  manifest_df = manifest_ufficial,
  filter_samples = all_exclude_samples,
  tr_type = "CpG"  # Solo per CpG
)
results_nonpara_cpg_trs <- perform_nonparametric_analysis(full_table_cpg_complete,  group_col = "Status")
}

results_nonpara_trs <- perform_nonparametric_analysis(full_table_complete,  group_col = "Status")

#PLOTS
plot_fisher_volcano(fisher,Row_Name= "region" ,save_path = file.path(output_dir,"fisher_volcano.png"))

# Manhattan plot  
plot_fisher_manhattan(fisher,Row_Name="region" ,save_path = file.path(output_dir,"fisher_manhattan.png"))

# Dot plot
plot_fisher_dot(fisher, Row_Name = "region", save_path = file.path(output_dir ,"fisher_dot.png"))


plot_results <- plot_wilcoxon_kruskal_flexible(
    test_results = results_nonpara_trs,
    count_data = full_table_complete,
    group_col = "Status",
    test_type = "wilcoxon",
    plot_option = "all",
    save_dir = file.path(output_dir,"plots")
)


# ============================================
# Save results
# ============================================

write.table(rates, 
            file.path(output_dir, "rate_summary.tsv"),
            sep = "\t", row.names = FALSE)

write.table(fisher,
            file.path(output_dir, "fisher_results.tsv"),
            sep = "\t", row.names = FALSE)

write.table(results_nonpara_trs,
            file.path(output_dir, "wilcoxon_per_region.tsv"),
            sep = "\t", row.names = FALSE)

write.table(results_nonpara_trs[results_nonpara_trs$metric == "Genome_wide", ],
            file.path(output_dir, "wilcoxon_genome_wide.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

if (include_cpg) {
write.table(fisher_cpg,
            file.path(output_dir, "fisher_cpg_results.tsv"),
            sep = "\t", row.names = FALSE)


write.table(results_nonpara_cpg_trs, file.path(output_dir, "wilcoxon_cpg_per_region.tsv"),sep ="\t",row.names = FALSE)
write.table(results_nonpara_cpg_trs[results_nonpara_cpg_trs$metric == "Genome_wide", ],file.path(output_dir, "wilcoxon_cpg_genome_wide.tsv"),sep = "\t", row.names = FALSE, quote = FALSE) 
write.table(rates_cpg,file.path(output_dir, "rate_confrontation_cpg.tsv"),sep ="\t",row.names = FALSE)   
}
cat("Genomic burden analysis complete!\n")
cat("All plots and results saved in:", output_dir, "\n")

done_file <- file.path(output_dir, "done")
file.create(done_file)
