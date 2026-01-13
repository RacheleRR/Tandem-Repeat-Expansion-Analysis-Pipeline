#!/usr/bin/env Rscript
#===========================================
# 6_prepare_data
#
# PURPOSE: Prepare and annotate EHDn results for downstream analysis
# 
# INPUT:
#   - args[1]: Annotated EHDn results (ehdn_DBSCAN_annotated.tsv)
#   - args[2]: Output directory path
#   - args[3]: Manifest mode ("complete" or "filtered")
#   - args[4]: Filtered manifest file path
#   - args[5]: Complete manifest file path (optional, for mode="complete")
#
# OUTPUT:
#   - ehdn_modified.tsv: Filtered and annotated tandem repeats
#   - manifest_ufficial.tsv: Unified sample manifest with status
#   - ehdn_results_annotated.tsv: Original data with additional annotations
#
# AUTHOR: Rachele Rubiu 
# DATE: 21/12/2025
# VERSION: 4.0
#===========================================

args <- commandArgs(trailingOnly = TRUE)
cat("Parameters:\n")
print(args)
cat("\n")

if (length(args) < 4) {
  stop("Usage: 6_prepare_data.r <input> <output_dir> <manifest_mode> <manifest> [manifest_complete]")
}

input <- args[1]
output_dir <- args[2]
manifest_mode <- args[3]
manifest <- args[4]
manifest_complete <- ifelse(length(args) >= 5, args[5], "")

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

#===========================================
# Load functions & libraries
#===========================================
source("scripts/fundamental_functions.r")
suppressPackageStartupMessages({
library(dplyr)
library(tidyr)
library(readr)
})
#===========================================
# Load data 
#===========================================
ehdn_results_annotated <- read.delim(file.path(input))

manifest_filtered <- read.delim(
  file.path(manifest),
  header = FALSE
)

if (manifest_mode == "complete") {
  if (manifest_complete == "" || !file.exists(manifest_complete)) {
    stop("manifest_mode='complete' but manifest_complete not provided or does not exist")
  }

  manifest_complete_df <- read.delim(manifest_complete)

  manifest_ufficial <- create_manifest_ufficial(
    manifest_filtered,
    manifest_complete_df,
    status_column_complete = "Status"
  )
} else {
  manifest_ufficial <- create_manifest_ufficial(manifest_filtered)
}

#===========================================
# Preprocess and annotate data
#===========================================
groups_info <- extract_groups_from_manifest(manifest_ufficial)
group_names <- groups_info$groups 
group_counts <- groups_info$counts

ehdn_results_annotated$outlier_label <- sapply(ehdn_results_annotated$outliers, check_outlier_label, manifest_df = manifest_ufficial)
ehdn_results_annotated$outlier_label2<- sapply(ehdn_results_annotated$outliers, get_outlier_labels, manifest_df = manifest_ufficial)

ehdn_results_annotated <- ehdn_results_annotated  %>%
    mutate(count = sapply(strsplit(outliers, ";"), length))%>%   mutate(CpG = ifelse(grepl("CG|GC|CGC|CGG|CCG|GCC|GGC|GCG", motif), 1, 0))

ehdn_modified <- filter_and_annotate(ehdn_results_annotated)


write.table(ehdn_modified, file.path(output_dir, "ehdn_modified.tsv"), sep="\t", row.names=FALSE, quote=FALSE)
write.table(manifest_ufficial, file.path(output_dir, "manifest_ufficial.tsv"), sep="\t", row.names=FALSE, quote=FALSE)
write.table(ehdn_results_annotated, file.path(output_dir, "ehdn_results_annotated.tsv"), sep="\t", row.names=FALSE, quote=FALSE)
write_session_info(output_dir)

#===========================================
#LOGS
#===========================================
cat("\n=== DATA LOADING ===\n")
cat("Input file:", input, "\n")
cat("  Rows:", nrow(ehdn_results_annotated), "\n")
cat("  Columns:", ncol(ehdn_results_annotated), "\n")

cat("\n=== MANIFEST ===\n")
cat("Manifest filter:", manifest_mode, "\n")
cat("  Rows:", nrow(manifest_ufficial), "\n")
cat("  Columns:", ncol(manifest_ufficial), "\n")
cat("Groups:\n")
print(group_counts)
cat("\n")

cat("All dataframes saved in:", output_dir, "\n")
