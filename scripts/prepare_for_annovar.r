#!/usr/bin/env Rscript

# =============================================================================
# ExpansionHunterDenovo Pipeline
# =============================================================================
# prepare tsv for annovar annotation 
# Usage: ./prepare_for_annovar.R <dbscan_input> <output_dir>
# =============================================================================
suppressPackageStartupMessages(library(dplyr))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: prepare_for_annovar.R <dbscan_input> <output_dir>")
}

input_file <- args[1]
output_dir <- args[2]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Load DBSCAN output
ehdn_results <- read.delim(input_file, stringsAsFactors = FALSE)

# Reorder / rename for ANNOVAR
ehdn_results_reorder <- ehdn_results %>%
  rename(contig = chr) %>%
  select(contig, start, end, everything())

# Write to output
output_file <- file.path(output_dir, "ehdn_DBSCAN_reorder.tsv")
write.table(ehdn_results_reorder, output_file, sep = "\t", row.names = FALSE, quote = FALSE)

message("Prepared DBSCAN file for ANNOVAR:")
message(output_file)
