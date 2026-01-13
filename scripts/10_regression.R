#!/usr/bin/env Rscript
# ============================================
# 10_regression_analysis.R
# PURPOSE: Regression analysis of TR features as predictors of group membership
# 
# INPUT:
#   - args[1]: Input directory (contains data_prepare outputs)
#   - args[2]: Output directory path
#   - args[3]: Privacy filter ("private" or "all")
#   - args[4]: Purity filter ("pure" or "mixed")
#   - args[5]: Group order (comma-separated or "auto")
#   - args[6]: Include CpG as predictor (TRUE/FALSE)
#
# OUTPUT:
#   - Binary regression results (2 groups): Odds ratios per predictor
#   - Multinomial regression results (>2 groups): Pairwise comparisons
#   - Odds ratio plots for visualization
#   - done file: Completion marker
#
# STATISTICS:
#   - Binary logistic regression (2 groups): glm with binomial family
#   - Multinomial logistic regression (>2 groups): nnet::multinom
#   - Pairwise comparisons via emmeans (>2 groups)
#
# PREDICTORS:
#   - Genomic regions (intergenic, intronic, etc.)
#   - CpG motif status (if include_cpg=TRUE)
#   - All regions as binary features
#
# FILTERS:
#   - Privacy: "private" (count=1) vs "all" TRs
#   - Purity: "pure" (no mixed labels) vs "mixed" (all labels)
#
# AUTHOR: Rachele R. Rubiu
# DATE: 12/12/2025
# VERSION: 5.0
# ============================================

args <- commandArgs(trailingOnly = TRUE)
cat("Parameters:\n")
print(args)
cat("\n")

input_dir  <- args[1]
output_dir <- args[2]
privacy    <- args[3]
purity     <- args[4]
group_order_arg <- args[5] 
include_cpg <- as.logical(args[6])

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================
# Load libraries & functions 
# ============================================
source("scripts/fundamental_functions.r")
source("scripts/statistics_functions.r")
source("scripts/plots_functions.r")

suppressPackageStartupMessages({
library(nnet)
library(broom)
library(dplyr)
library(emmeans)
})

# ============================================
# Load and prepare data
# ============================================

ehdn_modified <- read.delim(
  file.path(input_dir, "ehdn_modified.tsv")
)

manifest_ufficial <- read.delim(
  file.path(input_dir, "manifest_ufficial.tsv")
)

if (group_order_arg == "auto" || group_order_arg == "") {
  group_order <- sort(unique(manifest_ufficial$Status))
} else {
  group_order <- strsplit(group_order_arg, ",")[[1]]
}

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

groups_info <- extract_groups_from_manifest(manifest_ufficial)
group_names <- groups_info$groups 
group_counts <- groups_info$counts
n_groups     <- length(group_names)


tr_data <- build_outcome_variable(
  data = tr_data,
  group_col = "outlier_label",
  group_order = group_order
)

feature_out <- build_tr_features(tr_data)
tr_data_enriched <- feature_out$data
available_regions <- feature_out$regions

predictors <- available_regions

# Add CpG to predictors if requested
if (include_cpg) {
  predictors <- c("CpG", predictors)
}

# ===============================
# Binary regression (2 groups)
# ===============================
if (n_groups == 2) {

  binary_out <- run_binary_analysis(
    data = tr_data_enriched,
    outcome = "outcome",
    predictors = predictors
  )

  # Save full results table
  write.table(
    binary_out$results,
    file.path(output_dir, "binary_regression_results.tsv"),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )

  # Plot (auto-saved by plot_glm)
  plot <- plot_glm(
    binary_out$results,
    title = sprintf(
      "Binary Regression ORs by predictor - %s - %s",
      privacy,
      purity
    )
  )

  # CpG-only model
  if (include_cpg) {

    cpg_only_binary_out <- run_binary_analysis(
      data = tr_data_enriched,
      outcome = "outcome",
      predictors = "CpG"
    )

    write.table(
      cpg_only_binary_out$results,
      file.path(output_dir, "binary_regression_CpG_only_results.tsv"),
      sep = "\t",
      row.names = FALSE,
      quote = FALSE
    )

    plot_cpg <- plot_glm(
      cpg_only_binary_out$results,
      title = sprintf(
        "Binary Regression ORs for CpG - %s - %s",
        privacy,
        purity
      )
    )
  }
}

# ===============================
# Multinomial regression (>2 groups)
# ===============================
if (n_groups > 2) {

  multinom_out <- run_multinomial_analysis(
    data = tr_data_enriched,
    outcome = "outcome",
    predictors = predictors,
    ref = group_order[1]
  )

  write.table(
    multinom_out$pairwise_per_predictor,
    file.path(output_dir, "multinomial_regression_results.tsv"),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )

  plot <- plot_multinom_or(
    multinom_out$pairwise_per_predictor,
    title = sprintf(
      "Multinomial Regression ORs by predictor - %s - %s",
      privacy,
      purity
    )
  )

  if (include_cpg) {

    cpg_only_multinom_out <- run_multinomial_analysis(
      data = tr_data_enriched,
      outcome = "outcome",
      predictors = "CpG",
      ref = group_order[1]
    )

    write.table(
      cpg_only_multinom_out$pairwise_per_predictor,
      file.path(output_dir, "multinomial_regression_CpG_only_results.tsv"),
      sep = "\t",
      row.names = FALSE,
      quote = FALSE
    )

    plot_cpg <- plot_multinom_or(
      cpg_only_multinom_out$pairwise_per_predictor,
      title = sprintf(
        "Multinomial Regression ORs by CpG - %s - %s",
        privacy,
        purity
      )
    )
  }
}

done_file <- file.path(output_dir, "done")
file.create(done_file)

cat("Regression analysis complete!")
cat("All plots and results saved in:", output_dir, "\n")
