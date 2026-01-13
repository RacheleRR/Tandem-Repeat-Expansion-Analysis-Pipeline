#!/usr/bin/env Rscript
#======================
#!ENRICHMENT
# PURPOSE: Over-representation analysis (ORA) using gProfiler2 for gene enrichment
# 
# INPUT:
#   - args[1]: Input directory (contains data_prepare outputs)
#   - args[2]: Output directory path
#   - args[3]: Privacy filter ("private" or "all")
#   - args[4]: Purity filter ("pure" or "mixed")
#   - args[5]: GEM per group (TRUE: separate analyses, FALSE: multi-query)
#   - args[6]: Minimum term size (genes per term)
#   - args[7]: Maximum term size (genes per term)
#   - args[8]: GMT file path for Cytoscape networks
#
# OUTPUT:
#   - gProfiler2 results: TSV tables with enrichment statistics
#   - Interactive HTML and static PNG plots
#   - GEM files for Cytoscape EnrichmentMap
#   - Cytoscape network files (if Cytoscape running)
#   - Version information for reproducibility
#   - done file: Completion marker
#
# STATISTICS:
#   - Hypergeometric test for gene set enrichment
#
# ANALYSIS MODES:
#   - Single-query: Separate analysis for each group
#   - Multi-query: Combined analysis across all groups
#
# FILTERS:
#   - Privacy: "private" (count=1) vs "all" TRs
#   - Purity: "pure" (no mixed labels) vs "mixed" (all labels)
#
# AUTHOR: Rachele R. Rubiu
# DATE: 12/12/2025
# VERSION: 5.0
#======================
args <- commandArgs(trailingOnly = TRUE)
cat("Parameters:\n")
print(args)
cat("\n")

if (length(args) < 8) {
  stop("Expected 8 arguments: input_dir output_dir privacy purity gem_per_group min_term_size max_term_size gmt_file")
}

input_dir  <- args[1]
output_dir <- args[2]
privacy    <- args[3]
purity     <- args[4]
gem_per_group <- tolower(args[5]) == "true"
min_term_size  <- as.integer(args[6])
max_term_size  <- as.integer(args[7])
gmt_file <- args[8]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================
# load libraries and functions 
# ============================================
suppressPackageStartupMessages({
library(clusterProfiler)
library(org.Hs.eg.db)
library(gprofiler2)
library(dplyr)
library(ggplot2)
library(DOSE)
library(enrichplot)
library(RCy3)
library(htmlwidgets)
})

source("scripts/fundamental_functions.r")
source("scripts/plots_functions.r")
source("scripts/statistics_functions.r")

# ============================================
# Load data and prepare 
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

groups_info <- extract_groups_from_manifest(manifest_ufficial)
group_names <- groups_info$groups


# ============================================
# Enrichment: single-query or multi-query
# ============================================
if(gem_per_group) {
  # Single-query per group
  for (group in group_names) {
    message("Processing group: ", group)
    
    group_dir <- file.path(output_dir, group)
    dir.create(group_dir, recursive = TRUE, showWarnings = FALSE)
    
    group_genes <- unique(tr_data$Gene[tr_data$outlier_label == group])
    group_genes <- group_genes[!is.na(group_genes)]
    
    if(length(group_genes) == 0) {
      message("Skipping group ", group, ": no genes")
      next
    }
    
    res <- run_enrichment(group_genes)
    if(is.null(res) || nrow(res$result) == 0) {
      message("No enrichment results for group ", group)
      next
    }
    
    res$result <- res$result %>%
      filter(term_size >= min_term_size & term_size <= max_term_size)
    
    if(nrow(res$result) == 0) {
      message("No terms pass term size filter for group ", group)
      next
    }
    
    # Save results and plots
    save_gprofiler2_results(res$result, group, group_dir)
    
    # Generate plots
    png(file.path(group_dir, paste0("gprofiler_plot_", group, ".png")), width = 1200, height = 800)
    gostplot(res, capped = TRUE, interactive = FALSE)
    dev.off()
    
    htmlwidgets::saveWidget(
      gostplot(res, capped = TRUE, interactive = TRUE),
      file.path(group_dir, paste0("gprofiler_plot_", group, ".html"))
    )
    
    highlight_terms <- res$result %>%
      group_by(source) %>%
      slice_max(order_by = p_value, n = 10, with_ties = FALSE) %>%
      ungroup()
    
    publish_gosttable(
      res,
      highlight_terms = highlight_terms,
      use_colors = TRUE,
      show_columns = c("source","term_name","term_size"),
      filename = file.path(group_dir, paste0("gprofiler_table_", group, ".pdf"))
    )
    
    # Create GEM file
    gem_file <- file.path(group_dir, paste0("gProfiler_gem_", group, ".txt"))
    create_gem_file(res, gem_file)
  }
  
} else {
  # Multi-query enrichment
  message("Running multi-query enrichment for all groups...")
  
  multi_query <- lapply(group_names, function(g){
    genes <- unique(tr_data$Gene[tr_data$outlier_label == g])
    genes[!is.na(genes)]
  })
  names(multi_query) <- group_names
  
  res <- run_enrichment(multi_query)
  
  if(!is.null(res) && nrow(res$result) > 0) {
    res$result <- res$result %>%
      filter(term_size >= min_term_size & term_size <= max_term_size)
    
    # Save combined table
    save_gprofiler2_results(res$result, "multiquery", output_dir)
    
    # Generate one combined plot
    png(file.path(output_dir, "gprofiler_multiquery_plot.png"), width = 1600, height = 900)
    gostplot(res, capped = TRUE, interactive = FALSE)
    dev.off()
    
    htmlwidgets::saveWidget(
      gostplot(res, capped = TRUE, interactive = TRUE),
      file.path(output_dir, "gprofiler_multiquery_plot.html")
    )
    
    highlight_terms <- res$result %>%
      group_by(source) %>%
      slice_max(order_by = p_value, n = 20, with_ties = FALSE) %>%
      ungroup()
    
    publish_gosttable(
      res,
      highlight_terms = highlight_terms,
      use_colors = TRUE,
      show_columns = c("source","term_name","term_size"),
      filename = file.path(output_dir, "gprofiler_multiquery_table.pdf")
    )
    
    # Create GEM files per query
    gem <- res$result[, c("query", "term_id", "term_name", "p_value", "intersection")]
    colnames(gem) <- c("query","GO.ID","Description","p.Val","Genes")
    gem$FDR <- gem$p.Val
    gem$Phenotype <- "+1"
    
    gem %>% group_by(query) %>%
      group_walk(~ write.table(
        data.frame(.x[,c("GO.ID","Description","p.Val","FDR","Phenotype","Genes")]),
        file = file.path(output_dir, paste0("gProfiler_", unique(.y$query), "_gem.txt")),
        sep = "\t", quote = FALSE, row.names = FALSE
      ))
  }
}

# ============================================
# Cytoscape: build EnrichmentMap
# ============================================
if (!requireNamespace("RCy3", quietly = TRUE)) {
  message("RCy3 package not installed. Skipping Cytoscape networks.")
  message("Install with: BiocManager::install('RCy3')")
} else if (!cytoscapePing()) {
  message("Cytoscape not running. Skipping network generation.")
  message("To generate networks: Start Cytoscape and re-run this step.")
} else {
  tryCatch({
    if (gem_per_group) {
      for(group in group_names){
        group_dir <- file.path(output_dir, group)
        build_enrichmentmap_network(group_dir, group, gmt_file)
      }
    } else {
      build_enrichmentmap_network(output_dir, "combined", gmt_file)
    }
  }, error = function(e) {
    warning("Cytoscape network generation failed: ", e$message)
    warning("Networks will not be created, but analysis will continue.")
  })
}

# Save g:Profiler version info for reproducibility
version_info <- get_version_info(organism = "hsapiens")
write.table(as.data.frame(version_info), file = file.path(output_dir, "gProfiler_version_info.txt"), 
            sep = "\t", quote = FALSE, row.names = FALSE)

done_file <- file.path(output_dir, "done")
file.create(done_file)


cat("ORA network creation complete!")
cat("All plots and results saved in:", output_dir, "\n")