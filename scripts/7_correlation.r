#!/usr/bin/env Rscript

# ============================================
# 7_correlation.r
#
# PURPOSE: Analyze correlation between TR density and genomic features
# 
# INPUT:
#   - args[1]: Input directory (contains data_prepare outputs)
#   - args[2]: Output directory path
#   - args[3]: Privacy filter ("private" or "all")
#   - args[4]: Resource directory (genomic feature files)
#   - args[5]: Extra comparisons (comma-separated group pairs)
#
# OUTPUT:
#   - PNG plots: Odds ratio bar plots for genomic associations
#   - TSV files: Statistical results for each comparison
#   - done file: Completion marker
#
# STATISTICS:
#   - Logistic regression (binomial) for each genomic feature
#
# GENOMIC FEATURES:
#   - Fragile sites
#   - GC content
#   - PhastCons conservation scores
#   - phyloP conservation scores
#   - Known STR regions
#
# FILTERS:
#   - Privacy: "private" (count=1) vs "all" TRs
#
# AUTHOR: Rachele R. Rubiu
# DATE: 17/12/2025
# VERSION: 5.0
# ============================================
args <- commandArgs(trailingOnly = TRUE)

cat("Parameters:\n")
print(args)
cat("\n")

if (length(args) < 4) {
  stop("Usage: 7_correlation.r <input_dir> <output_dir> <privacy> <resource_dir> [extra_comparison]")
}

input_dir  <- args[1]
output_dir <- args[2]
privacy    <- args[3]
resource_dir <- args[4]
extra_comparison <- ifelse(length(args) >= 5, args[5], "")

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================
# Load functions
# ============================================

source("scripts/fundamental_functions.r")
source("scripts/statistics_functions.r")
source("scripts/plots_functions.r")

# ============================================
#Load libraries
# ============================================
suppressPackageStartupMessages({
   library(data.table)
    library(GenomicRanges)
    library(ggforce)
    library(ggplot2)
    library(ggrepel)
    library(cowplot)
    library(magrittr)
    library(phastCons100way.UCSC.hg38)
    library(ordinal)
    library(BSgenome.Hsapiens.UCSC.hg38)
    library(Biostrings)
})
# ============================================
# Load data
# ============================================

ehdn <- read.delim(
  file.path(input_dir, "ehdn_results_annotated.tsv")
)

manifest <- read.delim(
  file.path(input_dir, "manifest_ufficial.tsv")
)

# ============================================
#General data preparation 
# ============================================
gs <- phastCons100way.UCSC.hg38
bin.dt <- read.delim(file.path(resource_dir, "genome.bin.1k.tsv"), stringsAsFactors = F)
bin.dt <- bin.dt[bin.dt$seqnames %in% paste0("chr", c(1:22, "X", "Y")), ]
bin.g <- GRanges(bin.dt$seqnames, IRanges(bin.dt$start, bin.dt$end), "*")
tmp <- score(gs, bin.g)
bin.dt$PhastCons <- tmp

phylop <- read.delim(file.path(resource_dir,"subset_phylop_1kb_bins.tsv"), header = T)
phylop.g <- GRanges(phylop$seqnames, IRanges(phylop$start, phylop$end), "*")
olap <- data.frame(findOverlaps(bin.g, phylop.g))
olap$score <- phylop$phylop.max[olap$subjectHits]
olap <- aggregate(score ~ queryHits, olap, mean)
bin.dt$phylop <- NA
bin.dt$phylop[olap$queryHits] <- olap$score

fragile <- read.delim(file.path(resource_dir,"fagile_sites_hg38.tsv"),  stringsAsFactors = F)
fragile$V1 <- sub("chry", "chrY", fragile$V1)
fragile$V1 <- sub("chrx", "chrX", fragile$V1)
fragile <- fragile[fragile$V2 != "-", ]
fragile.g <- GRanges(fragile$V1, IRanges(as.numeric(fragile$V2), as.numeric(fragile$V3)), "*")
olap <- data.frame(findOverlaps(bin.g, fragile.g))
bin.dt$fragile <- 0
bin.dt$fragile[unique(olap$queryHits)] <- 1

known.str <- read.delim(file.path(resource_dir,"UCSC_simple_repeats_hg38.period_lte20.txt"), stringsAsFactors = F, header = F)
known.str <- data.frame(known.str) ### 1031708
known.str.g <- GRanges(known.str$V1, IRanges(known.str$V2, known.str$V3), "*")
olap <- data.frame(findOverlaps(bin.g, known.str.g))
bin.dt$str <- 0
bin.dt$str[unique(olap$queryHits)] <- 1

bin.dt <- bin.dt[bin.dt$seqnames %in% paste0("chr", c(1:22, "X", "Y")), ]
features <- c("fragile", "GC", "PhastCons", "phylop")
bin.dt[is.na(bin.dt)] <- 0

# ============================================
# Prepare data 
# ============================================
filtered_data <- filter_by_privacy(ehdn, privacy)

groups_info <- extract_groups_from_manifest(manifest)
group_names <- groups_info$groups



# ============================================
#! Comparisons 
# ============================================
# Add TR counts dynamically
for(g in group_names) {
    group_data <- filtered_data %>% filter(outlier_label == g)
    bin.dt <- add_tr_counts_to_bins(bin.dt, bin.g, group_data, g)
}
#Analysis 1 
type_data_list <- lapply(group_names, function(g) bin.dt[[g]])
names(type_data_list) <- group_names
result_1 <- run_comparison_analysis_dynamic(bin.dt, type_data_list, features = features)

# Result 1: all groups
title1 <- sprintf("Genomic Domain Associations: %s", paste(names(type_data_list), collapse = " vs "))
plot1 <- create_correlation_plot_dynamic(result_1, title = title1)

#Analysis 2 
type_data_list <- lapply(group_names, function(g) bin.dt[[g]])
names(type_data_list) <- group_names
type_data_list$KnownSTRs <- bin.dt$str 
result_2 <- run_comparison_analysis_dynamic(bin.dt, type_data_list, features = features)

# Result 2: all groups vs KnownSTRs
title2 <- sprintf("Genomic Domain Associations: %s", paste(names(type_data_list), collapse = " vs "))
plot2 <- create_correlation_plot_dynamic(result_2, title = title2)

# #Analysis 3 
if (!is.null(extra_comparison) && extra_comparison != "") {

  pairs <- strsplit(extra_comparison, ",")[[1]]
  comparisons <- lapply(pairs, function(p) strsplit(p, "-")[[1]])

  selected_groups <- unique(unlist(comparisons))

  # SAFETY: only keep groups that exist
  selected_groups <- selected_groups[selected_groups %in% colnames(bin.dt)]

  if (length(selected_groups) < 2) {
    message("Skipping analysis 3: <2 valid groups")
  } else {

    type_data_list <- lapply(selected_groups, function(g) bin.dt[[g]])
    names(type_data_list) <- selected_groups

    result_3 <- run_comparison_analysis_dynamic(
      bin.dt,
      type_data_list,
      features = features
    )
      
    title3 <- sprintf("Genomic Domain Associations: %s", paste(names(type_data_list), collapse = " vs "))
    plot3 <- create_correlation_plot_dynamic(result_3, title = title3)

  }
}


#!SAVE all 

# Example: saving plot1
fname1 <- title_to_filename(title1)
ggsave(file.path(output_dir, paste0("plot_", fname1, ".png")), plot1, width = 10, height = 6, dpi = 300)
write.table(result_1, file.path(output_dir, paste0("plot_", fname1, "_results.tsv")), sep = "\t", row.names = FALSE)

# Example: saving plot2
fname2 <- title_to_filename(title2)
ggsave(file.path(output_dir, paste0("plot_", fname2, ".png")), plot2, width = 10, height = 6, dpi = 300)
write.table(result_2, file.path(output_dir, paste0("plot_", fname2, "_results.tsv")), sep = "\t", row.names = FALSE)

# Optional plot3 if it exists
if (exists("plot3") && exists("result_3")) {
  fname3 <- title_to_filename(title3)
  ggsave(file.path(output_dir, paste0("plot_", fname3, ".png")), plot3, width = 10, height = 6, dpi = 300)
  write.table(result_3, file.path(output_dir, paste0("plot_", fname3, "_results.tsv")), sep = "\t", row.names = FALSE)
}

cat("All plots and results saved in:", output_dir, "\n")

done_file <- file.path(output_dir, "done")
file.create(done_file)
