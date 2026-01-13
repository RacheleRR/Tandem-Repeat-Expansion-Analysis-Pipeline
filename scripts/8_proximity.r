#!/usr/bin/env Rscript
#===========================================
# 8_proximity_analyis.r
# PURPOSE: Analyze proximity of TRs to transcriptional start sites and splice junctions
# 
# INPUT:
#   - args[1]: Input directory (contains data_prepare outputs)
#   - args[2]: Output directory path
#   - args[3]: Privacy filter ("private" or "all")
#   - args[4]: Comparison argument ("auto" or specific group pairs)
#   - args[5]: TSS window size (bp)
#   - args[6]: Splice junction window size (bp)
#   - args[7]: Resource directory (genome annotations)
#   - args[8]: Combined count directory (dbscan results)
#
# OUTPUT:
#   - TSS distance analysis: Plots and statistical results
#   - SJ distance analysis: Plots and statistical results
#   - Summary statistics tables
#   - done file: Completion marker
#
# STATISTICS:
#   - Distance calculations to nearest TSS/splice junction
#   - Wilcoxon rank-sum tests for distance comparisons
#   - Density distribution analysis
#
# FILTERS:
#   - Privacy: "private" (count=1) vs "all" TRs
#
# AUTHOR: Rachele R. Rubiu
# DATE: 17/12/2025
# VERSION: 5.0
#===========================================
args <- commandArgs(trailingOnly = TRUE)
cat("Parameters:\n")
print(args)
cat("\n")

if (length(args) < 8) {
  stop("Usage: 8_proximity.r  <input_dir> <output_dir> <privacy> <comparison_arg> <tss_window> <sj_window> <resource_dir> <combined_count_dir>")
}
input_dir  <- args[1]
output_dir <- args[2]
privacy    <- args[3]
comparison_arg <- args[4]
tss_window <- as.numeric(args[5])
sj_window <- as.numeric(args[6])
resource_dir <- args[7]
combined_count_dir <- args[8]
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

#===========================================
# Load funcctions 
#===========================================

source("scripts/fundamental_functions.r")
source("scripts/statistics_functions.r")
source("scripts/plots_functions.r")

#===========================================
# Load libraries
#===========================================
suppressPackageStartupMessages({
library(dplyr)
library(GenomicRanges)
library(ggplot2)
})

#===========================================
#load main data 
#===========================================

ehdn_results_annotated <- read.delim(
  file.path(input_dir, "ehdn_results_annotated.tsv")
)

manifest_ufficial <- read.delim(
  file.path(input_dir, "manifest_ufficial.tsv")
)

#===========================================
#prepare main data 
#===========================================

filtered_data <- filter_by_privacy(ehdn_results_annotated, privacy)
groups_info <- extract_groups_from_manifest(manifest_ufficial)
group_names <- groups_info$groups 

comparison_groups <- c(group_names, "KnownSTRs", "all_expansions")

# Define comparisons dynamically (all pairs)
if (comparison_arg == "auto" || is.na(comparison_arg) || comparison_arg == "") {
  comparisons <- combn(comparison_groups, 2, simplify = FALSE)
} else {
  pairs <- strsplit(comparison_arg, ",")[[1]]
  comparisons <- lapply(pairs, function(p) strsplit(p, "-")[[1]])
}

#===========================================
#! TSS proximity (1st type of analysis)
#===========================================
#load data
        refflat <- read.delim(file.path(resource_dir, "refFlat.txt"), stringsAsFactors = F, header = F)
        known.expansion <- read.delim(file.path(resource_dir, "UCSC_simple_repeats_hg38.period_lte20.txt"), stringsAsFactors = F, header = F)
        detected.expansion <- read.delim(file.path(combined_count_dir, "EHDn_DBSCAN.combinedCounts.bed"), header = F)

# Create a list for distance analysis
distance_list <- list()

# Add dynamic groups
for (g in group_names) {
  group_data <- filtered_data %>% filter(outlier_label == g)
  gR  <- prepare_granges(group_data, reduce = TRUE)$gr
  distance_list[[g]] <- calculate_tss_distance(gR, refflat, tss_window)
}

# Add the fixed special groups
known_strs <- prepare_granges(known.expansion,reduce =TRUE)$gr
all_expansions <- prepare_granges(detected.expansion,reduce=TRUE)$gr
distance_list$KnownSTRs <- calculate_tss_distance(known_strs, refflat, tss_window)
distance_list$all_expansions <- calculate_tss_distance(all_expansions, refflat, tss_window)

# Add privacy prefix if private
for (name in names(distance_list)) {
  if (privacy == "private") {
    distance_list[[name]]$rarity <- paste0("Private_", name)
  } else {
    distance_list[[name]]$rarity <- name
  }
  distance_list[[name]]$analysis <- "TSS"
}

# # Define comparisons dynamically (all pairs)
# comparison_groups <- c(group_names, "KnownSTRs", "all_expansions")
# comparisons <- combn(comparison_groups, 2, simplify = FALSE)
# #!should i just make it so that it needs to be specified in the configuration ?? 
# #!it either does the automatic comaprisons or it does the one we tell him in the configuration 

# Perform tests
tss_results <- perform_distance_tests(distance_list, comparisons)

# Combine for plotting
tss_combined <- do.call(rbind, distance_list)

# Prefix based on privacy
privacy_label <- ifelse(privacy == "private", "Private ", "")

# TSS title
tss_title <- sprintf("%sTSS Distance Analysis: %s", 
                     privacy_label, 
                     paste(comparison_groups, collapse = " vs "))

# Plot
tss_plot <- plot_distance_density(tss_combined, analysis_type ="TSS", analysis_name= tss_title)

fname <- title_to_filename(tss_title)

write.table(
  tss_results,
  file.path(output_dir, paste0(fname, "_results.tsv")),
  sep = "\t",
  row.names = FALSE
)

ggsave(
  file.path(output_dir, paste0(fname, ".png")),
  tss_plot,
  width = 10, height = 6, dpi = 300
)

cat("TSS analysis results and plot saved:", file.path(output_dir, fname), "\n")


#===========================================
#! splicing junctions proximity (2nd type of analysis)
#===========================================
#load data
  refflat <- read.delim(file.path(resource_dir, "refFlat.txt"), stringsAsFactors = F, header = F)
  known.expansion <- read.delim(file.path(resource_dir, "UCSC_simple_repeats_hg38.period_lte20.txt"), stringsAsFactors = F, header = F)
  introns <- read.delim(file.path(resource_dir, "hg38_intron_refFlat 1.txt"), stringsAsFactors = F)
  appris <- rbind(read.delim(file.path(resource_dir, "appris_data.principal.refseq108.hg38.txt"), stringsAsFactors = F, header = F),
                        read.delim(file.path(resource_dir, "appris_data.principal.refseq109.hg38.txt"), stringsAsFactors = F, header = F))
  exons <- read.delim(file.path(resource_dir, "hg38_exon_refFlat.txt"), stringsAsFactors = F)
  detected.expansion <- read.delim(file.path(combined_count_dir, "EHDn_DBSCAN.combinedCounts.bed"), header = F)


# Create a list for distance analysis
distance_list_SJ <- list()

# Add dynamic groups
for (g in group_names) {
  group_data <- filtered_data %>% filter(outlier_label == g)
  sj  <- prepare_for_splicing(group_data, max_size = sj_window)
  gR  <- prepare_granges(sj,reduce=FALSE)$gr
  distance_list_SJ[[g]] <- getsplicedistance(gR, introns, sj, exons)
}

appris$V3 <- sapply(sapply(appris$V3, strsplit, "\\."), "[", 1)
introns <- introns[introns$isoform %in% appris$V3, ]

all_expansions_sj <- prepare_for_splicing(detected.expansion, max_size = sj_window)
known_strs_sj <- prepare_for_splicing(known.expansion, max_size = sj_window)

known_strs_prep <- prepare_granges(known_strs_sj, reduce = TRUE)
known_strs_gr <- known_strs_prep$gr
known_strs_sj <- known_strs_prep$df  # Use reduced data frame
all_expansions_gr <- prepare_granges(all_expansions_sj, reduce = FALSE)$gr

distance_list_SJ$KnownSTRs <- getsplicedistance(known_strs_gr,  introns, known_strs_sj, exons) 
distance_list_SJ$all_expansions <- getsplicedistance(all_expansions_gr,  introns, all_expansions_sj, exons)

# Add privacy prefix if private
for (name in names(distance_list_SJ)) {
  if (privacy == "private") {
    distance_list_SJ[[name]]$rarity <- paste0("Private_", name)
  } else {
    distance_list_SJ[[name]]$rarity <- name
  }
  distance_list_SJ[[name]]$analysis <- "SJ"
}


#!should i just make it so that it needs to be specified in the configuration ?? 
#!it either does the automatic comaprisons or it does the one we tell him in the configuration 
# Perform tests
SJ_results <- perform_distance_tests(distance_list_SJ, comparisons)

# Combine for plotting
SJ_combined <- do.call(rbind, distance_list_SJ)

# Create dynamic plot title
SJ_title <- sprintf("SJ Distance Analysis: %s", paste(comparison_groups, collapse = " vs "))


# SJ title
SJ_title <- sprintf("%sSJ Distance Analysis: %s", 
                    privacy_label, 
                    paste(comparison_groups, collapse = " vs "))

# Plot
SJ_plot <- plot_distance_density(SJ_combined, analysis_type ="SJ" ,analysis_name= SJ_title)

fname <- title_to_filename(SJ_title)

write.table(
  SJ_results,
  file.path(output_dir, paste0(fname, "_results.tsv")),
  sep = "\t",
  row.names = FALSE
)

ggsave(
  file.path(output_dir, paste0(fname, ".png")),
  SJ_plot,
  width = 10, height = 6, dpi = 300
)

cat("SJ analysis results and plot saved:", file.path(output_dir, fname), "\n")

#===========================================
#! Create summary statistics
#===========================================
tss_summary <- summarize_distances(tss_combined)
sj_summary <- summarize_distances(SJ_combined)

# Save summaries
write.table(tss_summary,
            file.path(output_dir, sprintf("%s_TSS_summary.tsv", privacy)),
            sep = "\t", row.names = FALSE)

write.table(sj_summary,
            file.path(output_dir, sprintf("%s_splicing_summary.tsv", privacy)),
            sep = "\t", row.names = FALSE)


cat("Proximity analysis complete!\n")
cat(sprintf("Results saved to: %s\n", output_dir))
 

done_file <- file.path(output_dir, "done")
file.create(done_file)
