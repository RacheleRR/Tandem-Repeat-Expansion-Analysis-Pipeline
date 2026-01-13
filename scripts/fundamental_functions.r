# ============================================================================
# FILE: fundamental_functions.r
# PURPOSE: Core utility functions for data processing and analysis
# 
# CONTENTS:
#   - Data Preparation:
#     * create_manifest_ufficial: Unify manifest formats
#     * extract_groups_from_manifest: Extract group information
#     * filter_and_annotate: Process and annotate TR data
#     * load_genesets: Load psychiatric gene sets
#     * check_outlier_label, get_outlier_labels: Label assignment
#     * filter_by_privacy: Privacy-based filtering
#   
#   - Statistical Utilities:
#     * add_group_counts: Group count calculations
#     * add_tr_counts_to_bins: Add TRs counts to genomic bins 
#   
#   - Genomic Operations:
#     * prepare_granges: GRanges object creation
#     * calculate_tss_distance: TSS proximity calculations
#     * prepare_for_splicing, getsplicedistance: Splice junction analysis
#   
#   - Burden Anlsysis:
#     * create_fisher_count_tables, count_unique_samples: Fisher test preparation
#     * create_wilcoxon_kruskal_tables,create_complete_sample_table_simple: Non-parametric test preparation
#   
#   - Regression Analysis:
#     * build_tr_features: Feature engineering for regression
#     * build_outcome_variable: Outcome variable preparation
#     * tidy_logistic: prepare & clean model 
#   
#   - Reporting:
#     * create_gene_lists_excel: Excel report generation
#     * create_gene_list_summary: Summary of genes
#     * create_count_tables: summary of TRs 
#   
#   - ORA:
#     * create_gem_file: create file for cytoscape enrichmentmap 
#
#   - General Utilities:
#     * write_session_info: Session information logging
#     * title_to_filename: Safe filename generation
#     * safe_test: Error-wrapped test execution
#     * make_data_name: create name 
#
#
# AUTHOR: Rachele RR
# DATE: 12/12/2025
# VERSION: 3
# ============================================================================

# libries 
suppressPackageStartupMessages({
library(dplyr)
library(tidyr)
library(stringr)
library(gprofiler2)
library(openxlsx)
library(readr)
library(GenomicRanges)
library(openxlsx)
})

#! check if libaries are installed 
cran_packages <- c(
  "dplyr",
  "tidyr",
  "readr",
  "data.table",
  "ggforce",
  "ggplot2",
  "ggrepel",
  "cowplot",
  "magrittr",
  "ordinal",
  "nnet",
  "broom",
  "emmeans",
  "stringr",
  "openxlsx",
  "RColorBrewer",
  "gridExtra",
  "scales",
  "purrr",
  "tibble",
  "htmlwidgets",
  "gprofiler2",
  "RCy3"
)

bioc_packages <- c(
  "GenomicRanges",
  "Biostrings",
  "BSgenome.Hsapiens.UCSC.hg38",
  "phastCons100way.UCSC.hg38",
  "clusterProfiler",
  "org.Hs.eg.db",
  "DOSE",
  "enrichplot"
)

installed <- rownames(installed.packages())

for (pkg in cran_packages) {
  if (!pkg %in% installed) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

for (pkg in bioc_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  }
}


#!creat ufficial manifest 
create_manifest_ufficial <- function(
    manifest_filtered,
    manifest_complete = NULL,
    status_column_complete = "Status"
) {
    # -----------------------------
    # Case 1: manifest_complete is NOT provided
    # -----------------------------
    if (is.null(manifest_complete)) {
        # manifest_filtered has V1 and V2 by your pipeline definition
        mf <- manifest_filtered
        if (ncol(mf) < 2)
            stop("manifest_filtered must have at least two columns (sample ID and status).")

        colnames(mf)[1:2] <- c("sample_id", "Status")

        manifest_ufficial <- mf %>%
            dplyr::select(sample_id, Status) %>%
            dplyr::mutate(
                sample_id = trimws(as.character(sample_id)),
                Status = trimws(as.character(Status))
            ) %>%
            dplyr::filter(sample_id != "")

        return(manifest_ufficial)
    }

    # -----------------------------
    # Case 2: manifest_complete IS provided
    # -----------------------------
    # It must contain sample_id and the chosen Status column
    if (!"sample_id" %in% names(manifest_complete))
        stop("manifest_complete must contain column 'sample_id'.")

    if (!status_column_complete %in% names(manifest_complete))
        stop(paste0("Column '", status_column_complete, "' not found in manifest_complete."))

    # Rename and retain only these two columns
    mc <- manifest_complete %>%
        dplyr::select(sample_id, Status = !!rlang::sym(status_column_complete))

    # FILTER: keep only sample IDs that exist also in manifest_filtered
    filtered_ids <- unique(trimws(as.character(manifest_filtered[[1]])))
    mc <- mc %>%
        dplyr::mutate(sample_id = trimws(as.character(sample_id)),
                      Status    = trimws(as.character(Status))) %>%
        dplyr::filter(sample_id %in% filtered_ids)

    return(mc)
}

#! get groups 
extract_groups_from_manifest <- function(manifest_ufficial) {

    groups <- sort(unique(manifest_ufficial$Status))

    counts <- manifest_ufficial %>%
        dplyr::count(Status, name = "n")

    counts_list <- setNames(as.list(counts$n), counts$Status)

    return(list(
        groups = groups,
        counts = counts_list
    ))
}


#! label

    # Function to check the label for the outlier values
    check_outlier_label <- function(outlier_value, manifest_df) {
        outlier_values <- strsplit(outlier_value, ";")[[1]]
        labels <- sapply(outlier_values, function(val) {
            if (val %in% manifest_df$sample_id) {
                return(manifest_df$Status[manifest_df$sample_id == val])
            } else {
                return(NA)
            }
        })
        unique_labels <- unique(labels)
        if (length(unique_labels) > 1) {
            return("mixed")
        } else if (length(unique_labels) == 1) {
            return(unique_labels)
        } else {
            return(NA)
        }
    }

    # Function to get the labels (case or control) for the outlier values
    get_outlier_labels <- function(outlier_value, manifest_df) {
        outlier_values <- strsplit(outlier_value, ";")[[1]]
        labels <- sapply(outlier_values, function(val) {
            if (val %in% manifest_df$sample_id) {
                return(manifest_df$Status[manifest_df$sample_id == val])
            } else {
                return(NA)
            }
        })
        return(paste(labels, collapse = ", "))
    }


#!genesets 
#! load genesets 
load_genesets <- function(geneset_dir) {



  SCHEMA <- read_csv(file.path(geneset_dir, "SCHEMA.csv"),show_col_types = FALSE)
  GWAS_120 <- read.delim(file.path(geneset_dir, "GWAS_120.csv"))
  BipEx_Bipolar <- read_csv(file.path(geneset_dir, "Bipolar_Disorder.csv"),show_col_types = FALSE)
  brain_gene_consensus_filtered_consensus_no_pitular <-
    read.delim(file.path(geneset_dir, "brain_gene_consensus_filtered_consensus_no_pitular.tsv"))
  brain_gene_consensus_ntm_consensus_no_pitular <-
    read.delim(file.path(geneset_dir, "brain_gene_consensus_ntm_consensus_no_pitular.tsv"))


  # Convert gene names to a standard format for easier comparison
  # BipEx_Bipolar
  convert_col <- gconvert(query = BipEx_Bipolar$Gene, organism = "hsapiens", target = "ENSG", mthreshold = Inf, filter_na = TRUE)
  BipEx_Bipolar <- merge(BipEx_Bipolar, convert_col[, c("target", "name")], by.x = "Gene", by.y = "target", all.x = TRUE)
  BipEx_Bipolar <- BipEx_Bipolar %>% select(Gene, name, everything())

  # SCHEMA
  convert_col <- gconvert(query = SCHEMA$Gene, organism = "hsapiens", target = "ENSG", mthreshold = Inf, filter_na = TRUE)
  SCHEMA <- merge(SCHEMA, convert_col[, c("target", "name")], by.x = "Gene", by.y = "target", all.x = TRUE)
  SCHEMA <- SCHEMA %>% select(Gene, name, everything())

  # General clean-up of data
  # Remove unnecessary columns, rename columns, and clean up gene names
  BipEx_Bipolar <- BipEx_Bipolar %>% subset(select = -Gene) %>% dplyr::rename('Gene' = name) %>% mutate(Gene = gsub(' ', '', Gene)) %>% setNames(make.names(colnames(.))) 
  BipEx_Bipolar_p_val_PTV <- BipEx_Bipolar[!is.na(BipEx_Bipolar$PTV.Fisher.p.val) & BipEx_Bipolar$PTV.Fisher.p.val <= 0.05,]
  BipEx_Bipolar_p_val_Missense <- BipEx_Bipolar[!is.na(BipEx_Bipolar$Damaging.Missense.Fisher.p.val) & BipEx_Bipolar$Damaging.Missense.Fisher.p.val <= 0.05,]

  BipEx_Bipolar_p_val_Missense <- BipEx_Bipolar_p_val_Missense %>%filter(!is.na(Damaging.Missense.Fisher.odds.ratio), Damaging.Missense.Fisher.odds.ratio > 1)
  BipEx_Bipolar_p_val_PTV <- BipEx_Bipolar_p_val_PTV %>%filter(!is.na(PTV.Fisher.odds.ratio), PTV.Fisher.odds.ratio > 1)
  SCHEMA <- SCHEMA %>% subset(select = -Gene) %>% dplyr::rename('Gene' = name) %>% mutate(Gene = gsub(' ', '', Gene)) %>% setNames(make.names(colnames(.)))  
  SCHEMA_pVal <- SCHEMA[SCHEMA$P.meta <= 0.05, ]
  SCHEMA_qVal <- SCHEMA[SCHEMA$Q.meta <= 0.05, ]
  SCHEMA_pVal <- SCHEMA_pVal[ (!is.na(SCHEMA_pVal$OR..Class.I.) & SCHEMA_pVal$OR..Class.I. > 1) | (!is.na(SCHEMA_pVal$OR..Class.II.) & SCHEMA_pVal$OR..Class.II. > 1),]
  SCHEMA_qVal <- SCHEMA_qVal[ (!is.na(SCHEMA_qVal$OR..Class.I.) & SCHEMA_qVal$OR..Class.I. > 1) | (!is.na(SCHEMA_qVal$OR..Class.II.) & SCHEMA_qVal$OR..Class.II. >1),]
  GWAS_120 <- GWAS_120 %>% dplyr::rename('Gene' = GENE_name)

  list(
    genes_schema_pval = unique(SCHEMA_pVal$Gene),
    genes_schema_qval = unique(SCHEMA_qVal$Gene),
    genes_bipolar = unique(c(BipEx_Bipolar_p_val_PTV$Gene,
                             BipEx_Bipolar_p_val_Missense$Gene)),
    genes_gwas = unique(GWAS_120$Gene),
    genes_brain = unique(brain_gene_consensus_filtered_consensus_no_pitular$Gene.name),
    genes_brain_ntpm = unique(brain_gene_consensus_ntm_consensus_no_pitular$Gene.name)
  )

}

# load custom genesets 
load_custom_genesets <- function(custom_geneset_dir, 
                                 gene_column = "Gene",
                                 file_pattern = "\\.(csv|tsv|txt)$") {

  cat("\n=== Loading custom gene sets ===\n")
  cat("Directory:", custom_geneset_dir, "\n")

  if (!dir.exists(custom_geneset_dir)) {
    stop("Custom gene set directory does not exist: ", custom_geneset_dir)
  }

  gene_set_files <- list.files(
    path = custom_geneset_dir,
    pattern = file_pattern,
    full.names = TRUE
  )

  if (length(gene_set_files) == 0) {
    stop("No gene set files found in: ", custom_geneset_dir)
  }

  gene_sets <- list()

  for (file in gene_set_files) {

    gene_set_name <- tools::file_path_sans_ext(basename(file))
    message("Reading gene set: ", gene_set_name)

    df <- tryCatch({
      if (grepl("\\.csv$", file, ignore.case = TRUE)) {
        readr::read_csv(file, show_col_types = FALSE)
      } else {
        read.delim(file, stringsAsFactors = FALSE)
      }
    }, error = function(e) {
      warning("Failed to read file: ", file)
      return(NULL)
    })

    if (is.null(df)) next

    if (!gene_column %in% colnames(df)) {
      warning("Gene column '", gene_column, "' not found in: ", file)
      next
    }

    genes <- unique(na.omit(df[[gene_column]]))

    if (length(genes) == 0) {
      warning("No valid genes found in: ", file)
      next
    }

    gene_sets[[gene_set_name]] <- genes
    message("Loaded ", length(genes), " genes")
  }

  if (length(gene_sets) == 0) {
    stop("No valid gene sets loaded from: ", custom_geneset_dir)
  }

  cat("=== Custom gene sets loaded ===\n\n")
  return(gene_sets)
}

#select genesets 
select_genesets <- function(
  geneset_mode = c("basic", "combined", "different"),
  basic_genesets,
  supplementary_genesets = NULL
) {

  geneset_mode <- match.arg(geneset_mode)

  cat("\n=== Gene set selection ===\n")
  cat("Mode:", geneset_mode, "\n")

  if (geneset_mode == "basic") {
    cat("Using basic genesets only\n")
    return(basic_genesets)
  }

  if (is.null(supplementary_genesets) || length(supplementary_genesets) == 0) {
    stop("Supplementary genesets required for mode: ", geneset_mode)
  }

  if (geneset_mode == "different") {
    cat("Using supplementary genesets only\n")
    return(supplementary_genesets)
  }

  if (geneset_mode == "combined") {
    cat("Combining basic and supplementary genesets\n")

    combined <- c(basic_genesets, supplementary_genesets)

    # Optional safety: warn on name clashes
    dup_names <- intersect(names(basic_genesets), names(supplementary_genesets))
    if (length(dup_names) > 0) {
      warning("Duplicate geneset names overridden: ",
              paste(dup_names, collapse = ", "))
    }

    return(combined)
  }
}




#! filter 

# Function to filter and annotate the data frame
filter_and_annotate <- function(df) {
    df <- df %>%
        mutate(tr_id = row_number())%>%
        dplyr::rename("Gene" = gene) %>%
        mutate(Gene = gsub('\\([^\\)]+\\)', '', Gene)) %>%
        separate_rows(Gene, sep = ',') %>%
        mutate(
            motif_length = nchar(motif),
            expansion_length = end - start,
            Expansion_Double_Motif = expansion_length >= 2 * motif_length,
            T_count = str_count(motif, 'T'),
            C_count = str_count(motif, 'C'),
            G_count = str_count(motif, 'G'),
            A_count = str_count(motif, 'A'),
            type = ifelse(nchar(motif) <= 6, 'STR', 'VNTR')
        ) %>%
        mutate(
            CG_percentage = ((C_count + G_count) / motif_length) * 100,
            AT_percentage = ((A_count + T_count) / motif_length) * 100
        ) %>%
        mutate(
            T_percentage = (T_count / motif_length) * 100,
            C_percentage = (C_count / motif_length) * 100,
            G_percentage = (G_count / motif_length) * 100,
            A_percentage = (A_count / motif_length) * 100
        ) 

    return(df)
}

#save guard
safe_test <- function(expr, label = "") {
  tryCatch(
    expr,
    error = function(e) {
      message("Skipping ", label, ": ", e$message)
      return(NULL)
    }
  )
}


#session infos 
write_session_info <- function(outdir) {
  sink(file.path(outdir, "session_info.txt"))
  cat("Date:", Sys.time(), "\n\n")
  print(sessionInfo())
  sink()
}


#distance  functions for step 8 


#' Filter TRs by privacy status
#'
#' @param data Annotated TR data
#' @param privacy_type "private" or "all"
#' @return Filtered data frame
filter_by_privacy <- function(data, privacy_type = "all") {
  if (privacy_type == "private") {
    data <- data %>% filter(count == 1)
  }
  return(data)
}

#' Add TR counts to genomic bins
#' @param bin_df Data frame of genomic bins
#' @param bin_gr GRanges object of genomic bins
#' @param tr_data Data frame of TRs
#' @param col_name Name of the new column to add
#' @return Updated data frame of genomic bins
add_tr_counts_to_bins <- function(bin_df, bin_gr, tr_data, col_name) {
  tr.g <- GRanges(tr_data$contig, IRanges(tr_data$start, tr_data$end), "*")
  olap <- data.frame(findOverlaps(bin_gr, tr.g))
  olap <- aggregate(subjectHits ~ queryHits, olap, length)
  bin_df[[col_name]] <- 0
  bin_df[[col_name]][olap$queryHits] <- olap$subjectHits
  return(bin_df)
}

# Helper: convert title to safe file name
title_to_filename <- function(title) {
  fname <- gsub("[: ]+", "_", title)       # replace spaces or colons with underscore
  fname <- gsub("[^A-Za-z0-9_]", "", fname)  # remove any other special chars
  tolower(fname)                            # optional: make lowercase
  return(fname)
}



#' Clean and prepare GRanges objects for distance analysis
#'
#' @param data Data frame with columns: contig, start, end
#' @param chr_filter Vector of chromosomes to include (e.g., paste0("chr", 1:22))
#' @param reduce Logical, whether to reduce overlapping ranges
#' @return GRanges object


prepare_granges <- function(data, reduce = FALSE) {
  colnames(data) <- c("chr", "start", "end")
  chr_filter = paste0("chr", 1:22)
  data <- data[data$chr %in% chr_filter, ]
  # Create GRanges
  gr <- GRanges(data$chr, IRanges(data$start, data$end), "*")
  
  # Reduce if requested (only for known.expansion in original)
  if (reduce) {
    gr <- GenomicRanges::reduce(gr)
    # Convert back to data frame for consistency
    data <- as.data.frame(gr)[, c("seqnames", "start", "end")]
    colnames(data) <- c("chr", "start", "end")
    gr <- GRanges(data$chr, IRanges(data$start, data$end), "*")
  }
  
  return(list(gr = gr, df = data))
}
#' Get distance to nearest TSS
#'
#' @param tr_gr GRanges object of tandem repeats
#' @param refflat Data frame with gene annotation
#' @param tss_window Size of window around TSS (bp)
#' @return Data frame with distance measurements
calculate_tss_distance <- function(tr_gr, refflat, tss_window = 10000) {
  # Prepare TSS windows
  refflat$tis_start <- ifelse(refflat$V4 == "+", 
                               refflat$V5 - tss_window, 
                               refflat$V6 - tss_window)
  refflat$tis_end <- ifelse(refflat$V4 == "+", 
                             refflat$V5 + tss_window, 
                             refflat$V6 + tss_window)
  
  tss_gr <- GRanges(refflat$V3, IRanges(refflat$tis_start, refflat$tis_end))
  
  # Find overlaps
  overlaps <- data.frame(findOverlaps(tr_gr, tss_gr))
  overlaps$strand <- refflat$V4[overlaps$subjectHits]
  overlaps$tss <- ifelse(overlaps$strand == "+", 
                         refflat$tis_start[overlaps$subjectHits] + tss_window,
                         refflat$tis_end[overlaps$subjectHits] - tss_window)
  
  # Calculate distances
  tr_data <- as.data.frame(tr_gr)
  overlaps$expansion_start <- tr_data$start[overlaps$queryHits]
  overlaps$expansion_end <- tr_data$end[overlaps$queryHits]
  overlaps$mid_point <- (overlaps$expansion_start + overlaps$expansion_end) / 2
  
  overlaps$distance <- overlaps$tss - overlaps$mid_point
  overlaps$distance <- ifelse(overlaps$strand == "+", 
                              overlaps$distance, 
                              -overlaps$distance)
  
  # Keep nearest TSS for each TR
  overlaps <- overlaps[order(abs(overlaps$distance)), ]
  overlaps <- overlaps[!duplicated(overlaps$queryHits), ]
  
  return(overlaps)
}


# PREPARE FUNCTION - MATCH ORIGINAL LOGIC
prepare_for_splicing <- function(data, max_size = 10000) {
  colnames(data) <- c("chr", "start", "end")
  chr_filter = paste0("chr", 1:22)
  data <- data[data$chr %in% chr_filter, ]

  colnames(data) <- c("chr", "start", "end")
  # Add size calculation
  data$size <- data$end - data$start
  data <- data[data$size < max_size, ]
  
  # Keep unique genomic positions
  data <- unique(data[, c("chr", "start", "end")])
  
  return(data)
}

#' Calculate splicing distances for rare expansions
#' @param rare.expansion.gr GRanges object of rare expansions
#' @param introns Data frame of intron annotations
#' @param rare.expansion Data frame of rare expansions
#' @param exons Data frame of exon annotations
#' @return Data frame with splicing distance measurements   
getsplicedistance <- function(rare.expansion.g, introns, rare.expansion, exons){
        introns.g <- GRanges(introns$chr, IRanges(introns$start, introns$end), "*")
        rare.olap <- data.frame(findOverlaps(rare.expansion.g, introns.g))
        rare.olap$isoform <- introns$isoform[rare.olap$subjectHits]
        dup <- rare.olap$queryHits[which(duplicated(rare.olap[, c("queryHits", "isoform")]))]
        
        
        exons.g <- GRanges(exons$chr, IRanges(exons$start, exons$end), "*")
        rare.exon <-  data.frame(findOverlaps(rare.expansion.g, exons.g))
        rare.exon$sizeOlap <- width(pintersect(rare.expansion.g[rare.exon$queryHits], 
                                                exons.g[rare.exon$subjectHits]))
        rare.exon$sizeExon <- width(rare.expansion.g[rare.exon$queryHits])
        rare.exon <- rare.exon[rare.exon$sizeOlap == rare.exon$sizeExon, ]
        
        dup <- union(dup, rare.exon$queryHits)
        
        rare.olap <- data.frame(findOverlaps(rare.expansion.g, introns.g))
        rare.olap <- rare.olap[!rare.olap$queryHits %in% rare.olap$queryHits[dup], ]
        rare.olap$strand <- introns$strand[rare.olap$subjectHits]
        rare.olap$chr <- rare.expansion$chr[rare.olap$queryHits]
        rare.olap$start.expansion <- rare.expansion$start[rare.olap$queryHits]
        rare.olap$end.expansion <- rare.expansion$end[rare.olap$queryHits]
        rare.olap$start.intron <- introns$start[rare.olap$subjectHits]
        rare.olap$end.intron <- introns$end[rare.olap$subjectHits]
        rare.olap$mid.point <- round((rare.olap$start.expansion + rare.olap$end.expansion) / 2)
        
        rare.olap$donor.distance <- ifelse(rare.olap$strand == "+", rare.olap$mid.point - rare.olap$start.intron, 
                                            rare.olap$end.intron - rare.olap$mid.point)
        rare.olap$acceptor.distance <- ifelse(rare.olap$strand == "+", rare.olap$end.intron - rare.olap$mid.point, 
                                                rare.olap$mid.point - rare.olap$start.intron)
        
        
        rare.olap$nearest.site <- ifelse(abs(rare.olap$donor.distance) < abs(rare.olap$acceptor.distance), "Donor", "Acceptor")
        rare.olap$distance <- ifelse(rare.olap$nearest.site == "Donor", rare.olap$donor.distance, rare.olap$acceptor.distance) %>% abs()
        
        rare.olap <- rare.olap[order(abs(rare.olap$distance)), ]
        rare.olap <- rare.olap[!duplicated(rare.olap$queryHits), ]
        
        rare.olap <- rare.olap[abs(rare.olap$distance) < 10000, ]

        rare.olap$distance[which(rare.olap$nearest.site == "Acceptor")] <-
            -rare.olap$distance[which(rare.olap$nearest.site == "Acceptor")]
        return(rare.olap)
        }



add_group_counts <- function(data,group_names) {
  
  # Count from outlier_label2 (main groups)
  for (g in group_names) {
    count_col_name <- paste0("count_", g)
    data[[count_col_name]] <- str_count(data$outlier_label2, g)
  }
  
  # Add CpG column
  data$CpG <- ifelse(grepl("CG|GC|CGC|CGG|CCG|GCC|GGC|GCG", data$motif), 1, 0)
  
  return(data)
}


#' Count unique samples per region and group
#'
#' @param data TR data
#' @param manifest Manifest with sample-group mapping
#' @param sample_col Column with sample IDs in TR data (default: "outliers")
#' @param group_col Column with groups in manifest (default: "status")
#' @return Data frame with unique counts, column names are group values
count_unique_samples <- function(data, manifest, 
                                 sample_col = "outliers", 
                                 group_col = "Status") {
  
  result <- data %>%
    separate_rows(!!sym(sample_col), sep = ";") %>%
    mutate(sample_id = trimws(!!sym(sample_col))) %>%
    left_join(manifest, by = c("sample_id" = "sample_id")) %>%
    group_by(region, !!sym(group_col)) %>%
    summarise(unique_samples = n_distinct(sample_id), .groups = "drop") %>%
    pivot_wider(
      names_from = !!sym(group_col),
      values_from = unique_samples,
      values_fill = list(unique_samples = 0)
    )
  
  return(result)
}


#for wilcoxon table 
create_complete_sample_table_simple <- function(data, manifest_df, 
                                                filter_samples = NULL,
                                                tr_type = "all") {
  

  
  # 1. Prepare filter samples
  if (is.null(filter_samples)) filter_samples <- character(0)
  
  # 2. Filter out unwanted samples immediately
  manifest_df <- manifest_df %>%
    filter(!sample_id %in% filter_samples)
  
  # 3. Prepare data based on type
  if (tr_type == "CpG") {
    data <- data %>% filter(CpG == 1)
    count_col_name <- "cpg_count"
  } else {
    count_col_name <- "tr_count"
  }
  
  # 4. Unique regions
  all_regions <- unique(data$region)
  
  # 5. Count TR per sample-region
  tr_counts <- data %>%
    separate_rows(outliers, sep = ";") %>%
    mutate(sample_id = trimws(outliers)) %>%
    left_join(manifest_df, by = "sample_id") %>%
    distinct(repeatID, sample_id, region, .keep_all = TRUE) %>%
    group_by(sample_id, region, Status) %>%
    summarise(!!count_col_name := n(), .groups = "drop")
  
  # 6. Complete table
  complete_table <- manifest_df %>%
    tidyr::expand(sample_id, region = all_regions) %>%
    left_join(manifest_df, by = "sample_id") %>%
    left_join(tr_counts, by = c("sample_id", "region", "Status")) %>%
    mutate(!!count_col_name := ifelse(is.na(!!sym(count_col_name)), 0, !!sym(count_col_name)))
  
  # 7. Pivot regions as columns
  pivoted_table <- complete_table %>%
    pivot_wider(
      id_cols = c(sample_id, Status),
      names_from = region,
      values_from = !!sym(count_col_name),
      values_fill = 0
    )
  
  # 8. Total count
  pivoted_table <- pivoted_table %>%
    rowwise() %>%
    mutate(Genome_wide = sum(c_across(-c(sample_id, Status)))) %>%
    ungroup()
  
  return(pivoted_table)
}


#FOR GENESET BURDEN 
#!Create tabel for fisher 
create_count_dataframes_for_FISHER <- function(df, 
                                              data_name, 
                                              gene_lists, 
                                              private = FALSE,
                                              filter_samples = NULL,
                                              group_names =NULL
                                              ) {
  
  if (is.null(filter_samples)) filter_samples <- character(0)

  if (private) {
    df <- df %>% filter(count == 1)  # Only filter if private=TRUE
  }
  
  # Create unique TR identifier
  df <- df %>% mutate(tr_identifier = paste(contig,start,end, motif, sep = "_"))

  tr_annotations <- df %>%
    group_by(tr_identifier) %>%
    group_modify(~ add_gene_set_annotations(.x, gene_lists)) %>%
    ungroup()


  # Create expanded dataframe for counting individuals
  expanded_df <- df %>%
    separate_rows(outliers, sep = ";") %>%
    mutate(
      sample_id = trimws(outliers),
      group = outlier_label
    ) %>%
    filter(group %in% group_names) %>%
    filter(!sample_id %in% filter_samples) %>%
    distinct(tr_identifier, sample_id, .keep_all = TRUE)
  
  # Create variant dataframe for counting variants
  variant_df <- df %>%
    mutate(group = outlier_label) %>%
    filter(group %in% group_names) %>%
    distinct(tr_identifier, .keep_all = TRUE)

  # Add TR-level annotations to expanded_df
  expanded_df <- expanded_df %>%
    left_join(tr_annotations, by = "tr_identifier")
    # Add TR-level annotations to variant_df
  variant_df <- variant_df %>%
    left_join(tr_annotations, by = "tr_identifier")

  # Generate counts for all metrics
  count_results <- generate_counts(expanded_df, variant_df, group_names, gene_lists)
  
  # Create the final result dataframe
  result_df <- create_result_dataframe_simple(count_results, data_name, group_names, gene_lists)
  
  # Replace NA with 0
  result_df[is.na(result_df)] <- 0
  
  return(result_df)
}

# Add gene set annotations dynamically
add_gene_set_annotations <- function(df, gene_lists) {
    out <- lapply(gene_lists, function(glist) as.integer(any(df$Gene %in% glist)))
    tibble::as_tibble(out)
}


# Generate all count metrics
generate_counts <- function(expanded_df, variant_df, group_names, gene_lists) {
  counts <- list()
  
  counts$unique_individuals_with_variants <- expanded_df %>%
    group_by(group) %>%
    summarize(unique_individuals = n_distinct(sample_id), .groups = "drop")
  
  # Gene set specific counts
  for (gene_set_name in names(gene_lists)) {
    # Individual counts for this gene set
    counts[[paste0("unique_individuals_with_", gene_set_name, "_variants")]] <- expanded_df %>%
      filter(.data[[gene_set_name]] == 1) %>%
      group_by(group) %>%
      summarize(!!paste0("unique_individuals_", gene_set_name) := n_distinct(sample_id), .groups = "drop")
  }

  
  # Ensure all groups are represented (fill missing with zeros)
  for (count_name in names(counts)) {
    counts[[count_name]] <- counts[[count_name]] %>%
      complete(group = group_names, fill = list(0))
  }
  
  return(counts)
}

# Create the final result dataframe
create_result_dataframe_simple <- function(counts, data_name, group_names, gene_lists) {
  # Define row names based on available gene sets
  row_names <- c(
    paste("Number of unique individuals with", data_name, "variants")
  )
  
  # Add gene set specific rows
  for (gene_set_name in names(gene_lists)) {
    row_names <- c(
      row_names,
      paste("Number of unique individuals with", gene_set_name, data_name, "variants")
    )
  }
  
  # Create the result dataframe
  result_df <- data.frame(row_names = row_names)
  
  # Add columns for each group
  for (name in group_names) {
    group_data <- c(
      # Basic counts
      counts$unique_individuals_with_variants[[2]][counts$unique_individuals_with_variants$group == name],
      
      # Gene set specific counts
      unlist(lapply(names(gene_lists), function(gs) {
        c(
          counts[[paste0("unique_individuals_with_", gs, "_variants")]][[2]][counts[[paste0("unique_individuals_with_", gs, "_variants")]]$group == name]
        )
      }))
    )
    
    result_df[[name]] <- group_data
  }
  
  return(result_df)
}

# Batch processing function for multiple datasets
  
create_fisher_count_tables <- function(processed_data, data_name,
                                      gene_lists, 
                                      private = FALSE,filter_sample=NULL,group_names=NULL) {
  
  cat("=== Creating Fisher Count Tables ===\n")
  cat("Data name:", data_name, "\n")
  cat("Private variants only:", private, "\n")
  cat("Gene sets:", paste(names(gene_lists), collapse = ", "), "\n\n")
  
    
    result_tables <- create_count_dataframes_for_FISHER(
      df = processed_data,
      data_name = data_name,
      gene_lists = gene_lists,
      private = private,
      filter_samples = filter_sample,
      group_names = group_names
    )
    
    cat("done\n")
  
  
  cat("\n=== Fisher Count Tables Complete ===\n\n")
  return(result_tables)
}

#! Count table for Wilcoxon or KRUSKAL tests
create_count_dataframes_for_Wilcoxon_KRUSKAL <- function(df, 
                                                        name, 
                                                        gene_lists,
                                                        private = FALSE,
                                                        purity_filter = "pure",  # "pure", "all", or "both"
                                                        manifest_data = NULL,
                                                        filter_samples=NULL,
                                                        group_names= NULL) {
  
  # Validate purity_filter parameter
  purity_filter <- tolower(purity_filter)
  valid_filters <- c("pure", "all", "both")
  if (!purity_filter %in% valid_filters) {
    stop("purity_filter must be one of: ", paste(valid_filters, collapse = ", "))
  }
  
  cat("Using purity filter:", purity_filter, "\n")

    # 1. Prepare filter samples
  if (is.null(filter_samples)) filter_samples <- character(0)
  

 # colnames(manifest_data) <- c("Status", "sample_id", "group")
  # Validate and prepare manifest data
  if (is.null(manifest_data)) {
    stop("manifest_data must be provided")
  }

    manifest_data <- manifest_data %>%
    filter(!sample_id %in% filter_samples)
  
  # Prepare manifest based on group type
  # manifest_prepared <- prepare_manifest_for_analysis(manifest_data, group_type)
  
  if (private) {
    df <- df %>% filter(count == 1)  # Only filter if private=TRUE
  }
  
  # # Define groups based on group_type
  # group_config <- get_group_config(group_type)
  # valid_groups <- group_config$valid_groups

  # Create unique TR identifier
  df <- df %>%
        mutate(tr_identifier = paste(contig,start,end, motif, sep = "_"))
    
  # First, create a TR-level dataframe with gene set annotations
  # This will mark a TR if ANY of its genes are in the gene sets
  tr_annotations <- df %>%
    group_by(tr_identifier) %>%
    group_modify(~ add_gene_set_annotations(.x, gene_lists)) %>%
    ungroup()

  # Create expanded dataframe (all samples)
  expanded_df <- df %>%
    separate_rows(outliers, sep = ";") %>%
    mutate(sample_id = trimws(outliers)) %>%
    left_join(
      manifest_data %>% distinct(sample_id, group), 
      by = "sample_id"
    ) %>%
    filter(!sample_id %in% filter_samples) %>%
    distinct(tr_identifier, sample_id, .keep_all = TRUE)
  
  # Create pure dataframe (only samples with clean group labels)
  expanded_df_pure <- df %>% 
    separate_rows(outliers, sep = ";") %>%
    mutate(
      sample_id = trimws(outliers),
      group = outlier_label # Use pre-computed labels
    ) %>% 
    filter(group %in% group_names) %>%
    filter(!sample_id %in% filter_samples) %>%
    distinct(tr_identifier, sample_id, .keep_all = TRUE)
  
    # Add TR-level annotations to both dataframes
    expanded_df <- expanded_df %>%
        left_join(tr_annotations, by = "tr_identifier")
    
    expanded_df_pure <- expanded_df_pure %>%
        left_join(tr_annotations, by = "tr_identifier")

  
  # Count variants per individual for all samples
  result_df <- expanded_df %>%
    group_by(sample_id, group) %>%
    summarise(
      Number_of_Variants = n(),
      across(
        names(gene_lists),
        ~ sum(.x, na.rm = TRUE),
        .names = "Number_of_{.col}_Variants"
      ),
      .groups = "drop"
    )
  
  # Count variants per individual for pure samples only
  result_df_pure <- expanded_df_pure %>%  
    group_by(sample_id, group) %>%
    summarise(
      Number_of_Variants_pure = n(),
      across(
        names(gene_lists),
        ~ sum(.x, na.rm = TRUE),
        .names = "Number_of_{.col}_Variants_pure"
      ),
      .groups = "drop"
    )

  # Select which result to return based on purity_filter
  if (purity_filter == "pure") {
    cat("Selected: Pure samples only\n")
    final_result <- result_df_pure
    
  } else if (purity_filter == "all") {
    cat("Selected: All samples (pure + not pure)\n")
    final_result <- result_df
    
  } else { # purity_filter == "both"
    cat("Selected: Both pure and all samples - returning merged results\n")
    # Merge both results
    merged_results <- result_df %>%
      full_join(result_df_pure, by = "sample_id") %>%
      mutate(
        group = coalesce(group.x, group.y)  # Combine group columns
      ) %>%
      select(-group.x, -group.y) %>%
      distinct()
    
    final_result <- merged_results
  }
  
  return(final_result)
}  



# # Prepare manifest for analysis
# prepare_manifest_for_analysis <- function(manifest, group_type) {
#   manifest_clean <- manifest
  
#   # Standardize column names
#   if ("Case_control" %in% names(manifest_clean)) {
#     manifest_clean <- manifest_clean %>% rename(group = Case_control)
#   }
#   if ("Status" %in% names(manifest_clean)) {
#     manifest_clean <- manifest_clean %>% rename(group = Status)
#   }
#   if ("Sequencing_number" %in% names(manifest_clean)) {
#     manifest_clean <- manifest_clean %>% rename(sample_id = Sequencing_number)
#   }
  
#   # Recode group labels based on analysis type
#   if (group_type == "case_control") {
#     manifest_clean <- manifest_clean %>%
#       mutate(
#         group = case_when(
#           group %in% c("FEP_SCZ", "FEP_BD", "Converter") ~ "Case",
#           group == "Non_Converter" ~ "Control",
#           TRUE ~ group
#         )
#       )
# #   } else if (group_type == "scz_bd_converter") {
# #     manifest_clean <- manifest_clean %>%
# #       mutate(
# #         group = case_when(
# #           group == "FEP_SCZ" ~ "SCZ",
# #           group == "FEP_BD" ~ "BD",
# #           TRUE ~ group
# #         )
# #       )
#     }
  
#   return(manifest_clean)
# }


# Prepare complete dataset for analysis
prepare_complete_dataset <- function(manifest_data, count_results) {
  
  # # Prepare manifest
  # manifest_prepared <- prepare_manifest_for_analysis(manifest_data, group_type)
  
  # # Define groups based on group_type
  # group_config <- get_group_config(group_type)
  # valid_groups <- group_config$valid_groups
  
  # # Filter manifest to only include valid groups
  # manifest_filtered <- manifest_prepared %>%
  #   filter(group %in% valid_groups)
  
  # Merge with count results and fill missing values with 0
  complete_df <- manifest_data %>%
    left_join(count_results, by = c("sample_id", "group")) %>%
    mutate(across(starts_with("Number_of"), ~ ifelse(is.na(.), 0, .)))

  
  return(complete_df)
}



# Batch processing function for multiple datasets
create_wilcoxon_kruskal_tables <- function(processed_data,dataset_name,
                                          gene_lists, 
                                          manifest_data,
                                          private = FALSE,
                                          purity_filter = "pure",
                                          filter_sample = NULL,
                                          group_names=NULL
                                          ) {
  
  cat("=== Creating Wilcoxon/Kruskal Count Tables ===\n")
  cat("Private variants only:", private, "\n")
  cat("Gene sets:", paste(names(gene_lists), collapse = ", "), "\n\n")
  
  
  #change the name of the manifest ufficial column from status to group 
  manifest_data <- manifest_data %>% dplyr::rename(group = Status) 
    
    # Create count tables
    count_tables <- create_count_dataframes_for_Wilcoxon_KRUSKAL(
      df = processed_data,
      name = dataset_name,
      gene_lists = gene_lists,
      private = private,
      purity_filter = purity_filter,
      manifest_data = manifest_data,
      filter_samples = filter_sample,
      group_names = group_names

    )
    
    # Create complete datasets
    complete_dataset <- prepare_complete_dataset(
      manifest_data = manifest_data,
      count_results = count_tables
    )
    
    cat("done\n")
  
  
  cat("\n=== Wilcoxon/Kruskal Tables Complete ===\n\n")
  
  return(
    complete_dataset = complete_dataset
  )
}

make_data_name <- function(base = "ehdn", privacy, purity) {
  parts <- c(base)
  
  if (privacy == "private") {
    parts <- c(parts, "private")
  } else {
    parts <- c(parts, "all")
  }
  
  if (purity == "pure") {
    parts <- c(parts, "pure")
  } else {
    parts <- c(parts, "mixed")
  }
  
  paste(parts, collapse = "_")
}



#! for regression
build_tr_features <- function(tr_data,
                              region_col = "region",
                              motif_col  = "motif") {

  stopifnot(region_col %in% colnames(tr_data))
  stopifnot(motif_col  %in% colnames(tr_data))

  df <- tr_data


  cpg_regex <- "CG|GC|CGC|CGG|CCG|GCC|GGC|GCG"

  df <- df %>%
    mutate(
      CpG = as.integer(grepl(cpg_regex, .data[[motif_col]]))
    )

  df <- df %>%
    mutate(
      region_clean = strsplit(as.character(.data[[region_col]]), ";")
    )

  all_regions <- sort(unique(unlist(df$region_clean)))

  for (r in all_regions) {
    df[[r]] <- as.integer(sapply(df$region_clean, function(x) r %in% x))
  }

  df <- df %>%
    mutate(
      intergenic = as.integer(vapply(region_clean, function(x) "intergenic" %in% x, logical(1))),
      genic = 1L - intergenic
    )

  df <- df %>%
    dplyr::select(-region_clean)

  return(list(
    data = df,
    regions = all_regions
  ))
}


build_outcome_variable <- function(data,
                                   group_col,
                                   group_order) {

  stopifnot(group_col %in% colnames(data))

  # Keep only samples belonging to specified groups
  data <- data %>%
    filter(.data[[group_col]] %in% group_order)

  n_groups <- length(group_order)

  if (n_groups == 2) {

    # Binary outcome: 0 = reference, 1 = comparison
    data <- data %>%
      mutate(
        outcome = ifelse(.data[[group_col]] == group_order[2], 1L, 0L)
      )

  } else if (n_groups > 2) {

    # Multinomial outcome
    data <- data %>%
      mutate(
        outcome = factor(
          .data[[group_col]],
          levels = group_order
        )
      )

  } else {
    stop("group_order must contain at least two groups.")
  }

  return(data)
}


# tidy_logistic <- function(model, comparison = NA) {
#   broom::tidy(model, conf.int = TRUE) %>%
#     mutate(
#       odds_ratio = exp(estimate),
#       lower_OR   = exp(conf.low),
#       upper_OR   = exp(conf.high),
#       comparison = comparison
#     ) %>%
#     select(term, comparison, estimate, std.error,
#            odds_ratio, lower_OR, upper_OR, p.value)
# }

tidy_logistic <- function(model, comparison = NA) {
  # Get log-odds (with CI)
  log_odds <- broom::tidy(model, conf.int = TRUE)
  
  # Get odds ratios (exponentiated estimates with CI)
  odds <- broom::tidy(model, conf.int = TRUE, exponentiate = TRUE) %>%
    select(term, estimate, conf.low, conf.high)
  
  # Combine by term
  tab <- log_odds %>%
    left_join(odds, by = "term", suffix = c("_log", "_OR")) %>%
    mutate(
      comparison = comparison
    ) %>%
    dplyr::rename(
      predictor = term,
      log_odds = estimate_log,
      SE = std.error,
      odds_ratio = estimate_OR,
      lower_OR = conf.low_OR,
      upper_OR = conf.high_OR,
      p.value = p.value
    ) %>%
    select(predictor, comparison, log_odds, SE, odds_ratio, lower_OR, upper_OR, p.value)
  
  return(tab)
}


#14_ORA
# Function to create GEM file for Cytoscape EnrichmentMap with term size filter
# Create GEM file
create_gem_file <- function(gostres, filename){
  if(is.null(gostres) || is.null(gostres$result) || nrow(gostres$result) == 0) {
    message("No results for GEM file: ", filename)
    return(NULL)
  }
  
  gem_data <- gostres$result
  gem <- gem_data[, c("term_id","term_name","p_value","intersection")]
  colnames(gem) <- c("GO.ID","Description","p.Val","Genes")
  gem$FDR <- gem$p.Val
  gem$Phenotype <- "+1"
  gem <- gem[, c("GO.ID","Description","p.Val","FDR","Phenotype","Genes")]
  
  write.table(gem, file = filename, sep = "\t", quote = FALSE, row.names = FALSE)
  message("GEM file saved: ", filename, " (", nrow(gem), " terms)")
}


#for report 
#! GENERATE TABEL WITH ALL NUMBERS 
create_count_dataframes <- function(df, 
                                              data_name, 
                                              gene_lists, 
                                              group_type = "case_control",
                                              private = FALSE,
                                              filter_samples=NULL,
                                              group_names=NULL
                                              ) {
  
  if (private) {
    df <- df %>% filter(count == 1)  # Only filter if private=TRUE
  }

  if (is.null(filter_samples)) filter_samples <- character(0)

  


  # Create unique TR identifier
  df <- df %>% mutate(tr_identifier = paste(contig,start,end, motif, sep = "_"))


  tr_annotations <- df %>%
    group_by(tr_identifier) %>%
    group_modify(~ add_gene_set_annotations(.x, gene_lists)) %>%
    ungroup()


  # Create expanded dataframe for counting individuals
  expanded_df <- df %>%
    separate_rows(outliers, sep = ";") %>%
    mutate(
      sample_id = trimws(outliers),
      group = outlier_label
    ) %>%
    filter(group %in% group_names) %>%
    filter(!sample_id %in% filter_samples) %>%
    distinct(tr_identifier, sample_id, .keep_all = TRUE)
  
  # Create variant dataframe for counting variants
  variant_df <- df %>%
    mutate(group = outlier_label) %>%
    filter(group %in% group_names) %>%
    distinct(tr_identifier, .keep_all = TRUE)

  # Add TR-level annotations to expanded_df
  expanded_df <- expanded_df %>%
    left_join(tr_annotations, by = "tr_identifier")
    # Add TR-level annotations to variant_df
  variant_df <- variant_df %>%
    left_join(tr_annotations, by = "tr_identifier")

  # Generate counts for all metrics
  count_results <- generate_all_counts(expanded_df, variant_df, group_names, gene_lists)
  
  # Create the final result dataframe
  result_df <- create_result_dataframe(count_results, data_name, group_names, gene_lists)
  
  # Replace NA with 0
  result_df[is.na(result_df)] <- 0
  
  return(result_df)
}

# Add gene set annotations dynamically
add_gene_set_annotations <- function(df, gene_lists) {
    out <- lapply(gene_lists, function(glist) as.integer(any(df$Gene %in% glist)))
    tibble::as_tibble(out)
}


# Generate all count metrics
generate_all_counts <- function(expanded_df, variant_df, group_names, gene_lists) {
  counts <- list()
  
  # Basic variant counts
  counts$total_variants <- variant_df %>%
    group_by(group) %>%
    summarize(num_variants = n(), .groups = "drop")
  
  counts$individuals_with_variants <- expanded_df %>%
    group_by(group) %>%
    summarize(num_individuals = n(), .groups = "drop")
  
  counts$unique_individuals_with_variants <- expanded_df %>%
    group_by(group) %>%
    summarize(unique_individuals = n_distinct(sample_id), .groups = "drop")
  
  # Gene set specific counts
  for (gene_set_name in names(gene_lists)) {
    # Variant counts for this gene set
    counts[[paste0(gene_set_name, "_variants")]] <- variant_df %>%
      filter(.data[[gene_set_name]] == 1) %>%
      group_by(group) %>%
      summarize(!!paste0("num_", gene_set_name, "_variants") := n(), .groups = "drop")
    
    # Individual counts for this gene set
    counts[[paste0("individuals_with_", gene_set_name, "_variants")]] <- expanded_df %>%
      filter(.data[[gene_set_name]] == 1) %>%
      group_by(group) %>%
      summarize(!!paste0("num_individuals_", gene_set_name) := n(), .groups = "drop")
    
    counts[[paste0("unique_individuals_with_", gene_set_name, "_variants")]] <- expanded_df %>%
      filter(.data[[gene_set_name]] == 1) %>%
      group_by(group) %>%
      summarize(!!paste0("unique_individuals_", gene_set_name) := n_distinct(sample_id), .groups = "drop")
  }
  
  # Gene counts (unique genes)
  counts$genes <- expanded_df %>%
    group_by(group) %>%
    summarize(num_genes = n_distinct(Gene), .groups = "drop")
  
  for (gene_set_name in names(gene_lists)) {
    counts[[paste0(gene_set_name, "_genes")]] <- expanded_df %>%
      filter(.data[[gene_set_name]] == 1) %>%
      group_by(group) %>%
      summarize(!!paste0("num_", gene_set_name, "_genes") := n_distinct(Gene), .groups = "drop")
  }
  
  # Ensure all groups are represented (fill missing with zeros)
  for (count_name in names(counts)) {
    counts[[count_name]] <- counts[[count_name]] %>%
      complete(group = group_names, fill = list(0))
  }
  
  return(counts)
}

# Create the final result dataframe
create_result_dataframe <- function(counts, data_name, group_names, gene_lists) {
  # Define row names based on available gene sets
  row_names <- c(
    paste("Number of", data_name, "variants"),
    paste("Number of individuals with", data_name, "variants"),
    paste("Number of individuals unique with", data_name, "variants")
  )
  
  # Add gene set specific rows
  for (gene_set_name in names(gene_lists)) {
    row_names <- c(
      row_names,
      paste("Number of", gene_set_name, data_name, "variants"),
      paste("Number of individuals with", gene_set_name, data_name, "variants"),
      paste("Number of individuals unique with", gene_set_name, data_name, "variants")
    )
  }
  
  # Add gene count rows
  row_names <- c(
    row_names,
    paste("Number of genes in", data_name)
  )
  
  for (gene_set_name in names(gene_lists)) {
    row_names <- c(
      row_names,
      paste("Number of", gene_set_name, "genes in", data_name)
    )
  }
  
  # Create the result dataframe
  result_df <- data.frame(row_names = row_names)
  
  # Add columns for each group
  for (name in group_names) {
    group_data <- c(
      # Basic counts
      counts$total_variants[[2]][counts$total_variants$group == name],
      counts$individuals_with_variants[[2]][counts$individuals_with_variants$group == name],
      counts$unique_individuals_with_variants[[2]][counts$unique_individuals_with_variants$group == name],
      
      # Gene set specific counts
      unlist(lapply(names(gene_lists), function(gs) {
        c(
          counts[[paste0(gs, "_variants")]][[2]][counts[[paste0(gs, "_variants")]]$group == name],
          counts[[paste0("individuals_with_", gs, "_variants")]][[2]][counts[[paste0("individuals_with_", gs, "_variants")]]$group == name],
          counts[[paste0("unique_individuals_with_", gs, "_variants")]][[2]][counts[[paste0("unique_individuals_with_", gs, "_variants")]]$group == name]
        )
      })),
      
      # Gene counts
      counts$genes[[2]][counts$genes$group == name],
      
      # Gene set specific gene counts
      unlist(lapply(names(gene_lists), function(gs) {
        counts[[paste0(gs, "_genes")]][[2]][counts[[paste0(gs, "_genes")]]$group == name]
      }))
    )

    result_df[[name]] <- group_data
  }
  
  return(result_df)
}

# Batch processing function for multiple datasets

create_count_tables <- function(processed_data, data_name,
                                      gene_lists, 
                                      private = FALSE,
                                      filter_sample=NULL,
                                      group_names=NULL) {
  
  cat("=== Creating Count Tables ===\n")
  cat("Data_name:", data_name, "\n")
  cat("Private variants only:", private, "\n")
  cat("Gene sets:", paste(names(gene_lists), collapse = ", "), "\n\n")
  
    
    result_tables <- create_count_dataframes(
      df = processed_data,
      data_name = data_name,
      gene_lists = gene_lists,
      private = private,
      filter_samples = filter_sample,
      group_names = group_names 
    )
    
    cat("done\n")
  
  
  cat("\n=== Count Tables Complete ===\n\n")
  return(result_tables)
}



#' Create Gene Lists Excel Reports from a Single Dataframe
#'
#' This function processes a single dataframe to create Excel reports with gene lists filtered by various criteria.
#' It splits the dataframe by groups, filters by gene lists, and outputs multiple Excel files.
#'
#' @param df The input dataframe to process
#' @param dataset_name Name for the dataset (used in output filenames)
#' @param gene_lists A named list of gene vectors for filtering
#' @param output_prefix Prefix for output files (default: "group_datasets_")
#' @param group_column Column name to split by (default: "outlier_label")
#' @param gene_column Column name that stores gene names (default: "Gene")
#' @param dedup_columns Columns to consider for deduplication
#' @param selected_columns_gene_list Columns to select for gene list outputs
#' @param selected_columns Columns to select for basic outputs
#'
#' @return Invisible list containing processed datasets
#'
#' @import dplyr
#' @import openxlsx
#' @importFrom stats setNames
#' 
create_gene_lists_excel <- function(df,
                                    dataset_name = "dataset",
                                    gene_lists,
                                    output_prefix = "group_datasets_",
                                    group_column = "outlier_label",
                                    gene_column = "Gene",
                                    dedup_columns = c("contig", "start", "end", "repeatID", "Gene", "outliers"),
                                    selected_columns_gene_list = c("contig", "start", "end", "repeatID", "Gene", "region", 
                                                                   "motif", "outliers", "outlier_label", "CpG", "count",
                                                                   "type", "group"),
                                    selected_columns = c("contig", "start", "end", "repeatID", "Gene", "region", 
                                                         "motif", "outliers", "outlier_label", "CpG", "count")) {
  
  # Validate input
  if (!is.data.frame(df)) {
    stop("Input 'df' must be a dataframe")
  }
  
  if (!group_column %in% names(df)) {
    stop(sprintf("Group column '%s' not found in dataframe", group_column))
  }
  
  if (!gene_column %in% names(df)) {
    stop(sprintf("Gene column '%s' not found in dataframe", gene_column))
  }
  
  # --- Helper Functions ----------------------------------------------------
  
  # Split dataframe by group column
  split_by_group <- function(df, group_col = group_column) {
    # Remove any NA values in the group column
    df <- df %>% filter(!is.na(!!sym(group_col)))
    
    # Split by group
    groups <- unique(df[[group_col]])
    result <- list()
    
    for (group in groups) {
      group_df <- df %>% filter(!!sym(group_col) == group)
      # Deduplicate
      group_df <- distinct(group_df, across(all_of(dedup_columns)), .keep_all = TRUE)
      result[[as.character(group)]] <- group_df
    }
    
    return(result)
  }
  
  # Save unique genes for each group
  save_unique_genes_all <- function(split_data, prefix = output_prefix) {
    for (group in names(split_data)) {
      df_group <- split_data[[group]]
      if (nrow(df_group) > 0) {
        unique_df <- df_group %>% 
          distinct(across(all_of(gene_column)), .keep_all = TRUE) %>% 
          select(any_of(intersect(selected_columns, names(df_group))))
        
        file_name <- paste0(prefix, tolower(group), "_", dataset_name, ".tsv")
        write.table(unique_df, file = file_name, sep = "\t", row.names = FALSE, quote = FALSE)
        message(sprintf("Saved: %s (%d rows, %d unique genes)", 
                       file_name, nrow(unique_df), 
                       length(unique(unique_df[[gene_column]]))))
      }
    }
  }
  
  # Filter dataset by gene lists
  filter_by_gene_lists <- function(split_data, gene_lists, gene_col = gene_column) {
    results <- list()
    
    for (group in names(split_data)) {
      df_group <- split_data[[group]]
      
      if (nrow(df_group) > 0) {
        for (list_name in names(gene_lists)) {
          genes <- gene_lists[[list_name]]
          if (length(genes) > 0) {
            filtered <- df_group %>%
              filter(!!sym(gene_col) %in% genes) %>%
              mutate(type = list_name, group = group) %>% 
              select(any_of(intersect(selected_columns_gene_list, names(df_group))))
            
            if (nrow(filtered) > 0) {
              key <- paste(dataset_name, group, list_name, sep = "_")
              results[[key]] <- filtered
            }
          }
        }
      }
    }
    
    return(results)
  }
  
  # Write list of dataframes to Excel with safe sheet names
  write_list_to_excel <- function(data_list, filename) {
    wb <- createWorkbook()
    
    # Filter out empty dataframes
    data_list <- data_list[sapply(data_list, nrow) > 0]
    
    if (length(data_list) == 0) {
      message(sprintf("No data to write for %s", filename))
      return(NULL)
    }
    
    original_names <- names(data_list)
    
    # Generate safe sheet names: truncate to 31 chars and make unique
    safe_names <- substr(original_names, 1, 31)
    safe_names <- make.unique(safe_names, sep = "_")
    
    # Create mapping: original name → safe sheet name
    name_mapping <- data.frame(
      Original_Name = original_names,
      Sheet_Name = safe_names,
      Total_Rows = sapply(data_list, nrow),
      Unique_Genes = sapply(data_list, function(df) {
        if (nrow(df) > 0 && gene_column %in% names(df)) {
          n_distinct(df[[gene_column]])
        } else {
          0
        }
      }),
      stringsAsFactors = FALSE
    )
    
    # Write summary sheet
    addWorksheet(wb, "Summary")
    writeData(wb, "Summary", name_mapping)
    
    # Write each dataframe with safe sheet name
    for (i in seq_along(data_list)) {
      if (nrow(data_list[[i]]) > 0) {
        addWorksheet(wb, safe_names[i])
        writeData(wb, safe_names[i], data_list[[i]])
      }
    }
    
    saveWorkbook(wb, filename, overwrite = TRUE)
    message(sprintf("Saved Excel file: %s", filename))
    message(sprintf("\n"))
  }
  
  # Filter to unique genes
  filter_unique_genes <- function(df) {
    if (nrow(df) > 0 && gene_column %in% names(df)) {
      df %>% distinct(across(all_of(gene_column)), .keep_all = TRUE)
    } else {
      df
    }
  }
  
  # --- Main Processing -----------------------------------------------------
  message(sprintf("\n==== Excel reports ===\n" ))
  message(sprintf("Processing dataset: %s", dataset_name))
  message(sprintf("Total rows: %d", nrow(df)))
  message(sprintf("Unique groups in '%s': %s", group_column, 
                  paste(unique(df[[group_column]]), collapse = ", ")))
  
  # 1. Split dataframe by group
  split_data <- split_by_group(df)
  message(sprintf("Split into %d groups", length(split_data)))
  
  # 2. Save unique genes to TSV files
  save_unique_genes_all(split_data)
  
  # 3. Filter by gene lists
  filtered_results <- filter_by_gene_lists(split_data, gene_lists)
  message(sprintf("Created %d filtered datasets", length(filtered_results)))
  
  # 4. Create combined filtered dataset
  combined_filtered <- list()
  if (length(filtered_results) > 0) {
    combined_df <- bind_rows(filtered_results)
    if (nrow(combined_df) > 0) {
      combined_filtered[[paste0("combined_", dataset_name)]] <- combined_df
      message(sprintf("Combined filtered dataset: %d rows", nrow(combined_df)))
    }
  }
  
  # 5. Create unique versions
  unique_split_list <- lapply(split_data, filter_unique_genes)
  
  if (length(combined_filtered) > 0) {
    unique_combined_filtered <- lapply(combined_filtered, filter_unique_genes)
  } else {
    unique_combined_filtered <- list()
  }
  
  # 6. Write to Excel files
  if (length(split_data) > 0) {
    write_list_to_excel(split_data, sprintf("%s_grouped_datasets.xlsx", dataset_name))
  }
  
  if (length(unique_split_list) > 0) {
    write_list_to_excel(unique_split_list, sprintf("%s_unique_grouped_datasets.xlsx", dataset_name))
  }
  
  if (length(combined_filtered) > 0) {
    write_list_to_excel(combined_filtered, sprintf("%s_combined_filtered_by_gene_lists.xlsx", dataset_name))
  }
  
  if (length(unique_combined_filtered) > 0) {
    write_list_to_excel(unique_combined_filtered, sprintf("%s_unique_combined_filtered_by_gene_lists.xlsx", dataset_name))
  }
  
  # Return processed data invisibly
  invisible(list(
    split_data = split_data,
    filtered_results = filtered_results,
    combined_filtered = combined_filtered,
    unique_split_list = unique_split_list,
    unique_combined_filtered = unique_combined_filtered
  ))
}

#for genelists 
#' Create Detailed Gene List Summary
#'
#' Creates a comprehensive summary table showing which genes are in which gene lists,
#' their association with TRs, outliers, and other metadata.
#'
#' @param df The input dataframe containing TR information
#' @param gene_lists A named list of gene vectors for analysis
#' @param group_column Column name for outlier groups (default: "outlier_label")
#' @param gene_column Column name for gene names (default: "Gene")
#' @param sample_column Column name for sample IDs/outliers (default: "outliers")
#' @param region_column Column name for regions (default: "region")
#' @param motif_column Column name for motifs (default: "motif")
#' @param tr_id_column Column name for TR identifiers (default: "repeatID")
#' 
#' @return A list with detailed and wide format summaries
#' 
#' @import dplyr
#' @import tidyr
#'
create_gene_list_summary <- function(df,
                                     gene_lists,
                                     group_column = "outlier_label",
                                     gene_column = "Gene",
                                     sample_column = "outliers",
                                     region_column = "region",
                                     motif_column = "motif",
                                     tr_id_column = "repeatID") {
  
  # Validate input
  required_columns <- c(group_column, gene_column, sample_column, 
                       region_column, motif_column, tr_id_column)
  missing_cols <- setdiff(required_columns, names(df))
  
  if (length(missing_cols) > 0) {
    stop(sprintf("Missing required columns: %s", paste(missing_cols, collapse = ", ")))
  }
  
  # Ensure gene lists have unique names
  if (is.null(names(gene_lists))) {
    names(gene_lists) <- paste0("list_", seq_along(gene_lists))
  }
  
  # Function to split and count samples correctly
  split_and_count_samples <- function(sample_string) {
    if (is.na(sample_string) || sample_string == "") {
      return(list(samples = character(0), count = 0))
    }
    
    # Split by common delimiters: comma, semicolon, or space
    samples <- unlist(strsplit(sample_string, split = "[,; ]+"))
    
    # Remove empty strings
    samples <- samples[samples != ""]
    
    # Get unique samples
    unique_samples <- unique(samples)
    
    return(list(samples = unique_samples, count = length(unique_samples)))
  }
  
  # --- Main Processing -----------------------------------------------------
  
  message("\n===Creating detailed gene list summary...===\n")
  message(sprintf("Total rows in dataset: %d", nrow(df)))
  message(sprintf("Total unique genes: %d", n_distinct(df[[gene_column]])))
  message(sprintf("Gene lists to analyze: %s", paste(names(gene_lists), collapse = ", ")))
  
  # Get all unique groups
  groups <- unique(df[[group_column]])
  groups <- groups[!is.na(groups)]
  message(sprintf("Outlier groups found: %s", paste(groups, collapse = ", ")))
  
  all_details <- list()
  detail_counter <- 1
  
  # Process each group
  for (group in groups) {
    message(sprintf("  Processing group: %s", group))
    
    # Subset dataframe for this group
    df_group <- df %>% filter(!!sym(group_column) == group)
    
    if (nrow(df_group) == 0) next
    
    # Process each gene list for this group
    for (list_name in names(gene_lists)) {
      message(sprintf("    Analyzing gene list: %s", list_name))
      
      gene_list <- gene_lists[[list_name]]
      
      # Skip empty gene lists
      if (length(gene_list) == 0) next
      
      # Filter dataframe for genes in this list
      df_filtered <- df_group %>% 
        filter(!!sym(gene_column) %in% gene_list)
      
      if (nrow(df_filtered) == 0) next
      
      # Add gene list and process samples for each row
      df_with_list <- df_filtered %>%
        mutate(
          gene_list = list_name,
          # Process samples for each row
          sample_info = lapply(.[[sample_column]], split_and_count_samples)
        ) %>%
        rowwise() %>%
        mutate(
          sample_count = sample_info$count,
          all_samples = paste(sample_info$samples, collapse = "; ")
        ) %>%
        ungroup() %>%
        select(-sample_info)
      
      # Add to details list
      all_details[[detail_counter]] <- df_with_list
      detail_counter <- detail_counter + 1
    }
  }
  
  # Combine all details
  if (length(all_details) == 0) {
    message("No genes found in any gene list/group combination.")
    return(list(
      detailed_summary = data.frame(),
      wide_summary = data.frame()
    ))
  }
  
  detailed_summary <- bind_rows(all_details) %>%
    select(
      gene_list,
      outlier_label = !!sym(group_column),
      gene = !!sym(gene_column),
      tr_id = !!sym(tr_id_column),
      region = !!sym(region_column),
      motif = !!sym(motif_column),
      sample_count,
      all_samples
    ) %>%
    arrange(gene_list, outlier_label, gene, tr_id)
  
  # --- Create Wide Format Summary (aggregated) -----------------------------
  
  message("\nCreating wide format summary...")
  
  # Create aggregated wide format summary
  wide_summary <- detailed_summary %>%
    group_by(gene_list, outlier_label) %>%
    summarise(
      # All unique genes
      all_genes = paste(sort(unique(gene)), collapse = "; "),
      num_genes = n_distinct(gene),
      
      # All unique TRs
      all_tr_ids = paste(sort(unique(tr_id)), collapse = "; "),
      num_trs = n_distinct(tr_id),
      
      # All unique samples
      all_sample_ids = paste(unique(unlist(strsplit(all_samples, "; "))), collapse = "; "),
      num_samples = sum(sample_count),
      
      # All unique regions
      all_region_types = paste(sort(unique(region)), collapse = "; "),
      
      # All unique motifs
      all_motif_types = paste(sort(unique(motif)), collapse = "; "),
      
      .groups = 'drop'
    ) %>%
    arrange(gene_list, outlier_label)
  
  # --- Create Gene-List-Group Details (one row per gene) -------------------
  
  message("\nCreating gene-level summary...")
  
  gene_level_summary <- detailed_summary %>%
    group_by(gene_list, outlier_label, gene) %>%
    summarise(
      tr_count = n_distinct(tr_id),
      sample_count = sum(sample_count),
      tr_ids = paste(sort(unique(tr_id)), collapse = "; "),
      all_samples = paste(unique(unlist(strsplit(all_samples, "; "))), collapse = "; "),
      all_regions = paste(sort(unique(region)), collapse = "; "),
      all_motifs = paste(sort(unique(motif)), collapse = "; "),
      .groups = 'drop'
    ) %>%
    arrange(gene_list, outlier_label, desc(tr_count), gene)
  
  # --- Output Results ------------------------------------------------------
  
  message("\n=== SUMMARY STATISTICS ===")
  message(sprintf("Total TR rows in detailed summary: %d", nrow(detailed_summary)))
  message(sprintf("Total unique genes: %d", n_distinct(detailed_summary$gene)))
  message(sprintf("Total unique TRs: %d", n_distinct(detailed_summary$tr_id)))
  
  # Print a quick overview
  for (list_name in names(gene_lists)) {
    list_data <- detailed_summary %>% filter(gene_list == list_name)
    if (nrow(list_data) > 0) {
      message(sprintf("\nGene list: %s", list_name))
      for (group in groups) {
        group_data <- list_data %>% filter(outlier_label == group)
        if (nrow(group_data) > 0) {
          message(sprintf("  Group %s: %d genes, %d TRs, %d sample associations", 
                         group, n_distinct(group_data$gene),
                         n_distinct(group_data$tr_id), sum(group_data$sample_count)))
        }
      }
    }
  }
  
  # Write to Excel with three sheets
  if (requireNamespace("openxlsx", quietly = TRUE)) {

    
    wb <- createWorkbook()
    
    # 1. Detailed Summary sheet (one row per TR)
    addWorksheet(wb, "Detailed Summary (per TR)")
    writeData(wb, "Detailed Summary (per TR)", detailed_summary)
    setColWidths(wb, "Detailed Summary (per TR)", cols = 1:ncol(detailed_summary), 
                widths = c(15, 15, 20, 40, 20, 15, 10, 40))
    
    # 2. Gene-Level Summary (one row per gene)
    addWorksheet(wb, "Gene-List-Group Details")
    writeData(wb, "Gene-List-Group Details", gene_level_summary)
    setColWidths(wb, "Gene-List-Group Details", cols = 1:ncol(gene_level_summary), 
                widths = c(15, 15, 20, 10, 10, 40, 40, 30, 30))
    
    # 3. Wide Format Summary sheet (aggregated by gene list and group)
    addWorksheet(wb, "Wide Format Summary")
    writeData(wb, "Wide Format Summary", wide_summary)
    setColWidths(wb, "Wide Format Summary", cols = 1:ncol(wide_summary), 
                widths = c(15, 15, 50, 10, 50, 10, 50, 50, 50))
    
    # Save the workbook
    output_file <- file.path(output_dir, "gene_list_summary.xlsx")
    saveWorkbook(wb, output_file, overwrite = TRUE)
    message(sprintf("\nSummary saved to: %s", output_file))
    message("Sheets included:")
    message("  1. Detailed Summary (per TR) - One row per TR")
    message("  2. Gene-List-Group Details - One row per gene (aggregated)")
    message("  3. Wide Format Summary - Aggregated by gene list and group")
    message("\n")
  }
  
  # Return a list with all summaries
  return(list(
    detailed_summary = detailed_summary,      # One row per TR
    gene_level_summary = gene_level_summary,  # One row per gene
    wide_summary = wide_summary               # Aggregated by gene list and group
  ))
}
