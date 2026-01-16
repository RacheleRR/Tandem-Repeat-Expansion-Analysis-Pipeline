# ============================================================================
# FILE: statistics_functions.r
# PURPOSE: Statistical analysis functions for TR analysis pipeline
# 
# CONTENTS:
#   - Correlation Analysis:
#     * run_comparison_analysis_dynamic: Logistic regression for genomic features
#   
#   - Proximity Analysis:
#     * perform_distance_tests: Wilcoxon tests for distance comparisons
#     * summarize_distances: Descriptive statistics for distances
#   
#   - Burden Analysis:
#     * perform_fisher_analysis: Fisher's exact test with FDR correction
#     * rate_confront: Rate calculations per group
#     * perform_nonparametric_analysis: Wilcoxon/Kruskal-Wallis tests
#   
#   - Regression Analysis:
#     * run_binary_analysis: Binary logistic regression
#     * run_multinomial_analysis: Multinomial logistic regression
#   
#   - Enrichment Analysis:
#     * run_enrichment: gProfiler2 wrapper for ORA
#   
#   - Helper Functions:
#     * get_all_pairwise_comparisons: Generate comparison pairs
#     * perform_single_fisher_comparison: Core Fisher test
#     * perform_wilcoxon_tests, perform_kruskal_tests, perform_dunn_tests: Core nonParametric tests
#     * extract_type_from_row_name: Row name parsing
#
# STATISTICAL METHODS:
#   - Logistic regression (binomial)
#   - Fisher's exact test (2x2 contingency tables)
#   - Wilcoxon rank-sum test (Mann-Whitney U)
#   - Kruskal-Wallis test with Dunn's post-hoc
#   - FDR correction (Benjamini-Hochberg)
#   - Odds ratios with 95% confidence intervals
#
#
# AUTHOR: RRR
# DATE: 12/12/2025
# VERSION: 2.0
# ============================================================================


# load libries 
suppressPackageStartupMessages({
library(broom)
library(purrr)  # for map_dfr
library(tibble)
library(emmeans)
library(nnet)
})

#' Run logistic regression comparison analysis
#'
#' @param bin_data Data frame with genomic bins and features
#' @param type_data_list (a list of vector of counts) 
#' @param features Vector of feature names to test
#' @param adjust_method P-value adjustment method (default: "fdr")
#' @return Data frame with statistical results
run_comparison_analysis_dynamic <- function(
  bin_data,
  type_data_list,
  features = c("fragile", "GC", "PhastCons", "phylop"),
  adjust_method = "fdr"
) {
  cat("\n=== Dynamic comparison analysis ===\n")
  cat("Bins:", nrow(bin_data), "\n")
  cat("Types:", length(type_data_list), "\n")
  cat("Features:", paste(features, collapse = ", "), "\n\n")

  dt.out <- data.frame()

  for (type_name in names(type_data_list)) {
    cat("Type:", type_name, "\n")

    tmp_type <- as.factor(type_data_list[[type_name]] > 0)

    if (length(unique(tmp_type)) < 2) {
      message("Skipping type ", type_name, ": only one class present\n")
      next
    }
    
    bin_data$tmp_type <- tmp_type

    for (f in features) {

      if (!f %in% colnames(bin_data)) {
        message("Skipping feature ", f, ": not found\n")
        next
      }

      if (all(bin_data[[f]] == 0, na.rm = TRUE)) {
        message("Skipping feature ", f, ": all zeros\n")
        next
      }


      bin_data$feat <- bin_data[[f]]

      fit <- tryCatch(
        glm(tmp_type ~ feat, data = bin_data, family = binomial),
        error = function(e) NULL
      )
      if (is.null(fit)) next

      conf <- tryCatch(confint.default(fit), error = function(e) NULL)
      if (is.null(conf)) {
        cat("  confint failed for", f, "\n")
        next
      }

      coef <- coef(summary(fit))

      dt.out <- rbind(
        dt.out,
        data.frame(
          feature = f,
          Odds.ratio = exp(coef["feat", "Estimate"]),
          OR.lower = exp(conf["feat", 1]),
          OR.upper = exp(conf["feat", 2]),
          type = type_name,
          pvalue = coef["feat", "Pr(>|z|)"],
          stringsAsFactors = FALSE
        )
      )
    }
    cat("\n")
  }

  if (nrow(dt.out) == 0) {
    cat("No valid models fitted\n")
    return(data.frame())
  }

  dt.out$p_adj <- p.adjust(dt.out$pvalue, method = adjust_method)

  dt.out$stars <- cut(
    dt.out$p_adj,
    breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
    labels = c("***", "**", "*", ""),
    right = TRUE
  )

  feature_names <- c(
    fragile = "fragile sites",
    GC = "GC content",
    PhastCons = "PhastCons",
    phylop = "phyloP"
  )

  dt.out$feature_display <- factor(
    feature_names[dt.out$feature],
    levels = c("GC content", "phyloP", "PhastCons", "fragile sites")
  )

  cat("Total tests:", nrow(dt.out), "\n")
  cat("Significant (adj p < 0.05):",
      sum(dt.out$p_adj < 0.05, na.rm = TRUE), "\n")
  cat("=== Dynamic comparison complete ===\n\n")
  

  dt.out
}


#DISTANCE
perform_distance_tests <- function(distance_list, comparisons, 
                                   alternative = "less", adjust_method = "BH") {
  cat("\n=== Distance tests ===\n")
  cat("Comparisons:", length(comparisons), "\n")
  cat("Alternative:", alternative, "\n\n")

  results <- list()
  
  for (i in seq_along(comparisons)) {
    group1 <- comparisons[[i]][1]
    group2 <- comparisons[[i]][2]

    cat("Comparison", i, ":", group1, "vs", group2, "\n")
    
    if (!group1 %in% names(distance_list) || !group2 %in% names(distance_list)) {
      message("Skipping ", group1, " vs ", group2, ": group missing\n")
      next
    }

    d1 <- abs(distance_list[[group1]]$distance)
    d2 <- abs(distance_list[[group2]]$distance)

    cat("  n1 =", length(d1), "n2 =", length(d2), "\n")

    if (length(d1) < 5 || length(d2) < 5) {
      message("Skipping ", group1, " vs ", group2, ": insufficient data")
      next
    }

    # Perform Wilcoxon test
    test <- tryCatch(
            wilcox.test(d1, d2, alternative = alternative),
            error = function(e) {
            cat("  Wilcoxon failed:", e$message, "\n")
            NULL
      }
    )        

    if (is.null(test)) next
    
    results[[i]] <- data.frame(
      group1 = group1,
      group2 = group2,
      test = "Wilcoxon",
      alternative = alternative,
      W = test$statistic,
      p_value = test$p.value,
      stringsAsFactors = FALSE
    )
  }


  if (length(results) == 0) {
    cat("No valid distance tests\n")
    return(data.frame())
  }
  # Combine results
  results_df <- do.call(rbind, results)
  
  # Adjust p-values
  if (!is.null(adjust_method) && nrow(results_df) > 0) {
    results_df$p_adj <- p.adjust(results_df$p_value, method = adjust_method)
    
    # Add significance stars
    results_df$signif <- cut(
      results_df$p_adj,
      breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
      labels = c("***", "**", "*", ""),
      right = TRUE
    )
  }

  cat("Total tests:", nrow(results_df), "\n")
  cat("Significant (adj p < 0.05):",
      sum(results_df$p_adj < 0.05, na.rm = TRUE), "\n")
  cat("Distance tests complete\n\n")
  return(results_df)
}

#' Calculate descriptive statistics for distances
#'
#' @param distance_data Data frame with distance measurements
#' @param group_var Column name for grouping
#' @return Data frame with summary statistics
summarize_distances <- function(distance_data, group_var = "rarity") {
    cat("\n=== Distance summary ===\n")
    cat("Rows:", nrow(distance_data), "\n")
    cat("Grouping variable:", group_var, "\n")

    if (nrow(distance_data) == 0) {
        cat("Empty distance table\n\n")
        return(data.frame())
    }

    if (!group_var %in% colnames(distance_data)) {
        stop("Grouping variable not found: ", group_var)
    }

    summary_stats <- distance_data %>%
        group_by(!!sym(group_var)) %>%
        summarise(
        n = n(),
        mean_distance = mean(distance, na.rm = TRUE),
        median_distance = median(distance, na.rm = TRUE),
        min_distance = min(distance, na.rm = TRUE),
        max_distance = max(distance, na.rm = TRUE),
        q25 = quantile(distance, 0.25, na.rm = TRUE),
        q75 = quantile(distance, 0.75, na.rm = TRUE),
        .groups = "drop"
        )

  cat("Groups summarized:", nrow(summary_stats), "\n")
  cat("=== Distance summary complete ===\n\n")

  return(summary_stats)
}


#FISHER FOR GENESETS
#!---------------------------------------------------------------
#! FISHER test for GENOMIC BURDEN
#!---------------------------------------------------------------

#' Perform Fisher's exact tests with automatic handling of group numbers
#' 
#' @param count_table Dataframe with count data (first column = row names, other columns = group counts)
#' @param manifest_data Manifest dataframe with group information
#' @param group_type Type of grouping ("case_control", "scz_bd_converter", "fep_converter")
#' @param comparisons Specific comparisons to test (optional). If NULL, tests all possible comparisons.
#'                    Format: list(c("Group1", "Group2"), c("Group3", "Group4"))
#' @param fdr_method FDR correction method ("BH", "bonferroni", etc.)
#' @param return_format Return format: "combined" (all results together) or "separate" (list by comparison)
#' @return Either a combined dataframe or a list of dataframes with Fisher test results
perform_fisher_analysis <- function(
    count_table,  # CHANGED: count_tables -> count_table (singular)
    manifest_data,
    comparisons = NULL,
    fdr_method = "fdr",
    return_format = "combined",
    include_total_row = TRUE,
    group_counts
) {
    cat("\n")
    cat("=== Fisher's Exact Test Analysis ===\n")
    #cat("Group type:", group_counts, "\n\n")
    
    
    if (length(group_counts) < 2) {
        stop("Need at least 2 groups for Fisher test. Found: ", length(group_counts))
    }
    
    cat("Detected", length(group_counts), "groups:\n")
    for (group in names(group_counts)) {
        cat("  ", group, ":", group_counts[[group]], "individuals\n")
    }
    cat("\n")
    
    # Step 2: Determine comparisons to test
    if (is.null(comparisons)) {
        if (length(group_counts) == 2) {
            # For 2 groups, just compare them
            comparisons <- list(c(names(group_counts)[1], names(group_counts)[2]))
            cat("2 groups detected: Performing single Fisher test\n")
        } else {
            # For 3+ groups, do all pairwise comparisons
            comparisons <- get_all_pairwise_comparisons(names(group_counts))
            cat(length(group_counts), "groups detected: Performing pairwise Fisher tests\n")
            cat("Total comparisons:", length(comparisons), "\n")
        }
    }
    
    # ADD THIS CHECK FOR ALL CASES (BOTH 2 GROUPS AND 3+ GROUPS)
    if (!is.null(comparisons)) {
        # Remove any duplicate comparisons (A vs B is same as B vs A)
        # Convert each comparison to sorted character vector for uniqueness check
        comparisons_str <- sapply(comparisons, function(x) paste(sort(x), collapse = "::"))
        comparisons <- comparisons[!duplicated(comparisons_str)]
        
        # Remove any self-comparisons
        comparisons <- comparisons[sapply(comparisons, function(x) x[1] != x[2])]
    }
    
    if (length(comparisons) == 0) {
        stop("No valid comparisons to test")
    }



    # Step 3: Perform the tests
    all_results <- list()
    
    for (comp_idx in seq_along(comparisons)) {
        comp <- comparisons[[comp_idx]]
        group1 <- comp[1]
        group2 <- comp[2]
        
        cat("\nComparison ", comp_idx, "/", length(comparisons), ": ", 
            group1, " vs ", group2, "\n", sep = "")
        
        # Perform Fisher test for this comparison
        comp_results <- perform_single_fisher_comparison(
            count_table = count_table,  # CHANGED: Pass the dataframe directly
            group1 = group1,
            group2 = group2,
            total_group1 = group_counts[[group1]],
            total_group2 = group_counts[[group2]],
            fdr_method = fdr_method,
            include_total_row = include_total_row
        )
        
        if (!is.null(comp_results) && nrow(comp_results) > 0) {
            # Add comparison info
            comp_results$Comparison <- paste(group1, "vs", group2)
            all_results[[paste(group1, "vs", group2)]] <- comp_results
        }
    }
    
    # Step 4: Combine results if requested
    if (return_format == "combined" && length(all_results) > 0) {
        combined_results <- do.call(rbind, all_results)
        
        # Apply global FDR correction across ALL tests
        # Adjust within each comparison, then report both
        combined_results <- combined_results %>%
        group_by(Comparison) %>%
        mutate(p_adj_within_comparison = p.adjust(p_value, method = fdr_method)) %>%
        ungroup() %>%
        mutate(p_adj_global_all = p.adjust(p_value, method = fdr_method))    
        
        # Format p-values
        combined_results$p_value_formatted <- format.pval(combined_results$p_value, eps = 0.0001)
        combined_results$p_adj_within_comparison <- format.pval(combined_results$p_adj_within_comparison, eps = 0.0001)
        combined_results$p_adj_global_all_formatted <- format.pval(combined_results$p_adj_global_all, eps = 0.0001)
        
         # Extract clean Type from Row_Name
        # combined_results$Type <- extract_type_from_row_name(combined_results$Row_Name)
        # Clean row names
        # combined_results$Row_Name_Clean <- clean_fisher_row_names(combined_results$Row_Name)
        
        # Summary
        cat("\n")
        cat("\n=== Analysis Complete ===\n")
        cat("Total comparisons:", length(all_results), "\n")
        cat("Total tests:", nrow(combined_results), "\n")
        cat("Significant (p < 0.05):", sum(combined_results$p_value < 0.05, na.rm = TRUE), "\n")
        cat("Significant after global FDR:", sum(combined_results$p_adj_global_all < 0.05, na.rm = TRUE), "\n")
        cat("\n")
        return(combined_results)
        
    # } else if (return_format == "separate" && length(all_results) > 0) {
    #     # Clean up each separate result
    #     for (i in seq_along(all_results)) {
    #         all_results[[i]]$Type <- extract_type_from_row_name(all_results[[i]]$Row_Name)
            
    #         # Reorder columns
    #         all_results[[i]] <- all_results[[i]] %>%
    #             select(Type, Comparison, Group1, Group2, Group1_Count, Group2_Count,
    #                    Group1_Total, Group2_Total, OR, CI_low, CI_high, p_value,
    #                    p_adj_dataset, Row_Name)
    #     }

    #     # Return separate list
    #     cat("\n=== Analysis Complete ===\n")
    #     cat("Returning results for", length(all_results), "comparisons\n")
    #     return(all_results)
    } else {
        cat("\nNo valid results found.\n")
        return(NULL)
    }
}

#!---------------------------------------------------------------
#! CORE FISHER TEST FUNCTION (USED INTERNALLY)
#!---------------------------------------------------------------

#' Perform Fisher tests for a single comparison across all datasets
perform_single_fisher_comparison <- function(
    count_table,  # CHANGED: count_tables -> count_table (singular)
    group1,
    group2,
    total_group1,
    total_group2,
    fdr_method = "BH",
    include_total_row = TRUE
) {
    
    results <- data.frame()
    tests_performed <- 0
    
    # Determine which rows to process
    start_row <- ifelse(include_total_row, 1, 2)
    
    # Process each row
    for (i in start_row:nrow(count_table)) {
        row_name <- count_table[i, 1]
        
        # Get counts - use exact column names from the dataframe
        count1 <- as.numeric(count_table[i, group1])
        count2 <- as.numeric(count_table[i, group2])
        
        # Skip if NA
        if (any(is.na(c(count1, count2))) || length(count1) != 1 || length(count2) != 1) {
            next
        }
        
        # Skip if both counts are 0 (test would be meaningless)
        if (count1 == 0 && count2 == 0) {
            next
        }
        
        # Check counts don't exceed totals
        if (count1 > total_group1) {
            warning(paste("Count for", group1, "in row", i, "exceeds total:", 
                         count1, ">", total_group1))
            next
        }
        if (count2 > total_group2) {
            warning(paste("Count for", group2, "in row", i, "exceeds total:", 
                         count2, ">", total_group2))
            next
        }

        mat <- matrix(
            c(count1, total_group1 - count1,
            count2, total_group2 - count2),
            nrow = 2
            )

        if (any(rowSums(mat) == 0) || any(colSums(mat) == 0)) {
        next
        }
        
        # Perform Fisher test
        test_result <- tryCatch(
            fisher.test(mat),
            error = function(e) {
                warning(paste("Fisher test failed for row", i, ":", e$message))
                NULL
            }
        )

        if (is.null(test_result)) {
            next
        }
        
        # Determine row type
        row_type <- ifelse(i == 1, "Total_variants", "Gene_set_variants")
        
        # Create result row
        result_row <- data.frame(
            Row_Name = row_name,
            Group1 = group1,
            Group2 = group2,
            Group1_Count = count1,
            Group2_Count = count2,
            Group1_Total = total_group1,
            Group2_Total = total_group2,
            OR = round(test_result$estimate, 3),
            CI_low = round(test_result$conf.int[1], 3),
            CI_high = round(test_result$conf.int[2], 3),
            p_value = test_result$p.value,
            stringsAsFactors = FALSE
        )
        
        results <- rbind(results, result_row)
        tests_performed <- tests_performed + 1
    }
    
    # Apply FDR correction within this dataset for this comparison
    if (nrow(results) > 0) {
        results <- results %>%
            mutate(p_adj_dataset = p.adjust(p_value, method = fdr_method))
        
        cat(tests_performed, "tests performed\n")
        return(results)
    } else {
        cat("0 tests performed\n")
        return(NULL)
    }
}


#' Get all possible pairwise comparisons between groups
get_all_pairwise_comparisons <- function(groups) {
    if (length(groups) < 2) {
        return(list())
    }
    
    # Generate all unique pairs
    pairs <- combn(groups, 2, simplify = FALSE)
    return(pairs)
}


extract_type_from_row_name <- function(row_names) {
    # Remove the fixed prefix
    cleaned <- gsub("^Number of unique individuals with ", "", row_names)
    
    # Remove the last word before "variants" using sub() instead of gsub()
    # This pattern captures everything up to the last word, then replaces it
    cleaned <- sub("^(.*)\\s+\\S+\\s+variants$", "\\1 variants", cleaned)
    
    # Clean up
    cleaned <- trimws(cleaned)
    
    return(cleaned)
}

#!---------------------------------------------------------------
#! NON-PARAMETRIC TESTS FUNCTION (Wilcoxon/Kruskal)
#!---------------------------------------------------------------

#' Perform Wilcoxon or Kruskal-Wallis tests based on number of groups
#' 
#' @param df Dataframe with group column and numeric columns to test
#' @param group_col Name of group column (default: "group")
#' @param test_columns Columns to test (if NULL, uses all columns starting with "Number_")
#' @param fdr_method FDR correction method ("BH", "bonferroni", etc.)
#' @param pval_cutoff Significance cutoff for post-hoc tests (default: 0.05)
#' @param include_shapiro Include Shapiro-Wilk normality tests (default: TRUE)
#' @return List with test results and post-hoc results if applicable
perform_nonparametric_analysis <- function(
    df,
    group_col = "group",
    test_columns = NULL,
    fdr_method = "BH",
    pval_cutoff = 0.05,
    include_shapiro = TRUE
) {
    cat("=== Non-parametric Test Analysis ===\n\n")
    
    if (!group_col %in% names(df)) {
        stop("Group column '", group_col, "' not found in dataframe")
    }
    
    if (is.null(test_columns)) {
    # Take all columns except 'sample_id' and the group_col
    test_columns <- setdiff(names(df), c("sample_id", group_col))
    }
    if (length(test_columns) == 0) stop("No test columns found")

    cat("Testing", length(test_columns), "columns\n")
    
    groups <- unique(df[[group_col]])
    groups <- groups[!is.na(groups)]

    if (length(groups) < 2) {
        message("Skipping nonparametric test: <2 groups")
        return(data.frame())
    }
    
    cat("Detected", length(groups), "groups:", paste(groups, collapse = ", "), "\n\n")
    
    if (length(groups) == 2) {
        # ----------------------------------------------------------
        # Wilcoxon (2 groups)
        # ----------------------------------------------------------
        cat("2 groups detected: Performing Wilcoxon tests\n")
        
        results <- perform_wilcoxon_tests(
            df = df,
            group_col = group_col,
            test_columns = test_columns,
            groups = groups,
            include_shapiro = include_shapiro,
            fdr_method = fdr_method
        )
        
        if (!is.null(results) && nrow(results) > 0) {
            results <- results %>%
                mutate(
                    global_primary_p_adj   = p.adjust(primary_p,   method = fdr_method),
                    global_secondary_p_adj = p.adjust(secondary_p, method = fdr_method),
                    test_type = "wilcoxon",
                    posthoc   = FALSE
                )
        }
        
        return(results)
        
    } else if (length(groups) >= 3) {
        # ----------------------------------------------------------
        # Kruskal-Wallis (≥3 groups)
        # ----------------------------------------------------------
        cat(length(groups), "groups detected: Performing Kruskal-Wallis tests\n")
        
        kw <- perform_kruskal_tests(
            df = df,
            group_col = group_col,
            test_columns = test_columns,
            groups = groups,
            fdr_method = fdr_method
        )
        
        if (!is.null(kw) && nrow(kw) > 0) {
            kw <- kw %>%
                mutate(
                    global_kruskal_p_adj = p.adjust(kruskal_p, method = fdr_method),
                    test_type = "kruskal",
                    posthoc   = FALSE
                )
            
            significant <- kw %>%
                filter(global_kruskal_p_adj < pval_cutoff)
            
            cat("\nSignificant Kruskal-Wallis results (p_adj <", pval_cutoff, "):", 
                nrow(significant), "\n")
            
            dunn_res <- NULL
            if (nrow(significant) > 0) {
                cat("Running Dunn's post-hoc tests...\n")
                
                dunn_res <- perform_dunn_tests(
                    df         = df,
                    group_col  = group_col,
                    significant_metrics = significant$metric,
                    fdr_method = fdr_method
                )
                
                # standardize
                dunn_res <- dunn_res %>%
                    mutate(
                        test_type = "dunn_posthoc",
                        posthoc   = TRUE
                    )
            }
            
            # ------------------------------------------------------
            # Combine Kruskal + Dunn into a *single dataframe*
            # ------------------------------------------------------
            combined <- bind_rows(kw, dunn_res)
            return(combined)
        }
        
    } else {
        stop("Need at least 2 groups for analysis. Found: ", length(groups))
    }
}


#!---------------------------------------------------------------
#! HELPER FUNCTIONS
#!---------------------------------------------------------------

#' Perform Wilcoxon tests for 2 groups
perform_wilcoxon_tests <- function(df, group_col, test_columns, groups, include_shapiro = TRUE, fdr_method = "BH") {
    
    group1 <- groups[1]
    group2 <- groups[2]
    
    results <- map_dfr(test_columns, function(col) {
        # Initialize result row
        result_row <- tibble(
            metric = col,
            group1 = group1,
            group2 = group2,
            primary_test = NA_character_,
            primary_p = NA_real_,
            secondary_test = NA_character_,
            secondary_p = NA_real_,
            shapiro_group1_p = NA_real_,
            shapiro_group2_p = NA_real_,
            note = NA_character_
        )
        
        # Split data
        data1 <- df[[col]][df[[group_col]] == group1]
        data2 <- df[[col]][df[[group_col]] == group2]
        
        # Check for all-zero cases
        all_zero1 <- all(data1 == 0)
        all_zero2 <- all(data2 == 0)
        
        if (all_zero1 || all_zero2) {
            result_row$note <- case_when(
                all_zero1 && all_zero2 ~ "All zeros in both groups",
                all_zero1 ~ paste("All zeros in", group1, "group"),
                TRUE ~ paste("All zeros in", group2, "group")
            )
            return(result_row)
        }
        
        # Shapiro-Wilk normality tests (optional)
        if (include_shapiro) {
            safe_shapiro <- function(x) {
                if (length(unique(x)) < 3) return(list(p.value = NA))
                tryCatch(shapiro.test(x), error = function(e) list(p.value = NA))
            }
            
            shapiro1 <- safe_shapiro(data1)
            shapiro2 <- safe_shapiro(data2)
            result_row$shapiro_group1_p <- shapiro1$p.value
            result_row$shapiro_group2_p <- shapiro2$p.value
        }
        
        # Perform tests
        tryCatch({
            # Wilcoxon test
            w_res <- wilcox.test(df[[col]] ~ df[[group_col]]) %>% broom::tidy()
            
            # t-test for comparison
            t_res <- tryCatch({
                t.test(df[[col]] ~ df[[group_col]]) %>% broom::tidy()
            }, error = function(e) {
                list(p.value = NA)
            })
            
            # Determine primary test based on normality
            if (include_shapiro && 
                !is.na(result_row$shapiro_group1_p) && !is.na(result_row$shapiro_group2_p) &&
                result_row$shapiro_group1_p > 0.05 && result_row$shapiro_group2_p > 0.05 &&
                !is.na(t_res$p.value)) {
                # If normal, use t-test as primary
                result_row$primary_test <- "t-test"
                result_row$primary_p <- t_res$p.value
                result_row$secondary_test <- "wilcoxon"
                result_row$secondary_p <- w_res$p.value
            } else {
                # If not normal or t-test failed, use Wilcoxon as primary
                result_row$primary_test <- "wilcoxon"
                result_row$primary_p <- w_res$p.value
                result_row$secondary_test <- "t-test"
                result_row$secondary_p <- if (!is.na(t_res$p.value)) t_res$p.value else NA
            }
            
        }, error = function(e) {
            result_row$note <- paste("Test error:", e$message)
        })
        
        result_row
    })
    
    # Apply per-test FDR adjustments
    if (nrow(results) > 0) {
        results <- results %>%
            mutate(
                primary_p_adj = p.adjust(primary_p, method = fdr_method),
                secondary_p_adj = p.adjust(secondary_p, method = fdr_method)
            )
    }
    
    return(results)
}

#' Perform Kruskal-Wallis tests for multiple groups
perform_kruskal_tests <- function(df, group_col, test_columns, groups, fdr_method = "BH") {
    
    results <- map_dfr(test_columns, function(col) {
        # Initialize result row
        result_row <- tibble(
            metric = col,
            kruskal_statistic = NA_real_,
            kruskal_p = NA_real_,
            n_groups = length(groups),
            groups_tested = paste(groups, collapse = ", "),
            note = NA_character_
        )
        
        # Prepare data
        test_data <- df[c(group_col, col)]
        colnames(test_data) <- c("group", "value")
        test_data <- test_data %>% filter(!is.na(value))
        
        # Check for all-zero groups
        zero_groups <- test_data %>%
            group_by(group) %>%
            summarise(all_zero = all(value == 0)) %>%
            filter(all_zero) %>%
            pull(group)
        
        if (length(zero_groups) > 0) {
            result_row$note <- paste("All zeros in groups:", paste(zero_groups, collapse = ", "))
            return(result_row)
        }
        
        # Check if we have at least 2 groups with data
        groups_with_data <- unique(test_data$group)
        if (length(groups_with_data) < 2) {
            result_row$note <- paste("Insufficient groups with data:", paste(groups_with_data, collapse = ", "))
            return(result_row)
        }
        
        # Perform Kruskal-Wallis test
        tryCatch({
            kruskal_res <- kruskal.test(value ~ group, data = test_data)
            result_row$kruskal_statistic <- kruskal_res$statistic
            result_row$kruskal_p <- kruskal_res$p.value
        }, error = function(e) {
            result_row$note <- paste("Kruskal-Wallis error:", e$message)
        })
        
        result_row
    })
    
    return(results)
}

#' Perform Dunn's post-hoc tests for significant Kruskal-Wallis results
perform_dunn_tests <- function(df, group_col, significant_metrics, fdr_method = "BH") {
    
    all_dunn_results <- map_dfr(significant_metrics, function(col) {
        
        test_data <- df[c(group_col, col)]
        colnames(test_data) <- c("group", "value")
        test_data <- test_data %>% filter(!is.na(value))
        
        # Get unique groups in this metric
        groups <- unique(test_data$group)
        if (length(groups) < 2) return(NULL)
        
        # Perform Dunn's test
        tryCatch({
            if (!requireNamespace("dunn.test", quietly = TRUE)) {
                stop("Please install 'dunn.test' package: install.packages('dunn.test')")
            }
            
            dunn_res <- dunn.test::dunn.test(
                x = test_data$value,
                g = test_data$group,
                method = "none",
                list = FALSE
            )
            
            # Extract results
            comparisons <- strsplit(dunn_res$comparisons, " - ")
            
            tibble(
                metric = col,
                group1 = sapply(comparisons, function(x) x[1]),
                group2 = sapply(comparisons, function(x) x[2]),
                dunn_z = dunn_res$Z,
                dunn_p = dunn_res$P,
                comparison = paste(group1, "vs", group2)
            )
            
        }, error = function(e) {
            warning(paste("Dunn's test failed for", col, ":", e$message))
            NULL
        })
    })
    
    # Apply FDR adjustments
    if (!is.null(all_dunn_results) && nrow(all_dunn_results) > 0) {
        all_dunn_results <- all_dunn_results %>%
            group_by(metric) %>%
            mutate(
                dunn_p_adj_per_metric = p.adjust(dunn_p, method = fdr_method)
            ) %>%
            ungroup() %>%
            mutate(
                dunn_p_adj_global = p.adjust(dunn_p, method = fdr_method)
            )
    }

     return(all_dunn_results)
}

# rate confrontation 
rate_confront <- function(data, group_counts, id_column = "region") {
  cat("\n=== Rate confrontation ===\n")
  cat("Rows:", nrow(data), "\n")
  cat("ID column:", id_column, "\n")
  cat("Groups:\n")
  print(group_counts)
  cat("\n")

  # Create result dataframe starting with the ID column
  result <- data.frame(
    id = data[[id_column]],
    stringsAsFactors = FALSE
  )
  colnames(result)[1] <- id_column
  
  # Calculate rate for each group (vectorized - works on all rows at once)
#   for(group_name in names(group_counts)) {
#     result[[paste0(group_name, "_rate")]] <- data[[group_name]] / group_counts[[group_name]]
#   }

  for (group_name in names(group_counts)) {
    cat("Computing rate for group:", group_name,
        "(denominator =", group_counts[[group_name]], ")\n")
    result[[paste0(group_name, "_rate")]] <-
      data[[group_name]] / group_counts[[group_name]]
  }

  cat("Rate confrontation complete\n\n")

  return(result)
}

#REGRESSIOn
run_binary_analysis <- function(data, outcome, predictors) {
    
    if (length(unique(data[[outcome]])) < 2) {
        message("Binary outcome has only one class")
        return(list(results = data.frame()))
    }

    cat("\n=== Binary regression analysis ===\n")
    cat("Outcome:", outcome, "\n")
    cat("Predictors:", length(predictors), "\n")
    cat("Samples:", nrow(data), "\n\n")

    outcome_tab <- table(data[[outcome]])
    cat("Outcome distribution:\n")
    print(outcome_tab)
    cat("\n")

    cat("Predictor variability check:\n")
    for (p in predictors) {
        vals <- unique(data[[p]])
        cat("  ", p, "levels:", length(vals), "\n")
    }
    cat("\n")

    formula <- as.formula(
        paste(outcome, "~", paste(predictors, collapse = " + "))
    )

    cat("Fitting model:\n")
    cat(deparse(formula), "\n\n")

    model <- tryCatch(
        glm(formula, family = binomial, data = data),
        error = function(e) {
        message("Binary glm failed: ", e$message)
        NULL
        }
    )


    if (is.null(model)) {
        cat("Model fitting failed — returning empty results\n")
        return(list(model = NULL, results = data.frame()))
    }

    cat("Model fitted successfully\n")

    res <- tryCatch(
        tidy_logistic(model) %>%
        mutate(adj_p_val = p.adjust(p.value, method = "fdr")),
        error = function(e) {
        message("tidy_logistic failed: ", e$message)
        data.frame()
        }
    )
    cat("Returned coefficients:", nrow(res), "\n")
    cat("=== Binary regression complete ===\n\n")

    return(list(model = model, results = res))
}


run_multinomial_analysis <- function(data, outcome, predictors, ref = NULL, adjust_method = "fdr") {
   
    cat("\n=== Multinomial regression analysis ===\n")
    cat("Outcome:", outcome, "\n")
    cat("Predictors:", length(predictors), "\n")
    cat("Samples:", nrow(data), "\n")
    if (!is.null(ref)) cat("Reference level:", ref, "\n")
    cat("\n")

    outcome_tab <- table(data[[outcome]])
    cat("Outcome distribution:\n")
    print(outcome_tab)
    cat("\n")

    if (length(unique(data[[outcome]])) < 3) {
        message("Skipping multinomial: <3 outcome levels\n")
        return(list(
            model = NULL,
            single_vs_ref = data.frame(),
            pairwise = data.frame(),
            pairwise_per_predictor = data.frame()
        ))
    }


    
    # 1. Prepare outcome factor
    data[[outcome]] <- factor(data[[outcome]])
    if (!is.null(ref)) {
        data[[outcome]] <- relevel(data[[outcome]], ref = ref)
    }
    
    # 2. Build formula and fit model
    formula <- as.formula(
        paste(outcome, "~", paste(predictors, collapse = " + "))
    )

    cat("Fitting multinomial model:\n")
    cat(deparse(formula), "\n\n")
    
    model <- tryCatch(
        nnet::multinom(formula, data = data, Hess = TRUE, trace = FALSE),
        error = function(e) {
        message("Multinomial model failed: ", e$message)
        NULL
        }
    )


    if (is.null(model)) {
        cat("Model fitting failed — returning empty results\n")
        return(list(model = NULLL, results = data.frame()))
    }

    cat("Model fitted successfully\n")

    # 3. Single vs reference table
    single_table <- broom::tidy(model) %>%
        mutate(
            odds_ratio = exp(estimate),
            lower_OR = exp(estimate - 1.96 * std.error),
            upper_OR = exp(estimate + 1.96 * std.error),
            adj_p_val = p.adjust(2 * pnorm(-abs(estimate / std.error)), method = adjust_method)
        ) %>%
        dplyr::rename(comparison = y.level,
               predictor = term) %>%
        dplyr::select(comparison, predictor, estimate, std.error, odds_ratio, lower_OR, upper_OR, adj_p_val)

    # 4a. Pairwise contrasts marginal over all predictors
    emm_all <- tryCatch(emmeans(model, specs = "outcome"),
    error = function(e) {
        message("emmeans failed (global): ", e$message)
        NULL
    }
    )

    pairwise_all <- if (!is.null(emm_all)) {
        as.data.frame(pairs(emm_all, adjust = "none")) %>%    
        mutate(
            odds_ratio = exp(estimate),
            lower_OR = exp(estimate - 1.96 * SE),
            upper_OR = exp(estimate + 1.96 * SE)
        ) %>%
        dplyr:rename(
            comparison = contrast,
            estimate_log_odds = estimate,
            SE = SE,
            p.val = p.value
        ) %>%
        mutate(adj_p_val = p.adjust(p.val, method = adjust_method)) %>%
        dplyr::select(comparison, estimate_log_odds, SE, odds_ratio, lower_OR, upper_OR, p.val, adj_p_val)
        } else {
    data.frame()
    }
        
    # 4b. Pairwise contrasts per predictor
    pairwise_by_pred <- lapply(predictors, function(pred) {
        emm_pred <- tryCatch(emmeans(model, specs = "outcome", by = pred),error = function(e) NULL )
        
        if (is.null(emm_pred)) return(NULL)

        as.data.frame(pairs(emm_pred, adjust = "none")) %>%
            mutate(
                odds_ratio = exp(estimate),
                lower_OR = exp(estimate - 1.96 * SE),
                upper_OR = exp(estimate + 1.96 * SE),
                predictor = pred
            ) %>%
            dplyr::rename(
                comparison = contrast,
                estimate_log_odds = estimate,
                SE = SE,
                p.val = p.value
            ) %>%
            dplyr::select(predictor, comparison, estimate_log_odds, SE, odds_ratio, lower_OR, upper_OR, p.val)
    }) %>% bind_rows() %>%
        mutate(
            adj_p_val = p.adjust(p.val, method = adjust_method)  # global adjustment across all predictors
        )
    
    cat("Single-vs-reference rows:", nrow(single_table), "\n")
    cat("Pairwise (global) rows:", nrow(pairwise_all), "\n")
    cat("Pairwise per predictor rows:", nrow(pairwise_by_pred), "\n")

    cat("=== Multinomial regression complete ===\n\n")

    # 5. Return all three tables separately
    return(list(
        model = model,
        single_vs_ref = single_table,
        pairwise = pairwise_all,
        pairwise_per_predictor = pairwise_by_pred
    ))
}


# Enhanced enrichment function
run_enrichment <- function(genes) {
  cat("\n=== Gene set enrichment analysis ===\n")
  cat("Input genes:", length(genes), "\n")

  if (length(genes) < 5) {
    message("Skipping enrichment: too few genes")
    return(NULL)
  }
  # Run gprofiler2 enrichment
  gp <- tryCatch({
    cat("Running gProfiler...\n")
    gost(
      query = genes,
      organism = "hsapiens",
      sources = c("GO:BP", "GO:MF", "GO:CC"),
      correction_method = "fdr",
      significant = TRUE,
      evcodes = TRUE,
      domain_scope = "annotated"
    )
  }, error = function(e) {
    message("gProfiler failed: ", e$message)
    NULL
  })

  if (is.null(gp)) {
    cat("Enrichment failed — returning NULL\n")
    return(NULL)
  }

  cat("Enrichment successful\n")
  cat("Returned terms:", nrow(gp$result), "\n")
  cat("=== Enrichment complete ===\n\n")

  # Return the full gp object (or NULL if it failed)
  return(gp)
}

