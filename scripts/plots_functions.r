# ============================================================================
# FILE: plots_functions.r
# PURPOSE: Visualization functions for TR analysis pipeline
# 
# CONTENTS:
#   - Correlation Plots:
#     * create_correlation_plot_dynamic: Odds ratio bar plots
#   
#   - Proximity Plots:
#     * plot_distance_density: Density plots for TSS/SJ distances
#   
#   - Burden Analysis Plots:
#     * plot_fisher_volcano: Volcano plots for Fisher results
#     * plot_fisher_manhattan: Manhattan plots for Fisher results
#     * plot_fisher_dot: Dot plots for Fisher results
#     * plot_wilcoxon_kruskal_flexible: Boxplots for non-parametric tests
#   
#   - Regression Plots:
#     * plot_glm: Odds ratio plots for binary regression
#     * plot_multinom_or: Odds ratio plots for multinomial regression
#   
#   - Enrichment Plots:
#     * build_enrichmentmap_network: Cytoscape network creation
#     * save_gprofiler2_results: gProfiler result visualization
#   
#   - Summary Plots:
#     * compute_counts, compute_nucleotide_counts: Data aggregation
#     * plot_counts, plot_nucleotides: Basic bar plots
#   
#   - Utility Functions:
#     * save_multiple_plots_grid: Grid arrangement of multiple plots
#
# VISUALIZATION FEATURES:
#   - Consistent color palettes (RColorBrewer, Set3, Set1)
#   - Significance stars (*** p<0.001, ** p<0.01, * p<0.05)
#   - Error bars for confidence intervals
#   - Dynamic faceting for multiple comparisons
#   - Automatic saving with proper dimensions and DPI
#
#
# AUTHOR: RRRubiu
# DATE: 12/12/2025
# VERSION: 2.0
# ============================================================================

suppressPackageStartupMessages({
library(RColorBrewer)
library(ggplot2)
library(ggrepel)
library(gridExtra)
library(scales)
library(dplyr)
library(stringr)
})

# Correlation plot
create_correlation_plot_dynamic <- function(results_df, title) {
  cat("\n=== Correlation Plot ===\n")
  cat("Title:", title, "\n")

  # Check if results_df is empty
  if (is.null(results_df) || nrow(results_df) == 0) {
    warning("No data provided for correlation plot")
    return(NULL)
  }
  # Check for required columns
  required_cols <- c("type", "feature_display", "Odds.ratio", "OR.lower", "OR.upper", "stars")
  missing_cols <- setdiff(required_cols, names(results_df))
  if (length(missing_cols) > 0) {
    warning("Missing columns in results_df: ", paste(missing_cols, collapse = ", "))
    return(NULL)
  }

  present_types <- unique(results_df$type)
  n_types <- length(present_types)

  cat("Detected types:", paste(present_types, collapse = ", "), "\n")

  if (n_types == 0) {
    warning("No types found in results_df")
    return(NULL)
  }
  base_colors <- brewer.pal(8, "Set2")  # fixed palette of 8 colors
  color_palette <- setNames(rep(base_colors, length.out = n_types), present_types)
  
  results_df$feature_display <- factor(results_df$feature_display, levels = unique(results_df$feature_display))
  
  p <- ggplot(results_df, aes(x = feature_display, y = Odds.ratio, fill = type)) +
    geom_bar(stat = "identity", width = 0.6, position = position_dodge(width = 0.7), color = "black") +
    geom_errorbar(aes(ymin = OR.lower, ymax = OR.upper), position = position_dodge(width = 0.7), width = 0.2) +
    geom_text(aes(label = stars, y = OR.upper + 0.15), position = position_dodge(width = 0.7), vjust = 0, size = 5) +
    geom_hline(yintercept = 1, linetype = "dashed") +
    theme_classic() +
    ylab("Odds Ratio") +
    xlab("") +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
      axis.text.y = element_text(size = 14),
      axis.title.y = element_text(size = 14),
      plot.title = element_text(size = 16, hjust = 0.5),
      legend.position = "bottom"
    ) +
    scale_fill_manual(values = color_palette) +
    labs(title = title, fill = "Group / Type") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1)))
 
 cat("=== Correlation plot ready ===\n\n")
  return(p)
}



#proximty plots
plot_distance_density <- function(distance_data,
                                  analysis_type, 
                                  analysis_name,       # full title string
                                  x_limits = c(-5000, 5000),
                                  x_breaks = c(-5000, 0, 5000),
                                  y_breaks = c(0, 0.0002),
                                  color_palette = NULL) {

  cat("\n=== Distance Density Plot ===\n")
  cat("Analysis:", analysis_name, "\n")
  cat("Type:", analysis_type, "\n")

  # Check if data is empty
  if (is.null(distance_data) || nrow(distance_data) == 0) {
    warning("No distance data provided")
    return(NULL)
  }

  cat("Rows:", nrow(distance_data), "\n")

  # Check for required columns
  required_cols <- c("distance", "rarity")
  missing_cols <- setdiff(required_cols, names(distance_data))
  if (length(missing_cols) > 0) {
    warning("Missing columns in distance_data: ", paste(missing_cols, collapse = ", "))
    return(NULL)
  }

  # Extract unique rarity groups
  groups_present <- unique(distance_data$rarity)

  cat("Groups:", paste(groups_present, collapse = ", "), "\n")

  if (length(groups_present) == 0) {
    warning("No rarity groups found in distance_data")
    return(NULL)
  }

  # If user does NOT supply a palette, generate one dynamically
  if (is.null(color_palette)) {
    # Generate enough distinct colors
    generated_colors <- grDevices::rainbow(length(groups_present))
    names(generated_colors) <- groups_present
    color_palette <- generated_colors
  } else {
    # Only keep colors for groups present
    color_palette <- color_palette[names(color_palette) %in% groups_present]
    
    # Identify missing groups & auto-assign colors
    missing <- setdiff(groups_present, names(color_palette))
    if (length(missing) > 0) {
      extra <- grDevices::rainbow(length(missing))
      names(extra) <- missing
      color_palette <- c(color_palette, extra)
    }
  }
  
  # Decide whether the analysis is TSS or SJ based on the title
  if (grepl("TSS", analysis_type, ignore.case = TRUE)) {
    xlabel <- "Distance from TSS (bp)"
  } else if (grepl("SJ", analysis_type, ignore.case = TRUE)) {
    xlabel <- "Distance from splice junction (bp)"
  } else {
    xlabel <- "Distance (bp)"
  }
  
  # Build plot
  p <- ggplot(distance_data, aes(x = distance, y = after_stat(density), color = rarity)) +
    geom_density(alpha = 0.15, adjust = ifelse(analysis_type == "TSS", 1/10, 1/4)) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    theme_classic() +
    theme(
      panel.border = element_rect(fill = NA),
      legend.position = "bottom",
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 14),
      plot.title = element_text(size = 16, hjust = 0.5)
    ) +
    geom_hline(yintercept = 0, size = 1, color = "white") +
    coord_cartesian(xlim = x_limits) +
    scale_x_continuous(breaks = x_breaks) +
    scale_y_continuous(breaks = y_breaks) +
    scale_color_manual(values = color_palette, name = "Group") +
    labs(title = analysis_name, x = xlabel, y = "Density")

  cat("=== Distance plot ready ===\n\n")
  return(p)
}


#' Create volcano plot for Fisher test results
plot_fisher_volcano <- function(fisher_results, pval_threshold = 0.05, 
                                pval_column = "p_adj_global_all",Row_Name =  "Row_Name", save_path = NULL) {
    cat("\n=== Fisher Volcano Plot ===\n")
    # Check if data is empty
    if (is.null(fisher_results) || nrow(fisher_results) == 0) {
        warning("No Fisher test results provided\n")
        return(NULL)
    }

    cat("Rows:", nrow(fisher_results), "\n")
    cat("p-value column:", pval_column, "\n")
    cat("Threshold:", pval_threshold, "\n")

    # Check required columns
    required_cols <- c("OR", pval_column,Row_Name, "Comparison")
    missing_cols <- setdiff(required_cols, names(fisher_results))
    if (length(missing_cols) > 0) {
        warning("Fisher results missing columns: ", paste(missing_cols, collapse = ", "))
        return(NULL)
    }


    # Prepare data
    plot_data <- fisher_results %>%
        mutate(
            log_OR = log2(OR),
            log_p = -log10(!!sym(pval_column)),
            Significant = !!sym(pval_column) < pval_threshold
        )

    # Check if any p-values are valid
    if (all(is.na(plot_data$log_p))) {
        warning("All p-values are NA, cannot create volcano plot\n")
        return(NULL)
    }

    n_sig <- sum(plot_data$Significant, na.rm = TRUE)
    cat("Significant points:", n_sig, "\n")

 # Create plot
    p <- ggplot(plot_data, aes(x = log_OR, y = log_p)) +
        geom_point(aes(color = Comparison, shape = Significant), 
                   size = 3, alpha = 0.7) +
        geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
        geom_hline(yintercept = -log10(pval_threshold), 
                   color = "red", linetype = "dashed", alpha = 0.7) +
        labs(
            title = "Fisher Test Results - Volcano Plot",
            x = "log2(Odds Ratio)",
            y = paste0("-log10(", pval_column, ")"),
            color = "Group Comparison",
            shape = paste("p <", pval_threshold)
        ) +
        theme_bw() +
        theme(
            legend.position = "bottom",
            plot.title = element_text(hjust = 0.5, face = "bold")
        )
    
    # Add labels only if there are significant points
    if (any(plot_data$Significant, na.rm = TRUE)) {
        p <- p + ggrepel::geom_text_repel(
            data = plot_data %>% filter(Significant),
            aes(label = .data[[Row_Name]]),
            size = 3,
            max.overlaps = 20,
            box.padding = 0.5
        )
    }
    
    # Save if requested
    if (!is.null(save_path)) {
        tryCatch({
            ggsave(save_path, p, width = 12, height = 8)
        }, error = function(e) {
            warning("Failed to save plot: ", e$message)
        })
    }

    cat("=== Volcano plot complete ===\n\n")
    return(p)
}


#' Create Manhattan plot for Fisher test results
plot_fisher_manhattan <- function(fisher_results, pval_threshold = 0.05,
                                  pval_column = "p_adj_global_all",Row_Name="Row_Name" ,save_path = NULL) {

    cat("\n=== Fisher Manhattan Plot ===\n")

    # Check if data is empty
    if (is.null(fisher_results) || nrow(fisher_results) == 0) {
        warning("No Fisher test results provided\n")
        return(NULL)
    }

    cat("Rows:", nrow(fisher_results), "\n")
    cat("p-value column:", pval_column, "\n")
    cat("Threshold:", pval_threshold, "\n")

    # Check required columns
    required_cols <- c(pval_column, Row_Name)
    missing_cols <- setdiff(required_cols, names(fisher_results))
    if (length(missing_cols) > 0) {
        warning("Fisher results missing columns: ", paste(missing_cols, collapse = ", "))
        return(NULL)
    }
    
    
    # Prepare data
    plot_data <- fisher_results %>%
    mutate(
        log_p = -log10(.data[[pval_column]]),
        Significant = .data[[pval_column]] < pval_threshold,
        label_var = factor(.data[[Row_Name]], levels = unique(.data[[Row_Name]]))
    )

    # Check if any p-values are valid
    if (all(is.na(plot_data$log_p))) {
        warning("All p-values are NA, cannot create Manhattan plot")
        return(NULL)
    }

    n_sig <- sum(plot_data$Significant, na.rm = TRUE)
    cat("Significant tests:", n_sig, "\n")

    # Create plot
    p <- ggplot(plot_data, aes(x = label_var, y = log_p)) +
        geom_point(aes(color = Significant), size = 3) +
        geom_hline(yintercept = -log10(pval_threshold), 
                   color = "red", linetype = "dashed", alpha = 0.7) +
        scale_color_manual(values = c("FALSE" = "gray", "TRUE" = "blue")) +
        labs(
            title = "Fisher Test Results - Manhattan Plot",
            x = "Test",
            y = paste0("-log10(", pval_column, ")"),
            color = paste("Significant\n(p <", pval_threshold, ")")
        ) +
        theme_bw() +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1),
            plot.title = element_text(hjust = 0.5, face = "bold")
        )
    
    # Add labels only if there are significant points
    if (any(plot_data$Significant, na.rm = TRUE)) {
        p <- p + ggrepel::geom_text_repel(
            data = plot_data %>% filter(Significant),
            aes(label = label_var),
            size = 3
        )
    }
    
    # Save if requested
    if (!is.null(save_path)) {
        tryCatch({
            ggsave(save_path, p, width = 14, height = 8)
        }, error = function(e) {
            warning("Failed to save plot: ", e$message)
        })
    }

    cat("=== Manhattan plot complete ===\n\n")

    return(p)
}

#' Create dot plot for Fisher test results
plot_fisher_dot <- function(fisher_results, pval_threshold = 0.05,
                            pval_column = "p_adj_global_all",Row_Name ="Row_Name", save_path = NULL) {
    
    cat("\n=== Fisher Dot Plot ===\n")

    # Check if data is empty
    if (is.null(fisher_results) || nrow(fisher_results) == 0) {
        warning("No Fisher test results provided")
        return(NULL)
    }

    cat("Rows:", nrow(fisher_results), "\n")
    cat("p-value column:", pval_column, "\n")
    cat("Threshold:", pval_threshold, "\n")

    # Check required columns
    required_cols <- c(pval_column, Row_Name, "Comparison")
    missing_cols <- setdiff(required_cols, names(fisher_results))
    if (length(missing_cols) > 0) {
        warning("Fisher results missing columns: ", paste(missing_cols, collapse = ", "))
        return(NULL)
    }
    
    # Prepare data
    plot_data <- fisher_results %>%
    mutate(
        log_p = -log10(.data[[pval_column]]),
        Significant = .data[[pval_column]] < pval_threshold,
        label_var = factor(.data[[Row_Name]], levels = unique(.data[[Row_Name]]))
    )

    # Check if any p-values are valid
    if (all(is.na(plot_data$log_p))) {
        warning("All p-values are NA, cannot create dot plot")
        return(NULL)
    }

    n_sig <- sum(plot_data$Significant, na.rm = TRUE)
    cat("Significant tests:", n_sig, "\n")

    # Create plot
    p <- ggplot(plot_data, aes(x = log_p, y = label_var)) +
        geom_point(aes(color = Comparison, shape = Significant), size = 3) +
        geom_vline(xintercept = -log10(pval_threshold), 
                   linetype = "dashed", color = "red") +
        scale_shape_manual(values = c(16, 17)) +
        labs(
            title = "Fisher Test Results - Dot Plot",
            x = paste0("-log10(", pval_column, ")"),
            y = "Test",
            color = "Comparison",
            shape = paste("Significant\n(p <", pval_threshold, ")")
        ) +
        theme_bw() +
        theme(
            axis.text.y = element_text(size = 8),
            plot.title = element_text(hjust = 0.5, face = "bold")
        )
    
    # Add labels only if there are significant points
    if (any(plot_data$Significant, na.rm = TRUE)) {
        p <- p + ggrepel::geom_text_repel(
            data = plot_data %>% filter(Significant),
            aes(label = label_var),
            size = 3
        )
    }
    
    # Save if requested
    if (!is.null(save_path)) {
        tryCatch({
            ggsave(save_path, p, width = 10, height = max(6, nrow(plot_data) * 0.3))
        }, error = function(e) {
            warning("Failed to save plot: ", e$message)
        })
    }

    cat("=== Dot plot complete ===\n\n")
    return(p)
}



#! kruskal wilcoxon plots 
plot_wilcoxon_kruskal_flexible <- function(
    test_results, 
    count_data, 
    test_type = "wilcoxon",
    Row_Name = "metric",
    group_col = "group",
    plot_option = "significant",
    pval_threshold = 0.05,
    save_dir = NULL,
    ncol = 4,
    color_palette = NULL,
    jitter_alpha = 0.5,
    jitter_size = 0.8
) {
    
    # Validate plot_option
    valid_options <- c("significant", "all", "global_only")
    if (!plot_option %in% valid_options) {
        stop("plot_option must be one of: ", paste(valid_options, collapse = ", "))
    }

    # Check if data is empty
    if (is.null(test_results) || nrow(test_results) == 0) {
        warning("No test results provided")
        return(NULL)
    }
    
    if (is.null(count_data) || nrow(count_data) == 0) {
        warning("No count data provided")
        return(NULL)
    }

    # Determine p-value column
    pval_col <- if ("global_primary_p_adj" %in% names(test_results)) {
        "global_primary_p_adj"
    } else if ("primary_p_adj" %in% names(test_results)) {
        "primary_p_adj"
    } else if ("p_value" %in% names(test_results)) {
        "p_value"
    } else if ("dunn_p_adj_global"%in% names(test_results )){
        "dunn_p_adj_global"
    } else if ("global_kruskal_p_adj"%in% names(test_results)){
        "global_kruskal_p_adj"
    } else {
        stop("No p-value column found in test results")
        return(NULL)
    }
    
    # Determine metric column
    metric_col <- if ("metric" %in% names(test_results)) {
        "metric"
    } else if (Row_Name %in% names(test_results)) {
        Row_Name
    } else {
        stop("No metric/feature column found in test results")
        return(NULL)
    }
    
    # Get metrics based on option
    if (plot_option == "significant") {
        sig_metrics <- test_results %>%
            filter(!!sym(pval_col) < pval_threshold) %>%
            pull(!!sym(metric_col)) %>%
            unique()
        
        if (length(sig_metrics) == 0) {
            message("No significant results at p < ", pval_threshold, ". Returning NULL.")
            return(NULL)
        }
        selected_metrics <- sig_metrics
        plot_title <- paste("Significant", test_type, "Results (p <", pval_threshold, ")")
        
    } else if (plot_option == "all") {
        selected_metrics <- unique(test_results[[metric_col]])
        if (length(selected_metrics) == 0) {
            warning("No metrics found in test results")
            return(NULL)
        }
        plot_title <- paste("All", test_type, "Test Results")
        
    } else if (plot_option == "global_only") {
        global_metrics <- grep("^Number_of_Variants$", unique(test_results[[metric_col]]), value = TRUE)
        
        if (length(global_metrics) == 0) {
            message("No global 'Number_of_Variants' metric found")
            return(NULL)
        }
        selected_metrics <- global_metrics
        plot_title <- paste(test_type, "Test: Global Variants")
    }
    
    # Prepare annotation data
    annotation_data <- test_results %>%
        filter(!!sym(metric_col) %in% selected_metrics) %>%
        mutate(
            label = paste0(test_type, "\n",
                          ifelse(!!sym(pval_col) < 0.001, "p < 0.001", 
                                 paste0("p = ", format.pval(!!sym(pval_col), digits = 2)))),
            metric_clean = gsub("^Number_of_|_Variants$", "", !!sym(metric_col)),
            Significant = !!sym(pval_col) < pval_threshold
        )

    # Check if annotation data is empty
    if (nrow(annotation_data) == 0) {
        warning("No annotation data after filtering")
        return(NULL)
    }

    # Reshape count data
    plot_data <- count_data %>%
        select(sample_id, !!sym(group_col), any_of(selected_metrics)) %>%
        pivot_longer(
            cols = -c(sample_id, !!sym(group_col)),
            names_to = "metric",
            values_to = "value"
        ) %>%
        mutate(
            metric_clean = gsub("^Number_of_|_Variants$", "", metric)
        ) %>%
        filter(metric %in% selected_metrics)


    if (nrow(plot_data) == 0) {
        warning("No data to plot after filtering")
        return(NULL)
    }
    # Get unique groups and count
    unique_groups <- unique(plot_data[[group_col]])
    n_groups <- length(unique_groups)

    if (n_groups == 0) {
        warning("No groups found in plot data")
        return(NULL)
    }

    # Ensure group column is factor with proper ordering
    plot_data[[group_col]] <- factor(plot_data[[group_col]], 
                                     levels = unique(plot_data[[group_col]]))
    
    # Order metrics by p-value (most significant first)
    metric_order <- annotation_data %>%
        arrange(!!sym(pval_col)) %>%
        pull(metric_clean) %>%
        unique()
    
    plot_data$metric_clean <- factor(plot_data$metric_clean, levels = metric_order)
    annotation_data$metric_clean <- factor(annotation_data$metric_clean, levels = metric_order)
    
    # Create color palette based on number of groups
    if (is.null(color_palette)) {
        if (n_groups <= 12) {
            color_palette_name <- "Set3"
            colors <- RColorBrewer::brewer.pal(max(3, n_groups), color_palette_name)[1:n_groups]
        } else if (n_groups <= 20) {
            color_palette_name <- "Set1"
            colors <- colorRampPalette(RColorBrewer::brewer.pal(9, color_palette_name))(n_groups)
        } else {
            colors <- scales::hue_pal()(n_groups)
            color_palette_name <- "hue"
        }
    } else if (is.character(color_palette) && length(color_palette) == 1) {
        color_palette_name <- color_palette
        if (color_palette_name %in% row.names(RColorBrewer::brewer.pal.info)) {
            max_colors <- RColorBrewer::brewer.pal.info[color_palette_name, "maxcolors"]
            if (n_groups <= max_colors) {
                colors <- RColorBrewer::brewer.pal(max(3, n_groups), color_palette_name)[1:n_groups]
            } else {
                colors <- colorRampPalette(RColorBrewer::brewer.pal(
                    min(n_groups, max_colors), color_palette_name))(n_groups)
            }
        } else {
            colors <- scales::hue_pal()(n_groups)
            color_palette_name <- "hue"
        }
    } else {
        colors <- color_palette
        color_palette_name <- "custom"
    }
    
    if (length(colors) < n_groups) {
        warning("Not enough colors provided. Using default palette instead.")
        colors <- scales::hue_pal()(n_groups)
        color_palette_name <- "hue"
    }
    
    # Create summary plot - KEY FIX: map color to group for jitter points
    summary_plot <- ggplot(plot_data, aes(x = !!sym(group_col), y = value, 
                                          fill = !!sym(group_col),
                                          color = !!sym(group_col))) +
        geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.7) +
        geom_jitter(width = 0.2, height = 0, alpha = jitter_alpha, size = jitter_size) +
        facet_wrap(~ metric_clean, scales = "free_y", ncol = ncol) +
        scale_fill_manual(values = colors, name = group_col) +
        scale_color_manual(values = colors, name = group_col) +
        labs(
            title = plot_title,
            subtitle = if (plot_option %in% c("significant", "all")) 
                paste("Red labels: p <", pval_threshold, "(FDR)") else NULL,
            x = group_col,
            y = "Count"
        ) +
        theme_bw() +
        theme(
            strip.background = element_blank(),
            strip.text = element_text(face = "bold"),
            legend.position = "bottom",
            plot.subtitle = element_text(color = "red"),
            plot.title = element_text(hjust = 0.5, face = "bold"),
            axis.text.x = element_text(angle = 45, hjust = 1)
        ) +
        guides(fill = guide_legend(nrow = ceiling(n_groups/4), override.aes = list(size = 3)))
    
    # Add annotation only if there is data
    if (nrow(annotation_data) > 0) {
        summary_plot <- summary_plot + 
            geom_text(
                data = annotation_data,
                aes(x = Inf, y = Inf, label = label),
                inherit.aes = FALSE,
                vjust = 1.5, hjust = 1.1, size = 3,
                color = ifelse(annotation_data$Significant, "red", "black")
            )
    }

    # For global_only, adjust plot
    if (plot_option == "global_only") {
        summary_plot <- ggplot(plot_data, aes(x = !!sym(group_col), y = value, 
                                              fill = !!sym(group_col),
                                              color = !!sym(group_col))) +
            geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.7) +
            geom_jitter(width = 0.2, height = 0, alpha = jitter_alpha, size = 2) +
            scale_fill_manual(values = colors, name = group_col) +
            scale_color_manual(values = colors, name = group_col) +
            labs(
                title = plot_title,
                x = group_col,
                y = "Number of Variants"
            ) +
            theme_bw() +
            theme(
                plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
                legend.position = "bottom",
                axis.text = element_text(size = 12),
                axis.title = element_text(size = 13),
                axis.text.x = element_text(angle = 45, hjust = 1)
            ) +
            guides(fill = guide_legend(nrow = ceiling(n_groups/4), override.aes = list(size = 4)))
        
        # Add annotation if available
        if (nrow(annotation_data) > 0) {
            summary_plot <- summary_plot + 
                geom_text(
                    data = annotation_data,
                    aes(x = Inf, y = Inf, label = label),
                    inherit.aes = FALSE,
                    vjust = 1.5, hjust = 1.1, size = 5,
                    color = ifelse(annotation_data$Significant, "red", "black")
                )
        }
    }
    
    # Create individual plots for each metric
    individual_plots <- list()
    for (metric in selected_metrics) {
        metric_data <- plot_data %>% filter(metric == !!metric)
        metric_clean <- unique(gsub("^Number_of_|_Variants$", "", metric))
        metric_annotation <- annotation_data %>% filter(!!sym(metric_col) == !!metric)
        
        if (nrow(metric_data) == 0) {
            next
        }
        
        p <- ggplot(metric_data, aes(x = !!sym(group_col), y = value, 
                                     fill = !!sym(group_col),
                                     color = !!sym(group_col))) +
            geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.7) +
            geom_jitter(width = 0.2, height = 0, alpha = jitter_alpha, size = 1.5) +
            scale_fill_manual(values = colors, name = group_col) +
            scale_color_manual(values = colors, name = group_col) +
            labs(
                title = paste(test_type, "Test:", metric_clean),
                x = group_col,
                y = "Count"
            ) +
            theme_bw() +
            theme(
                plot.title = element_text(hjust = 0.5, face = "bold"),
                legend.position = "bottom",
                axis.text.x = element_text(angle = 45, hjust = 1)
            ) +
            guides(fill = guide_legend(nrow = ceiling(n_groups/4), override.aes = list(size = 3)))
        
        # Add p-value annotation if available
        if (nrow(metric_annotation) > 0) {
            pval_text <- ifelse(metric_annotation$Significant[1], 
                               paste0("FDR p = ", format.pval(metric_annotation[[pval_col]][1], digits = 2)),
                               paste0("p = ", format.pval(metric_annotation[[pval_col]][1], digits = 2)))
            
            p <- p + annotate(
                "text", x = -Inf, y = Inf,
                label = pval_text,
                hjust = -0.1, vjust = 1.5, size = 4,
                color = ifelse(metric_annotation$Significant[1], "red", "black")
            )
        }
        
        individual_plots[[metric_clean]] <- p
    }
    
    # Save if requested
   if (!is.null(save_dir) && length(selected_metrics) > 0) {
        tryCatch({
            if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
            
            if (plot_option == "significant") {
                filename <- paste0(test_type, "_SIGNIFICANT_boxplot.png")
            } else if (plot_option == "all") {
                filename <- paste0(test_type, "_ALL_metrics_boxplot.png")
            } else if (plot_option == "global_only") {
                filename <- paste0(test_type, "_GLOBAL_VARIANTS_boxplot.png")
            }
            
            if (plot_option == "global_only") {
                ggsave(file.path(save_dir, filename),
                       summary_plot, width = max(8, n_groups), height = 6)
            } else {
                ggsave(file.path(save_dir, filename),
                       summary_plot, 
                       width = min(20, ncol * 5),
                       height = ceiling(length(selected_metrics) / ncol) * 3.5 + 2)
            }
            
            if (plot_option != "global_only" && length(individual_plots) > 0) {
                save_multiple_plots_grid(
                    plots = individual_plots,
                    prefix = paste0(test_type, "_", plot_option, "_individual"),
                    output_dir = save_dir,
                    plots_per_page = 4
                )
            }
        }, error = function(e) {
            warning("Failed to save plots: ", e$message)
        })
    }
    
    # Print summary
    n_sig <- if (plot_option == "significant") {
        length(selected_metrics)
    } else if (plot_option == "all") {
        sum(annotation_data$Significant, na.rm = TRUE)
    } else {
        NA
    }
    
    cat("=== Plot Summary ===\n")
    cat("Plot option:", plot_option, "\n")
    cat("Test type:", test_type, "\n")
    cat("Metrics plotted:", length(selected_metrics), "\n")
    cat("Number of groups:", n_groups, "\n")
    cat("Color palette used:", color_palette_name, "\n")
    if (plot_option != "global_only") {
        cat("Significant metrics (p <", pval_threshold, "):", n_sig, "\n")
    }
    cat("\n")
    
    return(list(
        summary_plot = summary_plot,
        individual_plots = individual_plots,
        selected_metrics = selected_metrics,
        colors_used = colors,
        n_total = length(selected_metrics),
        n_significant = n_sig,
        n_groups = n_groups,
        plot_option = plot_option
    ))
}

# Utility function for saving multiple plots
save_multiple_plots_grid <- function(plots, prefix, output_dir, plots_per_page = 4) {
    if (length(plots) == 0) return()
    
    plot_chunks <- split(plots, ceiling(seq_along(plots) / plots_per_page))
    
    for (i in seq_along(plot_chunks)) {
        page_plots <- plot_chunks[[i]]
        
        grid_plot <- gridExtra::grid.arrange(
            grobs = page_plots,
            ncol = 2,
            top = paste(prefix, "- Page", i)
        )
        
        ggsave(
            file.path(output_dir, paste0(prefix, "_page_", i, ".png")),
            grid_plot,
            width = 14,
            height = 10
        )
    }
}


#regression plot 


plot_glm <- function(results_df, title = "Odds Ratio Plot", star_offset = 1.05) {
  
  cat("\n=== GLM Odds Ratio Plot ===\n")
  cat("Title:", title, "\n")

  # Check if data is empty
  if (is.null(results_df) || nrow(results_df) == 0) {
    warning("No GLM results provided\n")
    return(NULL)
  }

  cat("Rows:", nrow(results_df), "\n")

  # Check for required columns
  required_cols <- c("predictor", "odds_ratio", "lower_OR", "upper_OR", "adj_p_val")
  missing_cols <- setdiff(required_cols, names(results_df))
  if (length(missing_cols) > 0) {
    warning("Missing columns in GLM results: ", paste(missing_cols, collapse = ", "))
    return(NULL)
  }
  
  plot_data <- results_df %>% 
    filter(!grepl("^\\(Intercept\\)$|^intercept$", predictor, ignore.case = TRUE)) %>%
    filter(!is.na(odds_ratio)) %>%
    mutate(signif = case_when(
      adj_p_val < 0.001 ~ "***",
      adj_p_val < 0.01  ~ "**",
      adj_p_val < 0.05  ~ "*",
      TRUE ~ ""
    )) %>%
    mutate(star_y = upper_OR * star_offset)

  cat("Predictors plotted:", nrow(plot_data), "\n")
  cat("Significant predictors:", sum(plot_data$adj_p_val < 0.05), "\n")

  if (nrow(plot_data) == 0) {
    warning("No non-intercept terms found in results")
    return(NULL)
  }


  p <- ggplot(plot_data, aes(x = predictor, y = odds_ratio, fill = adj_p_val < 0.05)) +
    geom_col(stat = "identity", color = "black", width = 0.5, show.legend = FALSE) +
    geom_errorbar(aes(ymin = lower_OR, ymax = upper_OR),
                  position = position_dodge(width = 0.5), width = 0.2, size = 0.4) +
    scale_fill_manual(values = c("#4393C3", "#D6604D")) +
    geom_hline(yintercept = 1, linetype = "dashed") +
    theme_classic() +
    labs(title = title, y = "Odds Ratio", x = "")
  
  # Add significance stars only if there are any
  if (any(plot_data$signif != "")) {
    p <- p + geom_text(aes(y = star_y, label = signif),
              position = position_dodge(width = 0.5),
              vjust = 0, size = 5)
  }
  
  # Save plot automatically
  tryCatch({
    safe_title <- str_replace_all(title, "\\s+", "_")
    safe_title <- str_replace_all(safe_title, "[^A-Za-z0-9_]", "")
    ggsave(filename = file.path(output_dir, sprintf("%s.png", safe_title)), 
           plot = p, width = 8, height = 6, dpi = 300)
    message("Plot saved to: ", file.path(output_dir, sprintf("%s.png", safe_title)))
  }, error = function(e) {
    warning("Failed to save plot: ", e$message)
  })

  cat("=== GLM plot complete ===\n\n")
  return(p)
}


plot_multinom_or <- function(df, x_col = "predictor", y_col = "odds_ratio",
                             lower_col = "lower_OR", upper_col = "upper_OR",
                             comparison_col = "comparison", pval_col = "adj_p_val",
                             star_offset = 1.05,
                             title = "Odds Ratio Plot") {
  cat("\n=== Multinomial OR Plot ===\n")
  cat("Title:", title, "\n")

  # Check if data is empty
  if (is.null(df) || nrow(df) == 0) {
    warning("No multinomial regression results provided")
    return(NULL)
  }

  cat("Rows:", nrow(df), "\n")

  # Check for required columns
  required_cols <- c(x_col, y_col, lower_col, upper_col, comparison_col, pval_col)
  missing_cols <- setdiff(required_cols, names(df))
  if (length(missing_cols) > 0) {
    warning("Missing columns: ", paste(missing_cols, collapse = ", "))
    return(NULL)
  }



  df <- df %>%
    mutate(signif = case_when(
      !!sym(pval_col) < 0.001 ~ "***",
      !!sym(pval_col) < 0.01  ~ "**",
      !!sym(pval_col) < 0.05  ~ "*",
      TRUE ~ ""
    )) %>%
    filter(!grepl("^\\(Intercept\\)$|^intercept$", !!sym(x_col), ignore.case = TRUE)) %>%
    filter(!is.na(!!sym(y_col))) %>%
    mutate(star_y = !!sym(upper_col) * star_offset)

  cat("Terms plotted:", nrow(df), "\n")
  cat("Significant terms:", sum(df[[pval_col]] < 0.05), "\n")

  if (nrow(df) == 0) {
    warning("No non-intercept terms with valid odds ratios found")
    return(NULL)
  }

  p <- ggplot(df, aes_string(x = x_col, y = y_col, fill = comparison_col)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    geom_errorbar(aes_string(ymin = lower_col, ymax = upper_col),
                  position = position_dodge(width = 0.8), width = 0.2) +
    scale_y_log10() +
    theme_minimal() +
    labs(title = title, y = "Odds Ratio (log scale)", x = "Predictor / Group", fill = "Comparison")
  
  # Add significance stars only if there are any
  if (any(df$signif != "")) {
    p <- p + geom_text(aes(y = star_y, label = signif),
              position = position_dodge(width = 0.8),
              vjust = 0, size = 5)
  }
  
  # Save plot automatically
  tryCatch({
    safe_title <- str_replace_all(title, "\\s+", "_")
    safe_title <- str_replace_all(safe_title, "[^A-Za-z0-9_]", "")
    ggsave(filename = file.path(output_dir, sprintf("%s.png", safe_title)), 
           plot = p, width = 8, height = 6, dpi = 300)
    message("Plot saved to: ", file.path(output_dir, sprintf("%s.png", safe_title)))
  }, error = function(e) {
    warning("Failed to save plot: ", e$message)
  })

  cat("=== Multinomial plot complete ===\n\n")
  return(p)
}

#14ORA
build_enrichmentmap_network <- function(folder, name,gmt_file){
  folder <- normalizePath(folder, mustWork = TRUE)
  
  if(name == "combined") {
    # Find ALL GEM files for all groups
    gem_files <- list.files(folder, pattern = "gProfiler_.*_gem\\.txt$", full.names = TRUE)
    gem_files <- normalizePath(gem_files, mustWork = TRUE)
    if(length(gem_files) == 0) return(NULL)
    
    message("Building COMBINED network with ", length(gem_files), " groups")
    
    # Build the request body
    body <- list(
      analysisType = "generic",
      gmtFile = gmt_file,
      pvalue = 0.05,
      qvalue = 0.05,
      similaritycutoff = 0.375,
      coefficients = "COMBINED",
      combinedConstant = 0.5,
      edgeStrategy = "DISTINCT",
      runAutoAnnotate = FALSE
    )
    
    # Add each GEM file as a separate dataset
    for(i in seq_along(gem_files)) {
      body[[paste0("enrichmentsDataset", i)]] <- gem_files[i]
    }
    
    cyrestPOST("commands/enrichmentmap/build", body = body)
    
  } else {
    # Single group - use your existing code
    gem_file <- file.path(folder, paste0("gProfiler_gem_", name, ".txt"))
    gem_file <- normalizePath(gem_file, mustWork = TRUE)
    if(!file.exists(gem_file)) return(NULL)
    
    cyrestPOST(
      "commands/enrichmentmap/build",
      body = list(
        analysisType = "generic",
        gmtFile = gmt_file,
        enrichmentsDataset1 = gem_file,
        pvalue = 0.05,
        qvalue = 0.05,
        similaritycutoff = 0.375,
        coefficients = "COMBINED",
        combinedConstant = 0.5,
        edgeStrategy = "DISTINCT",
        runAutoAnnotate = FALSE
      )
    )
  }
  

commandsPOST(
    paste(
        "autoannotate annotate-clusterBoosted",
        "network=current",
        "edgeWeightColumn=EnrichmentMap::similarity_coefficient",
        "labelColumn=EnrichmentMap::GS_DESCR",
        "clusterAlgorithm=MCL",
        "maxWords=3",
        "minWordOccurrence=1",
        "adjacentWordBonus=8",
        "createSingletonClusters=true"
    )
)
  
cyrestPOST(
    "commands/layout/autoannotate-cose-cluster",
    body = list(
        network = "current",
        idealEdgeLength = "150",
        springStrength = "80",
        repulsionStrength = "100",
        gravityStrength = "70",
        compoundGravityStrength = "40",
        layoutQuality = "Proof"
    )
)


 exportImage(file.path(folder, paste0("network_", name, "_original.png")), type = "PNG", resolution = 300)
   saveSession(file.path(folder, paste0("network_", name, ".cys")))
    cyrestPOST(
      "commands/autoannotate/summary",
      body = list(
        network = "current",
        includeUnclustered = "true",
        createNewNetwork = "true"
      )
    )

 layoutNetwork("force-directed")

  exportImage(file.path(folder, paste0("network_", name, "_summary.png")), type = "PNG", resolution = 300)
  

  message("Finished Cytoscape network: ", name)
  
}


# Save gprofiler2 results with improved plots
save_gprofiler2_results <- function(gp_results,  group, output_dir) {
  # Flatten list columns
  gp_df <- gp_results
  gp_df[] <- lapply(gp_df, function(col) {
    if (is.list(col)) sapply(col, function(x) paste(x, collapse = ",")) else col
  })
  
  write.table(gp_df,
              file = file.path(output_dir, paste0(  group, "_gprofiler2.tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)
  
  # Enhanced barplot with better spacing
  if (nrow(gp_df) > 0) {
    top_terms <- gp_df %>%
      arrange(p_value) %>%
      slice_head(n = 20)
    
    pdf(file.path(output_dir, paste0( group, "_gprofiler2_barplot.pdf")), 
       width = 10, height = 8)
    
    p <- ggplot(top_terms, aes(x = reorder(term_name, -log10(p_value)), y = -log10(p_value))) +
      geom_bar(stat = "identity", fill = "steelblue") +
      coord_flip() +
      labs(title = paste( group, "_","gprofiler2"),
           x = "Term", y = "-log10(p-value)") +
      theme_minimal() +
      theme(
        axis.text.y = element_text(size = 10, margin = margin(r = 15)), # Increased right margin
        plot.margin = unit(c(1, 1, 1, 1), "cm"), # Overall plot margins
        panel.spacing = unit(1, "lines") # Space between facets if used
      )
    
    print(p)
    dev.off()
  }
}


#General plots 
compute_counts <- function(df, 
                           x_var, 
                           manifest_ufficial, 
                           cohort_col = "outlier_label", 
                           status_col = "Status",
                           top_n = 5) { 

  # Check if data is empty
  if (is.null(df) || nrow(df) == 0) {
    warning("No data provided for compute_counts")
    return(list(raw = NULL, normalized = NULL))
  }
  
  if (is.null(manifest_ufficial) || nrow(manifest_ufficial) == 0) {
    warning("No manifest data provided for compute_counts")
    return(list(raw = NULL, normalized = NULL))
  }

  # Add motif_length column only if x_var is "motif"
  if (x_var == "motif_length") {
    df <- df %>% mutate(motif_length = nchar(as.character(motif)))                       
  }

  # 1) raw counts per cohort
  raw_counts <- df %>%
    filter(!is.na(.data[[cohort_col]])) %>%
    count(cohort = .data[[cohort_col]], value = .data[[x_var]], name = "n_loci")

  if (nrow(raw_counts) == 0) {
    warning("No counts found after filtering")
    return(list(raw = NULL, normalized = NULL))
  }
  
  # Subset top N per cohort if requested
  if (!is.null(top_n) & x_var == "motif") {
    raw_counts <- raw_counts %>%
      group_by(cohort) %>%
      slice_max(order_by = n_loci, n = top_n, with_ties = FALSE) %>%
      ungroup()
  }

  # 2) number of samples per cohort
  sample_counts <- manifest_ufficial %>%
    count(cohort = .data[[status_col]], name = "n_samples")

  if (nrow(sample_counts) == 0) {
    warning("No sample counts found")
    return(list(raw = raw_counts, normalized = NULL))
  }

  # 3) normalized counts
  normalized_counts <- raw_counts %>%
    left_join(sample_counts, by = "cohort") %>%
    mutate(normalized_count = n_loci / n_samples)%>%
    filter(cohort != "mixed")

  list(
    raw = raw_counts,
    normalized = normalized_counts
  )
}


compute_nucleotide_counts <- function(df, manifest_ufficial, cohort_col = "outlier_label", status_col = "Status") {
  # Check if data is empty
  if (is.null(df) || nrow(df) == 0) {
    warning("No data provided for compute_nucleotide_counts")
    return(list(raw = NULL, normalized = NULL))
  }

  # Split motifs into nucleotides and count per cohort
  tryCatch({
    nucleotide_counts <- df %>%
      filter(!is.na(.data[[cohort_col]])) %>%
      mutate(nucleotides = strsplit(as.character(motif), "")) %>%
      unnest(nucleotides) %>%
      filter(nucleotides != "N") %>% # remove ambiguous bases
      count(cohort = .data[[cohort_col]], nucleotide = nucleotides, name = "n_loci")
    
    # Check if nucleotide_counts is empty
    if (nrow(nucleotide_counts) == 0) {
      warning("No nucleotide counts found after processing")
      return(list(raw = NULL, normalized = NULL))
    }

  # Get sample size per cohort
    sample_counts <- manifest_ufficial %>%
      count(cohort = .data[[status_col]], name = "n_samples")
    
    # Check if sample_counts is empty
    if (nrow(sample_counts) == 0) {
      warning("No sample counts found")
      return(list(raw = nucleotide_counts, normalized = NULL))
    }
  
  # Normalize by cohort size
  normalized_counts <- nucleotide_counts %>%
    left_join(sample_counts, by = "cohort") %>%
    mutate(normalized_count = n_loci / n_samples) %>%
    filter(cohort != "mixed")

 list(raw = nucleotide_counts, normalized = normalized_counts)
  }, error = function(e) {
    warning("Error processing nucleotide counts: ", e$message)
    return(list(raw = NULL, normalized = NULL))
  })
}

plot_counts <- function(df, y_var = "n_loci", title = NULL, x_label = NULL) {
  ggplot(df, aes(x = value, y = .data[[y_var]], fill = cohort)) +
    geom_col(position = "dodge") +
    labs(title = title, x = x_label, y = y_var) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          text = element_text(family = "Times", color = "black"))
}

plot_nucleotides <- function(df, y_var = "n_loci", title = NULL) {
  ggplot(df, aes(x = nucleotide, y = .data[[y_var]], fill = cohort)) +
    geom_col(position = "dodge") +
    labs(title = title, x = "Nucleotide", y = y_var) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
          text = element_text(family = "Times", color = "black"))
}
