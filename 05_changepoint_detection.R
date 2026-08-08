# ============================================================
# 05_CHANGEPOINT_DETECTION.R
# Multi-Method Change Point Analysis
# Combining PELT, Binary Segmentation, BIC, CUSUM, Bayesian,
# Moving-Average, Variance, and Visual heuristic methods
#
# Required local input (not included in this repository; see Data
# availability in README.md):
#   - Grouped_Dataset_Kainji_Jebba.xlsx, with a sheet matching the
#     pattern "WQV" or "Meteo", containing columns Site, Habitat,
#     Time_index, and the parameters listed below
# ============================================================

# Install missing packages if needed
required_packages <- c("strucchange", "changepoint", "bcp", "zoo",
                        "ggplot2", "readxl", "dplyr", "purrr", "tidyr", "gridExtra")
missing_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]

if(length(missing_packages) > 0) {
  cat("Installing missing packages:", paste(missing_packages, collapse = ", "), "\n")
  install.packages(missing_packages, dependencies = TRUE)
}

# Load required libraries
library(strucchange)
library(ggplot2)
library(readxl)
library(dplyr)
library(purrr)
library(tidyr)
library(gridExtra)
library(zoo)

# Try to load optional packages with fallbacks
changepoint_available <- requireNamespace("changepoint", quietly = TRUE)
bcp_available <- requireNamespace("bcp", quietly = TRUE)

if(changepoint_available) library(changepoint)
if(bcp_available) library(bcp)

# 1. Load and prepare data ----------------------------------------------------
file_path <- "Grouped_Dataset_Kainji_Jebba.xlsx"

# Check if file exists
if(!file.exists(file_path)) {
  stop("Excel file not found: ", file_path, "\nPlease ensure the file is in your working directory.")
}

# Get sheet name dynamically
target_sheet <- excel_sheets(file_path) |>
  grep(pattern = "WQV|Meteo", ignore.case = TRUE, value = TRUE) |>
  first()

if(length(target_sheet) == 0) {
  cat("Available sheets:\n")
  print(excel_sheets(file_path))
  stop("No sheet found matching 'WQV' or 'Meteo' pattern")
}

cat("Using sheet:", target_sheet, "\n")

# Read data
data <- read_excel(file_path, sheet = target_sheet) |>
  mutate(
    Date = case_when(
      Time_index == 1 ~ "Apr-2024",
      Time_index == 2 ~ "May-2024",
      Time_index == 3 ~ "Jun-2024",
      Time_index == 4 ~ "Jul-2024",
      Time_index == 5 ~ "Aug-2024",
      Time_index == 6 ~ "Sep-2024",
      Time_index == 7 ~ "Oct-2024",
      Time_index == 8 ~ "Nov-2024",
      Time_index == 9 ~ "Dec-2024",
      Time_index == 10 ~ "Jan-2025",
      Time_index == 11 ~ "Feb-2025",
      Time_index == 12 ~ "Mar-2025"
    ),
    Date = factor(Date, levels = c(
      "Apr-2024", "May-2024", "Jun-2024", "Jul-2024", "Aug-2024",
      "Sep-2024", "Oct-2024", "Nov-2024", "Dec-2024", "Jan-2025",
      "Feb-2025", "Mar-2025"
    )),
    Time_cont = as.numeric(Time_index)
  )

# 2. Multi-method change point detection --------------------------------------
detect_changes_multi <- function(y, param_name, site, habitat, significance = 0.05) {

  # Remove NA values
  y_clean <- y[!is.na(y)]
  n <- length(y_clean)

  if (n < 4) return(NULL)

  results <- list(
    param = param_name,
    site = site,
    habitat = habitat,
    n = n,
    methods = list()
  )

  # Method 1: BIC-based (strucchange)
  tryCatch({
    bp_model <- breakpoints(y_clean ~ 1, h = 0.15)
    bic_values <- BIC(bp_model)
    optimal_breaks <- which.min(bic_values) - 1

    if (optimal_breaks > 0) {
      bp_locations <- bp_model$breakpoints[1:optimal_breaks]
      results$methods$BIC <- list(
        detected = TRUE,
        breakpoints = bp_locations,
        n_breaks = optimal_breaks,
        bic = min(bic_values)
      )
    } else {
      results$methods$BIC <- list(detected = FALSE, bic = min(bic_values))
    }
  }, error = function(e) {
    results$methods$BIC <- list(detected = FALSE, error = e$message)
  })

  # Method 2: CUSUM Test (more sensitive to gradual changes)
  tryCatch({
    cusum_test <- efp(y_clean ~ 1, type = "OLS-CUSUM")
    cusum_pval <- sctest(cusum_test)$p.value

    results$methods$CUSUM <- list(
      detected = cusum_pval < significance,
      p_value = cusum_pval,
      test_statistic = cusum_test$process[length(cusum_test$process)]
    )
  }, error = function(e) {
    results$methods$CUSUM <- list(detected = FALSE, error = e$message)
  })

  # Method 3: PELT (Pruned Exact Linear Time) - changepoint package
  if(changepoint_available) {
    tryCatch({
      cpt_pelt <- changepoint::cpt.mean(y_clean, method = "PELT", minseglen = 2)
      pelt_cpts <- changepoint::cpts(cpt_pelt)

      results$methods$PELT <- list(
        detected = length(pelt_cpts) > 0,
        breakpoints = pelt_cpts,
        n_breaks = length(pelt_cpts),
        penalty = "BIC"
      )
    }, error = function(e) {
      results$methods$PELT <- list(detected = FALSE, error = e$message)
    })
  } else {
    results$methods$PELT <- list(detected = FALSE, error = "changepoint package not available")
  }

  # Method 4: Binary Segmentation
  if(changepoint_available) {
    tryCatch({
      cpt_binseg <- changepoint::cpt.mean(y_clean, method = "BinSeg", Q = 3, minseglen = 2)
      binseg_cpts <- changepoint::cpts(cpt_binseg)

      results$methods$BinSeg <- list(
        detected = length(binseg_cpts) > 0,
        breakpoints = binseg_cpts,
        n_breaks = length(binseg_cpts)
      )
    }, error = function(e) {
      results$methods$BinSeg <- list(detected = FALSE, error = e$message)
    })
  } else {
    results$methods$BinSeg <- list(detected = FALSE, error = "changepoint package not available")
  }

  # Method 5: Bayesian Change Point (bcp package)
  if(bcp_available) {
    tryCatch({
      bcp_result <- bcp::bcp(y_clean, mcmc = 1000, burnin = 100)
      prob_threshold <- 0.5
      change_probs <- bcp_result$prob.mean[-1]  # Remove first element
      bcp_cpts <- which(change_probs > prob_threshold)

      results$methods$Bayesian <- list(
        detected = length(bcp_cpts) > 0,
        breakpoints = bcp_cpts,
        n_breaks = length(bcp_cpts),
        max_prob = max(change_probs),
        mean_prob = mean(change_probs)
      )
    }, error = function(e) {
      results$methods$Bayesian <- list(detected = FALSE, error = e$message)
    })
  } else {
    results$methods$Bayesian <- list(detected = FALSE, error = "bcp package not available")
  }

  # Method 6: Visual/Statistical heuristic - largest jump
  tryCatch({
    diffs <- abs(diff(y_clean))
    if (length(diffs) > 0) {
      max_jump_idx <- which.max(diffs)
      jump_size <- diffs[max_jump_idx]
      relative_jump <- jump_size / sd(y_clean)

      # Heuristic: consider it a change if jump > 1.5 standard deviations
      results$methods$Visual <- list(
        detected = relative_jump > 1.5,
        breakpoints = max_jump_idx + 1,  # +1 because diff reduces length by 1
        jump_size = jump_size,
        relative_jump = relative_jump
      )
    } else {
      results$methods$Visual <- list(detected = FALSE)
    }
  }, error = function(e) {
    results$methods$Visual <- list(detected = FALSE, error = e$message)
  })

  # Method 7: Simple Moving Average Crossover
  tryCatch({
    if(n >= 6) {
      # Calculate short and long moving averages
      short_ma <- zoo::rollmean(y_clean, k = 2, fill = NA, align = "right")
      long_ma <- zoo::rollmean(y_clean, k = 4, fill = NA, align = "right")

      # Find crossover points
      crossovers <- which(diff(sign(short_ma - long_ma)) != 0)

      # Filter for significant crossovers
      if(length(crossovers) > 0) {
        # Calculate significance of crossovers
        crossover_strength <- sapply(crossovers, function(i) {
          if(i > 2 && i < (n-1)) {
            before_mean <- mean(y_clean[max(1, i-2):i])
            after_mean <- mean(y_clean[i:(min(n, i+2))])
            abs(after_mean - before_mean) / sd(y_clean)
          } else {
            0
          }
        })

        significant_crossovers <- crossovers[crossover_strength > 0.5]

        results$methods$MovingAverage <- list(
          detected = length(significant_crossovers) > 0,
          breakpoints = significant_crossovers,
          n_breaks = length(significant_crossovers),
          max_strength = ifelse(length(crossover_strength) > 0, max(crossover_strength), 0)
        )
      } else {
        results$methods$MovingAverage <- list(detected = FALSE)
      }
    } else {
      results$methods$MovingAverage <- list(detected = FALSE, error = "Insufficient data for moving average")
    }
  }, error = function(e) {
    results$methods$MovingAverage <- list(detected = FALSE, error = e$message)
  })

  # Method 8: Variance Change Detection
  tryCatch({
    if(n >= 6) {
      # Split data into segments and test variance changes
      split_points <- seq(3, n-2, by = 2)  # Possible split points
      variance_ratios <- sapply(split_points, function(i) {
        var1 <- var(y_clean[1:i])
        var2 <- var(y_clean[(i+1):n])
        max(var1, var2) / min(var1, var2)
      })

      # Find point with maximum variance ratio
      max_ratio_idx <- which.max(variance_ratios)
      max_ratio <- variance_ratios[max_ratio_idx]

      # Threshold for significant variance change
      results$methods$Variance <- list(
        detected = max_ratio > 4,  # 4:1 ratio threshold
        breakpoints = split_points[max_ratio_idx],
        variance_ratio = max_ratio
      )
    } else {
      results$methods$Variance <- list(detected = FALSE, error = "Insufficient data for variance test")
    }
  }, error = function(e) {
    results$methods$Variance <- list(detected = FALSE, error = e$message)
  })

  return(results)
}

# 3. Apply multi-method analysis ----------------------------------------------
parameters <- c("WQI", "DO", "EC", "Temp_Air", "Rainfall", "RH", "Temp_Water")
all_results <- list()

cat("=== MULTI-METHOD CHANGE POINT ANALYSIS ===\n")
cat("Available packages:\n")
cat("- strucchange: YES\n")
cat("- changepoint:", ifelse(changepoint_available, "YES", "NO"), "\n")
cat("- bcp:", ifelse(bcp_available, "YES", "NO"), "\n")
cat("- zoo: YES\n\n")

# Check available parameters in data
available_params <- intersect(parameters, names(data))
cat("Available parameters in data:", paste(available_params, collapse = ", "), "\n")

# Show data structure for debugging
cat("\nData structure:\n")
cat("Sites:", paste(unique(data$Site), collapse = ", "), "\n")
cat("Habitats:", paste(unique(data$Habitat), collapse = ", "), "\n")
cat("Total rows:", nrow(data), "\n\n")

analysis_count <- 0
successful_analyses <- 0

for (site in unique(data$Site)) {
  for (habitat in unique(data$Habitat)) {
    for (param in available_params) {

      analysis_count <- analysis_count + 1

      # Subset data
      df_sub <- data |>
        filter(Site == site, Habitat == habitat) |>
        arrange(Time_index)

      # Check if parameter has data
      non_na_count <- sum(!is.na(df_sub[[param]]))

      if (non_na_count < 4) {
        cat(sprintf("Skipping %s in %s (%s): insufficient data (%d non-NA values)\n",
                   param, habitat, site, non_na_count))
        next
      }

      # Apply multi-method detection
      cat(sprintf("Processing %s in %s (%s)...\n", param, habitat, site))

      result <- tryCatch({
        detect_changes_multi(df_sub[[param]], param, site, habitat)
      }, error = function(e) {
        cat(sprintf("Error processing %s in %s (%s): %s\n", param, habitat, site, e$message))
        return(NULL)
      })

      if (!is.null(result)) {
        result$data <- df_sub
        all_results[[length(all_results) + 1]] <- result
        successful_analyses <- successful_analyses + 1

        # Print summary
        methods_detected <- sapply(result$methods, function(x) {
          if(is.null(x$detected)) return(FALSE)
          return(x$detected)
        })

        n_methods <- sum(methods_detected, na.rm = TRUE)

        cat(sprintf("  -> %d/%d methods detected changes\n",
                   n_methods, length(methods_detected)))

        if (n_methods > 0) {
          detected_methods <- names(methods_detected)[methods_detected]
          cat(sprintf("  -> Methods: %s\n", paste(detected_methods, collapse = ", ")))
        }
      }
    }
  }
}

cat(sprintf("\n=== PROCESSING SUMMARY ===\n"))
cat(sprintf("Total analyses attempted: %d\n", analysis_count))
cat(sprintf("Successful analyses: %d\n", successful_analyses))
cat(sprintf("Results collected: %d\n", length(all_results)))

if(length(all_results) == 0) {
  stop("No successful analyses completed. Please check your data and parameters.")
}

# 4. Create comprehensive summary table ---------------------------------------
summary_table <- map_dfr(all_results, function(result) {

  # Create a safe function to extract method results
  extract_method_info <- function(method, method_name) {
    tryCatch({
      # Ensure all required fields exist
      detected <- ifelse(is.null(method$detected), FALSE, method$detected)
      n_breaks <- ifelse(is.null(method$n_breaks),
                        ifelse(detected, 1, 0),
                        method$n_breaks)

      # Handle breakpoints - ensure it's a character
      if(is.null(method$breakpoints) || length(method$breakpoints) == 0) {
        breakpoints_str <- NA_character_
      } else {
        breakpoints_str <- paste(method$breakpoints, collapse = ",")
      }

      # Create additional info string
      additional_info <- ""
      if(!is.null(method$p_value) && !is.na(method$p_value)) {
        additional_info <- sprintf("p=%.3f", method$p_value)
      } else if(!is.null(method$bic) && !is.na(method$bic)) {
        additional_info <- sprintf("BIC=%.1f", method$bic)
      } else if(!is.null(method$max_prob) && !is.na(method$max_prob)) {
        additional_info <- sprintf("prob=%.3f", method$max_prob)
      } else if(!is.null(method$relative_jump) && !is.na(method$relative_jump)) {
        additional_info <- sprintf("jump=%.2f", method$relative_jump)
      } else if(!is.null(method$variance_ratio) && !is.na(method$variance_ratio)) {
        additional_info <- sprintf("var_ratio=%.2f", method$variance_ratio)
      } else if(!is.null(method$max_strength) && !is.na(method$max_strength)) {
        additional_info <- sprintf("strength=%.2f", method$max_strength)
      } else if(!is.null(method$error)) {
        additional_info <- sprintf("ERROR: %s", method$error)
      }

      return(data.frame(
        method = method_name,
        detected = detected,
        n_breaks = n_breaks,
        breakpoints = breakpoints_str,
        additional_info = additional_info,
        stringsAsFactors = FALSE
      ))
    }, error = function(e) {
      # Return a safe default if something goes wrong
      return(data.frame(
        method = method_name,
        detected = FALSE,
        n_breaks = 0,
        breakpoints = NA_character_,
        additional_info = sprintf("PROCESSING_ERROR: %s", e$message),
        stringsAsFactors = FALSE
      ))
    })
  }

  # Process each method safely
  methods_list <- list()
  for(method_name in names(result$methods)) {
    methods_list[[method_name]] <- extract_method_info(result$methods[[method_name]], method_name)
  }

  # Combine all method results
  if(length(methods_list) > 0) {
    methods_df <- do.call(rbind, methods_list)

    # Add result metadata
    methods_df$parameter <- result$param
    methods_df$site <- result$site
    methods_df$habitat <- result$habitat
    methods_df$n_obs <- result$n

    return(methods_df)
  } else {
    # Return empty data frame with correct structure if no methods
    return(data.frame(
      method = character(0),
      detected = logical(0),
      n_breaks = numeric(0),
      breakpoints = character(0),
      additional_info = character(0),
      parameter = character(0),
      site = character(0),
      habitat = character(0),
      n_obs = numeric(0),
      stringsAsFactors = FALSE
    ))
  }
})

# 5. Create visualizations ---------------------------------------------------
if(nrow(summary_table) == 0) {
  cat("No data available for visualization - skipping plots\n")
} else {

  # Summary by method
  method_summary <- summary_table |>
    group_by(method) |>
    summarise(
      total_analyses = n(),
      detections = sum(detected, na.rm = TRUE),
      errors = sum(grepl("ERROR", additional_info)),
      detection_rate = ifelse(total_analyses > 0, detections / total_analyses * 100, 0),
      .groups = "drop"
    )

  p1 <- ggplot(method_summary, aes(x = reorder(method, detection_rate), y = detection_rate)) +
    geom_col(fill = "steelblue", alpha = 0.7) +
    geom_text(aes(label = sprintf("%.1f%%\n(%d/%d)", detection_rate, detections, total_analyses)),
              hjust = -0.1, size = 3) +
    coord_flip() +
    labs(
      title = "Change Point Detection Rate by Method",
      x = "Method",
      y = "Detection Rate (%)"
    ) +
    theme_minimal()

  # Summary by parameter
  param_summary <- summary_table |>
    group_by(parameter) |>
    summarise(
      total_analyses = n(),
      detections = sum(detected, na.rm = TRUE),
      detection_rate = ifelse(total_analyses > 0, detections / total_analyses * 100, 0),
      .groups = "drop"
    )

  p2 <- ggplot(param_summary, aes(x = reorder(parameter, detection_rate), y = detection_rate)) +
    geom_col(fill = "darkgreen", alpha = 0.7) +
    geom_text(aes(label = sprintf("%.1f%%", detection_rate)),
              hjust = -0.1, size = 3) +
    coord_flip() +
    labs(
      title = "Change Point Detection Rate by Parameter",
      x = "Parameter",
      y = "Detection Rate (%)"
    ) +
    theme_minimal()

  # Consensus analysis (cases where multiple methods agree)
  consensus_summary <- summary_table |>
    group_by(parameter, site, habitat) |>
    summarise(
      methods_detecting = sum(detected, na.rm = TRUE),
      total_methods = n(),
      consensus_strength = ifelse(total_methods > 0, methods_detecting / total_methods, 0),
      .groups = "drop"
    ) |>
    filter(methods_detecting >= 2)  # At least 2 methods agree

  if (nrow(consensus_summary) > 0) {
    p3 <- ggplot(consensus_summary, aes(x = parameter, y = consensus_strength,
                                       color = factor(methods_detecting))) +
      geom_jitter(size = 3, alpha = 0.7, width = 0.2) +
      facet_grid(site ~ habitat) +
      scale_y_continuous(labels = scales::percent) +
      labs(
        title = "Consensus Change Points (>=2 Methods)",
        x = "Parameter",
        y = "Consensus Strength",
        color = "# Methods\nDetecting"
      ) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  } else {
    p3 <- ggplot() +
      annotate("text", x = 0.5, y = 0.5, label = "No consensus detections\n(>=2 methods)",
               size = 6, color = "gray50") +
      theme_void()
  }
}

# 6. Save results and plots --------------------------------------------------
# Save detailed results
if(nrow(summary_table) > 0) {
  tryCatch({
    write.csv(summary_table, "multi_method_changepoint_results.csv", row.names = FALSE)
    cat("Saved: multi_method_changepoint_results.csv\n")
  }, error = function(e) {
    cat("Error saving summary_table:", e$message, "\n")
  })

  tryCatch({
    write.csv(method_summary, "method_performance_summary.csv", row.names = FALSE)
    cat("Saved: method_performance_summary.csv\n")
  }, error = function(e) {
    cat("Error saving method_summary:", e$message, "\n")
  })

  tryCatch({
    write.csv(param_summary, "parameter_changepoint_summary.csv", row.names = FALSE)
    cat("Saved: parameter_changepoint_summary.csv\n")
  }, error = function(e) {
    cat("Error saving param_summary:", e$message, "\n")
  })

  if (exists("consensus_summary") && nrow(consensus_summary) > 0) {
    tryCatch({
      write.csv(consensus_summary, "consensus_changepoints.csv", row.names = FALSE)
      cat("Saved: consensus_changepoints.csv\n")
    }, error = function(e) {
      cat("Error saving consensus_summary:", e$message, "\n")
    })
  }

  # Save plots
  tryCatch({
    ggsave("method_comparison.png", p1, width = 10, height = 6)
    cat("Saved: method_comparison.png\n")
  }, error = function(e) {
    cat("Error saving method_comparison.png:", e$message, "\n")
  })

  tryCatch({
    ggsave("parameter_comparison.png", p2, width = 10, height = 6)
    cat("Saved: parameter_comparison.png\n")
  }, error = function(e) {
    cat("Error saving parameter_comparison.png:", e$message, "\n")
  })

  tryCatch({
    ggsave("consensus_analysis.png", p3, width = 12, height = 8)
    cat("Saved: consensus_analysis.png\n")
  }, error = function(e) {
    cat("Error saving consensus_analysis.png:", e$message, "\n")
  })

  # Combined summary plot
  tryCatch({
    combined_plot <- grid.arrange(p1, p2, ncol = 2)
    ggsave("combined_summary.png", combined_plot, width = 16, height = 8)
    cat("Saved: combined_summary.png\n")
  }, error = function(e) {
    cat("Error saving combined_summary.png:", e$message, "\n")
  })

} else {
  cat("No data to save - skipping file outputs\n")
}

# 7. Print final summary ------------------------------------------------------
cat("\n=== FINAL SUMMARY ===\n")
cat("Total parameter-site-habitat combinations analyzed:",
    length(unique(paste(summary_table$parameter, summary_table$site, summary_table$habitat))), "\n")

cat("\nMethod Performance:\n")
print(method_summary)

cat("\nParameter Sensitivity:\n")
print(param_summary)

if (nrow(consensus_summary) > 0) {
  cat("\nConsensus Detections (>=2 methods):\n")
  print(consensus_summary)
} else {
  cat("\nNo consensus detections found (no cases where >=2 methods agree)\n")
}

cat("\nFiles saved:\n")
cat("- multi_method_changepoint_results.csv\n")
cat("- method_performance_summary.csv\n")
cat("- parameter_changepoint_summary.csv\n")
cat("- method_comparison.png\n")
cat("- parameter_comparison.png\n")
cat("- consensus_analysis.png\n")
cat("- combined_summary.png\n")

cat("\nAnalysis complete!\n")
