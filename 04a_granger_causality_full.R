# ============================================================
# 04a_GRANGER_CAUSALITY_FULL.R
# Complete Granger Causality Analysis - Full 480-Test Export
# Kainji-Jebba Reservoir Cascade Water Quality Study
#
# Required local input (not included in this repository; see Data
# availability in README.md):
#   - Grouped_Dataset_Kainji_Jebba.xlsx (sheet 1), with columns:
#     Site, Habitat, Month, Time_index, and the climate/water-quality
#     variables listed in Section 4 below
# ============================================================

# %% SECTION 1: LOAD PACKAGES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
required_packages <- c("readxl", "dplyr", "lmtest", "tseries")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
  library(pkg, character.only = TRUE)
}

# %% SECTION 2: LOAD DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data_path <- "Grouped_Dataset_Kainji_Jebba.xlsx"
df <- readxl::read_excel(data_path, sheet = 1)

df$Site    <- trimws(df$Site)
df$Habitat <- trimws(df$Habitat)
df$Month   <- trimws(df$Month)

cat("Data loaded:", nrow(df), "rows x", ncol(df), "columns\n")
cat("Sites:   ", paste(unique(df$Site),    collapse=", "), "\n")
cat("Habitats:", paste(unique(df$Habitat), collapse=", "), "\n")

# %% SECTION 3: DEFINE SEASON %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
wet_months <- c("Apr","May","Jun","Jul","Aug","Sep","Oct")
dry_months <- c("Jan","Feb","Mar","Nov","Dec")
df$Season  <- ifelse(df$Month %in% wet_months, "Wet", "Dry")
cat("Season counts:\n"); print(table(df$Season))

# %% SECTION 4: DEFINE ANALYSIS VARIABLES %%%%%%%%%%%%%%%%%%%%%
climate_vars <- c("Pressure", "Rainfall", "RH", "Temp_Air", "Windspeed")

# Temp_Water included as 8th WQ response variable
wq_vars      <- c("BOD", "DO", "EC", "NO3_N", "pH", "PO4_P", "Turbidity", "Temp_Water")

dams         <- c("Jebba", "Kainji")
habitats     <- c("Riverine", "Ecotonal", "Lacustrine")
seasons      <- c("Dry", "Wet")

cat("\nClimate drivers:", paste(climate_vars, collapse=", "), "\n")
cat("WQ parameters: ", paste(wq_vars,      collapse=", "), "\n")
cat("Expected tests:", length(dams) * length(habitats) * length(seasons) *
      length(climate_vars) * length(wq_vars), "\n\n")

# %% SECTION 5: HELPER FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# 5a. Stationarity via ADF + first differencing
make_stationary <- function(x) {
  x <- as.numeric(x)
  x <- x[!is.na(x)]
  if (length(x) < 5) return(list(series = x, differenced = FALSE))
  tryCatch({
    adf_result <- tseries::adf.test(x, alternative = "stationary")
    if (adf_result$p.value > 0.05) {
      return(list(series = diff(x), differenced = TRUE))
    }
    return(list(series = x, differenced = FALSE))
  }, error = function(e) {
    return(list(series = x, differenced = FALSE))
  })
}

# 5b. Select optimal lag by AIC (df-safe)
select_lag_aic <- function(cause, effect, max_lag = 3) {
  n        <- length(effect)
  best_lag <- 1
  best_aic <- Inf

  max_safe_lag <- floor((n - 2) / 3)
  upper        <- min(max_lag, max_safe_lag)

  if (upper < 1) return(1L)

  for (lag in 1:upper) {
    tryCatch({
      T_len  <- n - lag
      Y      <- effect[(lag + 1):n]
      X_cols <- list(intercept = rep(1, T_len))
      for (l in 1:lag) {
        X_cols[[paste0("e_lag", l)]] <- effect[(lag + 1 - l):(n - l)]
        X_cols[[paste0("c_lag", l)]] <- cause[(lag  + 1 - l):(n - l)]
      }
      X     <- do.call(cbind, X_cols)
      model <- lm(Y ~ X - 1)
      aic   <- AIC(model)
      if (is.finite(aic) && aic < best_aic) {
        best_aic <- aic
        best_lag <- lag
      }
    }, error = function(e) NULL)
  }
  return(best_lag)
}

# 5c. Run Granger test
run_granger <- function(cause, effect, lag) {
  tryCatch({
    series_df <- data.frame(effect = effect, cause = cause)
    gt        <- lmtest::grangertest(effect ~ cause,
                                     order = lag,
                                     data  = series_df)
    return(gt$`Pr(>F)`[2])
  }, error = function(e) NA_real_)
}

# %% SECTION 6: MAIN ANALYSIS LOOP %%%%%%%%%%%%%%%%%%%%%%%%%%%%
cat("Running Granger causality tests...\n")
results_list <- list()
row_idx      <- 0

for (dam in dams) {
  for (hab in habitats) {
    for (season in seasons) {

      stratum <- df %>%
        filter(Site == dam, Habitat == hab, Season == season) %>%
        arrange(Time_index)

      n_obs <- nrow(stratum)
      cat(sprintf("  %s | %s | %s   %d observations\n",
                  dam, hab, season, n_obs))

      for (cause_var in climate_vars) {
        for (effect_var in wq_vars) {

          row_idx    <- row_idx + 1
          cause_raw  <- as.numeric(stratum[[cause_var]])
          effect_raw <- as.numeric(stratum[[effect_var]])

          valid      <- complete.cases(cause_raw, effect_raw)
          cause_raw  <- cause_raw[valid]
          effect_raw <- effect_raw[valid]
          n_valid    <- length(cause_raw)

          if (n_valid < 5) {
            results_list[[row_idx]] <- data.frame(
              Site = dam, Habitat = hab, Season = season,
              Cause = cause_var, Effect = effect_var,
              N_obs = n_valid, Lag = NA_integer_,
              P_Value = NA_real_, Adj_P = NA_real_,
              Significant = NA, Notes = "Insufficient observations",
              stringsAsFactors = FALSE)
            next
          }

          cause_st  <- make_stationary(cause_raw)
          effect_st <- make_stationary(effect_raw)
          cause_s   <- cause_st$series
          effect_s  <- effect_st$series

          min_len  <- min(length(cause_s), length(effect_s))
          cause_s  <- tail(cause_s,  min_len)
          effect_s <- tail(effect_s, min_len)

          if (min_len < 5) {
            results_list[[row_idx]] <- data.frame(
              Site = dam, Habitat = hab, Season = season,
              Cause = cause_var, Effect = effect_var,
              N_obs = min_len, Lag = NA_integer_,
              P_Value = NA_real_, Adj_P = NA_real_,
              Significant = NA, Notes = "Too short after differencing",
              stringsAsFactors = FALSE)
            next
          }

          opt_lag <- select_lag_aic(cause_s, effect_s, max_lag = 3)
          p_val   <- run_granger(cause_s, effect_s, opt_lag)

          notes <- ifelse(
            cause_st$differenced | effect_st$differenced,
            "First-differenced", "Levels")

          results_list[[row_idx]] <- data.frame(
            Site = dam, Habitat = hab, Season = season,
            Cause = cause_var, Effect = effect_var,
            N_obs = min_len, Lag = opt_lag,
            P_Value = round(p_val, 6),
            Adj_P = NA_real_, Significant = NA,
            Notes = notes,
            stringsAsFactors = FALSE)
        }
      }
    }
  }
}

# %% SECTION 7: COMBINE & FDR CORRECT %%%%%%%%%%%%%%%%%%%%%%%%
results_df <- do.call(rbind, results_list)
cat(sprintf("\nTotal tests completed: %d\n", nrow(results_df)))

valid_mask <- !is.na(results_df$P_Value)
adj_p_vals <- p.adjust(results_df$P_Value[valid_mask], method = "BH")
results_df$Adj_P[valid_mask]       <- round(adj_p_vals, 6)
results_df$Significant[valid_mask] <- adj_p_vals < 0.05

# %% SECTION 8: SUMMARY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cat("\n=== GRANGER CAUSALITY RESULTS SUMMARY ===\n")
cat(sprintf("Total tests:               %d\n", nrow(results_df)))
cat(sprintf("Tests with valid p-values: %d\n", sum(valid_mask)))
cat(sprintf("Significant after FDR:     %d\n",
            sum(results_df$Significant == TRUE, na.rm = TRUE)))

near_sig <- results_df[!is.na(results_df$P_Value) & results_df$P_Value < 0.05, ]
cat(sprintf("Near-significant (p<0.05): %d\n", nrow(near_sig)))

near_sig_sorted <- near_sig[order(near_sig$P_Value), ]
print(near_sig_sorted[, c("Site","Habitat","Season","Cause",
                           "Effect","Lag","P_Value","Adj_P")])

# %% SECTION 9: EXPORT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
write.csv(results_df,
          "Granger_Full_480_Results.csv", row.names = FALSE)

write.csv(near_sig_sorted[, c("Site","Habitat","Season","Cause",
                               "Effect","Lag","P_Value","Adj_P","Notes")],
          "Granger_NearSig_Table3.csv", row.names = FALSE)

cat("\nGranger_Full_480_Results.csv saved (", nrow(results_df), "rows)\n")
cat("Granger_NearSig_Table3.csv   saved (", nrow(near_sig_sorted), "rows)\n")
cat("\nDone.\n")
