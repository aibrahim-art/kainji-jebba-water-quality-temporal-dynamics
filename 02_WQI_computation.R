# ===================================================================
# 02_WQI_COMPUTATION.R
# Comprehensive WQI Analysis Script
# Author: Ibrahim et al.
# Description:
#   - Reads in WQI dataset
#   - Computes summary statistics
#   - Creates monthly, spatial, and distribution plots
#   - Saves results as PNG plots and CSV
#
# Required local input (not included in this repository; see Data
# availability in README.md):
#   - WQI_data.csv, containing columns: Site, Habitat, Season, Year,
#     SimpleSeason, WQI_Class, WQI, Month, Longitude, Latitude
# ===================================================================

# -------------------- Load Packages --------------------
required_packages <- c("tidyverse", "ggplot2", "dplyr")
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)
lapply(required_packages, library, character.only = TRUE)

# -------------------- Load Data --------------------
wqi_data <- read.csv("WQI_data.csv")

# Ensure correct structure
str(wqi_data)

# -------------------- Summary Statistics --------------------
compute_summary_stats <- function(data) {
  stats <- data %>%
    group_by(Site, Habitat, Season, Year, SimpleSeason, WQI_Class) %>%
    summarise(
      Mean_WQI = mean(WQI, na.rm = TRUE),
      SD_WQI   = sd(WQI, na.rm = TRUE),
      Median_WQI = median(WQI, na.rm = TRUE),
      Min_WQI  = min(WQI, na.rm = TRUE),
      Max_WQI  = max(WQI, na.rm = TRUE),
      n = n(),
      .groups = "drop"
    )
  return(stats)
}

stats <- compute_summary_stats(wqi_data)

# -------------------- Plot Functions --------------------
create_comprehensive_plots <- function(data, stats) {

  # Monthly Trend Plot
  p1 <- ggplot(data, aes(x = Month, y = WQI, group = Site, color = Habitat)) +
    geom_line(alpha = 0.6) +
    geom_point(size = 2) +
    geom_smooth(aes(group = 1), method = "loess", se = FALSE, color = "black", linetype = "dashed") +
    labs(title = "Monthly WQI Trends", y = "WQI", x = "Month") +
    theme_minimal()

  # Spatial Plot
  p2 <- ggplot(data, aes(x = Longitude, y = Latitude, color = WQI, shape = Habitat)) +
    geom_point(size = 4) +
    scale_color_viridis_c() +
    labs(title = "Spatial WQI Distribution", x = "Longitude", y = "Latitude") +
    theme_minimal()

  # WQI Class Distribution Plot
  p3 <- ggplot(data, aes(x = WQI_Class, fill = Habitat)) +
    geom_bar(position = "dodge") +
    labs(title = "Distribution of WQI Classes by Habitat", x = "WQI Class", y = "Count") +
    theme_minimal()

  return(list(monthly = p1, spatial = p2, distribution = p3))
}

plots <- create_comprehensive_plots(wqi_data, stats)

# -------------------- Save Outputs --------------------
# Save plots
ggsave("WQI_Monthly_Trends.png", plots$monthly, width = 8, height = 5, dpi = 300)
ggsave("WQI_Spatial.png", plots$spatial, width = 8, height = 5, dpi = 300)
ggsave("WQI_Distribution.png", plots$distribution, width = 8, height = 5, dpi = 300)

# Save summary statistics
write.csv(stats, "WQI_Summary_Stats.csv", row.names = FALSE)

# -------------------- End of Script --------------------
cat("Analysis complete. Plots and summary stats saved in working directory:\n", getwd(), "\n")
