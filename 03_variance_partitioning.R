# ==============================================================================
# 03_VARIANCE_PARTITIONING.R
# Variance Partitioning Analysis
# Water Quality Drivers in Niger River Reservoir Systems
#
# Required local input (not included in this repository; see Data
# availability in README.md):
#   - water_quality_data.csv, with columns: Site, Habitat, Time_index, Month,
#     WQI, Temp_Air, RH, Rainfall, Windspeed, Pressure, Temp_Water, pH, DO,
#     BOD, NO3_N, PO4_P, Turbidity_NTU, Hardness, EC
# ==============================================================================

library(tidyverse)
library(vegan)
library(ggplot2)
library(viridis)

# ==============================================================================
# 1. LOAD AND PREPARE DATA
# ==============================================================================

wq_data <- read.csv("water_quality_data.csv", stringsAsFactors = FALSE) %>%
  mutate(
    Season = case_when(
      Month %in% c("Nov", "Dec", "Jan", "Feb", "Mar") ~ "Dry",
      TRUE ~ "Wet"
    ),
    Site_Habitat = paste(Site, Habitat, sep = "_")
  )

names(wq_data) <- gsub(" ", "", names(wq_data))

# Add land use data (ESA WorldCover 2021, 5-km buffer extraction;
# see 01_land_use_extraction.R for the extraction pipeline)
land_use_exact <- data.frame(
  Site_Habitat = c("Kainji_Riverine", "Kainji_Ecotonal", "Kainji_Locustrine",
                   "Jebba_Riverine", "Jebba_Ecotonal", "Jebba_Locustrine"),
  Site_Name = c("Garafini", "Awuru", "Yuna", "Gbajibo", "Juju Rock", "Kokodi"),
  Forest_pct = c(3.16, 8.45, 3.78, 2.17, 20.56, 19.86),
  Agriculture_pct = c(12.32, 46.12, 21.48, 17.99, 20.70, 11.08),
  Urban_pct = c(0.26, 0.28, 0.50, 0.25, 5.82, 0.15)
)

wq_data <- wq_data %>%
  left_join(land_use_exact, by = "Site_Habitat")

# ==============================================================================
# 2. PREPARE MATRICES FOR RDA
# ==============================================================================

# Response matrix (unstandardized for distance calculations)
Y_raw <- wq_data %>%
  select(DO, pH, BOD, Temp_Water, NO3_N, PO4_P, Turbidity_NTU, EC) %>%
  mutate(Turbidity_NTU = log10(Turbidity_NTU + 1))

# Standardized for RDA
Y_std <- Y_raw %>%
  mutate(across(everything(), ~scale(.)[,1]))

# Predictor matrices
X_spatial <- wq_data %>%
  mutate(
    Dam_Jebba = as.numeric(Site == "Jebba"),
    Habitat_Ecotonal = as.numeric(Habitat == "Ecotonal"),
    Habitat_Lacustrine = as.numeric(Habitat == "Locustrine")
  ) %>%
  select(Dam_Jebba, Habitat_Ecotonal, Habitat_Lacustrine)

X_temporal <- wq_data %>%
  mutate(
    Season_Wet = as.numeric(Season == "Wet"),
    Month_num = Time_index
  ) %>%
  select(Season_Wet, Month_num) %>%
  mutate(across(everything(), ~scale(.)[,1]))

X_anthropogenic <- wq_data %>%
  select(Agriculture_pct, Urban_pct, Forest_pct) %>%
  mutate(across(everything(), ~scale(.)[,1]))

# ==============================================================================
# 3. VARIANCE PARTITIONING
# ==============================================================================

set.seed(123)

# Pure effects using partial RDA
rda_spatial_pure <- rda(Y_std ~ Dam_Jebba + Habitat_Ecotonal + Habitat_Lacustrine +
                         Condition(Season_Wet + Month_num + Agriculture_pct +
                                   Urban_pct + Forest_pct),
                        data = cbind(X_spatial, X_temporal, X_anthropogenic))

rda_temporal_pure <- rda(Y_std ~ Season_Wet + Month_num +
                          Condition(Dam_Jebba + Habitat_Ecotonal + Habitat_Lacustrine +
                                    Agriculture_pct + Urban_pct + Forest_pct),
                         data = cbind(X_spatial, X_temporal, X_anthropogenic))

rda_anthropogenic_pure <- rda(Y_std ~ Agriculture_pct + Urban_pct + Forest_pct +
                               Condition(Dam_Jebba + Habitat_Ecotonal + Habitat_Lacustrine +
                                         Season_Wet + Month_num),
                              data = cbind(X_spatial, X_temporal, X_anthropogenic))

# Full model
rda_full <- rda(Y_std ~ ., data = cbind(X_spatial, X_temporal, X_anthropogenic))

# Significance testing
anova_spatial <- anova(rda_spatial_pure, permutations = 999)
anova_temporal <- anova(rda_temporal_pure, permutations = 999)
anova_anthropogenic <- anova(rda_anthropogenic_pure, permutations = 999)

# ==============================================================================
# 4. COMPILE RESULTS
# ==============================================================================

var_results <- tibble(
  Component = c("Temporal (Season + Month)",
                "Spatial (Dam + Habitat)",
                "Anthropogenic (Land Use)",
                "Residual"),

  Variance_pct = c(
    RsquareAdj(rda_temporal_pure)$adj.r.squared * 100,
    RsquareAdj(rda_spatial_pure)$adj.r.squared * 100,
    RsquareAdj(rda_anthropogenic_pure)$adj.r.squared * 100,
    (1 - RsquareAdj(rda_full)$adj.r.squared) * 100
  ),

  F_value = c(
    anova_temporal$F[1],
    anova_spatial$F[1],
    anova_anthropogenic$F[1],
    NA
  ),

  p_value = c(
    anova_temporal$`Pr(>F)`[1],
    anova_spatial$`Pr(>F)`[1],
    anova_anthropogenic$`Pr(>F)`[1],
    NA
  ),

  Significance = case_when(
    is.na(p_value) ~ "",
    p_value < 0.001 ~ "***",
    p_value < 0.01 ~ "**",
    p_value < 0.05 ~ "*",
    TRUE ~ "ns"
  )
)

# ==============================================================================
# 5. PERMANOVA
# ==============================================================================

# Option 1: Euclidean distance (equivalent to RDA)
dist_euclidean <- vegdist(Y_std, method = "euclidean")

permanova_euclidean <- adonis2(
  dist_euclidean ~ Dam_Jebba + Habitat_Ecotonal + Habitat_Lacustrine +
    Season_Wet + Month_num + Agriculture_pct + Urban_pct + Forest_pct,
  data = cbind(X_spatial, X_temporal, X_anthropogenic),
  permutations = 999,
  by = "margin"
)

# Option 2: Hellinger transformation + Bray-Curtis
Y_hellinger <- decostand(Y_raw, method = "hellinger")
dist_hellinger <- vegdist(Y_hellinger, method = "bray")

permanova_hellinger <- adonis2(
  dist_hellinger ~ Dam_Jebba + Habitat_Ecotonal + Habitat_Lacustrine +
    Season_Wet + Month_num + Agriculture_pct + Urban_pct + Forest_pct,
  data = cbind(X_spatial, X_temporal, X_anthropogenic),
  permutations = 999,
  by = "margin"
)

# ==============================================================================
# 6. NMDS
# ==============================================================================

set.seed(123)

# Use Hellinger transformation
nmds_hellinger <- metaMDS(Y_hellinger, distance = "bray", k = 2, trymax = 100)

# Extract scores
nmds_scores <- as.data.frame(scores(nmds_hellinger, display = "sites")) %>%
  bind_cols(wq_data %>% select(Site, Habitat, Season, Month))

# ==============================================================================
# 7. REPORTING
# ==============================================================================

cat("\n===============================================================\n")
cat("               VARIANCE PARTITIONING RESULTS                  \n")
cat("===============================================================\n")
print(var_results, width = 100)
cat("\n")

# Format for manuscript
cat("MANUSCRIPT TEXT:\n")
cat("---------------------------------------------------------------\n")
for(i in 1:3) {
  cat(sprintf("%s: %.1f%% (F = %.2f, p = %.3f) %s\n",
              var_results$Component[i],
              var_results$Variance_pct[i],
              var_results$F_value[i],
              var_results$p_value[i],
              var_results$Significance[i]))
}
cat("\n")

cat("PERMANOVA (Euclidean distance):\n")
cat("---------------------------------------------------------------\n")
print(permanova_euclidean)
cat("\n")

cat("PERMANOVA (Bray-Curtis on Hellinger):\n")
cat("---------------------------------------------------------------\n")
print(permanova_hellinger)
cat("\n")

cat(sprintf("NMDS Stress: %.3f ", nmds_hellinger$stress))
if(nmds_hellinger$stress < 0.05) {
  cat("(excellent)\n")
} else if(nmds_hellinger$stress < 0.1) {
  cat("(good)\n")
} else if(nmds_hellinger$stress < 0.2) {
  cat("(acceptable)\n")
} else {
  cat("(poor)\n")
}
cat("\n")

# ==============================================================================
# 8. PUBLICATION FIGURES
# ==============================================================================

# Variance partitioning figure
p_variance <- ggplot(var_results,
                     aes(x = Variance_pct,
                         y = reorder(Component, Variance_pct),
                         fill = Component)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = sprintf("%.1f%%", Variance_pct)),
            hjust = -0.1, size = 5.5, fontface = "bold") +
  geom_text(aes(label = Significance,
                x = Variance_pct + max(Variance_pct, na.rm = TRUE)*0.06),
            hjust = 0, size = 7) +
  scale_fill_viridis_d(option = "plasma", begin = 0.2, end = 0.9) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(
    title = "Drivers of Water Quality Variation",
    subtitle = "Niger River Reservoir Systems (Variance Partitioning)",
    x = "Variance Explained (%)",
    y = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    panel.grid.major.y = element_blank(),
    plot.title = element_text(face = "bold", size = 16),
    axis.text = element_text(size = 12)
  )

ggsave("FIG_VariancePartitioning.png", p_variance,
       width = 10, height = 6, dpi = 300)

# NMDS plot with environmental fitting
env_fit <- envfit(nmds_hellinger, cbind(X_spatial, X_temporal, X_anthropogenic),
                  permutations = 999)

# Extract significant vectors
env_scores <- as.data.frame(scores(env_fit, display = "vectors")) %>%
  rownames_to_column("Variable") %>%
  mutate(
    pval = env_fit$vectors$pvals,
    r2 = env_fit$vectors$r,
    Significant = pval < 0.05
  ) %>%
  filter(Significant)

p_nmds <- ggplot(nmds_scores, aes(x = NMDS1, y = NMDS2,
                                   color = Site, shape = Habitat)) +
  geom_point(size = 4, alpha = 0.7) +
  stat_ellipse(aes(group = Site), level = 0.95, linetype = 2) +
  geom_segment(data = env_scores,
               aes(x = 0, y = 0, xend = NMDS1*0.5, yend = NMDS2*0.5),
               arrow = arrow(length = unit(0.3, "cm")),
               color = "black", inherit.aes = FALSE, linewidth = 0.8) +
  geom_text(data = env_scores,
            aes(x = NMDS1*0.55, y = NMDS2*0.55, label = Variable),
            color = "black", size = 3.5, fontface = "bold",
            inherit.aes = FALSE) +
  scale_color_viridis_d(option = "viridis", end = 0.8) +
  scale_shape_manual(values = c(16, 17, 15)) +
  labs(
    title = "NMDS Ordination of Water Quality",
    subtitle = sprintf("Stress = %.3f (Hellinger-transformed, Bray-Curtis)",
                      nmds_hellinger$stress),
    color = "Reservoir",
    shape = "Habitat Type"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold")
  )

ggsave("FIG_NMDS.png", p_nmds,
       width = 12, height = 7, dpi = 300)

# Save results
write.csv(var_results, "TABLE_VariancePartitioning_Final.csv", row.names = FALSE)

cat("\nAnalysis complete. Figures and table saved.\n")
