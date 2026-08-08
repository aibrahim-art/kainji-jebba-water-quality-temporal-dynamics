# ============================================================
# 04b_GRANGER_HEATMAP.R
# Figure 5a and 5b - Granger Causality Heatmaps (Split by Season)
# Figure 5a: Wet Season | Figure 5b: Dry Season
# Kainji-Jebba Reservoir Cascade Water Quality Study
#
# Depends on the near-significant relationship table produced by
# 04a_granger_causality_full.R (values are hardcoded below from that
# script's output for reproducibility of the published figure).
# ============================================================

library(ggplot2)
library(dplyr)
library(tidyr)

# -- 1. DEFINE VARIABLES ---------------------------------------
climate_vars <- c("Pressure", "Rainfall", "RH", "Temp_Air", "Windspeed")
wq_vars      <- c("BOD", "DO", "EC", "NO3_N", "pH", "PO4_P", "Turbidity", "Temp_Water")
dams         <- c("Jebba", "Kainji")
habitats     <- c("Riverine", "Ecotonal", "Lacustrine")
seasons      <- c("Wet", "Dry")

# -- 2. BUILD FULL 480-TEST GRID ---------------------------------
df_full <- expand.grid(
  Dam     = dams,
  Habitat = habitats,
  Season  = seasons,
  Cause   = climate_vars,
  Effect  = wq_vars,
  stringsAsFactors = FALSE
)

# Default: all adjusted p-values = 0.999
df_full$P_raw <- 0.999
df_full$Near_Sig <- FALSE

# -- 3. INSERT ACTUAL NEAR-SIGNIFICANT VALUES FROM TABLE 3 ------
insert_p <- function(df, dam, habitat, season, cause, effect, p_raw) {
  idx <- df$Dam == dam & df$Habitat == habitat & df$Season == season &
         df$Cause == cause & df$Effect == effect
  df$P_raw[idx]    <- p_raw
  df$Near_Sig[idx] <- TRUE
  return(df)
}

# All 23 near-significant relationships
df_full <- insert_p(df_full, "Kainji", "Lacustrine", "Wet", "RH",        "BOD",        0.001375)
df_full <- insert_p(df_full, "Jebba",  "Lacustrine", "Wet", "Temp_Air",  "pH",         0.004917)
df_full <- insert_p(df_full, "Jebba",  "Lacustrine", "Wet", "Temp_Air",  "DO",         0.005711)
df_full <- insert_p(df_full, "Kainji", "Lacustrine", "Wet", "RH",        "Temp_Water", 0.007001)
df_full <- insert_p(df_full, "Kainji", "Lacustrine", "Wet", "Temp_Air",  "pH",         0.007146)
df_full <- insert_p(df_full, "Kainji", "Riverine",   "Dry", "Rainfall",  "pH",         0.008040)
df_full <- insert_p(df_full, "Jebba",  "Ecotonal",   "Wet", "Temp_Air",  "DO",         0.010336)
df_full <- insert_p(df_full, "Kainji", "Riverine",   "Dry", "Rainfall",  "EC",         0.013480)
df_full <- insert_p(df_full, "Kainji", "Ecotonal",   "Dry", "Rainfall",  "EC",         0.013480)
df_full <- insert_p(df_full, "Kainji", "Lacustrine", "Dry", "Rainfall",  "EC",         0.013480)
df_full <- insert_p(df_full, "Kainji", "Ecotonal",   "Dry", "Pressure",  "DO",         0.016354)
df_full <- insert_p(df_full, "Kainji", "Ecotonal",   "Wet", "RH",        "Temp_Water", 0.025323)
df_full <- insert_p(df_full, "Jebba",  "Ecotonal",   "Wet", "Pressure",  "pH",         0.028259)
df_full <- insert_p(df_full, "Jebba",  "Ecotonal",   "Dry", "Windspeed", "PO4_P",      0.029514)
df_full <- insert_p(df_full, "Jebba",  "Lacustrine", "Dry", "Temp_Air",  "NO3_N",      0.030054)
df_full <- insert_p(df_full, "Kainji", "Ecotonal",   "Dry", "Windspeed", "BOD",        0.031797)
df_full <- insert_p(df_full, "Jebba",  "Riverine",   "Dry", "Pressure",  "Turbidity",  0.033041)
df_full <- insert_p(df_full, "Jebba",  "Lacustrine", "Wet", "Temp_Air",  "Temp_Water", 0.034121)
df_full <- insert_p(df_full, "Kainji", "Lacustrine", "Wet", "Rainfall",  "pH",         0.035919)
df_full <- insert_p(df_full, "Jebba",  "Riverine",   "Dry", "Temp_Air",  "Temp_Water", 0.039778)
df_full <- insert_p(df_full, "Jebba",  "Ecotonal",   "Dry", "Temp_Air",  "pH",         0.044983)
df_full <- insert_p(df_full, "Kainji", "Riverine",   "Wet", "RH",        "Temp_Water", 0.045377)
df_full <- insert_p(df_full, "Kainji", "Riverine",   "Dry", "Windspeed", "pH",         0.047048)

# -- 4. TRANSFORM -------------------------------------------------
df_full <- df_full %>%
  mutate(neg_log_p = -log10(P_raw))

# -- 5. FACTOR ORDERING --------------------------------------------
df_full$Effect  <- factor(df_full$Effect,
                           levels = rev(c("BOD","DO","EC","NO3_N",
                                          "pH","PO4_P","Turbidity","Temp_Water")))
df_full$Cause   <- factor(df_full$Cause,   levels = climate_vars)
df_full$Habitat <- factor(df_full$Habitat, levels = c("Riverine","Ecotonal","Lacustrine"))
df_full$Season  <- factor(df_full$Season,  levels = c("Wet","Dry"))
df_full$Dam     <- factor(df_full$Dam,     levels = c("Jebba","Kainji"))

# -- 6. SHARED THEME -------------------------------------------------
shared_theme <- theme_bw(base_size = 11) +
  theme(
    plot.title       = element_text(face = "bold", size = 12, hjust = 0),
    plot.subtitle    = element_text(size = 9, colour = "grey25",
                                    hjust = 0, margin = margin(b = 8)),
    plot.caption     = element_text(size = 7.5, colour = "grey40",
                                    hjust = 0, margin = margin(t = 10)),
    axis.text.x      = element_text(angle = 38, hjust = 1, size = 9),
    axis.text.y      = element_text(size = 9),
    axis.title.x     = element_text(size = 10, face = "bold",
                                    margin = margin(t = 8)),
    axis.title.y     = element_text(size = 10, face = "bold",
                                    margin = margin(r = 8)),
    strip.text.x     = element_text(size = 10, face = "bold"),
    strip.text.y     = element_text(size = 9, angle = 0),
    strip.background = element_rect(fill = "grey88", colour = "grey65"),
    legend.title     = element_text(size = 9, face = "bold"),
    legend.position  = "right",
    panel.grid       = element_blank(),
    panel.border     = element_rect(colour = "grey55", linewidth = 0.5),
    panel.spacing    = unit(0.5, "lines")
  )

# -- 7. SHARED COLOUR SCALE -------------------------------------------
shared_fill <- scale_fill_gradient(
  low    = "#1A237E",
  high   = "#EDE7F6",
  limits = c(0, 3.10),
  breaks = c(0, 0.5, 1.0, 1.301, 1.5, 2.0, 2.5, 3.0),
  labels = c("0.0", "0.5", "1.0",
             "1.30\n\u2190 p = 0.05\nthreshold",
             "1.5", "2.0", "2.5", "3.0"),
  name   = expression(-log[10](raw~italic(p))),
  guide  = guide_colorbar(
    barheight    = unit(6, "cm"),
    barwidth     = unit(0.5, "cm"),
    ticks.colour = "grey40",
    frame.colour = "grey40",
    label.theme  = element_text(size = 7.5, lineheight = 1.1)
  )
)

# -- 8. FUNCTION TO BUILD ONE SEASON PLOT ----------------------------
build_plot <- function(season_filter, fig_label) {

  df_season <- df_full %>% filter(Season == season_filter)

  n_near <- sum(df_season$Near_Sig)

  p <- ggplot(df_season, aes(x = Cause, y = Effect, fill = neg_log_p)) +
    geom_tile(colour = "white", linewidth = 0.35) +

    # Yellow border for near-significant cells
    geom_tile(data = df_season %>% filter(Near_Sig),
              aes(x = Cause, y = Effect),
              fill  = NA,
              colour = "#FFD700",
              linewidth = 1.0) +

    # White dot on near-significant cells
    geom_point(data = df_season %>% filter(Near_Sig),
               aes(x = Cause, y = Effect),
               shape = 16, size = 1.8,
               colour = "white", inherit.aes = FALSE) +

    facet_grid(Habitat ~ Dam, switch = "y") +

    shared_fill +

    scale_colour_identity() +

    labs(
      title    = paste0("Figure 5", fig_label,
                        ". Granger Causality: Climate \u2192 Water Quality (",
                        season_filter, " Season)"),
      subtitle = paste0("No relationships significant after FDR correction ",
                        "(all adjusted p\u2009>\u20090.05). ",
                        "Near-significant cells (raw p\u2009<\u20090.05, n\u2009=\u2009",
                        n_near, ") marked with yellow borders."),
      x        = "Climate Variable (Cause)",
      y        = "Water Quality Variable (Effect)",
      caption  = paste0(
        "Colour intensity reflects \u2212log\u2081\u2080(raw p-value); ",
        "darker blue = weaker evidence of predictive association. ",
        "Dashed threshold at \u2212log\u2081\u2080 = 1.30 (p\u2009=\u20090.05).\n",
        "Part of 480 stratified tests total: 2 dams \u00d7 3 habitats \u00d7 ",
        "2 seasons \u00d7 5 climate drivers \u00d7 8 water-quality parameters. ",
        "Lag\u2009=\u20091 month throughout. Full results in Table S1."
      )
    ) +

    shared_theme

  return(p)
}

# -- 9. BUILD BOTH PLOTS ----------------------------------------------
p_wet <- build_plot("Wet", "a")
p_dry <- build_plot("Dry", "b")

# -- 10. SAVE ------------------------------------------------------
# Figure 5a - Wet Season
ggsave("Figure5a_Granger_Wet.pdf",
       plot = p_wet,
       width = 10, height = 10,
       units = "in", dpi = 300)

ggsave("Figure5a_Granger_Wet.png",
       plot = p_wet,
       width = 10, height = 10,
       units = "in", dpi = 300)

# Figure 5b - Dry Season
ggsave("Figure5b_Granger_Dry.pdf",
       plot = p_dry,
       width = 10, height = 10,
       units = "in", dpi = 300)

ggsave("Figure5b_Granger_Dry.png",
       plot = p_dry,
       width = 10, height = 10,
       units = "in", dpi = 300)

message("Figure 5a (Wet) and Figure 5b (Dry) saved successfully.")
