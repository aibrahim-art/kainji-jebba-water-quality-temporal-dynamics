# kainji-jebba-water-quality-temporal-dynamics

R scripts for variance partitioning, Granger causality, and changepoint detection of water-quality dynamics in the Kainji–Jebba reservoir cascade, Nigeria.

## Reference

Ibrahim, A., Sadiku, S. O. E., Robertson-Andersson, D. V., & Okpeku, M. (2026). Temporal dynamics outweigh spatial gradients in shaping water quality across a regulated Afrotropical reservoir cascade. *Environmental Monitoring and Assessment*, 198, 906. https://doi.org/10.1007/s10661-026-15681-8

This repository contains the R scripts used for data analysis in the above manuscript.

## Scripts

| Script | Purpose | Required local input |
|---|---|---|
| `01_land_use_extraction.R` | Land-cover extraction (forest, agriculture, urban cover) from ESA WorldCover 2021 within 5-km buffers around sampling sites, via Google Earth Engine (`rgee`); includes Moran's I spatial autocorrelation check | `sampling_site.gpkg`, `buffers_5km.gpkg` |
| `02_WQI_computation.R` | Computation of the modified NSF Water Quality Index (WQI) and summary statistics/plots by site, habitat, season, and month | `WQI_data.csv` |
| `03_variance_partitioning.R` | Partial redundancy analysis (RDA) and variance partitioning of water quality among spatial, temporal, and land-use predictor sets; PERMANOVA; NMDS ordination | `water_quality_data.csv` |
| `04a_granger_causality_full.R` | Full stratified Granger causality testing (480 tests: 5 climate drivers × 8 water-quality parameters × 12 strata), with stationarity checks, AIC-based lag selection, and Benjamini–Hochberg FDR correction | `Grouped_Dataset_Kainji_Jebba.xlsx` |
| `04b_granger_heatmap.R` | Generates the seasonally split Granger causality heatmaps (Fig. 5a–b) from the near-significant results of `04a_granger_causality_full.R` | *(values embedded from Table 3 of the manuscript)* |
| `05_changepoint_detection.R` | Multi-method changepoint detection (PELT, binary segmentation, CUSUM, Bayesian, variance-based, and moving-average methods) and consensus/sentinel-parameter analysis | `Grouped_Dataset_Kainji_Jebba.xlsx` |

All scripts were developed and tested in R version 4.5.1–4.5.3. None of the data files listed above are included in this repository — see **Data availability** below.

## Data availability

The water-quality, meteorological, and land-use datasets generated during this study are deposited in the Dryad Digital Repository: https://doi.org/10.5061/dryad.8pk0p2p4h

Scripts in this repository expect the corresponding data files (e.g., `water_quality_data.csv`, `Grouped_Dataset_Kainji_Jebba.xlsx`, `sampling_site.gpkg`, `buffers_5km.gpkg`) to be present locally and are not themselves bundled with the data — see the Dryad deposit above for the underlying dataset.

## Supplementary materials

Supplementary tables and figures referenced in the manuscript (Tables S1–S2, Figures S1a–S1b) are hosted by the journal and available via the article's online supplementary information: https://doi.org/10.1007/s10661-026-15681-8

## Requirements

Key R packages used across these scripts:

- `tidyverse`, `dplyr`, `tidyr`, `readxl`
- `vegan` (RDA, PERMANOVA, NMDS)
- `sf`, `rgee`, `spdep`, `ape` (spatial analysis / Earth Engine)
- `lmtest`, `tseries` (Granger causality, stationarity testing)
- `strucchange`, `changepoint`, `bcp`, `zoo` (changepoint detection)
- `ggplot2`, `viridis`, `gridExtra` (visualization)

## Contact

Ibrahim Abdullahi — a.ibrahim@futminna.edu.ng
