# kainji-jebba-water-quality-temporal-dynamics
R scripts for variance partitioning, Granger causality, and changepoint detection of water-quality dynamics in the Kainji–Jebba reservoir cascade, Nigeria
# Temporal Dynamics Outweigh Spatial Gradients in Afrotropical 
# Reservoir Water Quality — R Scripts

## Reference
Abdullahi, I., Sadiku, S.O.E., Robertson-Andersson, D.V. and Okpeku, M. 
(2026). Temporal dynamics outweigh spatial gradients in Afrotropical 
reservoir water quality: evidence from the Kainji–Jebba cascade, Nigeria. 
Environmental Monitoring and Assessment. [DOI to be inserted upon publication]

## Contents
This repository contains all R scripts used for data analysis in the 
above manuscript.

| Script | Purpose |
|--------|---------|
| 01_data_cleaning.R | Raw data import, quality checks, and preparation |
| 02_WQI_computation.R | NSF Water Quality Index calculation |
| 03_variance_partitioning.R | Partial RDA and variance partitioning |
| 04_granger_causality.R | Granger causality testing for hydroclimatic drivers |
| 05_changepoint_detection.R | PELT, binary segmentation, CUSUM, and Bayesian changepoint detection |

## Data
The underlying water-quality dataset is available from the corresponding 
author upon reasonable request.

## R version
All scripts were developed and tested in R version 4.6.0.

## Contact
Ibrahim Abdullahi
a.ibrahim@futminna.edu.ng
