# ============================================================================
# 01_LAND_USE_EXTRACTION.R
# LAND USE EXTRACTION FROM ESA WORLDCOVER 2021 USING R
# ============================================================================
# This script extracts agriculture, urban, and forest cover within 5km buffers
# around sampling sites using Google Earth Engine via rgee
#
# Required local inputs (not included in this repository; see Data availability
# in README.md):
#   - sampling_site.gpkg
#   - buffers_5km.gpkg
# ============================================================================

# Required packages
required_packages <- c("sf", "rgee", "tidyverse", "spdep", "ape", "viridis")

# Install missing packages
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

# Load libraries
library(sf)
library(rgee)
library(tidyverse)
library(spdep)
library(ape)
library(viridis)

# ============================================================================
# STEP 1: INITIALIZE GOOGLE EARTH ENGINE
# ============================================================================
# One-time setup (run once, then comment out):
# ee_install()
# ee_Initialize(user = 'your_email@example.com', drive = TRUE)

# Regular initialization:
ee_Initialize()

# ============================================================================
# STEP 2: LOAD SPATIAL DATA
# ============================================================================
cat("Loading spatial data...\n")

# Load sampling sites
sites <- st_read("sampling_site.gpkg")
cat("Loaded", nrow(sites), "sampling sites\n")

# Load 5km buffers
buffers <- st_read("buffers_5km.gpkg")
cat("Loaded", nrow(buffers), "5km buffers\n")

# Ensure CRS is WGS84 (EPSG:4326) for GEE
sites <- st_transform(sites, 4326)
buffers <- st_transform(buffers, 4326)

# ============================================================================
# STEP 3: DEFINE ESA WORLDCOVER LAND USE CLASSES
# ============================================================================
# ESA WorldCover 2021 classification:
# 10 = Tree cover (Forest)
# 20 = Shrubland
# 30 = Grassland
# 40 = Cropland (Agriculture)
# 50 = Built-up (Urban)
# 60 = Bare/sparse vegetation
# 70 = Snow and ice
# 80 = Permanent water bodies
# 90 = Herbaceous wetland
# 95 = Mangroves
# 100 = Moss and lichen

landuse_classes <- list(
  forest = 10,
  agriculture = 40,
  urban = 50
)

# ============================================================================
# STEP 4: EXTRACT LAND USE DATA FROM GOOGLE EARTH ENGINE
# ============================================================================
cat("\nExtracting land use data from Google Earth Engine...\n")

# Load ESA WorldCover 2021
esa_worldcover <- ee$ImageCollection("ESA/WorldCover/v200")$
  first()$
  select("Map")

# Function to calculate land use percentages for a single buffer
extract_landuse_buffer <- function(buffer_geom, buffer_id) {

  # Convert sf geometry to ee.Geometry
  ee_geom <- sf_as_ee(buffer_geom)

  # Calculate area of buffer (in hectares)
  buffer_area_ha <- ee_geom$area()$divide(10000)

  # Extract land use within buffer
  landuse_stats <- list()

  for (class_name in names(landuse_classes)) {
    class_value <- landuse_classes[[class_name]]

    # Create binary mask for this land use class
    class_mask <- esa_worldcover$eq(class_value)

    # Calculate area of this class within buffer
    class_area <- class_mask$
      multiply(ee$Image$pixelArea())$
      divide(10000)$  # Convert to hectares
      reduceRegion(
        reducer = ee$Reducer$sum(),
        geometry = ee_geom,
        scale = 10,  # ESA WorldCover is 10m resolution
        maxPixels = 1e13
      )

    # Get the value
    area_value <- class_area$get("Map")$getInfo()

    # Calculate percentage
    pct <- (area_value / buffer_area_ha$getInfo()) * 100

    landuse_stats[[paste0(class_name, "_ha")]] <- area_value
    landuse_stats[[paste0(class_name, "_pct")]] <- pct
  }

  return(as.data.frame(landuse_stats))
}

# Extract land use for all buffers
landuse_results <- list()

for (i in 1:nrow(buffers)) {
  cat("Processing buffer", i, "of", nrow(buffers), "\r")

  buffer_geom <- buffers[i, ]
  buffer_id <- if("id" %in% names(buffers)) buffers$id[i] else i

  tryCatch({
    stats <- extract_landuse_buffer(buffer_geom, buffer_id)
    stats$buffer_id <- buffer_id
    landuse_results[[i]] <- stats
  }, error = function(e) {
    cat("\nError processing buffer", i, ":", e$message, "\n")
    return(NULL)
  })
}

# Combine results
landuse_df <- bind_rows(landuse_results)

cat("\nLand use extraction complete!\n")

# ============================================================================
# STEP 5: JOIN RESULTS WITH SAMPLING SITES
# ============================================================================
# Join land use data with site information
sites_with_landuse <- sites %>%
  mutate(buffer_id = if("id" %in% names(.)) id else row_number()) %>%
  left_join(landuse_df, by = "buffer_id")

# Display summary
cat("\n=== LAND USE SUMMARY ===\n")
summary(sites_with_landuse %>%
          st_drop_geometry() %>%
          select(ends_with("_pct")))

# ============================================================================
# STEP 6: SPATIAL AUTOCORRELATION ANALYSIS (MORAN'S I)
# ============================================================================
cat("\n=== SPATIAL AUTOCORRELATION ANALYSIS ===\n")

# Create spatial weights matrix (k-nearest neighbors, k=4)
coords <- st_coordinates(st_centroid(sites_with_landuse))
knn_weights <- knearneigh(coords, k = 4)
nb <- knn2nb(knn_weights)
listw <- nb2listw(nb, style = "W")

# Calculate Moran's I for each land use type
morans_results <- data.frame(
  land_use = character(),
  morans_i = numeric(),
  p_value = numeric(),
  interpretation = character(),
  stringsAsFactors = FALSE
)

for (class_name in names(landuse_classes)) {
  pct_col <- paste0(class_name, "_pct")

  if (pct_col %in% names(sites_with_landuse)) {
    moran_test <- moran.test(
      sites_with_landuse[[pct_col]],
      listw,
      na.action = na.omit
    )

    interp <- ifelse(
      moran_test$p.value < 0.05,
      ifelse(moran_test$estimate[1] > 0,
             "Significant clustering",
             "Significant dispersion"),
      "No significant spatial pattern"
    )

    morans_results <- rbind(morans_results, data.frame(
      land_use = class_name,
      morans_i = round(moran_test$estimate[1], 4),
      p_value = round(moran_test$p.value, 4),
      interpretation = interp
    ))
  }
}

print(morans_results)

# Visualize Moran's I
ggplot(morans_results, aes(x = land_use, y = morans_i, fill = p_value < 0.05)) +
  geom_col() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = c("gray70", "coral2"),
                    labels = c("Not significant", "Significant (p < 0.05)")) +
  labs(title = "Spatial Autocorrelation of Land Use Variables",
       subtitle = "Moran's I statistic",
       x = "Land Use Type",
       y = "Moran's I",
       fill = "Significance") +
  theme_minimal()

# ============================================================================
# STEP 7: VISUALIZE RESULTS
# ============================================================================
# Create maps for each land use type
for (class_name in names(landuse_classes)) {
  pct_col <- paste0(class_name, "_pct")

  if (pct_col %in% names(sites_with_landuse)) {
    p <- ggplot() +
      geom_sf(data = buffers, fill = NA, color = "gray80", linewidth = 0.3) +
      geom_sf(data = sites_with_landuse,
              aes(color = .data[[pct_col]]),
              size = 3) +
      scale_color_viridis(option = "viridis", name = "% Cover") +
      labs(title = paste(str_to_title(class_name), "Cover within 5km Buffers"),
           subtitle = paste0("Mean: ",
                           round(mean(sites_with_landuse[[pct_col]], na.rm = TRUE), 2),
                           "%, Range: ",
                           round(min(sites_with_landuse[[pct_col]], na.rm = TRUE), 2),
                           "% - ",
                           round(max(sites_with_landuse[[pct_col]], na.rm = TRUE), 2),
                           "%")) +
      theme_minimal() +
      theme(legend.position = "right")

    print(p)

    # Save plot
    ggsave(paste0("map_", class_name, "_cover.png"),
           plot = p, width = 10, height = 8, dpi = 300)
  }
}

# ============================================================================
# STEP 8: EXPORT RESULTS
# ============================================================================
cat("\nExporting results...\n")

# Export as CSV
write_csv(
  sites_with_landuse %>% st_drop_geometry(),
  "landuse_extraction_results.csv"
)

# Export as shapefile (with spatial data)
st_write(
  sites_with_landuse,
  "landuse_extraction_results.gpkg",
  delete_dsn = TRUE
)

# Export Moran's I results
write_csv(morans_results, "spatial_autocorrelation_results.csv")

cat("\n=== ANALYSIS COMPLETE ===\n")
cat("Results saved to:\n")
cat("  - landuse_extraction_results.csv\n")
cat("  - landuse_extraction_results.gpkg\n")
cat("  - spatial_autocorrelation_results.csv\n")
cat("  - map_*_cover.png (visualization maps)\n")

# ============================================================================
# ALTERNATIVE: FASTER BATCH EXTRACTION (FOR MANY SITES)
# ============================================================================
# If you have many sites and GEE requests are slow, use this batch approach:

extract_landuse_batch <- function(buffers_sf) {

  # Convert all buffers to ee.FeatureCollection
  fc <- sf_as_ee(buffers_sf)

  # Function to extract stats for each feature
  extract_stats <- function(feature) {
    geom <- feature$geometry()

    # Calculate percentages for each class
    stats <- ee$Dictionary$fromLists(
      keys = ee$List(c("forest", "agriculture", "urban")),
      values = ee$List(c(10, 40, 50))$map(ee_utils_pyfunc(function(class_val) {
        mask <- esa_worldcover$eq(class_val)
        area <- mask$multiply(ee$Image$pixelArea())$
          reduceRegion(
            reducer = ee$Reducer$sum(),
            geometry = geom,
            scale = 10,
            maxPixels = 1e13
          )$get("Map")
        return(area)
      }))
    )

    return(feature$set(stats))
  }

  # Apply extraction to all features
  results_fc <- fc$map(extract_stats)

  # Convert back to data frame
  results_df <- ee_as_sf(results_fc) %>% st_drop_geometry()

  return(results_df)
}

# Uncomment to use batch extraction:
# landuse_batch <- extract_landuse_batch(buffers)
