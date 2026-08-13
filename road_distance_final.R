
# Load required libraries
library(sf)
library(terra)
library(sfnetworks)
library(tidygraph)
library(dplyr)
library(tidyr)

# Paths (adjust if different)
punjab_path  <- "noun/punjab_flood.tif"
central_path <- "noun/Sindh_new_flood.tif"

cat("Loading flood rasters...\n")
flood_punjab  <- rast(punjab_path)
flood_central <- rast(central_path)

# Combine them: flooded if either raster is flooded (pixel-wise max)
# Works for binary (0/1) or continuous -> you can adjust threshold later
cat("Combining flood rasters (pixel-wise max)...\n")
flood_stack <- c(flood_punjab, flood_central)
flood_combined <- app(flood_stack, fun = max, na.rm = TRUE)

# assign the raster used downstream
flood_raster <- flood_combined
cat("Combined flood raster ready\n")

# Optional: write combined raster for inspection
writeRaster(flood_raster, "noun/combined_flood.tif", overwrite=TRUE)

# -------------------------------------------------------------------------
# Build union bbox from both rasters (so we crop roads/health to entire affected area)
# -------------------------------------------------------------------------
# Note: use st_as_sfc(st_bbox(...)) to convert raster extents to simple features
bbox1 <- st_as_sfc(st_bbox(flood_punjab))
bbox2 <- st_as_sfc(st_bbox(flood_central))

# union the two extents and buffer to capture network effects
flood_union <- st_union(bbox1, bbox2)
# Add 50 km buffer (same as your original)
flood_bbox <- st_buffer(flood_union, 50000)

# -------------------------------------------------------------------------
# Load roads - crop to flood union bbox (whole-country roads cropped to union)
# -------------------------------------------------------------------------
roads_full <- st_read("noun/roads_all.shp")
roads <- roads_full %>% 
  st_crop(flood_bbox)
cat("Roads: Kept", nrow(roads), "of", nrow(roads_full), "segments in analysis area\n")

# -------------------------------------------------------------------------
# Load healthcare and districts, cropping/selecting by union bbox
# -------------------------------------------------------------------------
healthcare_full <- st_read("noun/healthcare_f.shp")
healthcare <- healthcare_full %>%
  st_crop(flood_bbox)
cat("Healthcare: Kept", nrow(healthcare), "of", nrow(healthcare_full), "facilities in analysis area\n")

districts_full <- st_read("noun/district.shp")

# Instead of selecting provinces manually, pick districts that intersect the union bbox
districts <- districts_full[ st_intersects(districts_full, flood_bbox, sparse = FALSE), ]
cat("Districts: Using", nrow(districts), "districts intersecting flood union area\n")



# ==============================================================================
# ROAD NETWORK PREPROCESSING (IMPORTANT!)
# ==============================================================================

valid_road_types <- c("primary", "secondary", "trunk", "motorway", "highway",
                      "primary_link", "secondary_link", "trunk_link", "motorway_link")  # include link roads

roads <- roads %>% filter(highway %in% valid_road_types)
cat("Filtered to valid road types:", nrow(roads), "segments remaining\n")
cat("Road type distribution:\n")
print(table(roads$highway))

# 2. Remove duplicate geometries (common in OSM data)
roads <- roads %>% distinct(geometry, .keep_all = TRUE)
cat("After removing duplicate geometries:", nrow(roads), "segments\n")

# 3. Fix invalid geometries
invalid_count <- sum(!st_is_valid(roads))
if(invalid_count > 0) {
  cat("Fixing", invalid_count, "invalid geometries...\n")
  roads <- st_make_valid(roads)
}

# 4. CRITICAL: Convert MULTILINESTRING to LINESTRING
# Your data has MultiLineString geometry which needs conversion
cat("Converting MultiLineString to LineString...\n")
original_count <- nrow(roads)
roads <- st_cast(roads, "LINESTRING")
cat("After conversion:", nrow(roads), "segments (was", original_count, ")\n")

# 5. Remove very short segments (< 5m) that can cause network issues
roads$length_m <- st_length(roads) %>% as.numeric()
short_segments <- sum(roads$length_m < 5)
if(short_segments > 0) {
  cat("Removing", short_segments, "very short segments (< 5m)...\n")
  roads <- roads %>% filter(length_m >= 5)
}

# 7. Add unique road ID and recalculate length after casting
roads$road_id <- 1:nrow(roads)
roads$length_m <- st_length(roads) %>% as.numeric()  # Recalculate after casting

cat("Final road network:", nrow(roads), "segments ready for analysis\n")
cat("Total road length:", round(sum(roads$length_m)/1000, 1), "km\n")


st_write(
  roads,
  "maps_f/roads_clean_nodup_lines.shp",
  delete_dsn = TRUE
)



# Verify CRS consistency
cat("\nVerifying coordinate systems...\n")
cat("Flood CRS:", crs(flood_raster, describe=TRUE)$name, "\n")
cat("Healthcare CRS:", st_crs(healthcare)$input, "\n")
cat("Roads CRS:", st_crs(roads)$input, "\n")
cat("Districts CRS:", st_crs(districts)$input, "\n")

# ==============================================================================
# 2. IDENTIFY FLOODED ROAD SEGMENTS (FAST METHOD - NO POLYGON CONVERSION)
# ==============================================================================

# Calculate road segment lengths
roads$length_m <- st_length(roads) %>% as.numeric()  # Already calculated in preprocessing

# METHOD 1: Sample points along each road and check flood raster
# This avoids the expensive raster-to-polygon conversion
cat("\nSampling roads to determine flood status...\n")

# Create sampling points along each road (every 20m to match raster resolution)
sampling_distance <- 20  # meters - adjust based on your raster resolution (20m for 20-pixel raster)

# Function to sample a single road segment
sample_road_flood <- function(road_geom, road_len, raster, sample_dist) {
  # Number of sample points
  n_samples <- ceiling(road_len / sample_dist)
  
  # Create points along the line at regular intervals
  sample_points <- st_line_sample(road_geom, n = n_samples)
  sample_points <- st_cast(sample_points, "POINT")
  
  # Extract raster values at these points
  flood_values <- terra::extract(raster, vect(sample_points), ID = FALSE)
  
  # Calculate percentage of points in flooded areas
  # Assuming flood raster: 1 = flooded, 0/NA = not flooded
  n_flooded <- sum(flood_values > 0, na.rm = TRUE)
  flood_pct <- (n_flooded / n_samples) * 100
  
  return(flood_pct)
}


# Apply to all road segments (with progress tracking)
# For faster processing with many roads, use parallel processing
cat("Processing", nrow(roads), "road segments...\n")

# OPTION A: Sequential (safer, easier to debug)
roads$flood_pct <- sapply(1:nrow(roads), function(i) {
  if(i %% 1000 == 0) cat("  Processed", i, "of", nrow(roads), "\n")
  sample_road_flood(roads$geometry[i], roads$length_m[i], 
                    flood_raster, sampling_distance)
})


cat("Sampling complete!\n")

# Mark roads as unusable if >= 70% flooded (or adjust threshold)
FLOOD_THRESHOLD <- 40  # percentage
roads$usable <- roads$flood_pct < FLOOD_THRESHOLD

cat("\nRoad Network Summary:\n")
cat("Total road segments:", nrow(roads), "\n")
cat("Unusable (>=", FLOOD_THRESHOLD, "% flooded):", sum(!roads$usable), "\n")
cat("Usable (<", FLOOD_THRESHOLD, "% flooded):", sum(roads$usable), "\n")


st_write(
  roads,
  "output/usable_roads_40.shp",
  delete_dsn = TRUE
)


write.csv(roads, 
          "output/usable_roads_40.csv", 
          row.names = FALSE)

library(sf)

# ==============================================================================
# 3. CREATE USABLE ROAD NETWORK
# ==============================================================================

# Filter to usable roads only
usable_roads <- roads %>% filter(usable)

# Create road network graph
road_network <- as_sfnetwork(usable_roads, directed = FALSE)

# ==============================================================================
# 4. DEFINE REFERENCE POINT(S) FOR DISTANCE CALCULATION
# ==============================================================================

# RECOMMENDED FOR PUNJAB PROVINCE: Use district headquarters as origins
# Adjust based on flood scenario

cat("\nSetting up origin points from district headquarters...\n")

# Use district centroids or actual HQ locations as origins
# Option A: Use centroids (simple)
origin_cities <- districts %>%
  st_centroid() %>%
  mutate(city = DISTRICT)  # Adjust column name as needed

# Option B: If you have actual district HQ coordinates, use those instead
# origin_cities <- districts %>%
#   # Add HQ point geometry here
#   mutate(city = district_name)

cat("Using", nrow(origin_cities), "district headquarters as origin points\n")


# ==============================================================================
# 5. CALCULATE NETWORK DISTANCES (OPTIMIZED VERSION - MUCH FASTER)
# ==============================================================================

cat("\nCalculating network distances (optimized method)...\n")

# Snap healthcare facilities to nearest road network node
healthcare$nearest_node <- st_nearest_feature(
  healthcare,
  road_network %>% activate("nodes") %>% st_as_sf()
)

# Snap origin cities to nearest road network nodes
origin_cities$nearest_node <- st_nearest_feature(
  origin_cities,
  road_network %>% activate("nodes") %>% st_as_sf()
)

# Pre-calculate distances from ALL origins to ALL nodes (one-time calculation)
cat("Pre-calculating distance matrix from all origins...\n")
all_nodes <- road_network %>% activate("nodes") %>% st_as_sf()

# Create distance matrix for each origin
distance_matrices <- list()

for(j in 1:nrow(origin_cities)) {
  cat("  Calculating distances from origin", j, "of", nrow(origin_cities), 
      "(", origin_cities$city[j], ")\n")
  
  distance_matrices[[j]] <- road_network %>%
    activate("nodes") %>%
    mutate(dist = node_distance_from(origin_cities$nearest_node[j], 
                                     weights = length_m)) %>%
    pull(dist)
}

# Now find minimum distance for each facility (FAST lookup)
cat("Finding nearest accessible origin for each facility...\n")

healthcare$network_distance_m <- sapply(1:nrow(healthcare), function(i) {
  if(i %% 1000 == 0) cat("  Processed", i, "of", nrow(healthcare), "\n")
  
  facility_node <- healthcare$nearest_node[i]
  
  # Get distances from this node to all origins
  distances_to_origins <- sapply(distance_matrices, function(dm) dm[facility_node])
  
  # Return minimum finite distance
  min_dist <- min(distances_to_origins[is.finite(distances_to_origins)])
  return(min_dist)
})

healthcare$nearest_origin <- sapply(1:nrow(healthcare), function(i) {
  facility_node <- healthcare$nearest_node[i]
  distances_to_origins <- sapply(distance_matrices, function(dm) dm[facility_node])
  
  if(all(!is.finite(distances_to_origins))) {
    return(NA)
  }
  
  origin_cities$city[which.min(distances_to_origins)]
})

cat("Distance calculation complete!\n")



# Convert to km
healthcare$network_distance_km <- healthcare$network_distance_m / 1000

# ==============================================================================
# 6. CLASSIFY ACCESSIBILITY AT 5, 10, 15 KM THRESHOLDS
# ==============================================================================

# Determine accessibility status for each threshold
healthcare <- healthcare %>%
  mutate(
    accessible_5km = network_distance_km <= 5 & is.finite(network_distance_km),
    accessible_10km = network_distance_km <= 10 & is.finite(network_distance_km),
    accessible_15km = network_distance_km <= 15 & is.finite(network_distance_km),
    # Categorize accessibility
    access_status_5km = case_when(
      !is.finite(network_distance_km) ~ "Disconnected",
      accessible_5km ~ "Accessible",
      TRUE ~ "Inaccessible"
    ),
    access_status_10km = case_when(
      !is.finite(network_distance_km) ~ "Disconnected",
      accessible_10km ~ "Accessible",
      TRUE ~ "Inaccessible"
    ),
    access_status_15km = case_when(
      !is.finite(network_distance_km) ~ "Disconnected",
      accessible_15km ~ "Accessible",
      TRUE ~ "Inaccessible"
    )
  )

# ==============================================================================
# 7. CREATE SUMMARY TABLE BY CATEGORY
# ==============================================================================


category_col <- "Category"  

# Create summary by category
summary_5km <- healthcare %>%
  st_drop_geometry() %>%
  group_by(!!sym(category_col), access_status_5km) %>%
  summarise(count = n(), .groups = "drop") %>%
  pivot_wider(names_from = access_status_5km, 
              values_from = count, 
              values_fill = 0) %>%
  mutate(distance_threshold = "5km")

summary_10km <- healthcare %>%
  st_drop_geometry() %>%
  group_by(!!sym(category_col), access_status_10km) %>%
  summarise(count = n(), .groups = "drop") %>%
  pivot_wider(names_from = access_status_10km, 
              values_from = count, 
              values_fill = 0) %>%
  mutate(distance_threshold = "10km")

summary_15km <- healthcare %>%
  st_drop_geometry() %>%
  group_by(!!sym(category_col), access_status_15km) %>%
  summarise(count = n(), .groups = "drop") %>%
  pivot_wider(names_from = access_status_15km, 
              values_from = count, 
              values_fill = 0) %>%
  mutate(distance_threshold = "15km")

# Combine all summaries
final_summary <- bind_rows(summary_5km, summary_10km, summary_15km) %>%
  select(distance_threshold, !!sym(category_col), everything())

# ==============================================================================
# 8. EXPORT RESULTS
# ==============================================================================
# Create output filenames
detailed_filename <- paste0("healthcare_accessibility_detailed_", 
                            FLOOD_SCENARIO, ".csv")
output_filename <- paste0("healthcare_accessibility_summary_", 
                          FLOOD_SCENARIO, ".csv")

st_write(healthcare, "output/healthcare_facilities.gpkg")

healthcare %>%
  st_drop_geometry() %>%
  select(!!sym(category_col),
         nearest_origin,
         network_distance_km,
         accessible_5km, accessible_10km, accessible_15km,
         access_status_5km, access_status_10km, access_status_15km) %>%
  write.csv("output/healthcare_accessibility_detailed.csv",
            row.names = FALSE)

write.csv(final_summary,
          "output/healthcare_accessibility_summary.csv",
          row.names = FALSE)



