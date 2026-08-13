# Flood Impact on Healthcare Facility Accessibility Analysis
# Analyzes which healthcare facilities become inaccessible via road network
# at 5km, 10km, and 15km distances when roads are flooded
getwd()
#setwd("..")


# Load required libraries
library(sf)
library(terra)
library(sfnetworks)
library(tidygraph)
library(dplyr)
library(tidyr)

# ==============================================================================
# 1. LOAD DATA
# ==============================================================================

# IMPORTANT: Set which flood scenario you're analyzing
FLOOD_SCENARIO <- "punjab"  # Options: "punjab" or "central_pak"

# Load appropriate flood extent raster
if(FLOOD_SCENARIO == "punjab") {
  flood_raster <- rast("noun/punjab_flood.tif")
  cat("Analyzing PUNJAB flood scenario\n")
} else if(FLOOD_SCENARIO == "central_pak") {
  flood_raster <- rast("noun/Sindh_new_flood.tif")
  cat("Analyzing CENTRAL PAKISTAN flood scenario\n")
}

# Get flood extent bbox with 50km buffer for cropping
flood_bbox <- st_bbox(flood_raster) %>% 
  st_as_sfc() %>%
  st_buffer(50000)  # 50km buffer to capture network effects

# Load roads - crop to flood area + buffer
roads_full <- st_read("noun/roads_all.shp")
roads <- roads_full %>% 
  st_crop(flood_bbox)
cat("Roads: Kept", nrow(roads), "of", nrow(roads_full), "segments in analysis area\n")

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

# 6. OPTIONAL: Merge connected segments to reduce complexity
# With ~28k segments, you probably DON'T need this unless it's slow
# Only uncomment if network analysis takes > 30 minutes
# cat("Merging connected segments by highway type...\n")
# roads <- roads %>%
#   group_by(highway) %>%
#   summarise(do_union = TRUE) %>%
#   st_cast("LINESTRING")
# cat("After merging:", nrow(roads), "segments\n")

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





# Load healthcare facilities - crop to flood area
healthcare_full <- st_read("noun/healthcare_f.shp")
healthcare <- healthcare_full %>%
  st_crop(flood_bbox)
cat("Healthcare: Kept", nrow(healthcare), "of", nrow(healthcare_full), "facilities in analysis area\n")

# Load districts - filter to affected provinces
districts_full <- st_read("noun/district.shp")

if(FLOOD_SCENARIO == "punjab") {
  # For Punjab floods, use Punjab district HQs as origins
  districts <- districts_full %>%
    filter(PROVINCE == "PUNJAB")  # Adjust column name as needed
} else if(FLOOD_SCENARIO == "central_pak") {
  # For Central Pak floods, use Punjab + Sindh + Balochistan districts
  districts <- districts_full %>%
    filter(PROVINCE %in% c("PUNJAB", "SINDH", "Balochistan"))
}
cat("Districts: Using", nrow(districts), "districts as origin points\n")

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


if (!dir.exists("output")) dir.create("output", recursive = TRUE)

st_write(roads,
         "output/usable_roads_40.gpkg",
         layer = "usable_roads_40",
         delete_dsn = TRUE)



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
# 5. CALCULATE NETWORK DISTANCES TO ALL HEALTHCARE FACILITIES
# ==============================================================================

# Snap healthcare facilities to nearest road network node
#healthcare$nearest_node <- st_nearest_feature(
 # healthcare,
#  road_network %>% activate("nodes") %>% st_as_sf()
#)

# Snap origin cities to nearest road network nodes
#origin_cities$nearest_node <- st_nearest_feature(
 # origin_cities,
#  road_network %>% activate("nodes") %>% st_as_sf()
#)

# Calculate shortest path distances from NEAREST origin to each facility
# This represents realistic access patterns (people go to nearest city)
#cat("\nCalculating network distances from nearest origin city...\n")

#healthcare$network_distance_m <- NA
#healthcare$nearest_origin <- NA

#for(i in 1:nrow(healthcare)) {
 # if(i %% 100 == 0) cat("  Processing facility", i, "of", nrow(healthcare), "\n")
  
  #facility_node <- healthcare$nearest_node[i]
  #min_distance <- Inf
  #nearest_city <- NA
  
  # Find minimum distance to any origin city
#  for(j in 1:nrow(origin_cities)) {
 #   tryCatch({
  #    path_dist <- road_network %>%
   #     activate("nodes") %>%
    #    mutate(dist = node_distance_from(origin_cities$nearest_node[j], 
     #                                    weights = length_m)) %>%
      #  st_as_sf() %>%
       # slice(facility_node) %>%
        #pull(dist)
      
    #  if(is.finite(path_dist) && path_dist < min_distance) {
     #   min_distance <- path_dist
      #  nearest_city <- origin_cities$city[j]
  #    }
    #}, error = function(e) {
    #  # Path doesn't exist from this origin
    #})
 # }
  
  #healthcare$network_distance_m[i] <- min_distance
  #healthcare$nearest_origin[i] <- nearest_city
#}

#cat("Distance calculation complete!\n")


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

st_write(healthcare,  # This has geometry + accessibility attributes
         paste0("healthcare_facilities_", FLOOD_SCENARIO, ".gpkg"),
         delete_dsn = TRUE)

# Export detailed facility-level data
healthcare %>%
  st_drop_geometry() %>%
  select(!!sym(category_col),
         nearest_origin,
         network_distance_km,
         accessible_5km, accessible_10km, accessible_15km,
         access_status_5km, access_status_10km, access_status_15km) %>%
  write.csv(detailed_filename, row.names = FALSE)

# Export summary by category
write.csv(final_summary, output_filename, row.names = FALSE)












