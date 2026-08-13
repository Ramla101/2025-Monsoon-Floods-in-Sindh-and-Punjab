library(sf)
library(dplyr)

# PARAMETERS
accessibility_distances <- c(5000, 10000, 15000)  # metres
isolation_threshold <- 0.8   # 80% flooded => isolated
major_classes <- c("Motorway", "Trunk", "Primary", "Secondary")  # your major classes

# Defensive checks (adjust field names if your healthcare table uses different column names)
if (!exists("healthcare") || nrow(healthcare) == 0) stop("healthcare missing or empty. Load the swat subset first.")
if (is.null(st_crs(healthcare))) stop("healthcare has no CRS. Set it first.")
if (st_crs(healthcare)$epsg != 32642) healthcare <- st_transform(healthcare, 32642)
if (is.null(st_crs(roads_master))) stop("roads_master missing CRS or not loaded.")
if (st_crs(roads_master)$epsg != 32642) roads_master <- st_transform(roads_master, 32642)

# Build flooded_road_ids if possible (fallback to spatial intersection if needed)
if (exists("roads_flood_with_attrs") && nrow(roads_flood_with_attrs) > 0 && "road_id" %in% names(roads_flood_with_attrs)) {
  flooded_road_ids <- unique(na.omit(roads_flood_with_attrs$road_id))
} else {
  if (!exists("flood_sf")) stop("No flood_sf available to compute flooded roads. Please load flood shapefile or provide roads_flood_with_attrs.")
  # subset roads_master to study area to speed up intersection
  study_area <- st_buffer(st_union(flood_sf), dist = 50000)
  roads_sub_all <- roads_master[lengths(st_intersects(roads_master, study_area, sparse = TRUE)) > 0, ]
  # ensure a road_id column exists for this temporary subset
  if (!"road_id" %in% names(roads_sub_all)) roads_sub_all$road_id <- seq_len(nrow(roads_sub_all))
  roads_buff <- st_buffer(roads_sub_all, dist = 10)
  ints <- st_intersects(roads_buff, flood_sf, sparse = TRUE)
  hit_idx <- which(lengths(ints) > 0)
  flooded_road_ids <- if (length(hit_idx) > 0) unique(roads_sub_all$road_id[hit_idx]) else integer(0)
  message("Fallback flooded_road_ids count: ", length(flooded_road_ids))
}

# OUTPUT containers
# overall summary (one row per radius)
healthcare_service_results <- data.frame(
  service_radius_km = numeric(),
  total_facilities = integer(),
  facilities_with_road_access = integer(),
  facilities_isolated = integer(),
  percent_isolated = numeric(),
  stringsAsFactors = FALSE
)

# per-category summary (one row per radius x category)
healthcare_by_category <- data.frame(
  service_radius_km = numeric(),
  category = character(),
  sub_cat = character(),
  gov_pvt = character(),
  total_facilities = integer(),
  facilities_with_road_access = integer(),
  facilities_isolated = integer(),
  percent_isolated = numeric(),
  stringsAsFactors = FALSE
)

# MAIN LOOP: compute isolation per facility, then summarize by category
for (d in accessibility_distances) {
  message("=== Processing service radius: ", d/1000, " km ===")
  # vectorized buffer of all healthcare points/features
  service_buffers <- st_buffer(healthcare, dist = d)
  # index roads intersecting each buffer
  roads_idx_per_fac <- st_intersects(service_buffers, roads_master, sparse = TRUE)
  n_fac <- length(roads_idx_per_fac)
  is_isolated <- logical(n_fac)
  has_any_road <- logical(n_fac)  # whether any road (major or not) exists in buffer
  
  # Progress bar
  pb <- txtProgressBar(min = 0, max = n_fac, style = 3)
  
  for (i in seq_len(n_fac)) {
    ix <- roads_idx_per_fac[[i]]
    if (length(ix) == 0) {
      # no roads in buffer -> isolated
      is_isolated[i] <- TRUE
      has_any_road[i] <- FALSE
      setTxtProgressBar(pb, i); next
    }
    has_any_road[i] <- TRUE
    subroads <- roads_master[ix, ]
    # only consider major roads
    major <- subroads %>% filter(road_class %in% major_classes)
    if (nrow(major) == 0) {
      # no major roads in buffer -> treat as isolated
      is_isolated[i] <- TRUE
      setTxtProgressBar(pb, i); next
    }
    total_major_len <- sum(as.numeric(st_length(major)))
    flooded_major_len <- 0
    
    # If we have a list of flooded road IDs and major has road_id column, use it
    if (length(flooded_road_ids) > 0 && "road_id" %in% names(major)) {
      flooded_idx_bool <- major$road_id %in% flooded_road_ids
      if (any(flooded_idx_bool)) flooded_major_len <- sum(as.numeric(st_length(major[flooded_idx_bool, ])))
    } else {
      # fallback: compute local intersection with flood_sf
      if (exists("flood_sf")) {
        flooded_local <- tryCatch({
          st_intersection(st_make_valid(major), st_make_valid(flood_sf))
        }, error = function(e) major[0, ])
        if (nrow(flooded_local) > 0) flooded_major_len <- sum(as.numeric(st_length(flooded_local)))
      }
    }
    
    # decide isolation
    if (total_major_len == 0) {
      is_isolated[i] <- TRUE
    } else {
      is_isolated[i] <- (flooded_major_len / total_major_len) > isolation_threshold
    }
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  # attach isolation results to healthcare attribute table (drop geometry for summarise)
  # Use your attribute names: here I assume 'Category', 'Sub_Cat', 'Gov_Pvt' exist
  attr_tbl <- healthcare %>%
    st_drop_geometry() %>%
    mutate(
      .facility_index = seq_len(n()),
      isolated = is_isolated,
      has_any_road = has_any_road
    )
  
  # overall counts for this radius
  total_fac <- nrow(attr_tbl)
  fac_iso <- sum(attr_tbl$isolated, na.rm = TRUE)
  fac_with_access <- sum(!attr_tbl$isolated, na.rm = TRUE)
  pct_iso <- round(100 * fac_iso / total_fac, 2)
  
  healthcare_service_results <- rbind(healthcare_service_results, data.frame(
    service_radius_km = d/1000,
    total_facilities = total_fac,
    facilities_with_road_access = fac_with_access,
    facilities_isolated = fac_iso,
    percent_isolated = pct_iso,
    stringsAsFactors = FALSE
  ))
  
  # ---- per-category summary ----
  # if your healthcare table uses different column names, replace 'Category', 'Sub_Cat', 'Gov_Pvt' accordingly
  per_cat <- attr_tbl %>%
    mutate(
      Category = ifelse("Category" %in% names(attr_tbl), as.character(Category), NA_character_),
      Sub_Cat = ifelse("Sub_Cat" %in% names(attr_tbl), as.character(Sub_Cat), NA_character_),
      Gov_Pvt = ifelse("Gov_Pvt" %in% names(attr_tbl), as.character(Gov_Pvt), NA_character_)
    ) %>%
    group_by(Category, Sub_Cat, Gov_Pvt) %>%
    summarise(
      total_facilities = n(),
      facilities_isolated = sum(isolated, na.rm = TRUE),
      facilities_with_road_access = total_facilities - facilities_isolated,
      percent_isolated = round(100 * facilities_isolated / total_facilities, 2),
      .groups = "drop"
    ) %>%
    mutate(service_radius_km = d/1000) %>%
    select(service_radius_km, Category, Sub_Cat, Gov_Pvt, total_facilities, facilities_with_road_access, facilities_isolated, percent_isolated)
  
  # align names for final combined table
  if (nrow(per_cat) > 0) {
    names(per_cat)[names(per_cat) == "Category"] <- "category"
    names(per_cat)[names(per_cat) == "Sub_Cat"] <- "sub_cat"
    names(per_cat)[names(per_cat) == "Gov_Pvt"] <- "gov_pvt"
    healthcare_by_category <- bind_rows(healthcare_by_category, per_cat)
  } else {
    # if no categories present, still record a row summarizing "unknown"
    healthcare_by_category <- bind_rows(healthcare_by_category, data.frame(
      service_radius_km = d/1000,
      category = NA_character_,
      sub_cat = NA_character_,
      gov_pvt = NA_character_,
      total_facilities = total_fac,
      facilities_with_road_access = fac_with_access,
      facilities_isolated = fac_iso,
      percent_isolated = pct_iso,
      stringsAsFactors = FALSE
    ))
  }
  
  message("Radius ", d/1000, " km: isolated = ", fac_iso, "/", total_fac, " (", pct_iso, "% )")
}

# Save outputs
write.csv(healthcare_service_results, "healthcare_service_results_swat_filtered_overall.csv", row.names = FALSE)
write.csv(healthcare_by_category, "healthcare_service_results_swat_by_category.csv", row.names = FALSE)

message("Saved healthcare_service_results_swat_filtered_overall.csv")
message("Saved healthcare_service_results_swat_by_category.csv")

# Optional: print first few rows of the by-category table
print(head(healthcare_by_category))
