library(sf); library(dplyr); library(units)
library(tidyr)


# read back saved flooded pieces and master roads (UTM)
roads_master <- st_read("roads/clean.shp", quiet = TRUE) %>% st_set_crs(4326) %>% st_transform(32642)
roads_master$road_id <- seq_len(nrow(roads_master))
# if you already have roads_sub in memory, use that instead

roads_flooded <- st_read("roads_flooded_swat_fixed.gpkg", quiet = TRUE) %>% st_transform(32642)

# ensure flooded pieces have road_id (they should)
if (!"road_id" %in% names(roads_flooded)) {
  warning("roads_flooded missing road_id — attempting spatial join to recover")
  roads_flooded <- st_join(roads_flooded, roads_master %>% select(road_id, highway), join = st_intersects, left = TRUE)
}

# recompute affected lengths (meters)
roads_flooded$affected_m <- as.numeric(st_length(st_transform(roads_flooded, 32642)))

# Recreate road_class on master (you said only Motorway, Trunk, Primary, Secondary exist)
roads_master <- roads_master %>% mutate(
  road_class = case_when(
    highway %in% c("motorway","motorway_link") ~ "Motorway",
    highway %in% c("trunk","trunk_link")       ~ "Trunk",
    highway %in% c("primary","primary_link")   ~ "Primary",
    highway %in% c("secondary","secondary_link") ~ "Secondary",
    TRUE ~ "Other"
  ),
  len_m = as.numeric(st_length(.))
)

# total length by class (km)
total_length_by_type <- roads_master %>%
  st_drop_geometry() %>%
  group_by(road_class) %>%
  summarise(total_length_km = sum(len_m, na.rm = TRUE)/1000) %>%
  arrange(desc(total_length_km))

# affected length by class (join back)
# map flooded segments to their road_class using spatial join to master (if needed)
roads_flooded <- st_join(roads_flooded, roads_master %>% select(road_id, road_class), join = st_intersects, left = TRUE)

affected_length_by_type <- roads_flooded %>%
  st_drop_geometry() %>%
  group_by(road_class) %>%
  summarise(affected_length_km = sum(affected_m, na.rm = TRUE)/1000)

road_impact_summary <- total_length_by_type %>%
  left_join(affected_length_by_type, by = "road_class") %>%
  mutate(affected_length_km = replace_na(affected_length_km, 0),
         percent_affected = round(100 * affected_length_km / total_length_km, 2)) %>%
  arrange(desc(affected_length_km))

print(road_impact_summary)
write.csv(road_impact_summary, "road_network_impact_summary_swat.csv", row.names = FALSE)
message("Saved road_network_impact_summary_swat.csv")
