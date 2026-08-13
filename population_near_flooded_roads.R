library(terra); library(sf)

#3) Population near flooded roads (rasterise flooded-road buffer & sum population rasters)

# flood roads buffer (in UTM), then transform to flood raster CRS if needed
roads_flooded_utm <- st_transform(roads_flooded, 32642)
flooded_roads_buffer <- st_buffer(roads_flooded_utm, dist = 500)  # road_exposure_dist
# transform to raster CRS (your flood raster may already be in UTM)
flood_raster <- rast("results/swat/flood_rasterized.tif")
flood_raster_crs <- crs(flood_raster)

# convert to SpatVector and rasterize
flooded_vect <- vect(st_transform(flooded_roads_buffer, st_crs(flood_raster)$wkt))
flooded_rast <- rasterize(flooded_vect, flood_raster, field = 1, background = 0)

# process population rasters (folder)
pop_folder <- "pak_agesex_structures_2025_CN_100m_R2025A_v1"
pop_files <- list.files(pop_folder, pattern = "\\.tif$", full.names = TRUE)
road_exposure_results <- data.frame(category = character(), pop_near_flooded_roads = numeric(), percent_of_total = numeric(), stringsAsFactors = FALSE)

for (fp in pop_files) {
  pr <- rast(fp)
  if (!compareGeom(pr, flooded_rast, stopOnError = FALSE)) {
    pr <- resample(pr, flooded_rast, method = "near")
  }
  total_pop <- global(pr, "sum", na.rm = TRUE)[1,1]
  pop_near <- global(pr * flooded_rast, "sum", na.rm = TRUE)[1,1]
  road_exposure_results <- rbind(road_exposure_results, data.frame(
    category = basename(fp),
    pop_near_flooded_roads = round(pop_near),
    percent_of_total = round(pop_near / total_pop * 100, 2)
  ))
}
print(road_exposure_results)
write.csv(road_exposure_results, "population_near_flooded_roads_swat.csv", row.names = FALSE)
message("Saved population_near_flooded_roads_swat.csv")
