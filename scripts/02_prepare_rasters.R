# scripts/02_prepare_rasters.R

# RStudio Configuration:
# R version: R 4.4.0 (Geospatial packages)  
# Number of cores: 4 (up to 32 available)   
# Account: saarman-np  
# Partition: saarman-shared-np (allows multiple simultaneous jobs)  
# Memory per job: 100G (cluster limit: 1000G total; avoid exceeding half)

# load only required packages
library(raster)
library(sf)
library(terra)
library(KernSmooth)
library(maps)
library(units)
library(future.apply)

# setwd from where you save the .R file
setwd("scripts")

# define data directory
data_dir <- "/uufs/chpc.utah.edu/common/home/saarman-group1/uganda-tsetse-LG/data"

# define coordinate reference system
crs_geo <- 4326     # EPSG code for WGS84

# Read template raster for extent and projection
pixels_raster <- raster(file.path(data_dir, "raw/slope_1KMmedian_MERIT_UgandaClip.tif"))

########################################
# Sampling Kernel Density
########################################

# Read pairwise CSE data and extract unique coordinates
pairs <- read.csv("../input/Gff_11loci_68sites_cse.csv", header = TRUE)
coords_df <- unique(rbind(
  data.frame(long = pairs$long1, lat = pairs$lat1),
  data.frame(long = pairs$long2, lat = pairs$lat2)
))

# Convert to sf and reproject to equal-distance projection
points_sf <- st_as_sf(coords_df, coords = c("long", "lat"), crs = crs_geo)
eqd_crs <- "+proj=eqdc +lat_0=39 +lon_0=-96 +lat_1=33 +lat_2=45 +datum=NAD83 +units=m"
points_m  <- st_transform(points_sf, crs = eqd_crs)
coords_m  <- st_coordinates(points_m)

# Build projected raster template at 1 km resolution
proj_temp <- terra::project(terra::rast(pixels_raster), eqd_crs, res = 1000)
ext_m     <- terra::ext(proj_temp)

# Compute 2D kernel density (20 km bandwidth)
est <- bkde2D(coords_m,
              bandwidth = c(20000, 20000),
              gridsize  = c(nrow(proj_temp), ncol(proj_temp)),
              range.x   = list(c(ext_m$xmin, ext_m$xmax), c(ext_m$ymin, ext_m$ymax)))

# Convert result to raster, project back to WGS84, and match grid
kd_raster <- raster(list(x = est$x1, y = est$x2, z = est$fhat),
                    crs = eqd_crs)
kd_wgs84  <- projectRaster(kd_raster, crs = crs(pixels_raster), res = res(pixels_raster))
sampling_kd <- resample(kd_wgs84, pixels_raster, method = "bilinear")
names(sampling_kd) <- "sampling_20km"
crs(sampling_kd) <- crs(pixels_raster)
plot(sampling_kd)

# Save output
writeRaster(sampling_kd,
            file.path(data_dir, "processed/sampling_kernel_density_20km.tif"),
            overwrite = TRUE)

########################################
# River Kernel Density
########################################

# Read river shapefile and ensure WGS84 projection
rivers <- st_read(file.path(data_dir, "raw/Uganda_rivers_shape.shp"))
rivers <- st_transform(rivers, crs = crs_geo)

# Set up parallel workers
options(future.rng.onMisuse = "ignore")  # suppress future RNG warnings
plan(multisession, workers = 4)          # 4 local cores, each gets a chunk

# Densify rivers using 1 km spacing in meters
river_dense_list <- future_lapply(st_geometry(rivers), function(geom) {
  st_segmentize(geom, dfMaxLength = set_units(1000, "m"))
})

# Rebuild full sf object from geometry list
rivers_dense <- st_sf(geometry = st_sfc(river_dense_list, crs = st_crs(rivers)))

# Reset plan
plan(sequential)

# Cast to POINT geometries and project to equal-distance CRS
eqd_crs <- "+proj=eqdc +lat_0=39 +lon_0=-96 +lat_1=33 +lat_2=45 +datum=NAD83 +units=m"
river_points <- st_cast(rivers_dense, "POINT")
river_points_eqd <- st_transform(river_points, crs = eqd_crs)
coords_river_m <- st_coordinates(river_points_eqd)

# --- 5 km bandwidth ---

# Use same raster template and extent as for sampling
est_riv <- bkde2D(coords_river_m,
                  bandwidth = c(5000, 5000),
                  gridsize  = c(nrow(proj_temp), ncol(proj_temp)),
                  range.x   = list(c(ext_m$xmin, ext_m$xmax), c(ext_m$ymin, ext_m$ymax)))

# Convert to raster, reproject, resample
kd_river <- raster(list(x = est_riv$x1, y = est_riv$x2, z = est_riv$fhat),
                   crs = eqd_crs)
kd_river_wgs84 <- projectRaster(kd_river, crs = crs(pixels_raster), res = res(pixels_raster))
river_kd <- resample(kd_river_wgs84, pixels_raster, method = "bilinear")
names(river_kd) <- "river_5km"
crs(river_kd) <- crs(pixels_raster)
plot(river_kd)

# Save output
writeRaster(river_kd,
            file.path(data_dir, "processed/river_kernel_density_5km.tif"),
            overwrite = TRUE)

# --- 2 km bandwidth ---

# Use same raster template and extent as for sampling
est_riv <- bkde2D(coords_river_m,
                  bandwidth = c(2000, 2000),
                  gridsize  = c(nrow(proj_temp), ncol(proj_temp)),
                  range.x   = list(c(ext_m$xmin, ext_m$xmax), c(ext_m$ymin, ext_m$ymax)))

# Convert to raster, reproject, resample
kd_river <- raster(list(x = est_riv$x1, y = est_riv$x2, z = est_riv$fhat),
                   crs = eqd_crs)
kd_river_wgs84 <- projectRaster(kd_river, crs = crs(pixels_raster), res = res(pixels_raster))
river_kd <- resample(kd_river_wgs84, pixels_raster, method = "bilinear")
names(river_kd) <- "river_2km"
crs(river_kd) <- crs(pixels_raster)
plot(river_kd)

# Save output
writeRaster(river_kd,
            file.path(data_dir, "processed/river_kernel_density_2km.tif"),
            overwrite = TRUE)

# --- 1 km bandwidth ---

# Use same raster template and extent as for sampling
est_riv <- bkde2D(coords_river_m,
                  bandwidth = c(1000, 1000),
                  gridsize  = c(nrow(proj_temp), ncol(proj_temp)),
                  range.x   = list(c(ext_m$xmin, ext_m$xmax), c(ext_m$ymin, ext_m$ymax)))

# Convert to raster, reproject, resample
kd_river <- raster(list(x = est_riv$x1, y = est_riv$x2, z = est_riv$fhat),
                   crs = eqd_crs)
kd_river_wgs84 <- projectRaster(kd_river, crs = crs(pixels_raster), res = res(pixels_raster))
river_kd <- resample(kd_river_wgs84, pixels_raster, method = "bilinear")
names(river_kd) <- "river_1km"
crs(river_kd) <- crs(pixels_raster)
plot(river_kd)

# Save output
writeRaster(river_kd,
            file.path(data_dir, "processed/river_kernel_density_1km.tif"),
            overwrite = TRUE)