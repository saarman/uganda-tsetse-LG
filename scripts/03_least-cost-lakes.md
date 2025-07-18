Least-cost paths avoiding lakes
================
Norah Saarman
2025-06-10

- [Setup](#setup)
- [Inputs](#inputs)
- [Outputs](#outputs)
- [1. Prepare your cost surface from lake
  layer](#1-prepare-your-cost-surface-from-lake-layer)
- [2. Build and save paths as SpatialLines for each
  pair](#2-build-and-save-paths-as-spatiallines-for-each-pair)

RStudio Configuration:  
- **R version:** R 4.4.0 (Geospatial packages)  
- **Number of cores:** 4 (up to 32 available)  
- **Account:** saarman-np  
- **Partition:** saarman-shared-np (allows multiple simultaneous jobs)  
- **Memory per job:** 100G (cluster limit: 1000G total; avoid exceeding
half)

# Setup

``` r
# load only required packages
library(raster)
library(gdistance)
library(sp)
library(sf)
library(foreach)
library(doParallel)
library(ggplot2)
library(rnaturalearth)
library(viridis)
library(units)
library(future.apply)

# base directories
data_dir  <- "/uufs/chpc.utah.edu/common/home/saarman-group1/uganda-tsetse-LG/data"
input_dir <- "../input"

# read the combined CSE + coords table
G.table <- read.csv(file.path(input_dir, "Gff_11loci_68sites_cse.csv"),
                    header = TRUE)

# simple mode helper
get_mode <- function(x) {
  ux <- unique(x[!is.na(x)])
  ux[ which.max(tabulate(match(x, ux))) ]
}

# define coordinate reference system
crs_geo <- 4326     # EPSG code for WGS84
```

# Inputs

- `../input/Gff_11loci_68sites_cse.csv` - Combined CSE table with
  coordinates (long1, lat1, long2, lat2) river+lake edge density
- `../data_dir/processed/sample_kernel_density_20km.tif` \# 20 km
  sampling density
- `../data_dir/processed/lake_binary.tif` \# binary lake mask (1 =
  water, 0 = land)

# Outputs

- `../data_dir/processed/LC_paths.shp` - Least cost paths spatialLines
  shapefile

# 1. Prepare your cost surface from lake layer

Assume lake_binary.tif has 1 = water, 0 = land

``` r
# read your lake mask (0 = land, 1 = lake)
lakes <- raster(file.path(data_dir, "processed", "lake_binary.tif"))

# make a cost/connectivity  surface
costRaster <- lakes
costRaster[lakes == 0] <- 10   # land cells cost/connectivity 10
costRaster[lakes == 1] <- 1   # lake cells cost/connectivity 1

# build the transition (minimize mean cost/connectivity)
tr <- transition(
  costRaster,
  transitionFunction = function(x) mean(x),
  directions         = 16
)

# correct to meters
trCorr <- geoCorrection(tr, type = "c")
```

# 2. Build and save paths as SpatialLines for each pair

\*\* NOTE:\*\* eval = FALSE so that skips on knit

``` r
# assume trCorr (cost transition) and G.table (with long1,lat1,long2,lat2,CSEdistance) exist

# use 4 cores (or however many you want)
plan(multisession, workers = 4)

# build least-cost SpatialLinesDataFrame for each pair
paths <- future_lapply(seq_len(nrow(G.table)), function(i) {
  from <- c(G.table$long1[i], G.table$lat1[i])
  to   <- c(G.table$long2[i], G.table$lat2[i])
  sl   <- shortestPath(trCorr, from, to, output = "SpatialLines")
  sp::proj4string(sl) <- sp::CRS("+proj=longlat +datum=WGS84 +no_defs")

  # build row-specific metadata
  metadata <- data.frame(
    Var1 = G.table$Var1[i],
    Var2 = G.table$Var2[i],
    id   = paste(G.table$Var1[i], G.table$Var2[i], sep = "_"),
    CSE  = G.table$CSEdistance[i],
    row.names = row.names(sl)
  )

  SpatialLinesDataFrame(sl, data = metadata)
})

# return to single-threaded mode
plan(sequential)

# combine into one object and convert to sf
lines_sldf <- do.call(rbind, paths)
lines_sf   <- st_as_sf(lines_sldf)
lines_sf   <- st_set_crs(lines_sf, 4326)

# write to shapefile with full metadata
st_write(lines_sf,
         dsn = file.path(data_dir, "processed", "LC_paths.shp"),
         delete_layer = TRUE)

lines_sf
```

Plot least cost paths with Uganda outline and lakes to check behavior

``` r
# assumes shape file of least-cost paths have already been created
lines_sf <- st_read(file.path(data_dir,"processed","LC_paths.shp"), quiet=TRUE)

# This was added only after completing LOPOCV...
# Filter out western outlier "50-KB" 
lines_sf <- lines_sf[
  lines_sf$Var1 != "50-KB" &
  lines_sf$Var2 != "50-KB",
]

# Read and transform lakes, then crop to Uganda
# Extract map extent
r_ext <- extent(lakes)
xlim <- c(r_ext@xmin, r_ext@xmax)
ylim <- c(r_ext@ymin, r_ext@ymax)

# Natural Earth background
uganda <- ne_countries(scale = "medium", continent = "Africa", returnclass = "sf") %>% st_transform(4326)
lakes_ne <- ne_download(scale = 10, type = "lakes", category = "physical", returnclass = "sf") %>% st_transform(4326)
```

    ## Reading layer `ne_10m_lakes' from data source `/tmp/RtmpeZa9c3/ne_10m_lakes.shp' using driver `ESRI Shapefile'
    ## Simple feature collection with 1355 features and 41 fields
    ## Geometry type: MULTIPOLYGON
    ## Dimension:     XY
    ## Bounding box:  xmin: -165.9656 ymin: -50.66967 xmax: 177.1544 ymax: 81.95521
    ## Geodetic CRS:  WGS 84

``` r
# Clean up invalid geometries
lakes_ne <- st_make_valid(lakes_ne)

# plot least‐cost paths
ggplot() +
  geom_sf(data = uganda, fill="white", color="black") +
  geom_sf(data = lines_sf,
          aes(color = CSE),
          size = 0.6) +
  geom_sf(data = lakes_ne,
          fill = "lightblue",
          color = NA,
          alpha = 0.5) +
  scale_color_viridis_c(
    name      = "CSE",
    option    = "magma",
    direction = -1
  ) +
  coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
  theme_void() +
  theme(
    legend.position = c(1.1, 0.5),
    legend.background = element_rect(fill = alpha("white", 0.6), color = NA)
  )
```

<figure>
<img
src="03_least-cost-lakes_files/figure-gfm/plot_least_cost_paths-1.png"
alt="Least‐cost paths colored by CSE" />
<figcaption aria-hidden="true">Least‐cost paths colored by
CSE</figcaption>
</figure>
