5b. RF model full – residuals (IBD by using LC lakes paths)
================
Norah Saarman
2026-03-09

- [Inputs](#inputs)
- [1. Prepare the data](#1-prepare-the-data)
- [2. IBD regression to find
  residuals](#2-ibd-regression-to-find-residuals)
  - [Plot IBD linear model, CSE vs km, CSE vs
    log10(km)](#plot-ibd-linear-model-cse-vs-km-cse-vs-log10km)
- [3. Run IBD-residuals rf model (with residuals as response, drop
  pix_dist for
  predictor)](#3-run-ibd-residuals-rf-model-with-residuals-as-response-drop-pix_dist-for-predictor)
  - [Full random forest model with IBD
    residuals](#full-random-forest-model-with-ibd-residuals)
    - [Load saved (residuals) model](#load-saved-residuals-model)
- [4. Project predicted values from full IBD-residuals
  model](#4-project-predicted-values-from-full-ibd-residuals-model)
  - [Build Projection](#build-projection)
  - [Plot predicted residuals](#plot-predicted-residuals)
- [5. Scale and plot predicted connectivity (residuals) and
  SDM](#5-scale-and-plot-predicted-connectivity-residuals-and-sdm)
  - [Scale 0-1, habitat suitability and inverse of predicted
    connectivity
    (residuals)](#scale-0-1-habitat-suitability-and-inverse-of-predicted-connectivity-residuals)
  - [Plot scaled predicted residuals and
    SDM](#plot-scaled-predicted-residuals-and-sdm)
- [6. Variable importance plots](#6-variable-importance-plots)
  - [Percent Improvement MSE](#percent-improvement-mse)
  - [Node Purity](#node-purity)

RStudio Configuration:  
- **R version:** R 4.4.0 (Geospatial packages)  
- **Number of cores:** 4 (up to 32 available)  
- **Account:** saarman-np  
- **Partition:** saarman-np (allows multiple simultaneous jobs)  
- **Memory per job:** 100G (cluster limit: 1000G total; avoid exceeding
half)  
\# Setup

``` r
# load only required packages
library(randomForest)
library(doParallel)
library(raster)
library(sf)
library(viridis)
library(dplyr)
library(terra)
library(sf)
library(classInt)
library(raster)
library(RColorBrewer)
library(ggplot2)
library(factoextra)   # for nice PCA plots
library(ggpubr)

# base directories
data_dir  <- "/uufs/chpc.utah.edu/common/home/saarman-group1/uganda-tsetse-LG/data"
input_dir <- "../input"
results_dir <- "/uufs/chpc.utah.edu/common/home/saarman-group1/uganda-tsetse-LG/results"

# read the combined CSE + coords table + pix_dist + Env variables
V.table <- read.csv(file.path(input_dir, "Gff_cse_envCostPaths.csv"),
                    header = TRUE)
# This was added only after completing LOPOCV...
# Filter out western outlier "50-KB" 
V.table <- V.table %>%
  filter(Var1 != "50-KB", Var2 != "50-KB")

# define coordinate reference system
crs_geo <- 4326     # EPSG code for WGS84

# simple mode helper
get_mode <- function(x) {
  ux <- unique(x[!is.na(x)])
  ux[ which.max(tabulate(match(x, ux))) ]
}

# setup running in parallel
cl <- makeCluster(4)
registerDoParallel(cl)
clusterExport(cl, "get_mode")
```

# Inputs

- `../input/Gff_cse_envCostPaths.csv` - Combined CSE table with
  coordinates (long1, lat1, long2, lat2), pix_dist = geographic distance
  in sum of pixels, and mean, median, mode of each Env parameter

# 1. Prepare the data

``` r
# Assign input, checking for any rows with NA
sum(!complete.cases(V.table))  # should return 0
```

    ## [1] 0

``` r
rf_data <- na.omit(V.table)    # should omit zero rows

# Confirm that CSEdistance is numeric
rf_data$CSEdistance <- as.numeric(rf_data$CSEdistance)

# Select variables: all predictors (mean, median, mode)  
predictor_vars <- c("pix_dist",                      # geo dist
  paste0("BIO", 1:7, "_mean"),                       # mean 
  paste0("BIO", 8:11, "S_mean"),                     # mean
  paste0("BIO", 12:15, "_mean"),                     # mean
  paste0("BIO", 16:19, "S_mean"),                    # mean
  "alt_mean", "slope_mean", "riv_3km_mean",          # mean
  "samp_20km_mean", "lakes_mean",                    # mean
  paste0("BIO", 1:7, "_median"),                     # median
  paste0("BIO", 8:11, "S_median"),                   # median
  paste0("BIO", 12:15, "_median"),                   # median
  paste0("BIO", 16:19, "S_median"),                  # median
  "alt_median", "slope_median", "riv_3km_median",    # median
  "samp_20km_median", "lakes_median",                # median
  paste0("BIO", 1:7, "_mode"),                       # mode
  paste0("BIO", 8:11, "S_mode"),                     # mode
  paste0("BIO", 12:15, "_mode"),                     # mode
  paste0("BIO", 16:19, "S_mode"),                    # mode
  "alt_mode", "slope_mode", "riv_3km_mode",          # mode
  "samp_20km_mode", "lakes_mode"                     # mode
)

# subset predictors that we want to use
rf_data <- rf_data[, c("CSEdistance", predictor_vars)]

g <- lm(rf_data$CSEdistance~rf_data$pix_dist)
plot(rf_data$pix_dist, rf_data$CSEdistance)
abline(g)
```

![](05b_RFmodel_residuals_files/figure-gfm/prep-1.png)<!-- -->

``` r
# Extract groups of variables by suffix
mean_vars   <- grep("_mean$", names(rf_data), value = TRUE)
median_vars <- grep("_median$", names(rf_data), value = TRUE)
mode_vars   <- grep("_mode$", names(rf_data), value = TRUE)
```

# 2. IBD regression to find residuals

``` r
# Load raw data
V.table_full <- read.csv(file.path(input_dir, "Gff_cse_envCostPaths.csv"))

# estimate mean sampling density
mean(V.table_full$samp_20km_mean, na.rm = TRUE)
```

    ## [1] 1.027064e-11

``` r
# Filter out western outlier "50-KB" 
V.table <- V.table_full %>%
  filter(Var1 != "50-KB", Var2 != "50-KB")

# Create unique ID after filtering
V.table$id <- paste(V.table$Var1, V.table$Var2, sep = "_")

# Define site list
sites <- sort(unique(c(V.table$Var1, V.table$Var2)))

# How many rows of data for each?
table(V.table$Pop1_cluster)
```

    ## 
    ## north south 
    ##   595   496

``` r
# How many unique sites?
length(sites)
```

    ## [1] 67

``` r
# Choose predictors for RF model (all but pix_dist)
predictor_vars <- c("BIO1_mean","BIO2_mean","BIO3_mean","BIO4_mean", "BIO5_mean","BIO6_mean","BIO7_mean", "BIO8S_mean", "BIO9S_mean","BIO10S_mean", "BIO11S_mean","BIO12_mean", "BIO13_mean","BIO14_mean","BIO15_mean","BIO16S_mean","BIO17S_mean", "BIO18S_mean","BIO19S_mean","slope_mean","alt_mean", "lakes_mean","riv_3km_mean", "samp_20km_mean")
# "pix_dist") # REMOVED

# Fit IBD model
lm_ibd <- lm(CSEdistance ~ pix_dist, data = V.table)
summary(lm_ibd)
```

    ## 
    ## Call:
    ## lm(formula = CSEdistance ~ pix_dist, data = V.table)
    ## 
    ## Residuals:
    ##       Min        1Q    Median        3Q       Max 
    ## -0.170128 -0.038628 -0.001543  0.036318  0.183730 
    ## 
    ## Coefficients:
    ##              Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept) 2.814e-01  3.255e-03   86.45   <2e-16 ***
    ## pix_dist    4.746e-04  1.151e-05   41.22   <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.05633 on 1089 degrees of freedom
    ## Multiple R-squared:  0.6094, Adjusted R-squared:  0.609 
    ## F-statistic:  1699 on 1 and 1089 DF,  p-value: < 2.2e-16

``` r
# Add residuals to predictors table (V.table)
V.table$resid_ibd <- resid(lm_ibd)

# Filter modeling-relevant columns of V.table
rf_mean_data <- V.table[, c("resid_ibd", predictor_vars)]

# Rename predictors by removing "_mean" for later projections
names(rf_mean_data) <- gsub("_mean$", "", names(rf_mean_data))
```

### Plot IBD linear model, CSE vs km, CSE vs log10(km)

``` r
# colors
colors <- c("north" = "#1f78b4", "south" = "#e66101")

# raw geo dist
lm_ibd <- lm(CSEdistance ~ pix_dist, data = V.table)
summary(lm_ibd)
```

    ## 
    ## Call:
    ## lm(formula = CSEdistance ~ pix_dist, data = V.table)
    ## 
    ## Residuals:
    ##       Min        1Q    Median        3Q       Max 
    ## -0.170128 -0.038628 -0.001543  0.036318  0.183730 
    ## 
    ## Coefficients:
    ##              Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept) 2.814e-01  3.255e-03   86.45   <2e-16 ***
    ## pix_dist    4.746e-04  1.151e-05   41.22   <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.05633 on 1089 degrees of freedom
    ## Multiple R-squared:  0.6094, Adjusted R-squared:  0.609 
    ## F-statistic:  1699 on 1 and 1089 DF,  p-value: < 2.2e-16

``` r
plot(V.table$pix_dist, V.table$CSEdistance, col = colors[ V.table$Pop1_cluster],pch = 19, xlab = "Geo. distance (km)", ylab = "Gen. distance (CSE)")
abline(lm_ibd)
```

![](05b_RFmodel_residuals_files/figure-gfm/unnamed-chunk-1-1.png)<!-- -->

``` r
# log10 geo dist
lm_ibd_log <- lm(CSEdistance ~ log10(pix_dist), data = V.table)
summary(lm_ibd_log)
```

    ## 
    ## Call:
    ## lm(formula = CSEdistance ~ log10(pix_dist), data = V.table)
    ## 
    ## Residuals:
    ##       Min        1Q    Median        3Q       Max 
    ## -0.163873 -0.043257 -0.006829  0.040182  0.234622 
    ## 
    ## Coefficients:
    ##                  Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)     -0.040200   0.011552   -3.48 0.000521 ***
    ## log10(pix_dist)  0.191883   0.005024   38.19  < 2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.05893 on 1089 degrees of freedom
    ## Multiple R-squared:  0.5726, Adjusted R-squared:  0.5722 
    ## F-statistic:  1459 on 1 and 1089 DF,  p-value: < 2.2e-16

``` r
plot(log10(V.table$pix_dist), V.table$CSEdistance, col = colors[ V.table$Pop1_cluster],pch = 19, xlab = "Log10 Geo. distance (km)", ylab = "Gen. distance (CSE)")
abline(lm_ibd_log)
```

![](05b_RFmodel_residuals_files/figure-gfm/unnamed-chunk-1-2.png)<!-- -->

# 3. Run IBD-residuals rf model (with residuals as response, drop pix_dist for predictor)

## Full random forest model with IBD residuals

``` r
# Tune mtry (number of variables tried at each split)
set.seed(92834)
rf_resid_tuned <- tuneRF(
  x = rf_mean_data[, -1],   # exclude response variable
  y = rf_mean_data$resid_ibd,
  ntreeTry = 500,
  stepFactor = 1.5,         # factor by which mtry is increased/decreased
  improve = 0.01,           # minimum improvement to continue search
  trace = TRUE,             # print progress
  plot = TRUE,              # plot OOB error vs mtry
  doBest = TRUE,             # return the model with lowest OOB error
  importance = TRUE
)
```

    ## mtry = 8  OOB error = 0.001222935 
    ## Searching left ...
    ## mtry = 6     OOB error = 0.001235964 
    ## -0.01065346 0.01 
    ## Searching right ...
    ## mtry = 12    OOB error = 0.001220842 
    ## 0.001711967 0.01

![](05b_RFmodel_residuals_files/figure-gfm/rf%20full%20residuals-1.png)<!-- -->

``` r
print(rf_resid_tuned)
```

    ## 
    ## Call:
    ##  randomForest(x = x, y = y, mtry = res[which.min(res[, 2]), 1],      importance = TRUE) 
    ##                Type of random forest: regression
    ##                      Number of trees: 500
    ## No. of variables tried at each split: 12
    ## 
    ##           Mean of squared residuals: 0.001232708
    ##                     % Var explained: 61.08

``` r
importance(rf_resid_tuned)
```

    ##            %IncMSE IncNodePurity
    ## BIO1      11.76271    0.10056247
    ## BIO2      20.00268    0.19410970
    ## BIO3      33.65062    0.51685192
    ## BIO4      16.42328    0.11451964
    ## BIO5      15.64219    0.10812344
    ## BIO6      15.36999    0.15328888
    ## BIO7      16.08151    0.11926059
    ## BIO8S     11.61654    0.09646008
    ## BIO9S     19.54143    0.14925619
    ## BIO10S    20.12315    0.10359776
    ## BIO11S    18.01946    0.11473141
    ## BIO12     17.17012    0.10550581
    ## BIO13     30.03003    0.20312245
    ## BIO14     20.52008    0.15939933
    ## BIO15     19.01378    0.12658511
    ## BIO16S    16.39778    0.14256920
    ## BIO17S    14.43686    0.06435874
    ## BIO18S    10.75113    0.05133021
    ## BIO19S    13.71012    0.05100297
    ## slope     23.58799    0.15269396
    ## alt       15.19737    0.09111589
    ## lakes     14.98634    0.10148212
    ## riv_3km   30.00681    0.13639393
    ## samp_20km 36.19544    0.24478694

``` r
varImpPlot(rf_resid_tuned)
```

![](05b_RFmodel_residuals_files/figure-gfm/rf%20full%20residuals-2.png)<!-- -->

``` r
# Save the tuned random forest model to disk
saveRDS(rf_resid_tuned, file = file.path(results_dir, "rf_residuals.rds"))
```

### Load saved (residuals) model

``` r
# load saved model
rf_residuals <- readRDS(file.path(results_dir, "rf_residuals.rds"))

#double check they look correct
print(rf_residuals)
```

    ## 
    ## Call:
    ##  randomForest(x = x, y = y, mtry = res[which.min(res[, 2]), 1],      importance = TRUE) 
    ##                Type of random forest: regression
    ##                      Number of trees: 500
    ## No. of variables tried at each split: 12
    ## 
    ##           Mean of squared residuals: 0.001232708
    ##                     % Var explained: 61.08

``` r
print(rf_residuals$importance)
```

    ##                %IncMSE IncNodePurity
    ## BIO1      0.0002484851    0.10056247
    ## BIO2      0.0004917423    0.19410970
    ## BIO3      0.0020390824    0.51685192
    ## BIO4      0.0003720112    0.11451964
    ## BIO5      0.0002608535    0.10812344
    ## BIO6      0.0004497172    0.15328888
    ## BIO7      0.0003463552    0.11926059
    ## BIO8S     0.0003049704    0.09646008
    ## BIO9S     0.0005704390    0.14925619
    ## BIO10S    0.0002972806    0.10359776
    ## BIO11S    0.0003859081    0.11473141
    ## BIO12     0.0002520567    0.10550581
    ## BIO13     0.0004797717    0.20312245
    ## BIO14     0.0008521505    0.15939933
    ## BIO15     0.0003629246    0.12658511
    ## BIO16S    0.0003539647    0.14256920
    ## BIO17S    0.0002175375    0.06435874
    ## BIO18S    0.0001599565    0.05133021
    ## BIO19S    0.0001636527    0.05100297
    ## slope     0.0003519104    0.15269396
    ## alt       0.0003478297    0.09111589
    ## lakes     0.0001856338    0.10148212
    ## riv_3km   0.0002491449    0.13639393
    ## samp_20km 0.0006070110    0.24478694

# 4. Project predicted values from full IBD-residuals model

## Build Projection

``` r
# Load env stack with named layers
env <- stack(file.path(data_dir, "processed", "env_stack.grd"))

# Neutralize sampling layer to average
env$samp_20km <- 1.027064e-11 #neutralize sampling bias

# Load rdf of final model
rf_predicted <- readRDS(file.path(results_dir, "rf_residuals.rds"))
rf_predicted
```

    ## 
    ## Call:
    ##  randomForest(x = x, y = y, mtry = res[which.min(res[, 2]), 1],      importance = TRUE) 
    ##                Type of random forest: regression
    ##                      Number of trees: 500
    ## No. of variables tried at each split: 12
    ## 
    ##           Mean of squared residuals: 0.001232708
    ##                     % Var explained: 61.08

``` r
prediction_raster <- predict(env, rf_predicted, type = "response")

# Write Prediction Raster to file
writeRaster(prediction_raster, file.path(results_dir,"fullRF_residuals.tif"), format = "GTiff", overwrite = TRUE)
```

## Plot predicted residuals

``` r
# Create base plot with viridis
plot(prediction_raster,
     col = viridis::magma(100),
     main = "Predicted Residuals",
     axes = FALSE,
     box = FALSE,
     legend.args = list(text = "IBD residuals", side = 3, line = 1, cex = 1))

# Overlay lakes in dark gray
lakes <- st_read(file.path(data_dir, "raw/ne_10m_lakes.shp"), quiet = TRUE)
lakes <- st_transform(lakes, crs = st_crs(prediction_raster))  # match CRS 
lakes <- st_make_valid(lakes) # fix geometries
r_ext <- st_as_sfc(st_bbox(prediction_raster)) # extent
st_crs(r_ext) <- st_crs(prediction_raster) # match CRS
lakes <- st_intersection(lakes, r_ext) # clip to extent
plot(st_geometry(lakes), col = "gray20", border = NA, add = TRUE)

# Overlay country outline
uganda <- rnaturalearth::ne_countries(continent = "Africa", scale = "medium", returnclass = "sf")
uganda <- st_intersection(uganda, r_ext) # clip to extent
plot(st_geometry(uganda), col = NA, border = "black", lwd = 1.2, add = TRUE)
```

![](05b_RFmodel_residuals_files/figure-gfm/plot-projection-1.png)<!-- -->

# 5. Scale and plot predicted connectivity (residuals) and SDM

## Scale 0-1, habitat suitability and inverse of predicted connectivity (residuals)

``` r
# Load raster layers
con_raster <- rast(file.path(results_dir, "fullRF_residuals.tif"))
fao <- rast(file.path(data_dir, "FAO_fuscipes_2001.tif"))
update <- rast(file.path(data_dir, "SDM_2018update.tif"))

# Match extent and resolution first
fao_crop <- crop(fao, update)
update_crop <- crop(update, fao_crop)
fao_resamp <- resample(fao_crop, update_crop)  # if needed to match resolution

# Combine
sdm_raw <- max(fao_resamp, update_crop, na.rm = TRUE)

# Crop to overlapping extent
sdm <- crop(sdm_raw, con_raster)
con <- crop(con_raster, sdm)

# Mask low-suitability areas
sdm[sdm <= 0.05] <- NA


# Rescale to 0–1
sdm_min <- global(sdm, "min", na.rm = TRUE)$min
sdm_max <- global(sdm, "max", na.rm = TRUE)$max
sdm <- (sdm - sdm_min) / (sdm_max - sdm_min)

# Mask to common suitable area
con <- mask(con, sdm)

# Rescale inverse of prediction to 0-1
con_min <- global(con, "min", na.rm = TRUE)$min
con_max <- global(con, "max", na.rm = TRUE)$max
con <- 1 - ((con - con_min) / (con_max - con_min))

# Convert back to raster for compatibility with bivariate.map function
sdm_r <- raster(sdm)
con_r <- raster(con)
```

## Plot scaled predicted residuals and SDM

``` r
# Plot Genetic Connectivity (inverse predicted values)
plot(con,
     col = rev(viridis::plasma(100)),  # high connectivity = dark
     main = "Genetic Connectivity (inverse predicted residuals)",
     axes = FALSE, box = FALSE,
     legend.args = list(text = "Connectivity", side = 2, line = 2.5, cex = 0.8))
plot(st_geometry(lakes), col = "black", border = NA, add = TRUE)
plot(st_geometry(uganda), border = "black", lwd = 0.25, add = TRUE)
```

![](05b_RFmodel_residuals_files/figure-gfm/plot-sdm-con-1.png)<!-- -->

``` r
# Plot Habitat Suitability
plot(sdm,
     col = viridis::viridis(100),  # high suitability = dark
     main = "Habitat Suitability",
     axes = FALSE, box = FALSE,
     legend.args = list(text = "Suitability", side = 2, line = 2.5, cex = 0.8))
plot(st_geometry(lakes), col = "black", border = NA, add = TRUE)
plot(st_geometry(uganda), border = "black", lwd = 0.25, add = TRUE)
```

![](05b_RFmodel_residuals_files/figure-gfm/plot-sdm-con-2.png)<!-- -->

``` r
# Plot with custom colors

# Custom palettes based on Bishop et al.
connectivity_colors <- colorRampPalette(c("#FFFF00", "#FFA500", "#FF4500", "#700E40", "#2E003E"))(100)
suitability_colors  <- colorRampPalette(c("white", "lightblue", "blue4"))(100)     # white → light blue → dark blue

# Plot Genetic Connectivity (inverse predicted) with custom colors
plot(con,
     col = connectivity_colors,
     main = "Genetic Connectivity (inverse predicted residuals)",
     axes = FALSE, box = FALSE,
     legend.args = list(text = "Connectivity", side = 2, line = 2.5, cex = 0.8))
plot(st_geometry(lakes), col = "black", border = NA, add = TRUE)
plot(st_geometry(uganda), border = "black", lwd = 0.25, add = TRUE)
```

![](05b_RFmodel_residuals_files/figure-gfm/plot-sdm-con-3.png)<!-- -->

``` r
# Plot Habitat Suitability with custom colors
plot(sdm,
     col = suitability_colors,
     main = "Habitat Suitability",
     axes = FALSE, box = FALSE,
     legend.args = list(text = "Suitability", side = 2, line = 2.5, cex = 0.8))
plot(st_geometry(lakes), col = "black", border = NA, add = TRUE)
plot(st_geometry(uganda), border = "black", lwd = .25, add = TRUE)
```

![](05b_RFmodel_residuals_files/figure-gfm/plot-sdm-con-4.png)<!-- -->

# 6. Variable importance plots

## Percent Improvement MSE

``` r
library(dplyr)
library(tibble)
library(ggplot2)
library(randomForest)

# Load full model (as opposed to LOPOCV later)
full_model <- rf_predicted 
full_imp <- importance(full_model, type = 1) %>%
  as.data.frame() %>%
  rownames_to_column("variable") %>%
  rename(IncMSE = `%IncMSE`) %>%
  mutate(model = "full")

# Define custom labels
label_map <- c(
  BIO1   = "Annual Mean Temperature (BIO1)",
  BIO2   = "Mean Diurnal Temp Range (BIO2)",
  BIO3   = "Isothermality (BIO3)",
  BIO4   = "Temperature Seasonality (BIO4)",
  BIO5   = "Max Temp of Warmest Month (BIO5)",
  BIO6   = "Min Temp of Coldest Month (BIO6)",
  BIO7   = "Temperature Annual Range (BIO7)",
  BIO8S  = "Mean Temp of Wettest Season (BIO8S)",
  BIO9S  = "Mean Temp of Driest Season (BIO9S)",
  BIO10S = "Mean Temp of Warmest Season (BIO10S)",
  BIO11S = "Mean Temp of Coldest Season (BIO11S)",
  BIO12  = "Annual Precipitation (BIO12)",
  BIO13  = "Precipitation of Wettest Month (BIO13)",
  BIO14  = "Precipitation of Driest Month (BIO14)",
  BIO15  = "Precipitation Seasonality (BIO15)",
  BIO16S = "Precipitation of Wettest Season (BIO16S)",
  BIO17S = "Precipitation of Driest Season (BIO17S)",
  BIO18S = "Precipitation of Warmest Season (BIO18S)",
  BIO19S = "Precipitation of Coldest Season (BIO19S)",
  slope  = "Slope",
  alt    = "Altitude",
  lakes  = "Lake Presence/Absence",
  riv_3km = "River Kernel Density (3 km bandwidth)",
  samp_20km = "Sampling Density (20 km bandwidth)",
  pix_dist = "Geographic Distance (km)"
)


# Order by full model's %IncMSE (top to bottom)
full_order <- full_imp %>%
  arrange(desc(IncMSE)) %>%
  pull(variable)

full_imp$variable <- factor(full_imp$variable, levels = rev(full_order))

# Plot
# pdf("../figures/VarImpPlot_residModel.pdf",width =6, height=6)
ggplot(full_imp, aes(x = variable, y = IncMSE)) +
  geom_point(data = filter(full_imp, model == "full"),
             color = "black", size = 3) +
  coord_flip() +
  scale_y_continuous(name = "%IncMSE") +
  scale_x_discrete(labels = label_map) +
  labs(x = NULL, title = "Variable Importance of IBD Residuals Model") +
  theme_minimal()
```

![](05b_RFmodel_residuals_files/figure-gfm/variable-imp-mse-1.png)<!-- -->

``` r
#dev.off()
```

## Node Purity

``` r
library(dplyr)
library(tibble)
library(ggplot2)
library(randomForest)

# Load full model
full_model <- rf_predicted 

full_imp <- importance(full_model, type = 2) %>%
  as.data.frame() %>%
  rownames_to_column("variable") %>%
  mutate(model = "full")

# Define custom labels
label_map <- c(
  BIO1   = "Annual Mean Temperature (BIO1)",
  BIO2   = "Mean Diurnal Temp Range (BIO2)",
  BIO3   = "Isothermality (BIO3)",
  BIO4   = "Temperature Seasonality (BIO4)",
  BIO5   = "Max Temp of Warmest Month (BIO5)",
  BIO6   = "Min Temp of Coldest Month (BIO6)",
  BIO7   = "Temperature Annual Range (BIO7)",
  BIO8S  = "Mean Temp of Wettest Season (BIO8S)",
  BIO9S  = "Mean Temp of Driest Season (BIO9S)",
  BIO10S = "Mean Temp of Warmest Season (BIO10S)",
  BIO11S = "Mean Temp of Coldest Season (BIO11S)",
  BIO12  = "Annual Precipitation (BIO12)",
  BIO13  = "Precipitation of Wettest Month (BIO13)",
  BIO14  = "Precipitation of Driest Month (BIO14)",
  BIO15  = "Precipitation Seasonality (BIO15)",
  BIO16S = "Precipitation of Wettest Season (BIO16S)",
  BIO17S = "Precipitation of Driest Season (BIO17S)",
  BIO18S = "Precipitation of Warmest Season (BIO18S)",
  BIO19S = "Precipitation of Coldest Season (BIO19S)",
  slope  = "Slope",
  alt    = "Altitude",
  lakes  = "Lake Presence/Absence",
  riv_3km = "River Kernel Density (3 km bandwidth)",
  samp_20km = "Sampling Density (20 km bandwidth)",
  pix_dist = "Geographic Distance (km)"
)

# Order variables by node purity
full_order <- full_imp %>%
  arrange(desc(IncNodePurity)) %>%
  pull(variable)

full_imp$variable <- factor(full_imp$variable, levels = rev(full_order))

# Plot
ggplot(full_imp, aes(x = variable, y = IncNodePurity)) +
  geom_point(color = "black", size = 3) +
  coord_flip() +
  scale_y_continuous(name = "Increase in Node Purity") +
  scale_x_discrete(labels = label_map) +
  labs(x = NULL, title = "Variable Importance of IBD Residuals Model") +
  theme_minimal()
```

![](05b_RFmodel_residuals_files/figure-gfm/variable-imp-nodepurity-1.png)<!-- -->
