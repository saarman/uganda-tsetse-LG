Spatial Evaluation of Alternative Response Formulations
================
Norah Saarman
2026-08-12

- [Overview of script](#overview-of-script)
  - [Load libraries](#load-libraries)
- [Raw CSE](#raw-cse)
  - [Step 1: Inputs](#step-1-inputs)
  - [Step 2: Check raster alignment and assign high cost to
    water](#step-2-check-raster-alignment-and-assign-high-cost-to-water)
  - [Step 3: Prepare pair table, site objects, and transition
    object](#step-3-prepare-pair-table-site-objects-and-transition-object)
  - [Step 4: Sanity checks](#step-4-sanity-checks)
  - [Step 5: Compute least-cost paths in parallel (once) and save
    output](#step-5-compute-least-cost-paths-in-parallel-once-and-save-output)
  - [Step 7: Extract predicted CSE values along LCPs and calculate mean
    and
    sum](#step-7-extract-predicted-cse-values-along-lcps-and-calculate-mean-and-sum)
  - [Step 8: Compute circuit-based effective resistance
    (CBER)](#step-8-compute-circuit-based-effective-resistance-cber)
  - [Step 9: Combine LCP and CBER
    results](#step-9-combine-lcp-and-cber-results)
  - [Step 10: Calculate evaluation metrics from
    eval_results](#step-10-calculate-evaluation-metrics-from-eval_results)
- [IBD residuals](#ibd-residuals)
  - [Step 1: Inputs](#step-1-inputs-1)
  - [Step 2: Check raster alignment and assign high cost to
    water](#step-2-check-raster-alignment-and-assign-high-cost-to-water-1)
  - [Step 3: Prepare pair table, site objects, and transition
    object](#step-3-prepare-pair-table-site-objects-and-transition-object-1)
  - [Step 4: Sanity checks](#step-4-sanity-checks-1)
  - [Step 5: Compute least-cost paths in parallel (once) and save
    output](#step-5-compute-least-cost-paths-in-parallel-once-and-save-output-1)
  - [Step 7: Extract predicted CSE values along LCPs and calculate mean
    and
    sum](#step-7-extract-predicted-cse-values-along-lcps-and-calculate-mean-and-sum-1)
  - [Step 8: Compute circuit-based effective resistance
    (CBER)](#step-8-compute-circuit-based-effective-resistance-cber-1)
  - [Step 9: Combine LCP and CBER
    results](#step-9-combine-lcp-and-cber-results-1)
  - [Step 10: Calculate evaluation metrics from
    eval_results](#step-10-calculate-evaluation-metrics-from-eval_results-1)
- [CSE per km](#cse-per-km)
  - [Step 1: Inputs](#step-1-inputs-2)
  - [Step 2: Check raster alignment and assign high cost to
    water](#step-2-check-raster-alignment-and-assign-high-cost-to-water-2)
  - [Step 3: Prepare pair table, site objects, and transition
    object](#step-3-prepare-pair-table-site-objects-and-transition-object-2)
  - [Step 4: Sanity checks](#step-4-sanity-checks-2)
  - [Step 5: Compute least-cost paths in parallel (once) and save
    output](#step-5-compute-least-cost-paths-in-parallel-once-and-save-output-2)
  - [Step 7: Extract predicted CSE values along LCPs and calculate mean
    and
    sum](#step-7-extract-predicted-cse-values-along-lcps-and-calculate-mean-and-sum-2)
  - [Step 8: Compute circuit-based effective resistance
    (CBER)](#step-8-compute-circuit-based-effective-resistance-cber-2)
  - [Step 9: Combine LCP and CBER
    results](#step-9-combine-lcp-and-cber-results-2)
  - [Step 10: Calculate evaluation metrics from
    eval_results](#step-10-calculate-evaluation-metrics-from-eval_results-2)
- [Screening results](#screening-results)
  - [Screening metrics](#screening-metrics)
  - [R-squared comparison](#r-squared-comparison)
  - [Correlation comparison](#correlation-comparison)
  - [Best-performing combination by
    R-squared](#best-performing-combination-by-r-squared)
  - [Observed versus predicted plots from existing saved RDS
    files](#observed-versus-predicted-plots-from-existing-saved-rds-files)

RStudio Configuration:  
- **R version:** R 4.4.0 (Geospatial packages)  
- **Number of cores:** 4 (up to 32 available)  
- **Account:** saarman-np  
- **Partition:** saarman-np (now auto allows multiple simultaneous
jobs)  
- **Memory per job:** 100G (cluster limit: 1000G total; avoid exceeding
half)

# Overview of script

**This script is not part of the final workflow.**

Instead, it evaluates alternative formulations of genetic distance and
alternative approaches for summarizing predicted resistance across the
landscape. This is a screening step using the full model only, with the
goal of selecting the response formulation and spatial summary to carry
forward into the final LOPOCV validation framework.

Response formulations compared:

- raw conditional squared Euclidean (CSE) genetic distance
- CSE residuals from a simple isolation-by-distance (IBD) model
- CSE standardized by geographic distance (CSE/km)

For each response formulation, spatial evaluation is performed across
the corresponding random forest resistance surface using three spatial
summaries:

1.  Mean predicted resistance across the least-cost path (LCP-mean)
2.  Sum of predicted resistance across the least-cost path (LCP-sum)
3.  Circuit-based effective resistance (CBER)

Metrics calculated for this screening step:

**Primary metrics:**

1.  Correlation - evaluates whether the spatial prediction recovers the
    relative ordering of observed genetic distances.
2.  R² (1 - SSE/SST) - calculated after linear calibration to correct
    for scale mismatch between spatial predictions and observed CSE.

**Secondary metrics:**

3.  RMSE - calculated after linear calibration
4.  MSE - calculated after linear calibration
5.  MAE - calculated after linear calibration

Correlation is calculated from the raw spatial predictions because it is
scale-invariant. R², RMSE, MSE, and MAE are calculated from calibrated
predictions and therefore evaluate predictive performance after
correcting for simple scale mismatch.

These screening metrics are used only to select the response formulation
and spatial summary carried forward into the final LOPOCV validation
framework and permutation-based significance testing. It is not part of
the final workflow.

## Load libraries

``` r
# load only required packages
library(doParallel)
library(foreach)
library(raster)
library(gdistance)
library(sf)
library(dplyr)
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)
```

# Raw CSE

### Step 1: Inputs

``` r
# Define Paths to directories
input_dir <- "../input"
data_dir  <- "/uufs/chpc.utah.edu/common/home/saarman-group1/uganda-tsetse-LG/data"
results_dir_big <- "/uufs/chpc.utah.edu/common/home/saarman-group1/uganda-tsetse-LG/results/"
results_dir <- "../results/"

# define coordinate reference system
crs_geo <- 4326     # EPSG code for WGS84

######################################
# What surface are you using?
surf_file  <- file.path(results_dir_big, "fullRF_CSE.tif")
########################################

lake_file  <- file.path(data_dir, "processed", "lake_binary.tif")

cse_file   <- file.path("..", "input", "Gff_cse_envCostPaths.csv")

# Read pairwise CSE table
cse_df <- read.csv(cse_file)
# This was added only after completing LOPOCV...
# Filter out western outlier "50-KB" 
cse_df <- cse_df %>%
  filter(Var1 != "50-KB", Var2 != "50-KB")


# Read projected raw CSE surface
cse_surface <- raster(surf_file)

# Read lake mask
lake_mask <- raster(lake_file)
```

### Step 2: Check raster alignment and assign high cost to water

``` r
if (!compareRaster(cse_surface, lake_mask,
                   extent = TRUE, rowcol = TRUE,
                   crs = TRUE, res = TRUE,
                   stopiffalse = FALSE)) {
  lake_mask <- projectRaster(lake_mask, cse_surface, method = "ngb")
}

# Inspect projected raw CSE values before altering water cells
summary(values(cse_surface))
```

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    ##  0.1877  0.3098  0.3323  0.3322  0.3608  0.4338

``` r
# Keep water traversible but very costly (max value from model)
cse_surface[lake_mask[] == 1] <- max(values(cse_surface))

# Check again after assigning lakes
summary(values(cse_surface))
```

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    ##  0.1922  0.3204  0.3386  0.3487  0.3728  0.4338

``` r
# Quick plot
plot(cse_surface, main = "Projected raw CSE surface with water = high cost")
```

![](/tmp/RtmpaVWtN9/knitted_mds/raw-check-alignment-1.png)<!-- -->

``` r
# Check column names of cse_df
names(cse_df)
```

    ##  [1] "Var1"             "Var2"             "CSEdistance"      "lat1"            
    ##  [5] "long1"            "lat2"             "long2"            "Pop1_cluster"    
    ##  [9] "Pop2_cluster"     "pix_dist"         "BIO1_mean"        "BIO2_mean"       
    ## [13] "BIO3_mean"        "BIO4_mean"        "BIO5_mean"        "BIO6_mean"       
    ## [17] "BIO7_mean"        "BIO8S_mean"       "BIO9S_mean"       "BIO10S_mean"     
    ## [21] "BIO11S_mean"      "BIO12_mean"       "BIO13_mean"       "BIO14_mean"      
    ## [25] "BIO15_mean"       "BIO16S_mean"      "BIO17S_mean"      "BIO18S_mean"     
    ## [29] "BIO19S_mean"      "alt_mean"         "slope_mean"       "riv_3km_mean"    
    ## [33] "samp_20km_mean"   "lakes_mean"       "BIO1_median"      "BIO2_median"     
    ## [37] "BIO3_median"      "BIO4_median"      "BIO5_median"      "BIO6_median"     
    ## [41] "BIO7_median"      "BIO8S_median"     "BIO9S_median"     "BIO10S_median"   
    ## [45] "BIO11S_median"    "BIO12_median"     "BIO13_median"     "BIO14_median"    
    ## [49] "BIO15_median"     "BIO16S_median"    "BIO17S_median"    "BIO18S_median"   
    ## [53] "BIO19S_median"    "alt_median"       "slope_median"     "riv_3km_median"  
    ## [57] "samp_20km_median" "lakes_median"     "BIO1_mode"        "BIO2_mode"       
    ## [61] "BIO3_mode"        "BIO4_mode"        "BIO5_mode"        "BIO6_mode"       
    ## [65] "BIO7_mode"        "BIO8S_mode"       "BIO9S_mode"       "BIO10S_mode"     
    ## [69] "BIO11S_mode"      "BIO12_mode"       "BIO13_mode"       "BIO14_mode"      
    ## [73] "BIO15_mode"       "BIO16S_mode"      "BIO17S_mode"      "BIO18S_mode"     
    ## [77] "BIO19S_mode"      "alt_mode"         "slope_mode"       "riv_3km_mode"    
    ## [81] "samp_20km_mode"   "lakes_mode"

### Step 3: Prepare pair table, site objects, and transition object

Diagonal movement is allowed in the 16-way transition matrix, and
conductance is corrected for geographic distance, with
geoCorrection(type = “c”), so diagonal steps are appropriately penalized
relative to orthogonal steps.

``` r
# keep only within-cluster comparisons, matching your earlier workflow
cse_df <- cse_df %>%
  filter(Pop1_cluster == Pop2_cluster) %>%
  mutate(
    id  = paste(Var1, Var2, sep = "_"),
    CSE = CSEdistance
  )

# unique site coordinates from both columns
sites_coords <- bind_rows(
  cse_df %>% dplyr::select(Site = Var1, lon = long1, lat = lat1),
  cse_df %>% dplyr::select(Site = Var2, lon = long2, lat = lat2)
) %>%
  distinct()

# convert to sf and then to SpatialPoints in raster CRS
sites_df <- st_as_sf(sites_coords, coords = c("lon", "lat"), crs = crs_geo) %>%
  st_transform(crs(cse_surface))

sites_sp <- as(sites_df, "Spatial")

# lookup index for site names
site_index <- setNames(seq_len(nrow(sites_df)), sites_coords$Site)

# pair table for evaluation
site_pairs <- cse_df %>%
  dplyr::select(Var1, Var2, id, CSE) %>%
  distinct()

# treat projected raw CSE as local resistance
resistance_rast <- cse_surface

# protect against zero or negative values before inversion
min_pos <- min(values(resistance_rast)[values(resistance_rast) > 0], na.rm = TRUE)
resistance_rast[resistance_rast <= 0] <- min_pos

# convert resistance to conductance for gdistance
conductance_rast <- 1 / resistance_rast

# transition object:
# directions = 16 allows diagonal/extended neighborhood movement
tr <- transition(conductance_rast, transitionFunction = mean, directions = 16)

# correction for least-cost paths
# geoCorrection(type = "c") corrects with longer diagonal
tr_lcp <- geoCorrection(tr, type = "c")

# correction for circuit / random-walk distances
# geoCorrection(type = "r") corrects for random walk
tr_cber <- geoCorrection(tr, type = "r")
```

### Step 4: Sanity checks

``` r
nrow(cse_df)
nrow(site_pairs)
nrow(sites_coords)

head(site_pairs)

plot(cse_surface, main = "Projected raw CSE surface")
plot(sites_sp, add = TRUE, pch = 16, cex = 0.5)
```

### Step 5: Compute least-cost paths in parallel (once) and save output

Test one path first:

``` r
# Test one least-cost path before parallel run
i <- site_index[site_pairs$Var1[100]]
j <- site_index[site_pairs$Var2[100]]

test_path <- shortestPath(tr_lcp, sites_sp[i, ], sites_sp[j, ], output = "SpatialLines")

plot(cse_surface, main = "Test least-cost path")
lines(test_path, col = "red", lwd = 2)
points(sites_sp[c(i, j), ], pch = 16)
```

Run full loop:

``` r
output_dir <- file.path(results_dir_big, "least_cost_paths")
dir.create(output_dir, showWarnings = FALSE)

# Set number of cores and register cluster
n_cores <- 16
cl <- makeCluster(n_cores)
registerDoParallel(cl)

paths_list <- foreach(k = 1:nrow(site_pairs),
                      .packages = c("gdistance", "sp", "sf")) %dopar% {

  i <- site_index[site_pairs$Var1[k]]
  j <- site_index[site_pairs$Var2[k]]

  path <- tryCatch({
    shortestPath(tr_lcp, sites_sp[i, ], sites_sp[j, ], output = "SpatialLines")
  }, error = function(e) NULL)

  if (!is.null(path)) {
    path_sf <- st_as_sf(path)
    path_sf$Var1 <- site_pairs$Var1[k]
    path_sf$Var2 <- site_pairs$Var2[k]
    path_sf$id   <- site_pairs$id[k]
    return(path_sf)
  } else {
    return(NULL)
  }
}

stopCluster(cl)

# drop NULLs and combine
paths_list <- paths_list[!sapply(paths_list, is.null)]
paths_sf <- do.call(rbind, paths_list)

# assign CRS and transform back to geographic coordinates for saving
st_crs(paths_sf) <- st_crs(cse_surface)
paths_sf <- st_transform(paths_sf, crs = st_crs(crs_geo))

# join observed CSE values
paths_sf <- left_join(paths_sf, site_pairs, by = "id")

# save paths
st_write(paths_sf,
         file.path(output_dir, "LC_paths_fullRF_rawCSE.shp"),
         delete_layer = TRUE)
```

### Step 7: Extract predicted CSE values along LCPs and calculate mean and sum

``` r
# if needed later, reload:
paths_sf <- st_read(file.path(results_dir_big, "least_cost_paths", "LC_paths_fullRF_rawCSE.shp"), quiet = TRUE)

# transform paths to raster CRS for extraction
paths_extract <- st_transform(paths_sf, crs = crs(cse_surface))

# extract raster values along each path
path_vals <- raster::extract(cse_surface, as(paths_extract, "Spatial"))

# summarize values along each least-cost path
lcp_summary <- data.frame(
  id = paths_extract$id,
  LCP_mean = sapply(path_vals, function(x) {
    if (is.null(x) || all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
  }),
  LCP_sum = sapply(path_vals, function(x) {
    if (is.null(x) || all(is.na(x))) NA_real_ else sum(x, na.rm = TRUE)
  })
)

# join back to pair table
LCP_results <- site_pairs %>%
  left_join(lcp_summary, by = "id")


summary(LCP_results$LCP_mean)
summary(LCP_results$LCP_sum)

# Save LCP output as text file in results directories
write.csv(LCP_results,file.path(results_dir, "spatial_eval_rawCSE_LCPpred.csv"), row.names = FALSE)

write.csv(LCP_results,file.path(results_dir_big, "spatial_eval_rawCSE_LCPpred.csv"), row.names = FALSE)
```

### Step 8: Compute circuit-based effective resistance (CBER)

Test one CBER calculation first:

``` r
# test one CBER calculation first
i <- site_index[site_pairs$Var1[1]]
j <- site_index[site_pairs$Var2[1]]

dmat_test <- commuteDistance(tr_cber, sites_sp[c(i, j), ])

as.numeric(dmat_test)[1] # returns the single number for pair
```

Run as a loop sequentially, also save output from steps 6-8 to file
again…

``` r
library(doParallel)
library(foreach)
library(gdistance)
library(sp)

# number of cores
n_cores <- 16

# on Linux/CHPC, FORK is usually more efficient for large objects like tr_cber
cl <- parallel::makeForkCluster(n_cores)
doParallel::registerDoParallel(cl)

# optional: quick test outside parallel
i <- site_index[site_pairs$Var1[1]]
j <- site_index[site_pairs$Var2[1]]

dmat_test <- commuteDistance(tr_cber, sites_sp[c(i, j), ])
as.numeric(dmat_test)[1]

# parallel loop across all 1091 pairs
cber_results <- foreach(
  k = seq_len(nrow(site_pairs)),
  .combine = rbind,
  .packages = c("gdistance", "sp")
) %dopar% {
  
  site1 <- site_pairs$Var1[k]
  site2 <- site_pairs$Var2[k]
  pair_id <- site_pairs$id[k]
  
  i <- site_index[site1]
  j <- site_index[site2]
  
  # compute commute distance for the pair
  dmat <- commuteDistance(tr_cber, sites_sp[c(i, j), ])
  
  # extract the single numeric value from the dist object
  cber_val <- as.numeric(dmat)[1]
  
  data.frame(
    Var1 = site1,
    Var2 = site2,
    id   = pair_id,
    CBER = cber_val
  )
}

# stop cluster
stopCluster(cl)

# inspect
head(cber_results)
str(cber_results)

# save output in results directories
write.csv(cber_results,file.path(results_dir, "spatial_eval_rawCSE_cber.csv"), row.names = FALSE)

write.csv(cber_results,file.path(results_dir_big, "spatial_eval_rawCSE_cber.csv"), row.names = FALSE)
```

### Step 9: Combine LCP and CBER results

``` r
# load and bind your existing evaluation data
library(dplyr)

# load the LCP sum and mean results you just computed
LCP_results <- read.csv(file.path(results_dir, "spatial_eval_rawCSE_LCPpred.csv"))

# load the CBER results you just computed
cber_results <- read.csv(file.path(results_dir, "spatial_eval_rawCSE_cber.csv"))

# bind CBER and LCP results
eval_results <- LCP_results %>%
  left_join(
    cber_results %>% select(id, CBER),
    by = "id"
  )

# Save results output as RDS for later ease of loading
saveRDS(
  eval_results,
  file.path(results_dir_big, "spatial_eval_rawCSE_predictions.rds")
)
```

### Step 10: Calculate evaluation metrics from eval_results

``` r
library(dplyr)
library(ggplot2)

# Start eval_df from eval_results
eval_df <- eval_results

# Fit linear calibration and return calibrated predictions
calibrate_pred <- function(obs, pred) {
  keep <- complete.cases(obs, pred)
  df <- data.frame(obs = obs[keep], pred = pred[keep])
  
  fit <- lm(obs ~ pred, data = df)
  
  out <- rep(NA_real_, length(obs))
  out[keep] <- predict(fit, newdata = data.frame(pred = pred[keep]))
  out
}

# Create calibrated predictions
eval_df$LCP_sum_cal  <- calibrate_pred(eval_df$CSE, eval_df$LCP_sum)
eval_df$LCP_mean_cal <- calibrate_pred(eval_df$CSE, eval_df$LCP_mean)
eval_df$CBER_cal     <- calibrate_pred(eval_df$CSE, eval_df$CBER)

# Plots: raw predicted vs observed, with fitted line
plot1 <- ggplot(eval_df, aes(x = LCP_sum, y = CSE)) +
  geom_point() +
  geom_smooth(method = "lm") +
  labs(
    x = "Predicted CSE (LCP sum, raw scale)",
    y = "Observed CSE"
  )
plot1

plot2 <- ggplot(eval_df, aes(x = LCP_mean, y = CSE)) +
  geom_point() +
  geom_smooth(method = "lm") +
  labs(
    x = "Predicted CSE (LCP mean, raw scale)",
    y = "Observed CSE"
  )
plot2

plot3 <- ggplot(eval_df, aes(x = CBER, y = CSE)) +
  geom_point() +
  geom_smooth(method = "lm") +
  labs(
    x = "Predicted CSE (CBER, raw scale)",
    y = "Observed CSE"
  )
plot3

# Calculate formal metrics
# Correlation is calculated on raw predictions because it is scale-invariant
# Error metrics and R2 are calculated on calibrated predictions
calc_metrics <- function(obs, pred_raw, pred_cal) {
  keep <- complete.cases(obs, pred_raw, pred_cal)
  obs      <- obs[keep]
  pred_raw <- pred_raw[keep]
  pred_cal <- pred_cal[keep]
  
  sse  <- sum((obs - pred_cal)^2)
  sst  <- sum((obs - mean(obs))^2)
  mse  <- mean((obs - pred_cal)^2)
  rmse <- sqrt(mse)
  mae  <- mean(abs(obs - pred_cal))
  rsq  <- 1 - (sse / sst)
  corr <- cor(obs, pred_raw)
  
  data.frame(
    n = length(obs),
    MSE = mse,
    RMSE = rmse,
    MAE = mae,
    RSQ = rsq,
    Correlation = corr
  )
}

metrics_df <- bind_rows(
  calc_metrics(eval_df$CSE, eval_df$LCP_mean, eval_df$LCP_mean_cal) %>%
    mutate(method = "LCP_mean"),
  
  calc_metrics(eval_df$CSE, eval_df$LCP_sum, eval_df$LCP_sum_cal) %>%
    mutate(method = "LCP_sum"),
  
  calc_metrics(eval_df$CSE, eval_df$CBER, eval_df$CBER_cal) %>%
    mutate(method = "CBER")
) %>%
  dplyr::select(method, everything())

# Print metrics
print(metrics_df)

# Save outputs
saveRDS(
  metrics_df,
  file.path(results_dir_big, "spatial_eval_rawCSE_eval_metrics.rds")
)

write.csv(
  metrics_df,
  file.path(results_dir_big, "spatial_eval_rawCSE_eval_metrics.csv"),
  row.names = FALSE
)

write.csv(
  metrics_df,
  file.path(results_dir, "spatial_eval_rawCSE_eval_metrics.csv"),
  row.names = FALSE
)
```

# IBD residuals

### Step 1: Inputs

``` r
# Define Paths to directories
input_dir <- "../input"
data_dir  <- "/uufs/chpc.utah.edu/common/home/saarman-group1/uganda-tsetse-LG/data"
results_dir_big <- "/uufs/chpc.utah.edu/common/home/saarman-group1/uganda-tsetse-LG/results/"
results_dir <- "../results/"

# define coordinate reference system
crs_geo <- 4326     # EPSG code for WGS84

######################################
# What surface are you using?
surf_file  <- file.path(results_dir_big, "fullRF_residuals.tif")
########################################

lake_file  <- file.path(data_dir, "processed", "lake_binary.tif")

cse_file   <- file.path("..", "input", "Gff_cse_envCostPaths.csv")

# Read pairwise CSE table
cse_df <- read.csv(cse_file)
# This was added only after completing LOPOCV...
# Filter out western outlier "50-KB" 
cse_df <- cse_df %>%
  filter(Var1 != "50-KB", Var2 != "50-KB")


# Read projected IBD residuals surface
cse_surface <- raster(surf_file)

# Read lake mask
lake_mask <- raster(lake_file)
```

### Step 2: Check raster alignment and assign high cost to water

``` r
if (!compareRaster(cse_surface, lake_mask,
                   extent = TRUE, rowcol = TRUE,
                   crs = TRUE, res = TRUE,
                   stopiffalse = FALSE)) {
  lake_mask <- projectRaster(lake_mask, cse_surface, method = "ngb")
}

# Inspect projected IBD residuals values before altering water cells
summary(values(cse_surface))
```

    ##     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
    ## -0.13626 -0.00128  0.01825  0.01800  0.04173  0.13413

``` r
# Keep water traversible but very costly (max value from model)
cse_surface[lake_mask[] == 1] <- max(values(cse_surface))

# Check again after assigning lakes
summary(values(cse_surface))
```

    ##      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
    ## -0.133161  0.005422  0.024874  0.036160  0.057104  0.134129

``` r
# Quick plot
plot(cse_surface, main = "Projected IBD residuals surface with water = high cost")
```

![](/tmp/RtmpaVWtN9/knitted_mds/resid-check-alignment-1.png)<!-- -->

``` r
# Check column names of cse_df
names(cse_df)
```

    ##  [1] "Var1"             "Var2"             "CSEdistance"      "lat1"            
    ##  [5] "long1"            "lat2"             "long2"            "Pop1_cluster"    
    ##  [9] "Pop2_cluster"     "pix_dist"         "BIO1_mean"        "BIO2_mean"       
    ## [13] "BIO3_mean"        "BIO4_mean"        "BIO5_mean"        "BIO6_mean"       
    ## [17] "BIO7_mean"        "BIO8S_mean"       "BIO9S_mean"       "BIO10S_mean"     
    ## [21] "BIO11S_mean"      "BIO12_mean"       "BIO13_mean"       "BIO14_mean"      
    ## [25] "BIO15_mean"       "BIO16S_mean"      "BIO17S_mean"      "BIO18S_mean"     
    ## [29] "BIO19S_mean"      "alt_mean"         "slope_mean"       "riv_3km_mean"    
    ## [33] "samp_20km_mean"   "lakes_mean"       "BIO1_median"      "BIO2_median"     
    ## [37] "BIO3_median"      "BIO4_median"      "BIO5_median"      "BIO6_median"     
    ## [41] "BIO7_median"      "BIO8S_median"     "BIO9S_median"     "BIO10S_median"   
    ## [45] "BIO11S_median"    "BIO12_median"     "BIO13_median"     "BIO14_median"    
    ## [49] "BIO15_median"     "BIO16S_median"    "BIO17S_median"    "BIO18S_median"   
    ## [53] "BIO19S_median"    "alt_median"       "slope_median"     "riv_3km_median"  
    ## [57] "samp_20km_median" "lakes_median"     "BIO1_mode"        "BIO2_mode"       
    ## [61] "BIO3_mode"        "BIO4_mode"        "BIO5_mode"        "BIO6_mode"       
    ## [65] "BIO7_mode"        "BIO8S_mode"       "BIO9S_mode"       "BIO10S_mode"     
    ## [69] "BIO11S_mode"      "BIO12_mode"       "BIO13_mode"       "BIO14_mode"      
    ## [73] "BIO15_mode"       "BIO16S_mode"      "BIO17S_mode"      "BIO18S_mode"     
    ## [77] "BIO19S_mode"      "alt_mode"         "slope_mode"       "riv_3km_mode"    
    ## [81] "samp_20km_mode"   "lakes_mode"

### Step 3: Prepare pair table, site objects, and transition object

Diagonal movement is allowed in the 16-way transition matrix, and
conductance is corrected for geographic distance, with
geoCorrection(type = “c”), so diagonal steps are appropriately penalized
relative to orthogonal steps.

``` r
# keep only within-cluster comparisons, matching your earlier workflow
cse_df <- cse_df %>%
  filter(Pop1_cluster == Pop2_cluster) %>%
  mutate(
    id  = paste(Var1, Var2, sep = "_"),
    CSE = CSEdistance
  )

# unique site coordinates from both columns
sites_coords <- bind_rows(
  cse_df %>% dplyr::select(Site = Var1, lon = long1, lat = lat1),
  cse_df %>% dplyr::select(Site = Var2, lon = long2, lat = lat2)
) %>%
  distinct()

# convert to sf and then to SpatialPoints in raster CRS
sites_df <- st_as_sf(sites_coords, coords = c("lon", "lat"), crs = crs_geo) %>%
  st_transform(crs(cse_surface))

sites_sp <- as(sites_df, "Spatial")

# lookup index for site names
site_index <- setNames(seq_len(nrow(sites_df)), sites_coords$Site)

# pair table for evaluation
site_pairs <- cse_df %>%
  dplyr::select(Var1, Var2, id, CSE) %>%
  distinct()

# treat projected IBD residuals as local resistance
resistance_rast <- cse_surface

# protect against zero or negative values before inversion
min_pos <- min(values(resistance_rast)[values(resistance_rast) > 0], na.rm = TRUE)
resistance_rast[resistance_rast <= 0] <- min_pos

# convert resistance to conductance for gdistance
conductance_rast <- 1 / resistance_rast

# transition object:
# directions = 16 allows diagonal/extended neighborhood movement
tr <- transition(conductance_rast, transitionFunction = mean, directions = 16)

# correction for least-cost paths
# geoCorrection(type = "c") corrects with longer diagonal
tr_lcp <- geoCorrection(tr, type = "c")

# correction for circuit / random-walk distances
# geoCorrection(type = "r") corrects for random walk
tr_cber <- geoCorrection(tr, type = "r")
```

### Step 4: Sanity checks

``` r
nrow(cse_df)
nrow(site_pairs)
nrow(sites_coords)

head(site_pairs)

plot(cse_surface, main = "Projected IBD residuals surface")
plot(sites_sp, add = TRUE, pch = 16, cex = 0.5)
```

### Step 5: Compute least-cost paths in parallel (once) and save output

Test one path first:

``` r
# Test one least-cost path before parallel run
i <- site_index[site_pairs$Var1[100]]
j <- site_index[site_pairs$Var2[100]]

test_path <- shortestPath(tr_lcp, sites_sp[i, ], sites_sp[j, ], output = "SpatialLines")

plot(cse_surface, main = "Test least-cost path")
lines(test_path, col = "red", lwd = 2)
points(sites_sp[c(i, j), ], pch = 16)
```

Run full loop:

``` r
output_dir <- file.path(results_dir_big, "least_cost_paths")
dir.create(output_dir, showWarnings = FALSE)

# Set number of cores and register cluster
n_cores <- 16
cl <- makeCluster(n_cores)
registerDoParallel(cl)

paths_list <- foreach(k = 1:nrow(site_pairs),
                      .packages = c("gdistance", "sp", "sf")) %dopar% {

  i <- site_index[site_pairs$Var1[k]]
  j <- site_index[site_pairs$Var2[k]]

  path <- tryCatch({
    shortestPath(tr_lcp, sites_sp[i, ], sites_sp[j, ], output = "SpatialLines")
  }, error = function(e) NULL)

  if (!is.null(path)) {
    path_sf <- st_as_sf(path)
    path_sf$Var1 <- site_pairs$Var1[k]
    path_sf$Var2 <- site_pairs$Var2[k]
    path_sf$id   <- site_pairs$id[k]
    return(path_sf)
  } else {
    return(NULL)
  }
}

stopCluster(cl)

# drop NULLs and combine
paths_list <- paths_list[!sapply(paths_list, is.null)]
paths_sf <- do.call(rbind, paths_list)

# assign CRS and transform back to geographic coordinates for saving
st_crs(paths_sf) <- st_crs(cse_surface)
paths_sf <- st_transform(paths_sf, crs = st_crs(crs_geo))

# join observed CSE values
paths_sf <- left_join(paths_sf, site_pairs, by = "id")

# save paths
st_write(paths_sf,
         file.path(output_dir, "LC_paths_fullRF_residuals.shp"),
         delete_layer = TRUE)
```

### Step 7: Extract predicted CSE values along LCPs and calculate mean and sum

``` r
# if needed later, reload:
paths_sf <- st_read(file.path(results_dir_big, "least_cost_paths", "LC_paths_fullRF_residuals.shp"), quiet = TRUE)

# transform paths to raster CRS for extraction
paths_extract <- st_transform(paths_sf, crs = crs(cse_surface))

# extract raster values along each path
path_vals <- raster::extract(cse_surface, as(paths_extract, "Spatial"))

# summarize values along each least-cost path
lcp_summary <- data.frame(
  id = paths_extract$id,
  LCP_mean = sapply(path_vals, function(x) {
    if (is.null(x) || all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
  }),
  LCP_sum = sapply(path_vals, function(x) {
    if (is.null(x) || all(is.na(x))) NA_real_ else sum(x, na.rm = TRUE)
  })
)

# join back to pair table
LCP_results <- site_pairs %>%
  left_join(lcp_summary, by = "id")

head(LCP_results)
summary(LCP_results$LCP_mean)
summary(LCP_results$LCP_sum)

# Save LCP output as text file in results directories
write.csv(LCP_results,file.path(results_dir, "spatial_eval_residuals_LCPpred.csv"), row.names = FALSE)

write.csv(LCP_results,file.path(results_dir_big, "spatial_eval_residuals_LCPpred.csv"), row.names = FALSE)
```

### Step 8: Compute circuit-based effective resistance (CBER)

Test one CBER calculation first:

``` r
# test one CBER calculation first
i <- site_index[site_pairs$Var1[1]]
j <- site_index[site_pairs$Var2[1]]

dmat_test <- commuteDistance(tr_cber, sites_sp[c(i, j), ])

as.numeric(dmat_test)[1] # returns the single number for pair
```

Run as a loop sequentially, also save output from steps 6-8 to file
again…

``` r
library(doParallel)
library(foreach)
library(gdistance)
library(sp)

# number of cores
n_cores <- 16

# on Linux/CHPC, FORK is usually more efficient for large objects like tr_cber
cl <- parallel::makeForkCluster(n_cores)
doParallel::registerDoParallel(cl)

# optional: quick test outside parallel
i <- site_index[site_pairs$Var1[1]]
j <- site_index[site_pairs$Var2[1]]

dmat_test <- commuteDistance(tr_cber, sites_sp[c(i, j), ])
as.numeric(dmat_test)[1]

# parallel loop across all 1091 pairs
cber_results <- foreach(
  k = seq_len(nrow(site_pairs)),
  .combine = rbind,
  .packages = c("gdistance", "sp")
) %dopar% {
  
  site1 <- site_pairs$Var1[k]
  site2 <- site_pairs$Var2[k]
  pair_id <- site_pairs$id[k]
  
  i <- site_index[site1]
  j <- site_index[site2]
  
  # compute commute distance for the pair
  dmat <- commuteDistance(tr_cber, sites_sp[c(i, j), ])
  
  # extract the single numeric value from the dist object
  cber_val <- as.numeric(dmat)[1]
  
  data.frame(
    Var1 = site1,
    Var2 = site2,
    id   = pair_id,
    CBER = cber_val
  )
}

# stop cluster
stopCluster(cl)

# inspect
head(cber_results)
str(cber_results)

# save output in results directories
write.csv(cber_results,file.path(results_dir, "spatial_eval_residuals_cber.csv"), row.names = FALSE)

write.csv(cber_results,file.path(results_dir_big, "spatial_eval_residuals_cber.csv"), row.names = FALSE)
```

### Step 9: Combine LCP and CBER results

``` r
# load and bind your existing evaluation data
library(dplyr)

# load the LCP sum and mean results you just computed
LCP_results <- read.csv(file.path(results_dir, "spatial_eval_residuals_LCPpred.csv"))

# bind CBER and LCP results
eval_results <- LCP_results %>%
  left_join(
    cber_results %>% select(id, CBER),
    by = "id"
  )

# Save results output as RDS for later ease of loading
saveRDS(
  eval_results,
  file.path(results_dir_big, "spatial_eval_residuals_predictions.rds")
)
```

### Step 10: Calculate evaluation metrics from eval_results

``` r
library(dplyr)
library(ggplot2)

# Start eval_df from eval_results
eval_df <- eval_results

# Fit linear calibration and return calibrated predictions
calibrate_pred <- function(obs, pred) {
  keep <- complete.cases(obs, pred)
  df <- data.frame(obs = obs[keep], pred = pred[keep])
  
  fit <- lm(obs ~ pred, data = df)
  
  out <- rep(NA_real_, length(obs))
  out[keep] <- predict(fit, newdata = data.frame(pred = pred[keep]))
  out
}

# Create calibrated predictions
eval_df$LCP_sum_cal  <- calibrate_pred(eval_df$CSE, eval_df$LCP_sum)
eval_df$LCP_mean_cal <- calibrate_pred(eval_df$CSE, eval_df$LCP_mean)
eval_df$CBER_cal     <- calibrate_pred(eval_df$CSE, eval_df$CBER)

# Plots: raw predicted vs observed, with fitted line
plot1 <- ggplot(eval_df, aes(x = LCP_sum, y = CSE)) +
  geom_point() +
  geom_smooth(method = "lm") +
  labs(
    x = "Predicted CSE (LCP sum, raw scale)",
    y = "Observed CSE"
  )
plot1

plot2 <- ggplot(eval_df, aes(x = LCP_mean, y = CSE)) +
  geom_point() +
  geom_smooth(method = "lm") +
  labs(
    x = "Predicted CSE (LCP mean, raw scale)",
    y = "Observed CSE"
  )
plot2

plot3 <- ggplot(eval_df, aes(x = CBER, y = CSE)) +
  geom_point() +
  geom_smooth(method = "lm") +
  labs(
    x = "Predicted CSE (CBER, raw scale)",
    y = "Observed CSE"
  )
plot3

# Calculate formal metrics
# Correlation is calculated on raw predictions because it is scale-invariant
# Error metrics and R2 are calculated on calibrated predictions
calc_metrics <- function(obs, pred_raw, pred_cal) {
  keep <- complete.cases(obs, pred_raw, pred_cal)
  obs      <- obs[keep]
  pred_raw <- pred_raw[keep]
  pred_cal <- pred_cal[keep]
  
  sse  <- sum((obs - pred_cal)^2)
  sst  <- sum((obs - mean(obs))^2)
  mse  <- mean((obs - pred_cal)^2)
  rmse <- sqrt(mse)
  mae  <- mean(abs(obs - pred_cal))
  rsq  <- 1 - (sse / sst)
  corr <- cor(obs, pred_raw)
  
  data.frame(
    n = length(obs),
    MSE = mse,
    RMSE = rmse,
    MAE = mae,
    RSQ = rsq,
    Correlation = corr
  )
}

metrics_df <- bind_rows(
  calc_metrics(eval_df$CSE, eval_df$LCP_mean, eval_df$LCP_mean_cal) %>%
    mutate(method = "LCP_mean"),
  
  calc_metrics(eval_df$CSE, eval_df$LCP_sum, eval_df$LCP_sum_cal) %>%
    mutate(method = "LCP_sum"),
  
  calc_metrics(eval_df$CSE, eval_df$CBER, eval_df$CBER_cal) %>%
    mutate(method = "CBER")
) %>%
  dplyr::select(method, everything())

# Print metrics
print(metrics_df)

# Save outputs
saveRDS(
  metrics_df,
  file.path(results_dir_big, "spatial_eval_residuals_eval_metrics.rds")
)

write.csv(
  metrics_df,
  file.path(results_dir_big, "spatial_eval_residuals_eval_metrics.csv"),
  row.names = FALSE
)

write.csv(
  metrics_df,
  file.path(results_dir, "spatial_eval_residuals_eval_metrics.csv"),
  row.names = FALSE
)
```

# CSE per km

### Step 1: Inputs

``` r
# Define Paths to directories
input_dir <- "../input"
data_dir  <- "/uufs/chpc.utah.edu/common/home/saarman-group1/uganda-tsetse-LG/data"
results_dir_big <- "/uufs/chpc.utah.edu/common/home/saarman-group1/uganda-tsetse-LG/results/"
results_dir <- "../results/"

# define coordinate reference system
crs_geo <- 4326     # EPSG code for WGS84

######################################
# What surface are you using?
surf_file  <- file.path(results_dir_big, "fullRF_rate.tif")
########################################

lake_file  <- file.path(data_dir, "processed", "lake_binary.tif")

cse_file   <- file.path("..", "input", "Gff_cse_envCostPaths.csv")

# Read pairwise CSE table
cse_df <- read.csv(cse_file)
# This was added only after completing LOPOCV...
# Filter out western outlier "50-KB" 
cse_df <- cse_df %>%
  filter(Var1 != "50-KB", Var2 != "50-KB")


# Read projected CSE-per-km surface
cse_surface <- raster(surf_file)

# Read lake mask
lake_mask <- raster(lake_file)
```

### Step 2: Check raster alignment and assign high cost to water

``` r
if (!compareRaster(cse_surface, lake_mask,
                   extent = TRUE, rowcol = TRUE,
                   crs = TRUE, res = TRUE,
                   stopiffalse = FALSE)) {
  lake_mask <- projectRaster(lake_mask, cse_surface, method = "ngb")
}

# Inspect projected CSE-per-km values before altering water cells
summary(values(cse_surface))
```

    ##     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
    ## 0.001359 0.005825 0.014867 0.015209 0.022207 0.042114

``` r
# Keep water traversible but very costly (max value from model)
cse_surface[lake_mask[] == 1] <- max(values(cse_surface))

# Check again after assigning lakes
summary(values(cse_surface))
```

    ##     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
    ## 0.001359 0.005825 0.014867 0.017235 0.028132 0.042114

``` r
# Quick plot
plot(cse_surface, main = "Projected CSE-per-km surface with water = high cost")
```

![](/tmp/RtmpaVWtN9/knitted_mds/rate-check-alignment-1.png)<!-- -->

``` r
# Check column names of cse_df
names(cse_df)
```

    ##  [1] "Var1"             "Var2"             "CSEdistance"      "lat1"            
    ##  [5] "long1"            "lat2"             "long2"            "Pop1_cluster"    
    ##  [9] "Pop2_cluster"     "pix_dist"         "BIO1_mean"        "BIO2_mean"       
    ## [13] "BIO3_mean"        "BIO4_mean"        "BIO5_mean"        "BIO6_mean"       
    ## [17] "BIO7_mean"        "BIO8S_mean"       "BIO9S_mean"       "BIO10S_mean"     
    ## [21] "BIO11S_mean"      "BIO12_mean"       "BIO13_mean"       "BIO14_mean"      
    ## [25] "BIO15_mean"       "BIO16S_mean"      "BIO17S_mean"      "BIO18S_mean"     
    ## [29] "BIO19S_mean"      "alt_mean"         "slope_mean"       "riv_3km_mean"    
    ## [33] "samp_20km_mean"   "lakes_mean"       "BIO1_median"      "BIO2_median"     
    ## [37] "BIO3_median"      "BIO4_median"      "BIO5_median"      "BIO6_median"     
    ## [41] "BIO7_median"      "BIO8S_median"     "BIO9S_median"     "BIO10S_median"   
    ## [45] "BIO11S_median"    "BIO12_median"     "BIO13_median"     "BIO14_median"    
    ## [49] "BIO15_median"     "BIO16S_median"    "BIO17S_median"    "BIO18S_median"   
    ## [53] "BIO19S_median"    "alt_median"       "slope_median"     "riv_3km_median"  
    ## [57] "samp_20km_median" "lakes_median"     "BIO1_mode"        "BIO2_mode"       
    ## [61] "BIO3_mode"        "BIO4_mode"        "BIO5_mode"        "BIO6_mode"       
    ## [65] "BIO7_mode"        "BIO8S_mode"       "BIO9S_mode"       "BIO10S_mode"     
    ## [69] "BIO11S_mode"      "BIO12_mode"       "BIO13_mode"       "BIO14_mode"      
    ## [73] "BIO15_mode"       "BIO16S_mode"      "BIO17S_mode"      "BIO18S_mode"     
    ## [77] "BIO19S_mode"      "alt_mode"         "slope_mode"       "riv_3km_mode"    
    ## [81] "samp_20km_mode"   "lakes_mode"

### Step 3: Prepare pair table, site objects, and transition object

Diagonal movement is allowed in the 16-way transition matrix, and
conductance is corrected for geographic distance, with
geoCorrection(type = “c”), so diagonal steps are appropriately penalized
relative to orthogonal steps.

``` r
# keep only within-cluster comparisons, matching your earlier workflow
cse_df <- cse_df %>%
  filter(Pop1_cluster == Pop2_cluster) %>%
  mutate(
    id  = paste(Var1, Var2, sep = "_"),
    CSE = CSEdistance
  )

# unique site coordinates from both columns
sites_coords <- bind_rows(
  cse_df %>% dplyr::select(Site = Var1, lon = long1, lat = lat1),
  cse_df %>% dplyr::select(Site = Var2, lon = long2, lat = lat2)
) %>%
  distinct()

# convert to sf and then to SpatialPoints in raster CRS
sites_df <- st_as_sf(sites_coords, coords = c("lon", "lat"), crs = crs_geo) %>%
  st_transform(crs(cse_surface))

sites_sp <- as(sites_df, "Spatial")

# lookup index for site names
site_index <- setNames(seq_len(nrow(sites_df)), sites_coords$Site)

# pair table for evaluation
site_pairs <- cse_df %>%
  dplyr::select(Var1, Var2, id, CSE) %>%
  distinct()

# treat projected CSE-per-km as local resistance
resistance_rast <- cse_surface

# protect against zero or negative values before inversion
min_pos <- min(values(resistance_rast)[values(resistance_rast) > 0], na.rm = TRUE)
resistance_rast[resistance_rast <= 0] <- min_pos

# convert resistance to conductance for gdistance
conductance_rast <- 1 / resistance_rast

# transition object:
# directions = 16 allows diagonal/extended neighborhood movement
tr <- transition(conductance_rast, transitionFunction = mean, directions = 16)

# correction for least-cost paths
# geoCorrection(type = "c") corrects with longer diagonal
tr_lcp <- geoCorrection(tr, type = "c")

# correction for circuit / random-walk distances
# geoCorrection(type = "r") corrects for random walk
tr_cber <- geoCorrection(tr, type = "r")
```

### Step 4: Sanity checks

``` r
nrow(cse_df)
nrow(site_pairs)
nrow(sites_coords)

head(site_pairs)

plot(cse_surface, main = "Projected CSE-per-km surface")
plot(sites_sp, add = TRUE, pch = 16, cex = 0.5)
```

### Step 5: Compute least-cost paths in parallel (once) and save output

Test one path first:

``` r
# Test one least-cost path before parallel run
i <- site_index[site_pairs$Var1[100]]
j <- site_index[site_pairs$Var2[100]]

test_path <- shortestPath(tr_lcp, sites_sp[i, ], sites_sp[j, ], output = "SpatialLines")

plot(cse_surface, main = "Test least-cost path")
lines(test_path, col = "red", lwd = 2)
points(sites_sp[c(i, j), ], pch = 16)
```

Run full loop:

``` r
output_dir <- file.path(results_dir_big, "least_cost_paths")
dir.create(output_dir, showWarnings = FALSE)

# Set number of cores and register cluster
n_cores <- 16
cl <- makeCluster(n_cores)
registerDoParallel(cl)

paths_list <- foreach(k = 1:nrow(site_pairs),
                      .packages = c("gdistance", "sp", "sf")) %dopar% {

  i <- site_index[site_pairs$Var1[k]]
  j <- site_index[site_pairs$Var2[k]]

  path <- tryCatch({
    shortestPath(tr_lcp, sites_sp[i, ], sites_sp[j, ], output = "SpatialLines")
  }, error = function(e) NULL)

  if (!is.null(path)) {
    path_sf <- st_as_sf(path)
    path_sf$Var1 <- site_pairs$Var1[k]
    path_sf$Var2 <- site_pairs$Var2[k]
    path_sf$id   <- site_pairs$id[k]
    return(path_sf)
  } else {
    return(NULL)
  }
}

stopCluster(cl)

# drop NULLs and combine
paths_list <- paths_list[!sapply(paths_list, is.null)]
paths_sf <- do.call(rbind, paths_list)

# assign CRS and transform back to geographic coordinates for saving
st_crs(paths_sf) <- st_crs(cse_surface)
paths_sf <- st_transform(paths_sf, crs = st_crs(crs_geo))

# join observed CSE values
paths_sf <- left_join(paths_sf, site_pairs, by = "id")

# save paths
st_write(paths_sf,
         file.path(output_dir, "LC_paths_fullRF_rate.shp"),
         delete_layer = TRUE)
```

### Step 7: Extract predicted CSE values along LCPs and calculate mean and sum

``` r
# if needed later, reload:
paths_sf <- st_read(file.path(results_dir_big, "least_cost_paths", "LC_paths_fullRF_rate.shp"), quiet = TRUE)

# transform paths to raster CRS for extraction
paths_extract <- st_transform(paths_sf, crs = crs(cse_surface))

# extract raster values along each path
path_vals <- raster::extract(cse_surface, as(paths_extract, "Spatial"))

# summarize values along each least-cost path
lcp_summary <- data.frame(
  id = paths_extract$id,
  LCP_mean = sapply(path_vals, function(x) {
    if (is.null(x) || all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
  }),
  LCP_sum = sapply(path_vals, function(x) {
    if (is.null(x) || all(is.na(x))) NA_real_ else sum(x, na.rm = TRUE)
  })
)

# join back to pair table
LCP_results <- site_pairs %>%
  left_join(lcp_summary, by = "id")

head(LCP_results)
summary(LCP_results$LCP_mean)
summary(LCP_results$LCP_sum)

# Save LCP output as text file in results directories
write.csv(LCP_results,file.path(results_dir, "spatial_eval_rate_LCPpred.csv"), row.names = FALSE)

write.csv(LCP_results,file.path(results_dir_big, "spatial_eval_rate_LCPpred.csv"), row.names = FALSE)
```

### Step 8: Compute circuit-based effective resistance (CBER)

Test one CBER calculation first:

``` r
# test one CBER calculation first
i <- site_index[site_pairs$Var1[1]]
j <- site_index[site_pairs$Var2[1]]

dmat_test <- commuteDistance(tr_cber, sites_sp[c(i, j), ])

as.numeric(dmat_test)[1] # returns the single number for pair
```

Run as a loop sequentially, also save output from steps 6-8 to file
again…

``` r
library(doParallel)
library(foreach)
library(gdistance)
library(sp)

# number of cores
n_cores <- 16

# on Linux/CHPC, FORK is usually more efficient for large objects like tr_cber
cl <- parallel::makeForkCluster(n_cores)
doParallel::registerDoParallel(cl)

# optional: quick test outside parallel
i <- site_index[site_pairs$Var1[1]]
j <- site_index[site_pairs$Var2[1]]

dmat_test <- commuteDistance(tr_cber, sites_sp[c(i, j), ])
as.numeric(dmat_test)[1]

# parallel loop across all 1091 pairs
cber_results <- foreach(
  k = seq_len(nrow(site_pairs)),
  .combine = rbind,
  .packages = c("gdistance", "sp")
) %dopar% {
  
  site1 <- site_pairs$Var1[k]
  site2 <- site_pairs$Var2[k]
  pair_id <- site_pairs$id[k]
  
  i <- site_index[site1]
  j <- site_index[site2]
  
  # compute commute distance for the pair
  dmat <- commuteDistance(tr_cber, sites_sp[c(i, j), ])
  
  # extract the single numeric value from the dist object
  cber_val <- as.numeric(dmat)[1]
  
  data.frame(
    Var1 = site1,
    Var2 = site2,
    id   = pair_id,
    CBER = cber_val
  )
}

# stop cluster
stopCluster(cl)

# inspect
head(cber_results)
str(cber_results)

# save output in results directories
write.csv(cber_results,file.path(results_dir, "spatial_eval_rate_cber.csv"), row.names = FALSE)

write.csv(cber_results,file.path(results_dir_big, "spatial_eval_rate_cber.csv"), row.names = FALSE)
```

### Step 9: Combine LCP and CBER results

``` r
# load and bind your existing evaluation data
library(dplyr)

# load the LCP sum and mean results you just computed
LCP_results <- read.csv(file.path(results_dir, "spatial_eval_rate_LCPpred.csv"))

# bind CBER and LCP results
eval_results <- LCP_results %>%
  left_join(
    cber_results %>% select(id, CBER),
    by = "id"
  )

# Save results output as RDS for later ease of loading
saveRDS(
  eval_results,
  file.path(results_dir_big, "spatial_eval_rate_predictions.rds")
)
```

### Step 10: Calculate evaluation metrics from eval_results

``` r
library(dplyr)
library(ggplot2)

# Start eval_df from eval_results
eval_df <- eval_results

# Fit linear calibration and return calibrated predictions
calibrate_pred <- function(obs, pred) {
  keep <- complete.cases(obs, pred)
  df <- data.frame(obs = obs[keep], pred = pred[keep])
  
  fit <- lm(obs ~ pred, data = df)
  
  out <- rep(NA_real_, length(obs))
  out[keep] <- predict(fit, newdata = data.frame(pred = pred[keep]))
  out
}

# Create calibrated predictions
eval_df$LCP_sum_cal  <- calibrate_pred(eval_df$CSE, eval_df$LCP_sum)
eval_df$LCP_mean_cal <- calibrate_pred(eval_df$CSE, eval_df$LCP_mean)
eval_df$CBER_cal     <- calibrate_pred(eval_df$CSE, eval_df$CBER)

# Plots: raw predicted vs observed, with fitted line
plot1 <- ggplot(eval_df, aes(x = LCP_sum, y = CSE)) +
  geom_point() +
  geom_smooth(method = "lm") +
  labs(
    x = "Predicted CSE (LCP sum, raw scale)",
    y = "Observed CSE"
  )
plot1

plot2 <- ggplot(eval_df, aes(x = LCP_mean, y = CSE)) +
  geom_point() +
  geom_smooth(method = "lm") +
  labs(
    x = "Predicted CSE (LCP mean, raw scale)",
    y = "Observed CSE"
  )
plot2

plot3 <- ggplot(eval_df, aes(x = CBER, y = CSE)) +
  geom_point() +
  geom_smooth(method = "lm") +
  labs(
    x = "Predicted CSE (CBER, raw scale)",
    y = "Observed CSE"
  )
plot3

# Calculate formal metrics
# Correlation is calculated on raw predictions because it is scale-invariant
# Error metrics and R2 are calculated on calibrated predictions
calc_metrics <- function(obs, pred_raw, pred_cal) {
  keep <- complete.cases(obs, pred_raw, pred_cal)
  obs      <- obs[keep]
  pred_raw <- pred_raw[keep]
  pred_cal <- pred_cal[keep]
  
  sse  <- sum((obs - pred_cal)^2)
  sst  <- sum((obs - mean(obs))^2)
  mse  <- mean((obs - pred_cal)^2)
  rmse <- sqrt(mse)
  mae  <- mean(abs(obs - pred_cal))
  rsq  <- 1 - (sse / sst)
  corr <- cor(obs, pred_raw)
  
  data.frame(
    n = length(obs),
    MSE = mse,
    RMSE = rmse,
    MAE = mae,
    RSQ = rsq,
    Correlation = corr
  )
}

metrics_df <- bind_rows(
  calc_metrics(eval_df$CSE, eval_df$LCP_mean, eval_df$LCP_mean_cal) %>%
    mutate(method = "LCP_mean"),
  
  calc_metrics(eval_df$CSE, eval_df$LCP_sum, eval_df$LCP_sum_cal) %>%
    mutate(method = "LCP_sum"),
  
  calc_metrics(eval_df$CSE, eval_df$CBER, eval_df$CBER_cal) %>%
    mutate(method = "CBER")
) %>%
  dplyr::select(method, everything())

# Print metrics
print(metrics_df)

# Save outputs
saveRDS(
  metrics_df,
  file.path(results_dir_big, "spatial_eval_rate_eval_metrics.rds")
)

write.csv(
  metrics_df,
  file.path(results_dir_big, "spatial_eval_rate_eval_metrics.csv"),
  row.names = FALSE
)

write.csv(
  metrics_df,
  file.path(results_dir, "spatial_eval_rate_eval_metrics.csv"),
  row.names = FALSE
)
```

# Screening results

This section reads the existing screening outputs for the three response
formulations and compares the spatial summaries. It does **not** save or
overwrite analysis outputs.

``` r
library(dplyr)
library(ggplot2)
library(knitr)

results_dir_big <- "/uufs/chpc.utah.edu/common/home/saarman-group1/uganda-tsetse-LG/results/"

result_sets <- list(
  `Raw CSE` = "rawCSE",
  `IBD residuals` = "residuals",
  `CSE per km` = "rate"
)

cols <- c(
  "Raw CSE" = "#FF2DA0",
  "IBD residuals" = "#1bc8ea",
  "CSE per km" = "grey60"
)

read_metrics <- function(label, stem) {
  f <- file.path(results_dir_big, paste0("spatial_eval_", stem, "_eval_metrics.csv"))
  if (!file.exists(f)) stop("Missing existing result file: ", f)
  read.csv(f) %>% mutate(response = label, .before = 1)
}

metrics_all <- bind_rows(Map(read_metrics, names(result_sets), unname(result_sets)))
metrics_all
```

    ##        response   method    n         MSE       RMSE        MAE         RSQ
    ## 1       Raw CSE LCP_mean 1091 0.007288414 0.08537221 0.06685447 0.101254278
    ## 2       Raw CSE  LCP_sum 1091 0.003114933 0.05581158 0.04540066 0.615892809
    ## 3       Raw CSE     CBER 1091 0.003508397 0.05923172 0.04633363 0.567374124
    ## 4 IBD residuals LCP_mean 1091 0.008004301 0.08946676 0.07061746 0.012977113
    ## 5 IBD residuals  LCP_sum 1091 0.008085445 0.08991910 0.07127313 0.002971155
    ## 6 IBD residuals     CBER 1091 0.007085320 0.08417434 0.06665123 0.126298189
    ## 7    CSE per km LCP_mean 1091 0.007946381 0.08914247 0.07077996 0.020119399
    ## 8    CSE per km  LCP_sum 1091 0.004766389 0.06903904 0.05503143 0.412249177
    ## 9    CSE per km     CBER 1091 0.007274920 0.08529314 0.06762456 0.102918244
    ##   Correlation
    ## 1   0.3182048
    ## 2   0.7847884
    ## 3   0.7532424
    ## 4   0.1139171
    ## 5  -0.0545083
    ## 6   0.3553846
    ## 7  -0.1418429
    ## 8   0.6420663
    ## 9   0.3208087

## Screening metrics

Spatial screening metrics for all response formulations and spatial
summaries.

``` r
metrics_all
```

    ##        response   method    n         MSE       RMSE        MAE         RSQ
    ## 1       Raw CSE LCP_mean 1091 0.007288414 0.08537221 0.06685447 0.101254278
    ## 2       Raw CSE  LCP_sum 1091 0.003114933 0.05581158 0.04540066 0.615892809
    ## 3       Raw CSE     CBER 1091 0.003508397 0.05923172 0.04633363 0.567374124
    ## 4 IBD residuals LCP_mean 1091 0.008004301 0.08946676 0.07061746 0.012977113
    ## 5 IBD residuals  LCP_sum 1091 0.008085445 0.08991910 0.07127313 0.002971155
    ## 6 IBD residuals     CBER 1091 0.007085320 0.08417434 0.06665123 0.126298189
    ## 7    CSE per km LCP_mean 1091 0.007946381 0.08914247 0.07077996 0.020119399
    ## 8    CSE per km  LCP_sum 1091 0.004766389 0.06903904 0.05503143 0.412249177
    ## 9    CSE per km     CBER 1091 0.007274920 0.08529314 0.06762456 0.102918244
    ##   Correlation
    ## 1   0.3182048
    ## 2   0.7847884
    ## 3   0.7532424
    ## 4   0.1139171
    ## 5  -0.0545083
    ## 6   0.3553846
    ## 7  -0.1418429
    ## 8   0.6420663
    ## 9   0.3208087

## R-squared comparison

``` r
ggplot(metrics_all, aes(x = method, y = RSQ, fill = response)) +
  geom_col(position = "dodge", alpha = 0.5) +
  scale_fill_manual(values = cols) +
  labs(x = "Spatial summary", y = expression(R^2), fill = "Response formulation") +
  theme_bw(base_size = 14) +
  theme(panel.grid.minor = element_blank())
```

![](/tmp/RtmpaVWtN9/knitted_mds/results-rsq-plot-1.png)<!-- -->

## Correlation comparison

``` r
ggplot(metrics_all, aes(x = method, y = Correlation, fill = response)) +
  geom_col(position = "dodge", alpha = 0.5) +
  scale_fill_manual(values = cols) +
  labs(x = "Spatial summary", y = "Correlation", fill = "Response formulation") +
  theme_bw(base_size = 14) +
  theme(panel.grid.minor = element_blank())
```

![](/tmp/RtmpaVWtN9/knitted_mds/results-correlation-plot-1.png)<!-- -->

## Best-performing combination by R-squared

``` r
best_rsq <- metrics_all %>% arrange(desc(RSQ)) 

best_rsq
```

    ##        response   method    n         MSE       RMSE        MAE         RSQ
    ## 1       Raw CSE  LCP_sum 1091 0.003114933 0.05581158 0.04540066 0.615892809
    ## 2       Raw CSE     CBER 1091 0.003508397 0.05923172 0.04633363 0.567374124
    ## 3    CSE per km  LCP_sum 1091 0.004766389 0.06903904 0.05503143 0.412249177
    ## 4 IBD residuals     CBER 1091 0.007085320 0.08417434 0.06665123 0.126298189
    ## 5    CSE per km     CBER 1091 0.007274920 0.08529314 0.06762456 0.102918244
    ## 6       Raw CSE LCP_mean 1091 0.007288414 0.08537221 0.06685447 0.101254278
    ## 7    CSE per km LCP_mean 1091 0.007946381 0.08914247 0.07077996 0.020119399
    ## 8 IBD residuals LCP_mean 1091 0.008004301 0.08946676 0.07061746 0.012977113
    ## 9 IBD residuals  LCP_sum 1091 0.008085445 0.08991910 0.07127313 0.002971155
    ##   Correlation
    ## 1   0.7847884
    ## 2   0.7532424
    ## 3   0.6420663
    ## 4   0.3553846
    ## 5   0.3208087
    ## 6   0.3182048
    ## 7  -0.1418429
    ## 8   0.1139171
    ## 9  -0.0545083

## Observed versus predicted plots from existing saved RDS files

``` r
for (label in names(result_sets)) {
  stem <- result_sets[[label]]
  f <- file.path(results_dir_big, paste0("spatial_eval_", stem, "_predictions.rds"))
  if (!file.exists(f)) {
    cat("\n\nExisting prediction RDS not found for ", label, ": ", f, "\n\n", sep = "")
    next
  }
  d <- readRDS(f)
  for (m in c("LCP_sum", "LCP_mean", "CBER")) {
    if (!all(c("CSE", m) %in% names(d))) next
    p <- ggplot(d, aes(x = .data[[m]], y = CSE)) +
      geom_point(color = cols[[label]], alpha = 0.5) +
      geom_smooth(method = "lm", se = TRUE,
                  color = cols[[label]], fill = cols[[label]], alpha = 0.2) +
      labs(title = paste(label, "-", m), x = paste(m, "prediction"), y = "Observed CSE") +
      theme_bw(base_size = 14) +
      theme(panel.grid.minor = element_blank())
    print(p)
  }
}
```

    ## `geom_smooth()` using formula = 'y ~ x'

![](/tmp/RtmpaVWtN9/knitted_mds/results-prediction-plots-1.png)<!-- -->

    ## `geom_smooth()` using formula = 'y ~ x'

![](/tmp/RtmpaVWtN9/knitted_mds/results-prediction-plots-2.png)<!-- -->

    ## `geom_smooth()` using formula = 'y ~ x'

![](/tmp/RtmpaVWtN9/knitted_mds/results-prediction-plots-3.png)<!-- -->

    ## `geom_smooth()` using formula = 'y ~ x'

![](/tmp/RtmpaVWtN9/knitted_mds/results-prediction-plots-4.png)<!-- -->

    ## `geom_smooth()` using formula = 'y ~ x'

![](/tmp/RtmpaVWtN9/knitted_mds/results-prediction-plots-5.png)<!-- -->

    ## `geom_smooth()` using formula = 'y ~ x'

![](/tmp/RtmpaVWtN9/knitted_mds/results-prediction-plots-6.png)<!-- -->

    ## `geom_smooth()` using formula = 'y ~ x'

![](/tmp/RtmpaVWtN9/knitted_mds/results-prediction-plots-7.png)<!-- -->

    ## `geom_smooth()` using formula = 'y ~ x'

![](/tmp/RtmpaVWtN9/knitted_mds/results-prediction-plots-8.png)<!-- -->

    ## `geom_smooth()` using formula = 'y ~ x'

![](/tmp/RtmpaVWtN9/knitted_mds/results-prediction-plots-9.png)<!-- -->
