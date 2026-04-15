LOPOCV: Random Forest (Internal) and Spatial Cross-Validation
================
Norah Saarman
2026-04-15

- [Setup](#setup)
  - [Overview of script](#overview-of-script)
  - [Inputs](#inputs)
  - [Outputs](#outputs)
  - [Libraries](#libraries)
  - [Define Paths to Directories](#define-paths-to-directories)
- [Permutation test for LOPOCV R²](#permutation-test-for-lopocv-r²)
  - [Chunk 1: LOPOCV on permuted
    models](#chunk-1-lopocv-on-permuted-models)
  - [Chunk 2: Calculate the metrics](#chunk-2-calculate-the-metrics)
  - [Chunk 3: Compare Distributions](#chunk-3-compare-distributions)
  - [Chunk 5: Visualize permutated vs observed LOPOCV
    side-by-side](#chunk-5-visualize-permutated-vs-observed-lopocv-side-by-side)
- [TODO: Spatial Eval, doesn’t seem worth too much
  energy](#todo-spatial-eval-doesnt-seem-worth-too-much-energy)
  - [TODO: Chunk 6: Spatial evaluation for permuted LOPOCV
    models](#todo-chunk-6-spatial-evaluation-for-permuted-lopocv-models)
  - [TODO: Chunk 7: Spatial permutation
    plots](#todo-chunk-7-spatial-permutation-plots)
  - [TODO: Chunk 8: Delete RDS when completely
    done](#todo-chunk-8-delete-rds-when-completely-done)

# Setup

RStudio Configuration:  
- **R version:** R 4.4.0 (Geospatial packages)  
- **Number of cores:** 16 (up to 32 available)  
- **Account:** saarman-np  
- **Partition:** saarman-np (allows multiple simultaneous jobs
automatically now)  
- **Memory per job:** 400G (cluster limit: 1000G total; avoid exceeding
half)

## Overview of script

Although model configuration and projection parameters were selected
based on full-model performance (RSQ), during LOPOCV, some held-out
folds exhibited limited variance in observed CSE, making fold-specific
RSQ unstable. We therefore report fold-level performance using Spearman
and Pearson correlation alongside RMSE and MAE, which provide more
robust summaries under these conditions.

(In contrast, model configuration and projection parameters were
selected based on full-model performance (RSQ) becaues in these cases,
variance in CSE was sufficiently large for stable estimation.)

## Inputs

- `../input/Gff_11loci_68sites_cse.csv` - Combined CSE table with
  coordinates (long1, lat1, long2, lat2)
- `../results_dir/fullRF_CSE_resistance.tif` - Final full model
  projected resistance surface
- `../results_dir/LC_paths_fullRF.shp"` -
- `../data_dir/processed/env_stack.grd` - Final prediction env stack
  with named layers (18 variables) env \<- stack(file.path(
- `../results_dir/lopocv/rf_model_01.rds` - 67 LOPOCV rf models leaving
  one point out

I also have the observed data completed as well:

non-spatial observed fold Pearsons: results/LOPOCV_summary.csv column =
Pearson

non-spatial pooled Pearson: results/LOPOCV_pooled_summary.csv column =
pooled_Pearson

spatial observed fold Pearsons:
results/spatial_LOPOCV_LCPsum_k1_summary.csv column = Pearson

spatial pooled Pearson:
results/spatial_LOPOCV_LCPsum_k1_pooled_summary.csv column =
pooled_Pearson

## Outputs

- `../results_dir/lopocv/` - rf_model_01.rds to rf_model_67.rds for each
  fold of LOPOCV  
- `../results_dir/lopocv_spatial_LCPsum_k/` - Spatial lopocv predictions
  and metrics for each fold 1-67
- `../results_dir/spatial_LOPOCV_LCPsum_k1_summary.csv` - Spatial lopocv
  evaluation metrics

## Libraries

``` r
# load only required packages
library(doParallel)
library(foreach)
library(raster)
library(gdistance)
library(sf)
library(dplyr)
library(randomForest)
library(readr)
library(ggplot2)
library(sp)
library(tidyverse)
library(rnaturalearth)
library(rnaturalearthdata)
library(reshape2)
```

## Define Paths to Directories

``` r
# Define paths
data_dir  <- "/uufs/chpc.utah.edu/common/home/saarman-group1/uganda-tsetse-LG/data"
input_dir <- "../input"
results_dir <- "/uufs/chpc.utah.edu/common/home/saarman-group1/uganda-tsetse-LG/results/"
output_dir <- file.path(results_dir, "lopocv")
perm_model_dir <- file.path(results_dir, "perm_rds_temp")
dir.create(perm_model_dir, recursive = TRUE, showWarnings = FALSE)
stopifnot(dir.exists(perm_model_dir))
stopifnot(file.access(perm_model_dir, 2) == 0)

# High-quality figure output (for manuscript)
fig_dir <- file.path(results_dir, "figures_pub")
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)

# define coordinate reference system
crs_geo <- 4326     # EPSG code for WGS84

# define ggplot2 extent
xlim <- c(28.6, 35.4)
ylim <- c(-1.500000 , 4.733333)

# Input: shape file of least-cost paths (already filtered to 67 sites and has CSEdistance)
lcp_sf <- st_read(file.path(results_dir, "LC_paths_fullRF.shp"))
```

    ## Reading layer `LC_paths_fullRF' from data source 
    ##   `/uufs/chpc.utah.edu/common/home/saarman-group1/uganda-tsetse-LG/results/LC_paths_fullRF.shp' 
    ##   using driver `ESRI Shapefile'
    ## Simple feature collection with 1026 features and 4 fields
    ## Geometry type: LINESTRING
    ## Dimension:     XY
    ## Bounding box:  xmin: 31.12083 ymin: -0.5958333 xmax: 34.5125 ymax: 3.695833
    ## Geodetic CRS:  WGS 84

``` r
st_crs(lcp_sf) <- crs_geo

# Build list of unique sites from Var1 and Var2
sites <- sort(unique(c(lcp_sf$Var1, lcp_sf$Var2)))

output_dir <- file.path(results_dir, "lopocv")
#dir.create(output_dir, showWarnings = FALSE)
scratch_dir <- "/scratch/local/u6036559"
```

# Permutation test for LOPOCV R²

GOAL: Create a 67 x 100 matrix where each cell is the Pearson’s
correlation for one LOPOCV fold (site) under a permuted response
(CSEdistance).

## Chunk 1: LOPOCV on permuted models

Run models and save rds

``` r
# Chunk 1: setup + data + helpers + run permuted LOPOCV models and save .rds files


# Load required packages
library(dplyr)
library(readr)
library(randomForest)
library(doParallel)
library(foreach)

# Load data
V.table_full <- read.csv(file.path(input_dir, "Gff_cse_envCostPaths.csv"))

# Filter out western outlier
V.table <- V.table_full %>%
  filter(Var1 != "50-KB", Var2 != "50-KB")

# Create unique pair ID
V.table$id <- paste(V.table$Var1, V.table$Var2, sep = "_")

# Define site list
sites <- sort(unique(c(V.table$Var1, V.table$Var2)))

# Predictors used in the original LOPOCV models
predictor_vars <- c(
  "BIO1_mean","BIO2_mean","BIO3_mean","BIO4_mean","BIO5_mean","BIO6_mean",
  "BIO7_mean","BIO8S_mean","BIO9S_mean","BIO10S_mean","BIO11S_mean",
  "BIO12_mean","BIO13_mean","BIO14_mean","BIO15_mean","BIO16S_mean",
  "BIO17S_mean","BIO18S_mean","BIO19S_mean","slope_mean","alt_mean",
  "lakes_mean","riv_3km_mean","samp_20km_mean","pix_dist"
)

# Modeling table
V.model <- V.table[, c("CSEdistance", predictor_vars)]

# Set number of permutations
n_perm <- 100

# Set up parallel backend
cl <- makeCluster(8)
registerDoParallel(cl)

# Run permutations and save models
foreach(
  p = 1:n_perm,
  .packages = c("randomForest", "dplyr")
) %dopar% {
  set.seed(500 + p)

  # Shuffle response column while preserving row order
  V.model_perm <- V.model
  V.model_perm$CSEdistance <- sample(V.model_perm$CSEdistance)

  for (i in seq_along(sites)) {
    site <- sites[i]

    # Identify test/train indices
    test_ids  <- V.table$id[V.table$Var1 == site | V.table$Var2 == site]
    test_idx  <- which(V.table$id %in% test_ids)
    train_idx <- which(!V.table$id %in% test_ids)

    train_df <- V.model_perm[train_idx, ]
    test_df  <- V.model_perm[test_idx, ]

    # Fit Random Forest on permuted response
    rf_model <- randomForest(
      CSEdistance ~ .,
      data = train_df,
      ntree = 500
    )

    # Save RF model for this permutation–fold pair
    saveRDS(
  rf_model,
  file.path(perm_model_dir, sprintf("rf_model_perm_%02d_fold_%02d.rds", p, i))
)
    # Save metrics file for this permutation–fold pair
    #metrics_file <- file.path(results_dir, sprintf("perm_progress_%03d.csv", p))                               # commented this out because it is pretty much a blank file with just the permutation number, some mistake i made...
#write.csv(data.frame(perm = p), metrics_file, row.names = FALSE)
  }

  NULL
}

# Stop cluster
stopCluster(cl)
```

## Chunk 2: Calculate the metrics

``` r
# Chunk 2: calculate Pearson null matrix from saved permutation models

# Helper for calibration
calibrate_pred <- function(obs_train, pred_train, pred_test) {
  keep <- complete.cases(obs_train, pred_train)
  df <- data.frame(obs = obs_train[keep], pred = pred_train[keep])
  fit <- lm(obs ~ pred, data = df)

  out <- rep(NA_real_, length(pred_test))
  keep_test <- complete.cases(pred_test)
  out[keep_test] <- predict(fit, newdata = data.frame(pred = pred_test[keep_test]))
  out
}

# Set number of permutations
n_perm <- 100

# Preallocate matrix for results
pearson_null <- matrix(NA_real_, nrow = length(sites), ncol = n_perm)
rownames(pearson_null) <- sites
colnames(pearson_null) <- paste0("perm_", seq_len(n_perm))

for (p in seq_len(n_perm)) {
  set.seed(500 + p)

  # Recreate permuted response exactly as in Chunk 1
  V.model_perm <- V.model
  V.model_perm$CSEdistance <- sample(V.model_perm$CSEdistance)

  pearson_perm <- numeric(length(sites))

  for (i in seq_along(sites)) {
    site <- sites[i]

    fold_model_file <- file.path(
  perm_model_dir,
  sprintf("rf_model_perm_%02d_fold_%02d.rds", p, i)
)

    if (!file.exists(fold_model_file)) {
      warning(sprintf("Missing model file: %s", fold_model_file))
      pearson_perm[i] <- NA_real_
      next
    }

    # Identify test/train indices
    test_ids  <- V.table$id[V.table$Var1 == site | V.table$Var2 == site]
    test_idx  <- which(V.table$id %in% test_ids)
    train_idx <- which(!V.table$id %in% test_ids)

    train_df <- V.model_perm[train_idx, ]
    test_df  <- V.model_perm[test_idx, ]

    rf_model <- readRDS(fold_model_file)

    train_obs <- train_df$CSEdistance
    test_obs  <- test_df$CSEdistance

    train_df <- train_df %>% dplyr::select(-CSEdistance)
    test_df  <- test_df %>% dplyr::select(-CSEdistance)

    stopifnot(!"CSEdistance" %in% names(train_df))
    stopifnot(!"CSEdistance" %in% names(test_df))

    stopifnot(identical(names(train_df), predictor_vars))
    stopifnot(identical(names(test_df), predictor_vars))
    stopifnot(identical(attr(rf_model$terms, "term.labels"), predictor_vars))

    pred_train <- predict(rf_model, newdata = train_df)
    pred_test  <- predict(rf_model, newdata = test_df)

    pred_test_cal <- calibrate_pred(
      obs_train = train_obs,
      pred_train = pred_train,
      pred_test = pred_test
    )

    keep <- complete.cases(test_obs, pred_test_cal)

    pearson_perm[i] <- cor(
      test_obs[keep],
      pred_test_cal[keep],
      method = "pearson"
    )
  }

  pearson_null[, p] <- pearson_perm
}

# Save matrix to file for later use
write.csv(
  pearson_null,
  file.path(results_dir, "pearson_LOPOCV_null_matrix.csv"),
  row.names = TRUE
)

print(pearson_null)
message("Permutation Pearson null matrix created and saved.")
```

## Chunk 3: Compare Distributions

Compare observed vs permuted LOPOCV Pearson:

``` r
# observed LOPOCV summary from the true data
metrics_all <- read.csv(file.path(results_dir, "LOPOCV_summary.csv"))

# permutation matrix
pearson_null <- read.csv(file.path(results_dir, "pearson_LOPOCV_null_matrix.csv"))
  
# make sure site order matches the permutation matrix
metrics_all <- metrics_all[match(rownames(pearson_null), metrics_all$site), ]

pearson_obs <- as.numeric(metrics_all$Pearson)

pearson_null_mat <- as.matrix(pearson_null)

# reorder rows to observed site order if needed
pearson_null_mat <- pearson_null_mat[match(metrics_all$site, rownames(pearson_null_mat)), ]

plot_df_pearson <- data.frame(
  pearson = c(pearson_obs, as.numeric(pearson_null_mat)),
  type = c(
    rep("Observed", length(pearson_obs)),
    rep("Permuted", length(as.numeric(pearson_null_mat)))
  )
)

ggplot(plot_df_pearson, aes(x = pearson, fill = type)) +
  geom_density(alpha = 0.5, color = NA) +
  coord_cartesian(xlim = c(-1, 1)) +
  theme_minimal() +
  labs(
    title = "Observed vs. Permuted Pearson Correlation (LOPOCV)",
    x = "Pearson Correlation",
    y = "Density"
  )
```

    ## Warning: Removed 6834 rows containing non-finite outside the scale range
    ## (`stat_density()`).

![](../figures/knitted_mds/plot-pearson-lopocv-1.png)<!-- --> \## Chunk
4: Empirical p-value GOAL: Calculate an empirical p-value comparing
observed vs permuted Pearson’s correlation.

``` r
# Define region based on SiteMajCluster
Gff <- read.csv(file.path(input_dir, "Gff_11loci_allsites_indinfo.txt"),
                header = TRUE, sep = "\t")
north_sites <- unique(Gff$SiteCode[Gff$SiteMajCluster == "north"])

# Load observed Pearson
metrics_all <- read.csv(file.path(results_dir, "LOPOCV_summary.csv"))
pearson_obs_df <- data.frame(
  site = metrics_all$site,
  pearson = metrics_all$Pearson
)
pearson_obs_df$cluster <- ifelse(pearson_obs_df$site %in% north_sites, "north", "south")
pearson_obs_df$type <- "Observed"

pearson_null_raw <- read.csv(
  file.path(results_dir, "pearson_LOPOCV_null_matrix.csv"),
  row.names = 1,
  check.names = FALSE
)
pearson_null_mat <- as.matrix(pearson_null_raw)

# align rows to observed site order
pearson_null_mat <- pearson_null_mat[match(metrics_all$site, rownames(pearson_null_mat)), ]

obs_mean_pearson <- mean(pearson_obs, na.rm = TRUE)
null_mean_pearson <- apply(pearson_null_mat, 2, mean, na.rm = TRUE)

# one-sided empirical p-value: observed mean Pearson greater than null
p_empirical <- (sum(null_mean_pearson >= obs_mean_pearson, na.rm = TRUE) + 1) /
  (sum(!is.na(null_mean_pearson)) + 1)

cat("Observed mean Pearson:", round(obs_mean_pearson, 4), "\n")
```

    ## Observed mean Pearson: NaN

``` r
cat("Null mean Pearson range:", round(range(null_mean_pearson, na.rm = TRUE), 4), "\n")
```

    ## Null mean Pearson range: -0.0763 0.0669

``` r
cat("Empirical one-sided p-value:", round(p_empirical, 4), "\n\n")
```

    ## Empirical one-sided p-value: 0.0099

``` r
obs_north <- pearson_obs_df$pearson[pearson_obs_df$cluster == "north"]
obs_south <- pearson_obs_df$pearson[pearson_obs_df$cluster == "south"]

null_north_mat <- as.matrix(pearson_null[rownames(pearson_null) %in% north_sites, ])
null_south_mat <- as.matrix(pearson_null[!rownames(pearson_null) %in% north_sites, ])

obs_mean_north <- mean(obs_north, na.rm = TRUE)
obs_mean_south <- mean(obs_south, na.rm = TRUE)

null_mean_north <- apply(null_north_mat, 2, mean, na.rm = TRUE)
null_mean_south <- apply(null_south_mat, 2, mean, na.rm = TRUE)

p_north <- (sum(null_mean_north >= obs_mean_north) + 1) / (length(null_mean_north) + 1)
p_south <- (sum(null_mean_south >= obs_mean_south) + 1) / (length(null_mean_south) + 1)

cat("Observed mean Pearson, north:", round(obs_mean_north, 4), "\n")
```

    ## Observed mean Pearson, north: 0.9122

``` r
cat("Empirical p-value, north:", round(p_north, 4), "\n\n")
```

    ## Empirical p-value, north: NA

``` r
cat("Observed mean Pearson, south:", round(obs_mean_south, 4), "\n")
```

    ## Observed mean Pearson, south: 0.9025

``` r
cat("Empirical p-value, south:", round(p_south, 4), "\n")
```

    ## Empirical p-value, south: NA

## Chunk 5: Visualize permutated vs observed LOPOCV side-by-side

``` r
# Define region based on SiteMajCluster
Gff <- read.csv(file.path(input_dir, "Gff_11loci_allsites_indinfo.txt"),
                header = TRUE, sep = "\t")
north_sites <- unique(Gff$SiteCode[Gff$SiteMajCluster == "north"])

# Load observed Pearson
metrics_all <- read.csv(file.path(results_dir, "LOPOCV_summary.csv"))
pearson_obs_df <- data.frame(
  site = metrics_all$site,
  pearson = metrics_all$Pearson
)
pearson_obs_df$cluster <- ifelse(pearson_obs_df$site %in% north_sites, "north", "south")
pearson_obs_df$type <- "Observed"

# Load permuted Pearson matrix
pearson_null <- read.csv(
  file.path(results_dir, "pearson_LOPOCV_null_matrix.csv"),
  row.names = 1,
  check.names = FALSE
)

pearson_null_df <- data.frame(
  site = rownames(pearson_null),
  cluster = ifelse(rownames(pearson_null) %in% north_sites, "north", "south"),
  pearson_null,
  check.names = FALSE
)

# Stack all permutation columns
pearson_null_long <- reshape2::melt(
  pearson_null_df,
  id.vars = c("site", "cluster"),
  variable.name = "perm",
  value.name = "pearson"
)
pearson_null_long$type <- "Permuted"

# Combine
pearson_plot_df <- rbind(
  pearson_obs_df[, c("site", "cluster", "pearson", "type")],
  pearson_null_long[, c("site", "cluster", "pearson", "type")]
)

# Run just once to save high quality png for publication
#png(file.path(fig_dir, "Fig_LOPOCV_obs_vs_permuted.png"), 
#       width = 7, height = 2.8, units = "in", res = 600) 

# Plot
ggplot(pearson_plot_df, aes(x = pearson, fill = interaction(type, cluster))) +
  geom_density(
    data = subset(pearson_plot_df, type == "Observed"),
    aes(y = after_stat(count)),
    alpha = 0.6,
    color = NA,
    position = "identity"
  ) +
  geom_density(
    data = subset(pearson_plot_df, type == "Permuted"),
    aes(y = after_stat(count / 100)),
    alpha = 0.3,
    color = NA,
    position = "identity"
  ) +
  scale_fill_manual(
    name = "Group",
    values = c(
      "Observed.north" = "#1f78b4",
      "Observed.south" = "#e66101",
      "Permuted.north" = "#2c3e50",
      "Permuted.south" = "#4e342e"
    ),
    labels = c(
      "Observed.north" = "North",
      "Observed.south" = "South",
      "Permuted.north" = "Permuted (North)",
      "Permuted.south" = "Permuted (South)"
    )
  ) +
  coord_cartesian(xlim = c(-1, 1)) +
  theme_minimal() +
  labs(
    title = "Observed vs. Permuted Pearson by cluster",
    x = "Pearson's r (LOPOCV)",
    y = "Density"
  )
```

![](../figures/knitted_mds/plot-side-by-side-1.png)<!-- -->

``` r
#dev.off()
```

# TODO: Spatial Eval, doesn’t seem worth too much energy

## TODO: Chunk 6: Spatial evaluation for permuted LOPOCV models

``` r
# where to save spatial permutation outputs
spatial_output_dir <- file.path(results_dir, "lopocv_spatial_perm_LCPsum_k1")
dir.create(spatial_output_dir, showWarnings = FALSE, recursive = TRUE)

# where to find saved permuted .rds files
perm_model_dir <- file.path(results_dir, "perm_rds_temp")

# use same fixed projection distance as working script
kdist <- 1

# helper
calibrate_pred <- function(obs_train, pred_train, pred_test) {
  keep <- complete.cases(obs_train, pred_train)
  df <- data.frame(obs = obs_train[keep], pred = pred_train[keep])
  fit <- lm(obs ~ pred, data = df)

  out <- rep(NA_real_, length(pred_test))
  keep_test <- complete.cases(pred_test)
  out[keep_test] <- predict(fit, newdata = data.frame(pred = pred_test[keep_test]))
  out
}

calc_metrics <- function(obs, pred) {
  keep <- complete.cases(obs, pred)

  obs  <- obs[keep]
  pred <- pred[keep]

  mse      <- mean((obs - pred)^2)
  rmse     <- sqrt(mse)
  mae      <- mean(abs(obs - pred))
  pearson  <- cor(obs, pred, method = "pearson")
  spearman <- cor(obs, pred, method = "spearman")

  data.frame(
    n = length(obs),
    MSE = mse,
    RMSE = rmse,
    MAE = mae,
    Pearson = pearson,
    Spearman = spearman
  )
}

# use one level of parallel only: over permutations
n_cores <- 8
cl <- makeCluster(n_cores)
registerDoParallel(cl)

n_perm <- 100

foreach(
  p = seq_len(n_perm),
  .packages = c("raster", "gdistance", "sf", "sp", "dplyr", "randomForest")
) %dopar% {

  perm_file <- file.path(
    spatial_output_dir,
    sprintf("spatial_perm_LCPsum_k1_%03d.csv", p)
  )

  if (file.exists(perm_file)) {
    return(NULL)
  }

  set.seed(500 + p)

  # recreate permuted response exactly as in chunk 1 and 2
  V.model_perm <- V.model
  V.model_perm$CSEdistance <- sample(V.model_perm$CSEdistance)

  V.table_perm <- V.table
  V.table_perm$CSEdistance <- V.model_perm$CSEdistance

  # pair table used in spatial evaluation
  site_pairs_perm <- V.table_perm %>%
    dplyr::select(id, Var1, Var2, CSE = CSEdistance) %>%
    distinct()

  perm_metrics <- vector("list", length(sites))

  for (fold_idx in seq_along(sites)) {

    site <- sites[fold_idx]

    fold_model_file <- file.path(
      perm_model_dir,
      sprintf("rf_model_perm_%02d_fold_%02d.rds", p, fold_idx)
    )

    if (!file.exists(fold_model_file)) {
      perm_metrics[[fold_idx]] <- data.frame(
        permutation = p,
        site = site,
        projection_km = kdist,
        method = "LCP_sum",
        n = NA_real_,
        MSE = NA_real_,
        RMSE = NA_real_,
        MAE = NA_real_,
        Pearson = NA_real_,
        Spearman = NA_real_
      )
      next
    }

    rf_model <- readRDS(fold_model_file)

    # Step 1: project surface at kdist = 1
    env_k <- env
    env_k[["pix_dist"]] <- setValues(env[[1]], kdist)

    rf_vars <- attr(rf_model$terms, "term.labels")
    env_k <- env_k[[rf_vars]]

    stopifnot(identical(names(env_k), rf_vars))

    cse_surface <- raster::predict(env_k, rf_model, type = "response")

    # align lake mask if needed
    if (!compareRaster(cse_surface, lake_mask,
                       extent = TRUE, rowcol = TRUE,
                       crs = TRUE, res = TRUE,
                       stopiffalse = FALSE)) {
      lake_mask_use <- projectRaster(lake_mask, cse_surface, method = "ngb")
    } else {
      lake_mask_use <- lake_mask
    }

    # assign lakes high cost
    cse_surface[lake_mask_use[] == 1] <- max(values(cse_surface), na.rm = TRUE)

    # Step 2: prepare sites in raster CRS
    sites_df <- st_as_sf(sites_coords, coords = c("lon", "lat"), crs = crs_geo) %>%
      st_transform(crs(cse_surface))

    sites_sp <- as(sites_df, "Spatial")
    site_index <- setNames(seq_len(nrow(sites_df)), sites_coords$Site)

    # Step 3: build transition object
    resistance_rast <- cse_surface
    min_pos <- min(values(resistance_rast)[values(resistance_rast) > 0], na.rm = TRUE)
    resistance_rast[resistance_rast <= 0] <- min_pos

    conductance_rast <- 1 / resistance_rast
    tr <- transition(conductance_rast, transitionFunction = mean, directions = 16)
    tr_lcp <- geoCorrection(tr, type = "c")

    # Step 4: recompute least-cost paths
    paths_list <- vector("list", nrow(site_pairs_perm))

    for (ii in seq_len(nrow(site_pairs_perm))) {
      i <- site_index[site_pairs_perm$Var1[ii]]
      j <- site_index[site_pairs_perm$Var2[ii]]

      path <- tryCatch(
        shortestPath(tr_lcp, sites_sp[i, ], sites_sp[j, ], output = "SpatialLines"),
        error = function(e) NULL
      )

      if (!is.null(path)) {
        path_sf <- st_as_sf(path)
        path_sf$id <- site_pairs_perm$id[ii]
        paths_list[[ii]] <- path_sf
      }
    }

    paths_list <- paths_list[!sapply(paths_list, is.null)]

    if (length(paths_list) == 0) {
      perm_metrics[[fold_idx]] <- data.frame(
        permutation = p,
        site = site,
        projection_km = kdist,
        method = "LCP_sum",
        n = NA_real_,
        MSE = NA_real_,
        RMSE = NA_real_,
        MAE = NA_real_,
        Pearson = NA_real_,
        Spearman = NA_real_
      )
      next
    }

    paths_sf <- do.call(rbind, paths_list)
    st_crs(paths_sf) <- st_crs(cse_surface)

    # extract raster values along each path
    path_vals <- raster::extract(cse_surface, as(paths_sf, "Spatial"))

    lcp_summary <- data.frame(
      id = paths_sf$id,
      LCP_sum = sapply(path_vals, function(x) {
        if (is.null(x) || all(is.na(x))) NA_real_ else sum(x, na.rm = TRUE)
      })
    )

    # Step 5: combine and calibrate
    eval_results <- site_pairs_perm %>%
      left_join(lcp_summary, by = "id")

    test_ids <- eval_results$id[eval_results$Var1 == site | eval_results$Var2 == site]
    train_df <- eval_results[!eval_results$id %in% test_ids, ]
    test_df  <- eval_results[eval_results$id %in% test_ids, ]

    test_df$LCP_sum_cal <- calibrate_pred(
      obs_train = train_df$CSE,
      pred_train = train_df$LCP_sum,
      pred_test = test_df$LCP_sum
    )

    metrics_k <- calc_metrics(
      obs = test_df$CSE,
      pred = test_df$LCP_sum_cal
    ) %>%
      mutate(
        permutation = p,
        site = site,
        projection_km = kdist,
        method = "LCP_sum"
      ) %>%
      dplyr::select(permutation, site, projection_km, method, everything())

    perm_metrics[[fold_idx]] <- metrics_k
  }

  perm_metrics_df <- bind_rows(perm_metrics)
  write.csv(perm_metrics_df, perm_file, row.names = FALSE)

  NULL
}

stopCluster(cl)

# combine permutation outputs into null matrices
metric_names <- c("Pearson", "Spearman", "RMSE", "MAE")

for (metric_name in metric_names) {
  metric_null <- matrix(NA_real_, nrow = length(sites), ncol = n_perm)
  rownames(metric_null) <- sites
  colnames(metric_null) <- paste0("perm_", seq_len(n_perm))

  for (p in seq_len(n_perm)) {
    perm_file <- file.path(
      spatial_output_dir,
      sprintf("spatial_perm_LCPsum_k1_%03d.csv", p)
    )

    if (!file.exists(perm_file)) next

    df <- read.csv(perm_file)
    df <- df[match(sites, df$site), ]

    metric_null[, p] <- df[[metric_name]]
  }

  write.csv(
    metric_null,
    file.path(
      results_dir,
      sprintf("%s_spatial_LOPOCV_LCPsum_k1_null_matrix.csv", tolower(metric_name))
    ),
    row.names = TRUE
  )
}
```

## TODO: Chunk 7: Spatial permutation plots

Use this only after the spatial permutation chunk has finished.

``` r
rsq_null_raw <- read.csv(
  file.path(results_dir, "rsq_spatial_LOPOCV_LCPsum_k1_null_matrix.csv"),
  row.names = 1
)

rsq_null_mat <- as.matrix(rsq_null_raw)
rsq_null_long <- as.numeric(as.vector(rsq_null_mat))

metrics_all <- read.csv(file.path(results_dir, "spatial_LOPOCV_LCPsum_k1_summary.csv"))
rsq_obs <- as.numeric(metrics_all$RSQ)

plot_df_rsq <- data.frame(
  rsq = c(rsq_obs, rsq_null_long),
  type = c(rep("Observed", length(rsq_obs)),
           rep("Permuted", length(rsq_null_long)))
)

ggplot(plot_df_rsq, aes(x = rsq, fill = type)) +
  geom_density(alpha = 0.5, color = NA) +
  scale_fill_manual(values = c("Observed" = "#1f78b4", "Permuted" = "gray70")) +
  xlim(-1, 1) +
  theme_minimal() +
  labs(title = "Observed vs. Permuted R² (Spatial LOPOCV)",
       x = "Test R²",
       y = "Density")
```

## TODO: Chunk 8: Delete RDS when completely done

``` r
unlink(perm_model_dir, recursive = TRUE)
```
