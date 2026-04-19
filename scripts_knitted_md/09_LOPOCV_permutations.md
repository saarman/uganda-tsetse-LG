LOPOCV with permutated data as Null
================
Norah Saarman
2026-04-20

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
  - [Chunk 3: Compare observed vs permuted LOPOCV
    Spearman](#chunk-3-compare-observed-vs-permuted-lopocv-spearman)
  - [Chunk 4: Empirical p-value](#chunk-4-empirical-p-value)
  - [Chunk 5: Visualize permutated vs observed LOPOCV
    side-by-side](#chunk-5-visualize-permutated-vs-observed-lopocv-side-by-side)
  - [Chunk 6: Spatial included, all
    side-by-side](#chunk-6-spatial-included-all-side-by-side)
  - [TODO: Chunk 7: Delete RDS when completely
    done](#todo-chunk-7-delete-rds-when-completely-done)

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

helpers:

``` r
# Chunk 1: setup + data + helpers 

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
```

1b run the permutations and save .rds files

``` r
#run permuted LOPOCV models and save .rds files

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
# Chunk 2: calculate Spearman null matrix from saved permutation models

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
spearman_null <- matrix(NA_real_, nrow = length(sites), ncol = n_perm)
rownames(spearman_null) <- sites
colnames(spearman_null) <- paste0("perm_", seq_len(n_perm))

for (p in seq_len(n_perm)) {
  set.seed(500 + p)

  # Recreate permuted response exactly as in Chunk 1
  V.model_perm <- V.model
  V.model_perm$CSEdistance <- sample(V.model_perm$CSEdistance)

  spearman_perm <- numeric(length(sites))

  for (i in seq_along(sites)) {
    site <- sites[i]

    fold_model_file <- file.path(
      perm_model_dir,
      sprintf("rf_model_perm_%02d_fold_%02d.rds", p, i)
    )

    if (!file.exists(fold_model_file)) {
      warning(sprintf("Missing model file: %s", fold_model_file))
      spearman_perm[i] <- NA_real_
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

    spearman_perm[i] <- cor(
      test_obs[keep],
      pred_test_cal[keep],
      method = "spearman"
    )
  }

  spearman_null[, p] <- spearman_perm
}

# Save matrix to file for later use
write.csv(
  spearman_null,
  file.path(results_dir, "spearman_LOPOCV_null_matrix.csv"),
  row.names = TRUE
)

print(spearman_null)
message("Permutation Spearman null matrix created and saved.")
```

## Chunk 3: Compare observed vs permuted LOPOCV Spearman

``` r
# observed LOPOCV summary from the true data
metrics_all <- read.csv(file.path(results_dir, "LOPOCV_summary.csv"))

# permutation matrix
spearman_null <- read.csv(
  file.path(results_dir, "spearman_LOPOCV_null_matrix.csv"),
  row.names = 1,
  check.names = FALSE
)

# make sure site order matches the permutation matrix
metrics_all <- metrics_all[match(rownames(spearman_null), metrics_all$site), ]

spearman_obs <- as.numeric(metrics_all$Spearman)

spearman_null_mat <- as.matrix(spearman_null)

# reorder rows to observed site order if needed
spearman_null_mat <- spearman_null_mat[match(metrics_all$site, rownames(spearman_null_mat)), ]

plot_df_spearman <- data.frame(
  spearman = c(spearman_obs, as.numeric(spearman_null_mat)),
  type = c(
    rep("Observed", length(spearman_obs)),
    rep("Permuted", length(as.numeric(spearman_null_mat)))
  )
)

ggplot(plot_df_spearman, aes(x = spearman, fill = type)) +
  geom_density(alpha = 0.5, color = NA) +
  coord_cartesian(xlim = c(-1, 1)) +
  theme_minimal() +
  labs(
    title = "Observed vs. Permuted Spearman Correlation (LOPOCV)",
    x = "Spearman Correlation",
    y = "Density"
  )
```

![](../figures/knitted_mds/plot-spearman-lopocv-1.png)<!-- -->

## Chunk 4: Empirical p-value

GOAL: Calculate an empirical p-value comparing observed vs permuted
Spearman correlation.

``` r
# Define region based on SiteMajCluster
Gff <- read.csv(file.path(input_dir, "Gff_11loci_allsites_indinfo.txt"),
                header = TRUE, sep = "\t")
north_sites <- unique(Gff$SiteCode[Gff$SiteMajCluster == "north"])

# Load observed Spearman
metrics_all <- read.csv(file.path(results_dir, "LOPOCV_summary.csv"))
spearman_obs_df <- data.frame(
  site = metrics_all$site,
  spearman = metrics_all$Spearman
)
spearman_obs_df$cluster <- ifelse(spearman_obs_df$site %in% north_sites, "north", "south")
spearman_obs_df$type <- "Observed"

spearman_null_raw <- read.csv(
  file.path(results_dir, "spearman_LOPOCV_null_matrix.csv"),
  row.names = 1,
  check.names = FALSE
)
spearman_null_mat <- as.matrix(spearman_null_raw)

# align rows to observed site order
spearman_null_mat <- spearman_null_mat[match(metrics_all$site, rownames(spearman_null_mat)), ]

obs_mean_spearman <- mean(spearman_obs_df$spearman, na.rm = TRUE)
null_mean_spearman <- apply(spearman_null_mat, 2, mean, na.rm = TRUE)

# one-sided empirical p-value: observed mean Spearman greater than null
p_empirical <- (sum(null_mean_spearman >= obs_mean_spearman, na.rm = TRUE) + 1) /
  (sum(!is.na(null_mean_spearman)) + 1)

cat("Observed mean Spearman:", round(obs_mean_spearman, 4), "\n")
```

    ## Observed mean Spearman: 0.873

``` r
cat("Null mean Spearman range:", round(range(null_mean_spearman, na.rm = TRUE), 4), "\n")
```

    ## Null mean Spearman range: -0.0817 0.067

``` r
cat("Empirical one-sided p-value:", round(p_empirical, 4), "\n\n")
```

    ## Empirical one-sided p-value: 0.0099

``` r
obs_north <- spearman_obs_df$spearman[spearman_obs_df$cluster == "north"]
obs_south <- spearman_obs_df$spearman[spearman_obs_df$cluster == "south"]

null_north_mat <- spearman_null_mat[rownames(spearman_null_mat) %in% north_sites, , drop = FALSE]
null_south_mat <- spearman_null_mat[!rownames(spearman_null_mat) %in% north_sites, , drop = FALSE]

obs_mean_north <- mean(obs_north, na.rm = TRUE)
obs_mean_south <- mean(obs_south, na.rm = TRUE)

null_mean_north <- apply(null_north_mat, 2, mean, na.rm = TRUE)
null_mean_south <- apply(null_south_mat, 2, mean, na.rm = TRUE)

p_north <- (sum(null_mean_north >= obs_mean_north, na.rm = TRUE) + 1) /
  (sum(!is.na(null_mean_north)) + 1)
p_south <- (sum(null_mean_south >= obs_mean_south, na.rm = TRUE) + 1) /
  (sum(!is.na(null_mean_south)) + 1)

cat("Observed mean Spearman, north:", round(obs_mean_north, 4), "\n")
```

    ## Observed mean Spearman, north: 0.869

``` r
cat("Empirical p-value, north:", round(p_north, 4), "\n\n")
```

    ## Empirical p-value, north: 0.0099

``` r
cat("Observed mean Spearman, south:", round(obs_mean_south, 4), "\n")
```

    ## Observed mean Spearman, south: 0.8774

``` r
cat("Empirical p-value, south:", round(p_south, 4), "\n")
```

    ## Empirical p-value, south: 0.0099

## Chunk 5: Visualize permutated vs observed LOPOCV side-by-side

``` r
# Define region based on SiteMajCluster
Gff <- read.csv(file.path(input_dir, "Gff_11loci_allsites_indinfo.txt"),
                header = TRUE, sep = "\t")
north_sites <- unique(Gff$SiteCode[Gff$SiteMajCluster == "north"])

# Load observed Spearman
metrics_all <- read.csv(file.path(results_dir, "LOPOCV_summary.csv"))
spearman_obs_df <- data.frame(
  site = metrics_all$site,
  spearman = metrics_all$Spearman
)
spearman_obs_df$cluster <- ifelse(spearman_obs_df$site %in% north_sites, "north", "south")
spearman_obs_df$type <- "Observed"

# Load permuted Spearman matrix
spearman_null <- read.csv(
  file.path(results_dir, "spearman_LOPOCV_null_matrix.csv"),
  row.names = 1,
  check.names = FALSE
)

spearman_null_df <- data.frame(
  site = rownames(spearman_null),
  cluster = ifelse(rownames(spearman_null) %in% north_sites, "north", "south"),
  spearman_null,
  check.names = FALSE
)

# Stack all permutation columns
spearman_null_long <- reshape2::melt(
  spearman_null_df,
  id.vars = c("site", "cluster"),
  variable.name = "perm",
  value.name = "spearman"
)
spearman_null_long$type <- "Permuted"

# Combine
spearman_plot_df <- rbind(
  spearman_obs_df[, c("site", "cluster", "spearman", "type")],
  spearman_null_long[, c("site", "cluster", "spearman", "type")]
)

# Plot
ggplot(spearman_plot_df, aes(x = spearman, fill = interaction(type, cluster))) +
  geom_density(
    data = subset(spearman_plot_df, type == "Observed"),
    aes(y = after_stat(count)),
    alpha = 0.6,
    color = NA,
    position = "identity"
  ) +
  geom_density(
    data = subset(spearman_plot_df, type == "Permuted"),
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
      "Permuted.north" = "#4A545E",
      "Permuted.south" = "#5A4F49"
      #"Permuted.north" = "#2c3e50",
      #"Permuted.south" = "#4e342e"
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
    title = "Observed vs. Permuted Spearman by cluster",
    x = "Spearman's r (LOPOCV)",
    y = "Density"
  )
```

![](../figures/knitted_mds/plot-side-by-side-1.png)<!-- -->

## Chunk 6: Spatial included, all side-by-side

Now with the Spatial Spearman’s R included… “north” = “\#2A4F6E” “south”
= “\#7A4A2A”

``` r
results_dir <- "/uufs/chpc.utah.edu/common/home/saarman-group1/uganda-tsetse-LG/results/"
input_dir <- "../input"
fig_dir <- file.path(results_dir, "figures_pub")

# Run just once to save high quality png for publication
png(file.path(fig_dir, "Fig_LOPOCV_obs_vs_permuted.png"),
    width = 7, height = 2.8, units = "in", res = 600)

# Define region based on SiteMajCluster
Gff <- read.csv(
  file.path(input_dir, "Gff_11loci_allsites_indinfo.txt"),
  header = TRUE, sep = "\t"
)
north_sites <- unique(Gff$SiteCode[Gff$SiteMajCluster == "north"])

# 1. Non-spatial LOPOCV Spearman
metrics_obs <- read.csv(file.path(results_dir, "LOPOCV_summary.csv"))

obs_df <- data.frame(
  site = metrics_obs$site,
  cluster = ifelse(metrics_obs$site %in% north_sites, "north", "south"),
  spearman = metrics_obs$Spearman,
  method = "Observed"
)

# 2. Null LOPOCV Spearman from permutations
spearman_null <- read.csv(
  file.path(results_dir, "spearman_LOPOCV_null_matrix.csv"),
  row.names = 1,
  check.names = FALSE
)

perm_df <- data.frame(
  site = rownames(spearman_null),
  cluster = ifelse(rownames(spearman_null) %in% north_sites, "north", "south"),
  spearman_null,
  check.names = FALSE
)

perm_long <- reshape2::melt(
  perm_df,
  id.vars = c("site", "cluster"),
  variable.name = "perm",
  value.name = "spearman"
) %>%
  mutate(method = "Permuted")

# 3. Spatial LOPOCV Spearman
metrics_spatial <- read.csv(
  file.path(results_dir, "spatial_LOPOCV_LCPsum_k1_summary.csv")
)

spatial_df <- data.frame(
  site = metrics_spatial$site,
  cluster = ifelse(metrics_spatial$site %in% north_sites, "north", "south"),
  spearman = metrics_spatial$Spearman,
  method = "Spatial"
)

# Combine all three
plot_df <- bind_rows(
  obs_df[, c("site", "cluster", "spearman", "method")],
  perm_long[, c("site", "cluster", "spearman", "method")],
  spatial_df[, c("site", "cluster", "spearman", "method")]
)

# Set method and cluster order for legend/display
plot_df$method <- factor(
  plot_df$method,
  levels = c("Observed", "Spatial", "Permuted")
)

plot_df$cluster <- factor(
  plot_df$cluster,
  levels = c("north", "south")
)

# Plot
ggplot(plot_df, aes(x = spearman, fill = interaction(method, cluster, sep = "."))) +

  geom_density(
    data = subset(plot_df, method == "Observed"),
    aes(y = after_stat(count)),
    alpha = 0.6,
    color = NA,
    position = "identity"
  ) +

  geom_density(
    data = subset(plot_df, method == "Spatial"),
    aes(y = after_stat(count)),
    alpha = 0.5,
    color = NA,
    position = "identity"
  ) +

  geom_density(
    data = subset(plot_df, method == "Permuted"),
    aes(y = after_stat(count / 100)),
    alpha = 0.3,
    color = NA,
    position = "identity"
  ) +

  scale_fill_manual(
    name = "Prediction type",
    breaks = c(
      "Observed.north", "Observed.south",
      "Spatial.north",  "Spatial.south",
      "Permuted.north", "Permuted.south"
    ),
    values = c(
      "Observed.north" = "#1f78b4",
      "Observed.south" = "#e66101",
      "Spatial.north"  = "#2A4F6E",
      "Spatial.south"  = "#7A4A2A",
      "Permuted.north" = "#4A545E",
      "Permuted.south" = "#5A4F49"
    ),
    labels = c(
      "Observed.north" = "RF-predicted North",
      "Observed.south" = "RF-predicted South",
      "Spatial.north"  = "LCP-predicted North",
      "Spatial.south"  = "LCP-predicted South",
      "Permuted.north" = "Null-predicted North",
      "Permuted.south" = "Null-predicted South"
    )
  ) +
  coord_cartesian(xlim = c(-1, 1)) +
  theme_minimal() +
  labs(
    title = "Fold-level LOPOCV Spearman by prediction type and cluster",
    x = "Spearman's r",
    y = "Count"
  )


dev.off()
```

    ## png 
    ##   2

## TODO: Chunk 7: Delete RDS when completely done

``` r
unlink(perm_model_dir, recursive = TRUE)
```
