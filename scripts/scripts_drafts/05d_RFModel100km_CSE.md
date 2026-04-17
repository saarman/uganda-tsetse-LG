5a. RF model local (\<100 km) Raw CSE (LC lakes paths)
================
Norah Saarman
2026-03-09

- [Inputs](#inputs)
- [1. Prepare the data](#1-prepare-the-data)
- [2. Build full Random Forest model](#2-build-full-random-forest-model)
- [3. Mean-only pruning](#3-mean-only-pruning)
  - [Compare mean versus median versus
    mode:](#compare-mean-versus-median-versus-mode)
  - [PCA of mean env predictor variables, patterns of
    autocorrelation?](#pca-of-mean-env-predictor-variables-patterns-of-autocorrelation)
  - [(Optional) Prune more variables after narrowing to mean
    only?](#optional-prune-more-variables-after-narrowing-to-mean-only)
  - [(Optional) Scaling CSE](#optional-scaling-cse)
- [4. Final full model - (Tune random forest with mean
  variables)](#4-final-full-model---tune-random-forest-with-mean-variables)
  - [Choose final predictor
    variables](#choose-final-predictor-variables)
  - [Full random forest model with raw
    CSE](#full-random-forest-model-with-raw-cse)
    - [Load saved (raw CSE) model](#load-saved-raw-cse-model)
  - [(Optional) Compare full and full tuned models (mean-only
    predictors)](#optional-compare-full-and-full-tuned-models-mean-only-predictors)
- [5. Project predicted values from full CSE
  model](#5-project-predicted-values-from-full-cse-model)
  - [Build Projection](#build-projection)
  - [Plot predicted CSE](#plot-predicted-cse)
- [6. Scale and plot predicted connectivity (CSE) and
  SDM](#6-scale-and-plot-predicted-connectivity-cse-and-sdm)
  - [Scale 0-1, habitat suitability and inverse of predicted
    connectivity
    (CSE)](#scale-0-1-habitat-suitability-and-inverse-of-predicted-connectivity-cse)
  - [Plot scaled predicted CSE and
    SDM](#plot-scaled-predicted-cse-and-sdm)
- [6. Variable importance plots](#6-variable-importance-plots)
  - [Percent Improvement MSE](#percent-improvement-mse)
  - [Node Purity](#node-purity)

RStudio Configuration:  
- **R version:** R 4.4.0 (Geospatial packages)  
- **Number of cores:** 4 (up to 32 available)  
- **Account:** saarman-np  
- **Partition:** saarman-np (allows multiple simultaneous jobs
automatically now)  
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

# This was added only after completing LOPOCV on full model...
# Filter out western outlier "50-KB" 
V.table <- V.table %>%
  filter(Var1 != "50-KB", Var2 != "50-KB")

# This is only for the RFModel100km runs...
# Filter out pairs with geographic distance >100 km
# based on results from Mantel Correlogram
V.table <- V.table %>%
  filter(pix_dist < 100)

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

![](05d_RFModel100km_CSE_files/figure-gfm/prep-1.png)<!-- -->

``` r
# Extract groups of variables by suffix
mean_vars   <- grep("_mean$", names(rf_data), value = TRUE)
median_vars <- grep("_median$", names(rf_data), value = TRUE)
mode_vars   <- grep("_mode$", names(rf_data), value = TRUE)
```

# 2. Build full Random Forest model

Note: Marked eval = FALSE to avoid re-running on knit

``` r
# Build full RF model
set.seed(1234)  # ensures reproducibility
rf_full <- randomForest(
  CSEdistance ~ .,
  data = rf_data,
  importance = TRUE,
  ntree = 500
)

print(rf_full)

importance(rf_full)
```

# 3. Mean-only pruning

## Compare mean versus median versus mode:

Note: Marked eval = FALSE to avoid re-running on knit

``` r
# Extract groups of variables by suffix
mean_vars   <- grep("_mean$", names(rf_data), value = TRUE)
median_vars <- grep("_median$", names(rf_data), value = TRUE)
mode_vars   <- grep("_mode$", names(rf_data), value = TRUE)

# Always include geographic distance
common_var <- "pix_dist"

# Build and run each model
set.seed(123438972)  # ensures reproducibility
rf_mean <- randomForest(CSEdistance ~ ., data = rf_data[, c("CSEdistance", common_var, mean_vars)], ntree = 500, importance = TRUE)
rf_median <- randomForest(CSEdistance ~ ., data = rf_data[, c("CSEdistance", common_var, median_vars)], ntree = 500, importance = TRUE)
rf_mode <- randomForest(CSEdistance ~ ., data = rf_data[, c("CSEdistance", common_var, mode_vars)], ntree = 500, importance = TRUE)

# Compare performance
c(mean = rf_mean$rsq[500] * 100,
  median = rf_median$rsq[500] * 100,
  mode = rf_mode$rsq[500] * 100)
```

Including mean of env variable along least cost paths performs the best,
adding median and mode does not greatly improve the model and increases
risks of over fitting…

## PCA of mean env predictor variables, patterns of autocorrelation?

``` r
# Run PCA on the subset of variables
pca_res <- prcomp(rf_data[mean_vars], scale. = TRUE)

# Quick summary of variance explained
summary(pca_res)
```

    ## Importance of components:
    ##                           PC1    PC2    PC3     PC4     PC5     PC6     PC7
    ## Standard deviation     3.7578 1.8596 1.5946 1.11662 0.96397 0.72845 0.63985
    ## Proportion of Variance 0.5884 0.1441 0.1059 0.05195 0.03872 0.02211 0.01706
    ## Cumulative Proportion  0.5884 0.7325 0.8384 0.89037 0.92908 0.95119 0.96825
    ##                            PC8    PC9    PC10    PC11   PC12    PC13    PC14
    ## Standard deviation     0.60247 0.4490 0.35693 0.16078 0.1385 0.09374 0.08560
    ## Proportion of Variance 0.01512 0.0084 0.00531 0.00108 0.0008 0.00037 0.00031
    ## Cumulative Proportion  0.98338 0.9918 0.99709 0.99816 0.9990 0.99933 0.99963
    ##                           PC15    PC16    PC17    PC18    PC19    PC20     PC21
    ## Standard deviation     0.06319 0.05087 0.03208 0.02532 0.01608 0.01420 0.008238
    ## Proportion of Variance 0.00017 0.00011 0.00004 0.00003 0.00001 0.00001 0.000000
    ## Cumulative Proportion  0.99980 0.99991 0.99995 0.99998 0.99999 1.00000 1.000000
    ##                            PC22     PC23      PC24
    ## Standard deviation     0.004843 0.000527 3.668e-07
    ## Proportion of Variance 0.000000 0.000000 0.000e+00
    ## Cumulative Proportion  1.000000 1.000000 1.000e+00

``` r
# Scree plot of variance explained
ve <- (pca_res$sdev^2) / sum(pca_res$sdev^2)
cumve <- cumsum(ve)
df_ve <- data.frame(PC = seq_along(ve), Var = ve*100, Cum = cumve*100)
ggplot(df_ve, aes(PC, Var)) +
  geom_line() + geom_point() +
  geom_text(aes(label = sprintf("%.1f", Var)), vjust = -0.6, size = 3) +
  labs(x = "Principal component", y = "Variance explained (%)") +
  theme_classic()
```

![](05d_RFModel100km_CSE_files/figure-gfm/unnamed-chunk-1-1.png)<!-- -->

``` r
# PCA biplot with individuals (rows) and variables (arrows)
fviz_pca_biplot(pca_res,
                repel = TRUE, # avoid text overlap
                col.var = "steelblue", # variables
                col.ind = "gray30"    # individuals
)
```

![](05d_RFModel100km_CSE_files/figure-gfm/unnamed-chunk-1-2.png)<!-- -->

``` r
# PCA variables plot (correlation circle)
load <- as.data.frame(pca_res$rotation[, 1:2])
load$var <- rownames(load)
# circle helper
circle <- data.frame(
  x = cos(seq(0, 2*pi, length.out = 200)),
  y = sin(seq(0, 2*pi, length.out = 200))
)

ggplot() +
  geom_path(data = circle, aes(x, y), linewidth = 0.3) +
  geom_segment(data = load, aes(x = 0, y = 0, xend = PC1, yend = PC2),
               arrow = arrow(length = unit(0.15, "cm")), linewidth = 0.3) +
  geom_text(data = load, aes(PC1, PC2, label = var),
            hjust = 0.5, vjust = -0.5, size = 3) +
  coord_equal(xlim = c(-1,1), ylim = c(-1,1)) +
  labs(x = "PC1", y = "PC2") +
  theme_classic()
```

![](05d_RFModel100km_CSE_files/figure-gfm/unnamed-chunk-1-3.png)<!-- -->

``` r
fviz_pca_var(pca_res,
             col.var = "contrib", # color by contribution
             gradient.cols = c("gray70", "steelblue", "darkblue"),
             repel = TRUE
)
```

![](05d_RFModel100km_CSE_files/figure-gfm/unnamed-chunk-1-4.png)<!-- -->

## (Optional) Prune more variables after narrowing to mean only?

``` r
# Plot variable importance
par(mar = c(5, 10, 2, 2))  # bottom, left, top, right
varImpPlot(rf_mean, main = "Mean Model Importance",cex = 0.6, pch = 19)

# Rank variables by %IncMSE (from tuned model)
imp <- importance(rf_mean)[, "%IncMSE"]
imp <- sort(imp, decreasing = TRUE)

# Multiple runs with N top predictors
# Store results
prune_results <- list()
n_list <- c(5:length(imp))

for (n in n_list) {
  top_vars <- names(imp)[1:n]
  formula_n <- as.formula(paste("CSEdistance ~", paste(top_vars, collapse = " + ")))
  
  set.seed(1234783645)
  rf_n <- randomForest(
    formula = formula_n,
    data = rf_data,
    ntree = 500,
    importance = TRUE
  )
  
  prune_results[[paste0("Top", n)]] <- rf_n
}

sapply(prune_results, function(mod) {
  c(OOB_MSE = mod$mse[500], VarExpl = mod$rsq[500] * 100)
})
```

% Variance Explained increases rapidly up to around 18 variables, after
which it plateaus.

OOB MSE decreases quickly early on, with minimal gains beyond the top
~18 predictors.

However, there are not too many more than 18 total, and many of the env
variables have similar node purity, so given the special nature of the
top 2 (geo dist and sampling density), I will leave all mean env
variables.

## (Optional) Scaling CSE

Note: I did also try scaling CSE before runnning RF as a different idea,
but performance of the full model did not improve, so no need to scale
CSE before modeling.

# 4. Final full model - (Tune random forest with mean variables)

## Choose final predictor variables

``` r
# estimate mean sampling density
mean(V.table$samp_20km_mean, na.rm = TRUE)
```

    ## [1] 1.550707e-11

``` r
# Define site list
sites <- sort(unique(c(V.table$Var1, V.table$Var2)))

# How many rows of data for each?
table(V.table$Pop1_cluster)
```

    ## 
    ## north south 
    ##   124    70

``` r
# How many unique sites?
length(sites)
```

    ## [1] 66

``` r
# Choose predictors for RF model (adjust names if necessary)
predictor_vars <- c("BIO1_mean","BIO2_mean","BIO3_mean","BIO4_mean", "BIO5_mean","BIO6_mean","BIO7_mean", "BIO8S_mean", "BIO9S_mean","BIO10S_mean", "BIO11S_mean","BIO12_mean", "BIO13_mean","BIO14_mean","BIO15_mean","BIO16S_mean","BIO17S_mean", "BIO18S_mean","BIO19S_mean","slope_mean","alt_mean", "lakes_mean","riv_3km_mean", "samp_20km_mean","pix_dist")

# Filter to modeling-relevant columns only
rf_mean_data <- V.table[, c("CSEdistance", predictor_vars)]

# Rename predictors by removing "_mean" for later projections
names(rf_mean_data) <- gsub("_mean$", "", names(rf_mean_data))
```

## Full random forest model with raw CSE

``` r
# Tune mtry (number of variables tried at each split)
set.seed(92834)
rf_mean_full_tuned <- tuneRF(
  x = rf_mean_data[, -1],   # exclude response variable
  y = rf_mean_data$CSEdistance,
  ntreeTry = 500,
  stepFactor = 1.5,         # factor by which mtry is increased/decreased
  improve = 0.01,           # minimum improvement to continue search
  trace = TRUE,             # print progress
  plot = TRUE,              # plot OOB error vs mtry
  doBest = TRUE,             # return the model with lowest OOB error
  importance = TRUE
)
```

    ## mtry = 8  OOB error = 0.002158435 
    ## Searching left ...
    ## mtry = 6     OOB error = 0.002137213 
    ## 0.009832112 0.01 
    ## Searching right ...
    ## mtry = 12    OOB error = 0.002156232 
    ## 0.001020703 0.01

![](05d_RFModel100km_CSE_files/figure-gfm/tune_run-1.png)<!-- -->

``` r
print(rf_mean_full_tuned)
```

    ## 
    ## Call:
    ##  randomForest(x = x, y = y, mtry = res[which.min(res[, 2]), 1],      importance = TRUE) 
    ##                Type of random forest: regression
    ##                      Number of trees: 500
    ## No. of variables tried at each split: 6
    ## 
    ##           Mean of squared residuals: 0.002116133
    ##                     % Var explained: 56.85

``` r
importance(rf_mean_full_tuned)
```

    ##             %IncMSE IncNodePurity
    ## BIO1       9.484682    0.03688015
    ## BIO2      13.940412    0.06798981
    ## BIO3      16.424771    0.07749905
    ## BIO4      11.699304    0.04654505
    ## BIO5       8.836165    0.02778484
    ## BIO6       7.398376    0.02767627
    ## BIO7      14.969535    0.06810512
    ## BIO8S     10.315028    0.02800927
    ## BIO9S      9.832812    0.03596521
    ## BIO10S    10.731887    0.02797769
    ## BIO11S    11.909853    0.04384236
    ## BIO12      7.043478    0.01922448
    ## BIO13      8.683840    0.02850986
    ## BIO14     10.325783    0.02326266
    ## BIO15      7.561399    0.02057907
    ## BIO16S     7.787524    0.03256321
    ## BIO17S     9.222669    0.02968561
    ## BIO18S     7.130051    0.02297461
    ## BIO19S     9.637463    0.02680013
    ## slope      6.599890    0.01726052
    ## alt        8.609335    0.03036913
    ## lakes      7.879204    0.01811064
    ## riv_3km    7.614750    0.02186096
    ## samp_20km 15.874590    0.09509855
    ## pix_dist  11.384204    0.04400029

``` r
varImpPlot(rf_mean_full_tuned)
```

![](05d_RFModel100km_CSE_files/figure-gfm/tune_run-2.png)<!-- -->

``` r
# Save the tuned random forest model to disk
saveRDS(rf_mean_full_tuned, file = file.path(results_dir, "rf_mean_full_tuned_100km.rds"))
```

FYI: Later, to load the model back into R:
`rf_mean_full_tuned <- readRDS(file.path(results_dir, "rf_mean_full_tuned_100km.rds"))`

### Load saved (raw CSE) model

``` r
# load saved model
rf_full <- readRDS(file.path(results_dir, "rf_mean_full_tuned_100km.rds"))

#double check they look correct
print(rf_full)
```

    ## 
    ## Call:
    ##  randomForest(x = x, y = y, mtry = res[which.min(res[, 2]), 1],      importance = TRUE) 
    ##                Type of random forest: regression
    ##                      Number of trees: 500
    ## No. of variables tried at each split: 6
    ## 
    ##           Mean of squared residuals: 0.002116133
    ##                     % Var explained: 56.85

``` r
print(rf_full$importance)
```

    ##                %IncMSE IncNodePurity
    ## BIO1      7.181217e-04    0.03688015
    ## BIO2      1.115098e-03    0.06798981
    ## BIO3      9.334336e-04    0.07749905
    ## BIO4      3.818416e-04    0.04654505
    ## BIO5      3.681166e-04    0.02778484
    ## BIO6      1.466886e-04    0.02767627
    ## BIO7      8.948434e-04    0.06810512
    ## BIO8S     2.442712e-04    0.02800927
    ## BIO9S     3.199100e-04    0.03596521
    ## BIO10S    2.407870e-04    0.02797769
    ## BIO11S    6.517909e-04    0.04384236
    ## BIO12     1.435967e-04    0.01922448
    ## BIO13     1.868064e-04    0.02850986
    ## BIO14     2.253158e-04    0.02326266
    ## BIO15     1.960444e-04    0.02057907
    ## BIO16S    3.818734e-04    0.03256321
    ## BIO17S    2.260867e-04    0.02968561
    ## BIO18S    3.494583e-04    0.02297461
    ## BIO19S    2.365338e-04    0.02680013
    ## slope     9.749325e-05    0.01726052
    ## alt       3.134921e-04    0.03036913
    ## lakes     1.905998e-04    0.01811064
    ## riv_3km   1.030855e-04    0.02186096
    ## samp_20km 9.038861e-04    0.09509855
    ## pix_dist  3.232680e-04    0.04400029

## (Optional) Compare full and full tuned models (mean-only predictors)

Note: eval = FALSE to avoid re-running on knit

``` r
# Build full RF model
set.seed(10981234)  # ensures reproducibility
rf_mean_full <- randomForest(
  CSEdistance ~ .,
  data = rf_mean_data,
  importance = TRUE,
  ntree = 500
)

print(rf_mean_full)
importance(rf_mean_full)

print(rf_mean_full_tuned)
importance(rf_mean_full_tuned)

data.frame(
  Model = c("Full (default mtry)", paste("Tuned (mtry = ",rf_mean_full_tuned$mtry,")")),
  MSE = c(rf_mean_full$mse[rf_mean_full$ntree], rf_mean_full_tuned$mse[rf_mean_full_tuned$ntree]),
  Rsq = c(rf_mean_full$rsq[rf_mean_full$ntree], rf_mean_full_tuned$rsq[rf_mean_full_tuned$ntree])
)

# pad names to trick varImpPlot
rownames(rf_mean_full$importance) <- paste0("  ", rownames(rf_mean_full$importance), "  ")
rownames(rf_mean_full_tuned$importance) <- paste0("  ", rownames(rf_mean_full_tuned$importance), "  ")

# plot with varImpPlot
par(mar = c(5, 30, 2,2))  # bottom, left, top, right
varImpPlot(rf_mean_full, main = "Full Model Importance",cex = 0.6, pch = 19)
varImpPlot(rf_mean_full_tuned, main = "Tuned Full Model Importance",cex = 0.6, pch = 19)
```

Results from optional comparison: The tuned model performs slightly
better, but the gain may not be meaningful… however, it does confirm
that the model is stable and that the mean-only predictors carry strong
signal.

Top 18 mean-based predictors retain nearly all the explanatory power of
the original full model with 42 predictors.

# 5. Project predicted values from full CSE model

## Build Projection

``` r
# Load env stack with named layers
env <- stack(file.path(data_dir, "processed", "env_stack.grd"))

# Neutralize sampling layer to average
env$samp_20km <- mean(V.table$samp_20km_mean, na.rm = TRUE) #neutralize sampling bias

# Load rdf of final model
rf_predicted <- readRDS(file.path(results_dir, "rf_mean_full_tuned_100km.rds"))
rf_predicted
```

    ## 
    ## Call:
    ##  randomForest(x = x, y = y, mtry = res[which.min(res[, 2]), 1],      importance = TRUE) 
    ##                Type of random forest: regression
    ##                      Number of trees: 500
    ## No. of variables tried at each split: 6
    ## 
    ##           Mean of squared residuals: 0.002116133
    ##                     % Var explained: 56.85

``` r
prediction_raster <- predict(env, rf_predicted, type = "response")

# Write Prediction Raster to file
writeRaster(prediction_raster, file.path(results_dir,"fullRF_CSE.tif"), format = "GTiff", overwrite = TRUE)
```

## Plot predicted CSE

``` r
# Create base plot with viridis
plot(prediction_raster,
     col = viridis::magma(100),
     main = "Predicted CSE Distance",
     axes = FALSE,
     box = FALSE,
     legend.args = list(text = "CSE", side = 3, line = 1, cex = 1))

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

![](05d_RFModel100km_CSE_files/figure-gfm/plot-projection-1.png)<!-- -->

# 6. Scale and plot predicted connectivity (CSE) and SDM

## Scale 0-1, habitat suitability and inverse of predicted connectivity (CSE)

``` r
# Load raster layers
con_raster <- rast(file.path(results_dir, "fullRF_CSE.tif"))
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

## Plot scaled predicted CSE and SDM

``` r
# Plot Genetic Connectivity (inverse predicted values)
plot(con,
     col = rev(viridis::plasma(100)),  # high connectivity = dark
     main = "Genetic Connectivity (inverse predicted CSE)",
     axes = FALSE, box = FALSE,
     legend.args = list(text = "Connectivity", side = 2, line = 2.5, cex = 0.8))
plot(st_geometry(lakes), col = "black", border = NA, add = TRUE)
plot(st_geometry(uganda), border = "black", lwd = 0.25, add = TRUE)
```

![](05d_RFModel100km_CSE_files/figure-gfm/plot-sdm-con-1.png)<!-- -->

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

![](05d_RFModel100km_CSE_files/figure-gfm/plot-sdm-con-2.png)<!-- -->

``` r
# Plot with custom colors

# Custom palettes based on Bishop et al.
connectivity_colors <- colorRampPalette(c("#FFFF00", "#FFA500", "#FF4500", "#700E40", "#2E003E"))(100)
suitability_colors  <- colorRampPalette(c("white", "lightblue", "blue4"))(100)     # white → light blue → dark blue

# Plot Genetic Connectivity (inverse predicted) with custom colors
plot(con,
     col = connectivity_colors,
     main = "Genetic Connectivity (inverse predicted CSE)",
     axes = FALSE, box = FALSE,
     legend.args = list(text = "Connectivity", side = 2, line = 2.5, cex = 0.8))
plot(st_geometry(lakes), col = "black", border = NA, add = TRUE)
plot(st_geometry(uganda), border = "black", lwd = 0.25, add = TRUE)
```

![](05d_RFModel100km_CSE_files/figure-gfm/plot-sdm-con-3.png)<!-- -->

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

![](05d_RFModel100km_CSE_files/figure-gfm/plot-sdm-con-4.png)<!-- -->

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
  labs(x = NULL, title = "Variable Importance of Raw CSE Model") +
  theme_minimal()
```

![](05d_RFModel100km_CSE_files/figure-gfm/variable-imp-mse-1.png)<!-- -->

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
  labs(x = NULL, title = "Variable Importance of Raw CSE Model") +
  theme_minimal()
```

![](05d_RFModel100km_CSE_files/figure-gfm/variable-imp-nodepurity-1.png)<!-- -->
