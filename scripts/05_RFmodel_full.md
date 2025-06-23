RF model full – LC lakes paths
================
Norah Saarman
2025-06-17

- [Setup](#setup)
- [Inputs](#inputs)
- [1. Prepare the data](#1-prepare-the-data)
- [2. Build full Random Forest model](#2-build-full-random-forest-model)
- [2. Prune variables?](#2-prune-variables)
  - [Compare mean versus median versus
    mode:](#compare-mean-versus-median-versus-mode)
  - [Prune more variables after narrowing to mean
    only?](#prune-more-variables-after-narrowing-to-mean-only)
  - [Top 18 predictors](#top-18-predictors)
- [3. Tune random forest with chosen
  variables](#3-tune-random-forest-with-chosen-variables)
  - [Compare full and full tuned models (top 18 mean-only
    predictors)](#compare-full-and-full-tuned-models-top-18-mean-only-predictors)
- [4. Projection of model](#4-projection-of-model)

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
library(randomForest)
library(doParallel)
library(raster)
library(sf)
library(viridis)
library(dplyr)

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
  \# Outputs  
- `../results_dir/rf_mean18_tuned.rds` - Full RF model outpu, readRDS()
  to load  
- `../results_dir/predicted_CSEdistance.tif` - Projection of full RF
  model, with pix_dist and samp_20km nuetralized by using uniform layers
  in projection

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
```

# 2. Build full Random Forest model

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
```

    ## 
    ## Call:
    ##  randomForest(formula = CSEdistance ~ ., data = rf_data, importance = TRUE,      ntree = 500) 
    ##                Type of random forest: regression
    ##                      Number of trees: 500
    ## No. of variables tried at each split: 24
    ## 
    ##           Mean of squared residuals: 0.001119129
    ##                     % Var explained: 86.2

``` r
importance(rf_full)
```

    ##                     %IncMSE IncNodePurity
    ## pix_dist         58.3012131  3.0935298531
    ## BIO1_mean         9.0020495  0.0301285131
    ## BIO2_mean         6.4098370  0.0473834223
    ## BIO3_mean        13.6026451  0.2721208278
    ## BIO4_mean        10.2607730  0.0613038622
    ## BIO5_mean        10.4818330  0.0291173353
    ## BIO6_mean        11.8583425  0.1575919110
    ## BIO7_mean        10.7247498  0.0484384046
    ## BIO8S_mean        9.1198968  0.0389941449
    ## BIO9S_mean        9.1545944  0.0757165381
    ## BIO10S_mean      10.5581937  0.0348492182
    ## BIO11S_mean      12.6650592  0.0999023193
    ## BIO12_mean       10.2105738  0.0420781559
    ## BIO13_mean       10.3447647  0.0557322008
    ## BIO14_mean        8.6877132  0.1224180146
    ## BIO15_mean        9.2166452  0.0942930637
    ## BIO16S_mean       9.3860539  0.0266651745
    ## BIO17S_mean       5.5162695  0.0419393147
    ## BIO18S_mean       5.1426282  0.0318266777
    ## BIO19S_mean       9.4523005  0.0313136617
    ## alt_mean          9.6483998  0.0391622442
    ## slope_mean       10.8249810  0.0391778902
    ## riv_3km_mean     13.2872838  0.0527211941
    ## samp_20km_mean   19.3440414  1.0712456972
    ## lakes_mean        7.4531090  0.0780686943
    ## BIO1_median       9.1525088  0.0410832893
    ## BIO2_median       9.3369590  0.0473558188
    ## BIO3_median      14.0929216  0.0892901654
    ## BIO4_median      10.1173720  0.0507536064
    ## BIO5_median       9.6663029  0.0497769331
    ## BIO6_median      12.7588988  0.2462962100
    ## BIO7_median      10.6655843  0.0495553419
    ## BIO8S_median      9.9223161  0.0337600533
    ## BIO9S_median      9.9244398  0.0804881296
    ## BIO10S_median     9.5509941  0.0373292693
    ## BIO11S_median    10.8215075  0.0844560970
    ## BIO12_median     10.1568111  0.0327408721
    ## BIO13_median     11.7315217  0.0719496005
    ## BIO14_median      4.7695670  0.0209310192
    ## BIO15_median     13.4686933  0.0706187410
    ## BIO16S_median     8.9510772  0.0425605141
    ## BIO17S_median    10.9025714  0.0276578250
    ## BIO18S_median     7.7717184  0.0224679158
    ## BIO19S_median     9.2931810  0.0334663036
    ## alt_median       11.5067482  0.0603982245
    ## slope_median     10.2245463  0.0471151689
    ## riv_3km_median    7.1750996  0.0521492669
    ## samp_20km_median 15.8059008  0.5984020935
    ## lakes_median     -0.8027164  0.0002072142
    ## BIO1_mode        11.2642008  0.0261941733
    ## BIO2_mode        10.0272409  0.0582292777
    ## BIO3_mode        13.2646204  0.0893638675
    ## BIO4_mode        11.8018707  0.0319873963
    ## BIO5_mode         8.5177924  0.0239607111
    ## BIO6_mode         7.8889344  0.0905625124
    ## BIO7_mode         9.7704289  0.0303038567
    ## BIO8S_mode        7.5655875  0.0213594257
    ## BIO9S_mode        9.0357805  0.0330665378
    ## BIO10S_mode      12.3764817  0.0458989549
    ## BIO11S_mode      10.8438312  0.0548054849
    ## BIO12_mode        8.8893019  0.0300983290
    ## BIO13_mode        9.4327509  0.0315015252
    ## BIO14_mode        9.4888215  0.0268860615
    ## BIO15_mode        9.1871791  0.1212897856
    ## BIO16S_mode       8.6535206  0.0262811613
    ## BIO17S_mode       7.9266488  0.0311551856
    ## BIO18S_mode       8.6996598  0.0221932434
    ## BIO19S_mode      10.2959906  0.0395754563
    ## alt_mode          9.2601975  0.0278142951
    ## slope_mode        8.3837670  0.0458055496
    ## riv_3km_mode      7.3381090  0.0527243831
    ## samp_20km_mode   11.6645843  0.2253255309
    ## lakes_mode       -1.0010015  0.0000351668

# 2. Prune variables?

## Compare mean versus median versus mode:

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

    ##     mean   median     mode 
    ## 85.52060 84.86890 84.51402

Including mean of env variable along least cost paths performs the best,
adding median and mode does not greatly improve the model and increases
risks of over fitting…

Could add in minimum, maximum, standard deviation, range, turnover
(analogous to slope) etc… later!

## Prune more variables after narrowing to mean only?

``` r
# Plot variable importance
par(mar = c(5, 10, 2, 2))  # bottom, left, top, right
varImpPlot(rf_mean, main = "Mean Model Importance",cex = 0.6, pch = 19)
```

![](05_RFmodel_full_files/figure-gfm/mean-1.png)<!-- -->

``` r
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

    ##                 Top5         Top6         Top7         Top8         Top9
    ## OOB_MSE  0.001407692  0.001338065  0.001307951  0.001236427  0.001206078
    ## VarExpl 82.641529917 83.500115170 83.871456681 84.753429692 85.127669122
    ##                Top10        Top11        Top12        Top13        Top14
    ## OOB_MSE  0.001205162  0.001188499  0.001177653  0.001183068  0.001183889
    ## VarExpl 85.138958253 85.344431928 85.478174323 85.411404269 85.401277243
    ##                Top15        Top16       Top17        Top18        Top19
    ## OOB_MSE  0.001179652  0.001176377  0.00118273  0.001167503  0.001174329
    ## VarExpl 85.453528639 85.493909572 85.41557043 85.603339051 85.519159798
    ##                Top20        Top21        Top22        Top23        Top24
    ## OOB_MSE  0.001165533  0.001152764  0.001158101  0.001156222  0.001152734
    ## VarExpl 85.627630199 85.785081355 85.719270871 85.742451545 85.785451021
    ##                Top25
    ## OOB_MSE  0.001165264
    ## VarExpl 85.630943128

% Variance Explained increases rapidly up to around 18 variables, after
which it plateaus.

OOB MSE decreases quickly early on, with minimal gains beyond the top
~18 predictors.

## Top 18 predictors

``` r
# Get variable importance
var_imp <- importance(rf_mean)[, "%IncMSE"]

# Sort and get names of top 18 predictors
top18_vars <- names(sort(var_imp, decreasing = TRUE))[1:18]

rf_top18_data <- rf_data[, c("CSEdistance", top18_vars)]
```

# 3. Tune random forest with chosen variables

``` r
# Subset data for the mean-only model
#rf_mean_data <- rf_data[, c("CSEdistance", common_var, mean_vars)]
rf_mean_data <- rf_data[, c("CSEdistance", top18_vars)]

# Rename predictors by removing "_mean" for later projections
names(rf_mean_data) <- gsub("_mean$", "", names(rf_mean_data))

# Build full RF model
set.seed(10981234)  # ensures reproducibility
rf_mean18 <- randomForest(
  CSEdistance ~ .,
  data = rf_mean_data,
  importance = TRUE,
  ntree = 500
)

print(rf_mean18)
```

    ## 
    ## Call:
    ##  randomForest(formula = CSEdistance ~ ., data = rf_mean_data,      importance = TRUE, ntree = 500) 
    ##                Type of random forest: regression
    ##                      Number of trees: 500
    ## No. of variables tried at each split: 6
    ## 
    ##           Mean of squared residuals: 0.001169491
    ##                     % Var explained: 85.58

``` r
importance(rf_mean18)
```

    ##            %IncMSE IncNodePurity
    ## pix_dist  71.65795     3.5783744
    ## samp_20km 26.61125     1.3569066
    ## BIO13     17.21835     0.2040798
    ## BIO8S     19.43221     0.1802590
    ## BIO3      21.65268     0.5142565
    ## riv_3km   18.67469     0.1963006
    ## BIO9S     17.19416     0.2310162
    ## BIO7      21.73728     0.2253818
    ## BIO5      17.93580     0.1972963
    ## BIO11S    21.52975     0.2382724
    ## slope     17.26434     0.1276901
    ## alt       19.67267     0.1843154
    ## BIO10S    17.90859     0.1263081
    ## BIO6      17.99176     0.3491119
    ## BIO14     16.41491     0.3824678
    ## BIO15     17.00173     0.3248504
    ## BIO1      17.83020     0.1828608
    ## BIO12     17.34726     0.1383382

``` r
# Tune mtry (number of variables tried at each split)
set.seed(92834567)
rf_mean18_tuned <- tuneRF(
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

    ## mtry = 6  OOB error = 0.001171661 
    ## Searching left ...
    ## mtry = 4     OOB error = 0.001200961 
    ## -0.02500702 0.01 
    ## Searching right ...
    ## mtry = 9     OOB error = 0.001179634 
    ## -0.006804841 0.01

![](05_RFmodel_full_files/figure-gfm/tune-1.png)<!-- -->

``` r
# Save the tuned random forest model to disk
saveRDS(rf_mean18_tuned, file = file.path(results_dir, "rf_mean18_tuned.rds"))

# Preserve as-is for projection
rf_final <- rf_mean18_tuned
```

FYI: Later, to load the model back into R:
`rf_mean18_tuned <- readRDS(file.path(results_dir, "rf_mean18_tuned.rds"))`

## Compare full and full tuned models (top 18 mean-only predictors)

``` r
print(rf_mean18)
```

    ## 
    ## Call:
    ##  randomForest(formula = CSEdistance ~ ., data = rf_mean_data,      importance = TRUE, ntree = 500) 
    ##                Type of random forest: regression
    ##                      Number of trees: 500
    ## No. of variables tried at each split: 6
    ## 
    ##           Mean of squared residuals: 0.001169491
    ##                     % Var explained: 85.58

``` r
print(rf_mean18_tuned)
```

    ## 
    ## Call:
    ##  randomForest(x = x, y = y, mtry = res[which.min(res[, 2]), 1],      importance = TRUE) 
    ##                Type of random forest: regression
    ##                      Number of trees: 500
    ## No. of variables tried at each split: 6
    ## 
    ##           Mean of squared residuals: 0.001168978
    ##                     % Var explained: 85.59

``` r
data.frame(
  Model = c("Full (default mtry)", paste("Tuned (mtry = ",rf_mean18_tuned$mtry,")")),
  MSE = c(rf_mean18$mse[rf_mean18$ntree], rf_mean18_tuned$mse[rf_mean18_tuned$ntree]),
  Rsq = c(rf_mean18$rsq[rf_mean18$ntree], rf_mean18_tuned$rsq[rf_mean18_tuned$ntree])
)
```

    ##                 Model         MSE       Rsq
    ## 1 Full (default mtry) 0.001169491 0.8557882
    ## 2  Tuned (mtry =  6 ) 0.001168978 0.8558514

``` r
# pad names to trick varImpPlot
rownames(rf_mean18$importance) <- paste0("  ", rownames(rf_mean18$importance), "  ")
rownames(rf_mean18_tuned$importance) <- paste0("  ", rownames(rf_mean18_tuned$importance), "  ")

# plot with varImpPlot
par(mar = c(5, 30, 2,2))  # bottom, left, top, right
varImpPlot(rf_mean18, main = "Full Model Importance",cex = 0.6, pch = 19)
```

![](05_RFmodel_full_files/figure-gfm/compare-1.png)<!-- -->

``` r
varImpPlot(rf_mean18_tuned, main = "Tuned Full Model Importance",cex = 0.6, pch = 19)
```

![](05_RFmodel_full_files/figure-gfm/compare-2.png)<!-- -->

The tuned model performs slightly better, but the gain may not be
meaningful… however, it does confirm that the model is stable and that
the mean-only predictors carry strong signal.

Top 18 mean-based predictors retain nearly all the explanatory power of
the original full model with 42 predictors.

# 4. Projection of model

``` r
# Load env stack
env <- stack(file.path(data_dir, "processed", "env_stack.grd"))

# Replace samp_20km with neutral constant raster, 
# retains sampling bias in the model but neutralize 
# sampling bias during projection
samp_uniform <- env[["samp_20km"]]
values(samp_uniform) <- mean(rf_data$samp_20km_mean, na.rm = TRUE)
env[["samp_20km"]] <- samp_uniform # Replace in raster stack

prediction_raster <- predict(env, rf_final, type = "response")

# Write Prediction Raster to file
#writeRaster(prediction_raster, file.path(results_dir,"predicted_CSEdistance.tif"), format = "GTiff", overwrite = TRUE)
```

``` r
# Create base plot with viridis
plot(prediction_raster,
     col = viridis::magma(100),
     main = "Predicted CSE Distance",
     axes = FALSE,
     box = FALSE,
     legend.args = list(text = "CSE Distance", side = 2, line = 2.5, cex = 0.8))

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

![](05_RFmodel_full_files/figure-gfm/plot-projection-1.png)<!-- -->
