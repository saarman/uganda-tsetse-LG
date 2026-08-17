LOPOCV: Random Forest (Internal) and Spatial Cross-Validation
================
Norah Saarman
2026-08-17

- [Overview of script](#overview-of-script)
- [Mean fold-level Spearman for full model (north vs
  south)](#mean-fold-level-spearman-for-full-model-north-vs-south)
- [Mean fold-level Spearman for distance-only model (north vs
  south)](#mean-fold-level-spearman-for-distance-only-model-north-vs-south)

RStudio Configuration:  
- **R version:** R 4.4.0 (Geospatial packages)  
- **Number of cores:** 4 (up to 32 available)  
- **Account:** saarman-np  
- **Partition:** saarman-np (now auto allows multiple simultaneous
jobs)  
- **Memory per job:** 100G (cluster limit: 1000G total; avoid exceeding
half)

## Overview of script

This script simply calculates the mean pooled out-of-fold R² and
Spearman’s r for north and south separately, using pre-exisiting results
files.

# Mean fold-level Spearman for full model (north vs south)

``` r
input_dir <- "../input"
data_dir  <- "/uufs/chpc.utah.edu/common/home/saarman-group1/uganda-tsetse-LG/data"
results_dir <- "/uufs/chpc.utah.edu/common/home/saarman-group1/uganda-tsetse-LG/results/"

# -----------------------------
# Read files
# -----------------------------

# Fold-level metrics
metrics_all <- read.csv(
  file.path(results_dir, "LOPOCV_summary_v2.csv")
)

# Site cluster assignments
indinfo <- read.delim("../input/Gff_11loci_allsites_indinfo.txt")

site_clusters <- unique(
  indinfo[c("SiteCode", "SiteMajCluster")]
)

names(site_clusters) <- c("site", "Cluster")

# Combine west with south
site_clusters$Cluster[site_clusters$Cluster == "west"] <- "south"


# -----------------------------
# 1. Mean fold-level Spearman
# -----------------------------

metrics_all$Cluster <- site_clusters$Cluster[
  match(metrics_all$site, site_clusters$site)
]

aggregate(
  Spearman ~ Cluster,
  data = metrics_all,
  FUN = mean,
  na.rm = TRUE
)
```

    ##   Cluster  Spearman
    ## 1   north 0.8689949
    ## 2   south 0.8773816

``` r
# -----------------------------
# 2. Pooled out-of-fold R2
# -----------------------------

# Recreate site order used for LOPOCV
V.table <- read.csv("../input/Gff_cse_envCostPaths.csv")
V.table <- V.table[
  V.table$Var1 != "50-KB" & V.table$Var2 != "50-KB",
]

sites <- sort(unique(c(V.table$Var1, V.table$Var2)))

# Read all held-out predictions and add held-out site
eval_all <- do.call(
  rbind,
  lapply(seq_along(sites), function(i) {

    x <- read.csv(
      file.path(
        results_dir,
        sprintf("lopocv_eval_fold_%02d_v2.csv", i)
      )
    )

    x$site <- sites[i]
    x
  })
)

# Add North/South cluster
eval_all$Cluster <- site_clusters$Cluster[
  match(eval_all$site, site_clusters$site)
]

# Function matching your pooled R2 calculation
pooled_R2 <- function(x) {

  keep <- complete.cases(x$CSE, x$pred_cal)

  obs  <- x$CSE[keep]
  pred <- x$pred_cal[keep]

  1 - sum((obs - pred)^2) /
      sum((obs - mean(obs))^2)
}

# Pooled R2 by cluster
sapply(
  split(eval_all, eval_all$Cluster),
  pooled_R2
)
```

    ##     north     south 
    ## 0.7700091 0.8028254

# Mean fold-level Spearman for distance-only model (north vs south)

``` r
metrics_geodist <- read.csv(
  file.path(results_dir, "geodist_LOPOCV_summary.csv")
)

metrics_geodist$Cluster <- site_clusters$Cluster[
  match(metrics_geodist$site, site_clusters$site)
]

aggregate(
  Spearman ~ Cluster,
  data = metrics_geodist,
  FUN = mean,
  na.rm = TRUE
)
```

    ##   Cluster  Spearman
    ## 1   north 0.7412128
    ## 2   south 0.6453284
