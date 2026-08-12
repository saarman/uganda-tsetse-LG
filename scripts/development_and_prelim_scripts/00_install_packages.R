
# Install required packages for “Uganda-tsetse-LG”

# If you do not yet have devtools available,
# install it so you can install packages from GitHub if needed:
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}

# (Only needed if a package is not on CRAN and must be installed from GitHub.)
#    Example: install the development version of adegenet from its GitHub repo:
#    devtools::install_github("thibautjombart/adegenet")

# --------------------------------------------------------------------
# Install from CRAN if possible
# --------------------------------------------------------------------

# 01_PCA_diversity_CSE.Rmd
# --------------------------------------------------------------------
install.packages("adegenet")         # read Genepop, handle “genind” objects, perform PCA, genind→genpop
install.packages("geosphere")        # compute geographic distances (distm, distHaversine) for site filtering
install.packages("igraph")           # build adjacency graphs for grouping nearby sites (<2 km clusters)
install.packages("scatterplot3d")    # create 3D PCA scatterplots
install.packages("ggplot2")          # create density plots and CSE‐colored maps
install.packages("dplyr")            # data‐frame manipulation (filter, mutate, merge, etc.)
install.packages("hierfstat")        # convert “genind” to hierfstat format and compute diversity stats (basic.stats)
install.packages("pegas")            # calculate allelic richness (rarefied) per population
install.packages("PopGenReport")     # generate standardized population‐genetics summary reports
install.packages("poppr")            # compute per‐population diversity metrics (poppr())
install.packages("sf")               # build and plot spatial LINESTRINGs for CSE paths
install.packages("rnaturalearth")    # download country boundary shapefiles (Uganda outline)
install.packages("rnaturalearthdata")# underlying spatial data required by rnaturalearth
install.packages("cowplot")          # combine multiple ggplot panels into one figure
install.packages("viridis")          # provide “magma” and other Viridis palettes for CSE mapping
install.packages("ade4")             # perform mantel.rtest of isolation by distance (IBD) 
install.packages("vegan")            # perform Mantel Correlogram for scale of IBD

# 02_prepare_rasters.R
# --------------------------------------------------------------------
install.packages("raster")        # for reading, projecting, and resampling raster layers
install.packages("sf")            # for reading and manipulating vector data (lake and river edges)
install.packages("terra")         # fast raster I/O and reprojection (project, ext, resample)
install.packages("KernSmooth")    # for 2D kernel density estimation (bkde2D)
install.packages("maps")          # for base map data and coordinate utilities
install.packages("units")         # for handling physical units (set_units)
install.packages("future.apply")  # for parallel apply functions (future_lapply)

# 03_least-cost-lakes.Rmd
# --------------------------------------------------------------------
install.packages("gdistance")       # build transition matrices and shortest paths (least-cost paths)
install.packages("sp")              # create and manipulate SpatialLines objects
install.packages("foreach")         # enable foreach-style iteration (parallel-safe)
install.packages("doParallel")      # register parallel backend for foreach
# (already included in 01/02): raster, sf, ggplot2, rnaturalearth, viridis, units, future.apply

# 04_extract_envvar.Rmd
# --------------------------------------------------------------------
# (already installed above): raster, gdistance, sp, sf, foreach, doParallel, future.apply
# no new packages beyond those already installed in steps 01–03

# 05_RFmodel_full.Rmd
# --------------------------------------------------------------------
# (already installed above): raster, sf, viridis, doParallel
install.packages("randomForest")    # build Random Forest and Multivariate RF models



