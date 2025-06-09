
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
install.packages("ade4")             # perform mantel.rtest 



