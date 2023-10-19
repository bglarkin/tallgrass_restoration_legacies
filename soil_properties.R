#' ---
#' title: "Soil Abiotic Properties"
#' author: "Beau Larkin\n"
#' date: "Last updated: `r format(Sys.time(), '%d %B, %Y')`"
#' output:
#'   github_document:
#'     toc: true
#'     toc_depth: 3
#'     df_print: paged
#'     fig_width: 8
#'     fig_height: 7
#' ---
#'
#' # Description
#' Soil nutrients were analyzed by Ward Lab, analysis details available in local files.
#' Soil organic matter is in percent determined by the loss-on-ignition method.
#' pH is in a log scale as is typical, and all the other minerals are in ppm. This may need
#' to be converted to mg/kg or a similar ratio. 
#' 
#' This script provides a quick overview of the soil abiotic property data and produces
#' products (e.g., ordination axis values) for use in downstream analysis. 
#' 
#' # Packages and libraries
packages_needed = c("tidyverse", "vegan", "pca3d", "GGally")
packages_installed = packages_needed %in% rownames(installed.packages())
#+ packages,message=FALSE
if (any(!packages_installed)) {
    install.packages(packages_needed[!packages_installed])
}
#+ libraries,message=FALSE
for (i in 1:length(packages_needed)) {
    library(packages_needed[i], character.only = T)
}
#' 
#' # Functions
#' Cleanplot PCA produces informative visualizations of PCA ordinations 
#' [(Borcard et al. 2018)](http://link.springer.com/10.1007/978-3-319-71404-2)
#+ cleanplot_pca
source(paste0(getwd(), "/supporting_files/cleanplot_pca.txt"))
#' 
#' # Data