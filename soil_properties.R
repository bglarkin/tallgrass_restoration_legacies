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
#' Soil nutrients were analyzed by [Ward Laboratories, Inc.](https://www.wardlab.com/services/soil-health-analysis/), 
#' analysis methods available in local files or at the link included here.
#' Soil organic matter is in percent determined by the loss-on-ignition method.
#' Soil pH is in a log scale as is typical, and all the other minerals are in parts per million. 
#' This may need to be converted to $mg*kg^{-1}$ or a similar ratio. 
#' 
#' This script provides a quick overview of the soil abiotic property data and produces
#' products (e.g., ordination axis values) for use in downstream analysis. 
#' 
#' # Packages and libraries
packages_needed = c("tidyverse", "vegan", "GGally")
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
sites <-
    read_csv(paste0(getwd(), "/clean_data/sites.csv"), show_col_types = FALSE) %>%
    mutate(
        field_type = factor(
            field_type,
            ordered = TRUE,
            levels = c("corn", "restored", "remnant")),
        yr_since = replace(yr_since, which(field_type == "remnant"), "+"),
        yr_since = replace(yr_since, which(field_type == "corn"), "-")) %>%
    select(-lat, -long, -yr_restore, -yr_rank) %>% 
    arrange(field_key)
# Remove soil data from sites in old fields (rows 26, 27)
soil <- read_csv(paste0(getwd(), "/clean_data/soil.csv"), show_col_types = FALSE)[-c(26:27), ]
#' 
#' # Results
#' ## PCA ordination
soil_z <- decostand(data.frame(soil[, -1], row.names = 1), "standardize")
soil_pca <- rda(soil_z)
soil_pca %>% summary(., display = NULL)
#+ soil_screeplot_fig,fig.align='center'
screeplot(soil_pca, bstick = TRUE)
#' Eigenvalues from the first five axes exceed the broken stick model.
#+ soil_cleanplot_fig,fig.align='center'
cleanplot.pca(soil_pca)
#' Abiotic properties that exceed the unit circle (Scaling 1 plot) exert more influence on the 
#' ordination of sites. These are $P,~NO_3,~SO_4,~Ca,~Mg,$ and $OM$. $Mn$ and $pH$ are close
#' enough to warrant further investigation.
#' 
#' Sites sort in somewhat predictable ways (Scaling 2 plot). Cornfields are associated with phosphorus,
#' nitrate, and sulfate. Many, but not all, remnants are associates with soil organic matter. 
#' Blue Mounds fields are associated with manganese, but since manganese isn't a very strong
#' element in this ordination, these fields may also be very low in soild organic matter, magnesium, 
#' or calcium.