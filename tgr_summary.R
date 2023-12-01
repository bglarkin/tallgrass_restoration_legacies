#' ---
#' title: "Tallgrass Prairie Restoration Legacies, Summary"
#' author: "Beau Larkin\n"
#' date: "Last updated: `r format(Sys.time(), '%d %B, %Y')`"
#' output:
#'  pdf_document:
#'      df_print: tibble
#'      fig_caption: yes
#'      highlight: tango
#'      toc: true
#'      toc_depth: 2
#' ---
#' 
#' # Packages and Libraries
packages_needed = c("tidyverse",
                    "png",
                    "knitr",
                    "conflicted")
packages_installed = packages_needed %in% rownames(installed.packages())
#+ packages,message=FALSE
if (any(!packages_installed)) {
    install.packages(packages_needed[!packages_installed])
}
#+ libraries,message=FALSE
for (i in 1:length(packages_needed)) {
    library(packages_needed[i], character.only = T)
}
#+ conflicts,message=FALSE
conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")
#' 
#' # Description
#' 
#+ plot
plot(1:10, 1:10)


#+ site_map,echo=FALSE,fig.align='center',out.width="40%"
include_graphics("site_locations_files/figure-gfm/site_map-1.png")
