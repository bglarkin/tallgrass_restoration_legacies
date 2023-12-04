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
#' # Data
#' For this summary, I'll pull as many objects as possible from existing files to reduce the number of 
#' interspersed code chunks. A few quick new analyses will be necessary, though. Data are loaded here. 
#+ data_sites
sites <-
    read_csv(paste0(getwd(), "/clean_data/sites.csv"), show_col_types = FALSE) %>%
    mutate(
        field_type = factor(
            field_type,
            ordered = TRUE,
            levels = c("corn", "restored", "remnant")),
        yr_since = replace(yr_since, which(field_type == "remnant"), NA),
        yr_since = replace(yr_since, which(field_type == "corn"), NA)) %>%
    select(-lat, -long, -yr_restore, -yr_rank) %>% 
    arrange(field_key)
#' 
#' # Methods
#' The survey followed an unbalanced complete block design. I have called blocks "regions" so far.
#' We collected samples and data from four regions, shown on the map below. 
#+ site_map,echo=FALSE,fig.align='center',out.width="80%"
include_graphics("site_locations_files/figure-gfm/site_map-1.png")
#' 
#' The data are unbalanced because there are more restored fields than corn or remnant. In all but one case,
#' only single corn and remnant fields were available in each region. This means that we only have
#' replication to separate field types when using the entire block design. 
#+ fields_regions_types_table
kable(table(sites$region, sites$field_type),
      format = "pandoc",
      caption = "Count of fields by type and region:\nBM = Blue Mounds, FG = Faville Grove,\nFL = Fermilab, LP = Lake Petite")
