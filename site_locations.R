#' ---
#' title: "Site locations, metadata, and climate"
#' author: "Beau Larkin\n"
#' date: "Last updated: `r format(Sys.time(), '%d %B, %Y')`"
#' output:
#'   github_document:
#'     toc: true
#'     toc_depth: 3
#'     df_print: paged
#' ---
#'
#' # Description
#' Site locations, maps, and metadata are produced here.
#' Climate data (precipitation) are downloaded and processed here.
#'
#' - Precipitation data was downloaded from [PRISM](https://prism.oregonstate.edu/) on 2022-01-06
#' following this [tutorial](https://rpubs.com/collnell/get_prism).
#' - Raw (downloaded) data is stored locally; site level data was processed into .csv files
#' and included in this repository.
#' - Precipitation data are 30-year normals from 1991-2020
#'
#' Note that in this script, messages and verbose outputs are often suppressed for brevity.
#'
#' # Package and library installation
packages_needed = c(
    "tidyverse",
    "colorspace",
    "ggmap",
    "raster",
    "sp",
    "prism",
    "conflicted",
    "ggbeeswarm",
    "knitr",
    "ggrepel"
)
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
{
    conflict_prefer("filter", "dplyr")
    conflict_prefer("select", "dplyr")
    conflict_prefer("extract", "raster")
}
#+ map_settings
source(paste0(getwd(), "/supporting_files/map_settings.R"))
#' 
#' # Data and ETL
#' ## Sites
sites <-
    read_csv(paste0(getwd(), "/clean_data/sites.csv"), show_col_types = FALSE) %>%
    glimpse()
#'
#' ## Climate data
#' Climate data was accessed and downloaded on 2022-01-06. The following script does not need to be run
#' again and is commented out to prevent errors and overwrites. Raw downloaded data and summary files were written to
#' the working directory and included in this repository.
#'
#' ### Normals
#' Download archive rasters
# prism_set_dl_dir(paste0(getwd(), "/prism/"))
# get_prism_normals(type = "ppt", resolution = "4km", annual = TRUE, keepZip = FALSE)
# prism_archive_ls()
# RS <- pd_stack(prism_archive_ls())
# proj4string(RS) <- CRS("+proj=longlat +ellps=WGS84 +no_defs")
# sites_spdf <-
#     SpatialPointsDataFrame(
#         coords = sites[, c("long", "lat")],
#         data = sites,
#         proj4string = CRS("+proj=longlat +ellps=WGS84 +no_defs")
#     )
#' Extract locations from archive rasters
# sites_ppt <- raster::extract(RS, sites_spdf, fun = mean, na.rm = TRUE, sp = TRUE)@data %>%
#     rename(ppt_mm = PRISM_ppt_30yr_normal_4kmM4_annual_bil)
#' Create data table `site_precip_normal.csv`
# write_csv(sites_ppt %>% select(field_key, ppt_mm), paste0(getwd(), "/clean_data/site_precip_normal.csv"))
#'
#' # Results
#' ## Site types
#' How many sites are in each field type?
kable(table(sites$region, sites$field_type),
      format = "pandoc",
      caption = "Count of fields by type and region:\nBM = Blue Mounds, FG = Faville Grove,\nFL = Fermilab, LP = Lake Petite")
#'
#' ## Regional map
#' https://docs.stadiamaps.com/themes/
#+ map_object,message=FALSE
map <- ggmap(get_stadiamap(
    bbox = c(
        left = -90.3,
        bottom = 41.5,
        right = -87.4,
        top = 43.4
    ),
    zoom = 9,
    maptype = c("stamen_terrain"),
    color = c("color")
))
#+ site_map,message=FALSE,fig.height=7,fig.width=7,fig.align='center'
map +
    geom_label(
        data = sites %>% group_by(region) %>% summarize(
            long_cen = mean(long),
            lat_cen = mean(lat),
            .groups = "drop"
        ),
        aes(x = long_cen, y = lat_cen, label = region),
        fill = "gray90",
        size = 8
    ) +
    theme_void()
#'
#' The map labels show centroids for each region: BM = Blue Mounds, FG = Faville Grove, FL = Fermilab, LP = Lake Petite.
#' https://docs.stadiamaps.com/themes/.
#' 
#' ## Blue Mounds map
sites %>% 
    filter(region == "BM") %>% 
    select(lat, long) %>% 
    map(. %>% range(.))
#+ map_object_bm,message=FALSE
map_bm <- ggmap(get_stadiamap(
    bbox = c(
        left = -89.98,
        bottom = 42.76,
        right = -89.50,
        top = 43.12
    ),
    zoom = 12,
    maptype = c("stamen_terrain"),
    color = c("color")
))
#+ site_map_bm,message=FALSE,fig.width=7,fig.height=7,fig.align='center',warning=FALSE
map_bm +
    geom_label_repel(data = sites %>% filter(region == "BM"),
        aes(x = long, y = lat, label = field_name),
        fill = "gray90",
        size = 7
    ) +
    geom_point(data = sites %>% filter(region == "BM"),
               aes(x = long, y = lat), fill = "red", size = 3, shape = 21) +
    theme_void()
#' 
#' ## Regional areas
#' The area of convex hulls around field in each region is shown in the table below:
kable(data.frame(
    region = c("BM", "FG", "FL", "LP"),
    area_ha = c("39,500", "12", "330", "14")
),
format = "pandoc")
#'
#' ## Precipitation normals
#' 30-year normals in 4km grid cells around each field are displayed in the table below. Some
#' fields occupy the same 4km grid cell. 
sites_ppt <-
    read_csv(paste0(getwd(), "/clean_data/site_precip_normal.csv"),
             show_col_types = FALSE) %>%
    left_join(sites, by = "field_key")
#+ precipitation_normals_table
kable(
    sites_ppt %>% 
        select(region, field_name, field_type, ppt_mm) %>% 
        mutate(ppt_mm = round(ppt_mm, 0)) %>% 
        arrange(region, -ppt_mm),
    format = "pandoc"
)
#' See 30-year precipitation normals plotted by region and field type in the figure below. Regions
#' and fields within regions do differ in precipitation, but the differences are a tiny proportion
#' of the total precipitation in any field. 
#+ normals_plot,fig.align='center',fig.width=7,fig.height=5
ggplot(sites_ppt, aes(x = region, y = ppt_mm)) +
    geom_beeswarm(
        aes(fill = factor(
            field_type,
            ordered = TRUE,
            levels = c("corn", "restored", "remnant")
        )),
        dodge.width = 0.2,
        shape = 21,
        size = 4
    ) +
    labs(x = "", y = "Precipitation (mm)") +
    scale_fill_discrete_qualitative(name = "field_type", palette = "Harmonic") +
    theme_bw()

