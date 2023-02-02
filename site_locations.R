#' ---
#' title: "Site locations, metadata, and climate"
#' author: "Beau Larkin\n"
#' date: "Last updated: `r format(Sys.time(), '%d %B, %Y')`"
#' output:
#'   github_document:
#'     toc: true
#'     toc_depth: 3
#'     df_print: paged
#'     fig_width: 7
#'     fig_height: 6
#' ---
#'
#' # Description
#' Site locations, maps, and metadata are included here
#' Other raster data, like climate data, are also be processed here
#' Precipitation data downloaded from [PRISM](https://prism.oregonstate.edu/) 2022-01-06, 
#' Following this [tutorial](https://rpubs.com/collnell/get_prism)
#' - 30-year normal from 1991-2020 was attempted
#' - Also monthly data from 2007-2016 were downloaded and converted into BIOCLIM variables


# Resources ------------------------
# ——————————————————————————————————
packages_needed = c(
    "tidyverse",
    "colorspace",
    "ggmap",
    "lubridate",
    "devtools",
    "raster",
    "sp",
    "prism",
    "conflicted",
    "ggbeeswarm",
    "knitr",
    "dismo",
    "vegan"
)
packages_installed = packages_needed %in% rownames(installed.packages())

if (any(! packages_installed))
    install.packages(packages_needed[! packages_installed])
for (i in 1:length(packages_needed)) {
    library(packages_needed[i], character.only = T)
} 

conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("extract", "raster")

# Cleanplot PCA (Borcard et al.)
source("/Volumes/blarkin@mpgcloud/Private/blarkin/ΩMiscellaneous/R_global/cleanplot_pca.txt")


# Data -----------------------------
# ——————————————————————————————————
sites <- read_csv(paste0(getwd(), "/clean_data/site.csv")) %>% 
    filter(site_type != "oldfield") %>% 
    glimpse()

# Climate data
# ——————————————————————————————————
# All data now stored in the external hard drive; summary tables written to the working directory
# ——————————————————————————————————
# Directory and command to download the normals
# prism_set_dl_dir(paste0(getwd(), "/prism_rasters/normals/ppt"))
# get_prism_normals(type = 'ppt', resolution = '4km', annual = TRUE, keepZip = FALSE) # data already available in working directory
# prism_archive_ls()
# RS <- pd_stack(prism_archive_ls())
proj4string(RS) <- CRS("+proj=longlat +ellps=WGS84 +no_defs")

sites_spdf <-
    SpatialPointsDataFrame(
        coords = sites[, c("long", "lat")],
        data = sites,
        proj4string = CRS("+proj=longlat +ellps=WGS84 +no_defs")
    )
sites_ppt <- raster::extract(RS, sites_spdf, fun = mean, na.rm = TRUE, sp = TRUE)@data %>% 
    rename(ppt_mm = PRISM_ppt_30yr_normal_4kmM3_annual_bil)
write_csv(sites_ppt %>% select(site_key, ppt_mm), paste0(getwd(), "/clean_data/site_precip_normal.csv"))
# ——————————————————————————————————

# All data now stored in the external hard drive; summary tables written to the working directory
# ——————————————————————————————————
# Directory and command to download the monthly data
# run once each for ppt, tmin, tmax
# prism_set_dl_dir(paste0(getwd(), "/prism_rasters/monthly/åtmax"))
# get_prism_monthlys(type = "tmax", years = 2007:2016, mon = 1:12, keepZip = FALSE)
# prism_archive_ls()
# RS <- pd_stack(prism_archive_ls())
proj4string(RS) <- CRS("+proj=longlat +ellps=WGS84 +no_defs")
sites_spdf <-
    SpatialPointsDataFrame(
        coords = sites[, c("long", "lat")],
        data = sites,
        proj4string = CRS("+proj=longlat +ellps=WGS84 +no_defs")
    )
sites_tmax_month <- raster::extract(RS, sites_spdf, fun = mean, na.rm = TRUE, sp = TRUE)@data


sites_ppt_month[, c(2,10,11)]
sites_tmin_month[, c(2,10,11)]
sites_tmax_month[, c(2,10,11)]

wrangle_clim <- function(data, value) {
    data %>% 
        select(site_name, starts_with("PRISM")) %>% 
        pivot_longer(starts_with("PRISM"), names_to = "var", values_to = value) %>% 
        separate(var, c(NA, NA, NA, NA, "date_str", NA), sep = "_")
}

ppt_month <- wrangle_clim(sites_ppt_month, "ppt_mm")
tmin_month <- wrangle_clim(sites_tmin_month, "tmin_C")
tmax_month <- wrangle_clim(sites_tmax_month, "tmax_C")

tgr_clim <- 
    ppt_month %>% 
    left_join(tmin_month) %>% 
    left_join(tmax_month) %>% 
    mutate(date = ym(date_str), year = year(date), month = month(date)) %>% 
    select(site_name, year, month, ppt_mm, tmin_C, tmax_C) %>% 
    glimpse()

write_csv(tgr_clim, paste0(getwd(), "/clean_data/clim.csv"))

yrs <- unique(tgr_clim$year)
sites <- unique(tgr_clim$site_name)
tgr_bioclim <- vector("list", length(yrs))
names(tgr_bioclim) <- yrs
for(i in 1:length(yrs)) {
    data <- 
        tgr_clim %>% 
        filter(year == yrs[i])
    site_list <- split(data, data$site_name)
    tgr_bioclim[[i]] <- 
        bind_rows(lapply(site_list, function(x){data.frame(biovars(x$ppt_mm, x$tmin_C, x$tmax_C))}), .id = "site")
}

bind_rows(tgr_bioclim, .id = "year") %>% 
    pivot_longer(starts_with("bio"), names_to = "biovar", values_to = "value") %>% 
    group_by(site, biovar) %>% 
    summarize(mean = mean(value)) %>% 
    pivot_wider(names_from = "biovar", values_from = mean) %>% 
    write_csv(paste0(getwd(), "/clean_data/bioclim.csv"))


# ——————————————————————————————————



# Results --------------------------
# ——————————————————————————————————
# Sites in categories
kable(table(sites$region, sites$site_type), format = "pandoc")

# Site maps and extent
map <- ggmap(
    get_stamenmap(
        bbox = c(
            left = -90.3,
            bottom = 41.5,
            right = -87.4,
            top = 43.4
        ),
        zoom = 9,
        maptype = c("toner-lite"),
        color = c("color")
    )
)
map +
    geom_polygon(
        data = sites %>% group_by(region) %>% slice(chull(long, lat)),
        aes(x = long, y = lat, group = region),
        fill = "red", alpha = 0.5, color = "red"
    ) +
    geom_label(
        data = sites %>% group_by(region) %>% summarize(long_cen = mean(long), lat_cen = mean(lat), .groups = "drop"),
        aes(x = long_cen, y = lat_cen, label = region),
        color = "red", size = 6
    ) +
    theme_void()

# Precipitation normals
ggplot(sites_ppt, aes(x = region, y = ppt_mm)) +
    geom_beeswarm(aes(fill = factor(site_type, ordered = TRUE, levels = c("corn", "restored", "remnant"))), dodge.width = 0.2, shape = 21, size = 4) +
    labs(x = "", y = "Precipitation (mm)") +
    scale_fill_discrete_qualitative(name = "site_type", palette = "Harmonic") +
    theme_bw()

# Bioclimate variables from 2007-2016
bioclim <- read_csv(paste0(getwd(), "/clean_data/bioclim.csv")) %>% glimpse
bioclim_pca <- rda(data.frame(bioclim, row.names = 1), scale = TRUE)
summary(bioclim_pca, display = NULL)
screeplot(bioclim_pca, bstick = TRUE)
cleanplot.pca(bioclim_pca)
# Component loadings (correlations)
data.frame(scores(bioclim_pca, choices = c(1,2), display = "species", scaling = 0)) %>% arrange(PC1)
# Site scores
data.frame(scores(bioclim_pca, choices = c(1,2), display = "sites", scaling = 1)) %>% 
    rownames_to_column() %>% 
    rename(site_name = rowname) %>% 
    write_csv(paste0(getwd(), "/clean_data/clim_axes.csv"))
