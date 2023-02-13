#' ---
#' title: "Testing for autocorrelation among sites"
#' author: "Beau Larkin\n"
#' date: "Last updated: `r format(Sys.time(), '%d %B, %Y')`"
#' output:
#'   github_document:
#'     toc: true
#'     toc_depth: 3
#'     df_print: paged
#'     fig_width: 5.5
#'     fig_height: 5
#' ---
#'
#' # Description
#' The sites sampled in this project lie in different soil types,
#' some at great distances from each other and some adjacent.
#' Before attempting to characterize the sites as a chronosequence,
#' or comparing results in corn, restored, and remnant sites,
#' we must learn how important spatial autocorrelation is in these data.
#'
#' As an initial test, the microbial abundances determined from the ITS gene
#' and clustered into 97% similar OTUs will be used. To make the mantel
#' test easier, the species abundances will first be averaged among
#' repeated measures at sites.
#'
#' Depending on the results of this test, other tests may be performed. It must
#' be checked whether Mantel tests are still appropriate for this kind of
#' inference.
#'
#' # Packages and libraries
packages_needed = c("tidyverse",
                    "vegan",
                    "geosphere")
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
#' # Data
#' Species abundance data and geographic locations must be loaded. ITS sequence data
#' clustered into 97% OTUs are requested first. Soil properties may also be used.
spe_its_otu <-
    data.frame(read_csv(
        paste0(getwd(), "/clean_data/spe_ITS_otu_siteSpeMatrix_avg.csv"),
        show_col_types = FALSE
    ),
    row.names = 1)
sites <-
    read_csv(paste0(getwd(), "/clean_data/site.csv"), show_col_types = FALSE) %>%
    mutate(field_type = factor(
        site_type,
        ordered = TRUE,
        levels = c("corn", "restored", "remnant")
    )) %>%
    filter(field_type != "oldfield")
site_locs <-
    data.frame(sites %>% select(site_key, long, lat),
               row.names = 1)
soil <-
    data.frame(
        read_csv(paste0(getwd(), "/clean_data/soil.csv"), show_col_types = FALSE) %>%
            left_join(sites %>% select(site_key, region, field_type), by = join_by(site_key)) %>%
            filter(region == "BM"),
        row.names = 1
    ) %>%
    select(-site_name,-region,-field_type)
#' Create distance matrices:
#'
#' - Percent difference for the microbial species abundance data
#' - Shortest ellipsoid distance for the location coordinates using
#' [`distGeo()`](https://www.rdocumentation.org/packages/geosphere/versions/1.5-18/topics/distGeo)
#' - Euclidean distance on column standardized soil chemical variable data
#'
dist_spe_its_otu <- vegdist(spe_its_otu, method = "bray")
dist_geo <- as.dist(distm(site_locs, fun = distGeo))
dist_soil <- dist(decostand(soil, "standardize"))
#'
#' # Results
#' Test 1 is with ITS and OTU parameters on the species data. Eight classes specified 
#' because more than this results in NA values due to single sites in classes. 
man_spe_geo_1 <-
    mantel(
        dist_spe_its_otu,
        dist_geo,
        method = "spearman",
        permutations = 9999,
        na.rm = TRUE
    )
man_spe_geo_1
man_spe_correlog <-
    mantel.correlog(
        dist_spe_its_otu,
        dist_geo,
        n.class = 8,
        cutoff = FALSE,
        r.type = "spearman",
        nperm = 9999,
        mult = "holm"
    )

#' Result of test 1
#' Global test fails to reject the null of autocorrelation, but barely
#' Correlation test shows that two distance classes are weakly correlated.
#' At smaller scales (< 5 km), correlation is apparent (not shown). This is
#' because only replicate plots (e.g., switchgrass fields at Fermi or
#' multiple restored fields at Lake Petite) fall within that distance bin,
#' so similarity is expected.
#' 
#' I think we need to be cautious with an expectation that these sites can
#' be compared generally, but maybe we can do it. It may be easier, more
#' conservative, and appropriate to consider just the Blue Mounds sites
#' for comparisons over time. Even Blue Mounds will be difficult to justify as
#' a true chronosequence, but as a case study it could add a worthwhile dimension to the
#' project without reaching too far.
#'
#' It may also be important to test for autocorrelation using soil data. Let's do that in the
#' next section, with only the Blue Mounds sites.
