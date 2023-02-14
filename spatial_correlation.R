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
#' Spearman's Rho is used to test significance of correlations. The Sturges equation is
#' used to determine the number of distance classes.
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
#' 
#' ITS-based species data
#+ spe_its_otu
spe <-
    data.frame(read_csv(
        paste0(getwd(), "/clean_data/spe_ITS_otu_siteSpeMatrix_avg.csv"),
        show_col_types = FALSE
    ),
    row.names = 1)
#' Site metadata and locations dataframe
#+ sites_table
sites <-
    read_csv(paste0(getwd(), "/clean_data/site.csv"), show_col_types = FALSE) %>%
    mutate(field_type = factor(
        site_type,
        ordered = TRUE,
        levels = c("corn", "restored", "remnant")
    )) %>%
    filter(field_type != "oldfield")
locs <-
    data.frame(sites %>% select(site_key, long, lat),
               row.names = 1)
#' Soil chemical properties, filtered to the Blue Mounds area and restored fields only
#+ soil_chem
soil <-
    data.frame(
        read_csv(paste0(getwd(), "/clean_data/soil.csv"), show_col_types = FALSE) %>%
            left_join(sites %>% select(site_key, region, field_type), by = join_by(site_key)) %>%
            filter(region == "BM", field_type == "restored"),
        row.names = 1
    ) %>%
    select(-site_name,-region,-field_type)
bm_resto <- rownames(soil)
#' Filter species and location data to the subset in Blue Mounds restored fields only. Also
#' remove species that are zero abundance in these samples. 
#+ spe_bm
spe_bm <- spe[bm_resto, -c(which(apply(spe[bm_resto, ], 2, sum) == 0))]
locs_bm <- locs[bm_resto, ]
#+ site_locs_bm
#' 
#' Create distance matrices:
#'
#' - Percent difference for the microbial species abundance data
#' - Shortest ellipsoid distance for the location coordinates using
#' [`distGeo()`](https://www.rdocumentation.org/packages/geosphere/versions/1.5-18/topics/distGeo)
#' - Euclidean distance on column standardized soil chemical variable data (standardized 
#' to mean=0 and unit variance)
#'
#+ dist_spe
dist_spe <- vegdist(spe, method = "bray")
dist_spe_bm  <- vegdist(spe_bm, method = "bray")
#+ dist_locs
dist_geo <- as.dist(distm(locs, fun = distGeo))
dist_geo_bm <- as.dist(distm(locs_bm, fun = distGeo))
#+ dist_soil
dist_soil <- dist(decostand(soil, "standardize"))
#'
#' # Results
#' ## Test 1
#' Test 1 is with ITS and OTU parameters on the species data from all regions and fields.  
#+ mantel_1
man_1 <-
    mantel(
        dist_spe,
        dist_geo,
        method = "spearman",
        permutations = 9999,
        na.rm = TRUE
    )
#+ mantel_1_summary,echo=FALSE
man_1
#+ mantel_1_cor
man_1_cor <-
    mantel.correlog(
        dist_spe,
        dist_geo,
        n.class = 0,
        cutoff = FALSE,
        r.type = "spearman",
        nperm = 9999,
        mult = "holm"
    )
print(man_1_cor)
#+ fig_man_1,fig.align='center'
plot(man_1_cor)
#' 
#' Global test fails to reject the null of autocorrelation, but barely. The correlation 
#' was weak (R=0.11) but marginally significant (p<0.05). 
#' At smaller scales (< 5 km), correlation is apparent (not shown). This is
#' because only replicate plots (e.g., switchgrass fields at Fermi or
#' multiple restored fields at Lake Petite) fall within that distance bin,
#' so similarity is expected.
#' 
#' In the correlogram, no significant correlations are found. 
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
#' 
#' ## Test 2
#' Actually two tests: location compared with species or soil chemical data to 
#' test for autocorrelation within the restored Blue Mounds fields.
#' 
#' ### Autocorrelation with species data
#+ mantel_2,message=FALSE
man_2 <-
    mantel(
        dist_spe_bm,
        dist_geo_bm,
        method = "spearman",
        na.rm = TRUE
    )
#+ mantel_2_summary,echo=FALSE
man_2
#' Only 5039 permutations are possible due to the small number of sites (n=7). The global 
#' test shows a non-significant negative correlation. No correlogram test is warranted. 
#' 
#' ### Autocorrelation with soil chemical data
#+ mantel_3,message=FALSE
man_3 <-
    mantel(
        dist_soil,
        dist_geo_bm,
        method = "spearman",
        na.rm = TRUE
    )
#+ mantel_3_summary
man_3
#' Only 5039 permutations are possible due to the small number of sites (n=7). The global 
#' test shows a non-significant negative correlation. No correlogram test is warranted. 
