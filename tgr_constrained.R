#' ---
#' title: "Constrained and summary analysis"
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
#' This script presents a constrained analysis on restored sites in the Blue Mounds area. 
#' Lack of replication prevents constraints from being significants across regions, although
#' this was tried previously (not shown). 
#' 
#' # Packages and libraries
packages_needed = c("tidyverse", "vegan", "conflicted", "knitr", "GGally", "ape", "ade4", "GGally")
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
#' # Functions
#' Reference the sidecar file `tgr_constrained_functions.R`.
source("supporting_files/tgr_constrained_functions.R")
#' 
#' # Data
#' ## Site metadata and experimental design
#+ sites
sites <-
    read_csv(paste0(getwd(), "/clean_data/sites.csv"), show_col_types = FALSE) %>%
    mutate(
        field_type = factor(
            field_type,
            ordered = TRUE,
            levels = c("corn", "restored", "remnant"))) %>%
    select(-lat, -long, -yr_restore, -yr_rank) %>% 
    arrange(field_key)
#' 
#' ## Plant data
#' ### Traits
#' Plant abundance data was only available at 16 sites (none at Fermi). These were translated into 
#' percent cover of traits in `plant.R`. Site metadata are joined to allow custom filtering of sites
#' later.
ptr <- read_csv(paste0(getwd(), "/clean_data/plant_trait_abund.csv"), show_col_types = FALSE) %>% 
    left_join(sites %>% select(starts_with("field"), region), by = join_by("field_name")) %>% 
    left_join(read_csv(paste0(getwd(), "/clean_data/spe_plant_abund.csv"), show_col_types = FALSE) %>% 
                  rename(field_name = SITE) %>% select(field_name, BARESOIL, LITTER), by = join_by(field_name)) %>% 
    select(field_name, field_type, region, everything(), -field_key)
#' Plant releve data was available from four Fermi sites, but survey data aren't correctable 
#' between Fermi and other sites. Also, translating counts of species to counts of traits isn't appropriate.
#' Trait data isn't included for the plant presence dataset. 
#' 
#' ### Plant communities
#' Both abundance and presence data are provided. 
#' Plant site-species matrices use field names for rows, site metadata will be joined. 
pspe <- list(
    ab = read_csv(paste0(getwd(), "/clean_data/spe_plant_abund.csv"), show_col_types = FALSE),
    pr = read_csv(paste0(getwd(), "/clean_data/spe_plant_presence.csv"), show_col_types = FALSE)
) %>% map(function(x) x %>% 
              rename(field_name = SITE) %>% 
              left_join(sites %>% select(starts_with("field"), region), by = join_by("field_name")) %>% 
              select(field_name, field_type, region, everything(), -field_key, -BARESOIL, -LITTER))
#' 
#' ## Environmental data
#' Use precipitation as proxy for soil moisture
rain = read_csv(paste0(getwd(), "/clean_data/site_precip_normal.csv"), show_col_types = FALSE)
#' Create subsets of soil environmental data to align with plant abundance or presence sites
soil <- read_csv(paste0(getwd(), "/clean_data/soil.csv"), show_col_types = FALSE) %>% 
    left_join(rain, by = join_by("field_key")) %>% 
    filter(field_key %in% sites$field_key) %>% 
    left_join(sites %>% select(starts_with("field"), region), by = join_by("field_name", "field_key")) %>% 
    select(field_name, field_type, region, everything(), -field_key)
#' 
#' ## Microbial data
#' ### Fungal communities
#' Fungal species matrices must have field names instead of field keys to align all datasets.
#' Create subsets of fungal species matrices to align with plant abundance or presence sites
#' Sites-species tables contain rarefied sequence abundances. This list includes
#' composition summarized in fields.
#' 
#' - Fermi switchgrass prairies are removed because they have no plant data associated. 
#' - CSV files were produced in [process_data.R](process_data.md)
fspe <- list(
    its = read_csv(paste0(getwd(), "/clean_data/spe_ITS_rfy.csv"), show_col_types = FALSE),
    amf = read_csv(paste0(getwd(), "/clean_data/spe_18S_rfy.csv"), show_col_types = FALSE)
) %>% map(function(x) x %>% 
              left_join(sites %>% select(starts_with("field"), region), by = join_by("field_key")) %>% 
              filter(!(field_name %in% c("FLRSP1", "FLRSP2", "FLRSP3"))) %>% 
              select(field_name, field_type, region, everything(), -field_key))
#' 
#' ### Species metadata
#' The OTUs and sequence abundances in these files matches the rarefied data in `spe$` above.
#' CSV files were produced in the [microbial diversity script](microbial_diversity.md).
fspe_meta <- list(
    its = read_csv(paste0(getwd(), "/clean_data/speTaxa_ITS_rfy.csv"), show_col_types = FALSE),
    amf =  read_csv(paste0(getwd(), "/clean_data/speTaxa_18S_rfy.csv"), show_col_types = FALSE)
)
#' 
#' ## Response data
#' ### Fungal biomass
fb <- read_csv(paste0(getwd(), "/clean_data/plfa.csv"), show_col_types = FALSE) %>% 
    left_join(sites %>% select(starts_with("field"), region), by = join_by("field_key", "field_name")) %>% 
    select(field_name, field_type, region, everything(), -field_key, -starts_with("fa_"))
#' 
#' ### Water stable aggregates
# Remove rows from old field sites (26 and 27)
wsa <- read_csv(paste0(getwd(), "/clean_data/wsa.csv"), show_col_types = FALSE)[-c(26:27), ] %>% 
    left_join(sites, by = "field_key")
#' 
#' # Analysis and Results
#' A great number of symmetric and asymmetric comparative and constrained analyses have been attempted
#' with these data. The best and simplest is a db-RDA, executed with the custom function
#' `dbrda_fun()`. The strategy is to forward select variables 
#' of experimental interest. These are the variables that we expect to change 
#' as a result of the restoration, or affect the microbial community directly. 
#' Agricultural nutrients, soil organic matter, and plant community/traits data are in this group. 
#' 
#' ## Plant community axes
#' In this analysis, both plant community and traits data can be used. Let's use `pcoa_fun()` to produce 
#' community axes for the abundance-based plant data. 
#+ plant_pcoa_ab
(pspe_pcoa_ab <- pcoa_fun(pspe$ab))
#' 
#' ## Microbial communities with constraints
#' db-RDA in function `dbrda-fun()`
#' 
#' ### Fungal community (ITS sequence abundance)
#' #### Blue Mounds with plant traits data
#+ dbrda_bm_tr_its,message=FALSE,warning=FALSE
(dbrda_bm_tr_its <- dbrda_fun(s = fspe$its, pspe_pcoa = "none"))[c(3, 4, 5, 2)]
#' The number of sites limits the permutation design, but forb was selected. We know that this variable 
#' is inversely related to C4 grass and years since restoration.  
#+ plot_bm_tr_its,fig.align='center'
plot_dbrda(site_sc = dbrda_bm_tr_its$plot_data$sites,
           site_bp = dbrda_bm_tr_its$plot_data$biplot)
write_csv(as_tibble(dbrda_bm_tr_its$plot_data$sites, rownames = "field_name"), "tgr_constrained_files/bm_tr_its_sitelocs.csv")
write_csv(as_tibble(dbrda_bm_tr_its$plot_data$biplot, rownames = "envvar"), "tgr_constrained_files/bm_tr_its_bp.csv")
#' 
#' #### Blue Mounds with plant community data
#+ dbrda_bm_ab_its,message=FALSE,warning=FALSE
(dbrda_bm_ab_its <- dbrda_fun(s = fspe$its, pspe_pcoa = pspe_pcoa_ab$site_vectors))[c(3, 4, 5, 2)]
#' The number of sites limits the permutation design, but years since was selected. 
#+ plot_bm_ab_its,fig.align='center'
plot_dbrda(site_sc = dbrda_bm_ab_its$plot_data$sites,
           site_bp = dbrda_bm_ab_its$plot_data$biplot)
write_csv(as_tibble(dbrda_bm_ab_its$plot_data$sites, rownames = "field_name"), "tgr_constrained_files/bm_ab_its_sitelocs.csv")
write_csv(as_tibble(dbrda_bm_ab_its$plot_data$biplot, rownames = "envvar"), "tgr_constrained_files/bm_ab_its_bp.csv")
#'
#' ### AMF community (18S sequence abundance)
#' #### Blue Mounds with plant traits data
#+ dbrda_bm_tr_amf,message=FALSE,warning=FALSE
(dbrda_bm_tr_amf <- dbrda_fun(s = fspe$amf, pspe_pcoa = "none"))[c(3, 4, 5, 2)]
#' The number of sites limits the permutation design, but C4 grass was selected. 
#+ plot_bm_tr_amf,fig.align='center'
plot_dbrda(site_sc = dbrda_bm_tr_amf$plot_data$sites,
           site_bp = dbrda_bm_tr_amf$plot_data$biplot)
write_csv(as_tibble(dbrda_bm_tr_amf$plot_data$sites, rownames = "field_name"), "tgr_constrained_files/bm_tr_amf_sitelocs.csv")
write_csv(as_tibble(dbrda_bm_tr_amf$plot_data$biplot, rownames = "envvar"), "tgr_constrained_files/bm_tr_amf_bp.csv")
#+ dbrda_bm_ab_amf,message=FALSE,warning=FALSE
(dbrda_bm_ab_amf <- dbrda_fun(s = fspe$amf, pspe_pcoa = pspe_pcoa_ab$site_vectors))[c(3, 4, 5, 2)]
#' The number of sites limits the permutation design, but years since was selected. 
#+ plot_bm_ab_amf,fig.align='center'
plot_dbrda(site_sc = dbrda_bm_ab_amf$plot_data$sites,
           site_bp = dbrda_bm_ab_amf$plot_data$biplot)
write_csv(as_tibble(dbrda_bm_ab_amf$plot_data$sites, rownames = "field_name"), "tgr_constrained_files/bm_ab_amf_sitelocs.csv")
write_csv(as_tibble(dbrda_bm_ab_amf$plot_data$biplot, rownames = "envvar"), "tgr_constrained_files/bm_ab_amf_bp.csv")
#'
#' Microbial communities align with years since restoration across regions and types (general fungi and amf).
#' 
#' ## Effects of soil and plants
#' We know that restored plant communities differ among fields, and that those differences change 
#' in a systematic way over time. Within our chronosequence sites at Blue Mounds, what is the relative 
#' effect of soil and plant data on soil microbes?
#' 
#' Variation partitioning will be used. Since forward selection failed to find a parsimonious number
#' of important explanatory variables, we'll just use the entire plant and soil datasets here to look for 
#' the relative contribution of these weak effects. 
#' 
#' With so many explanatory axes, the analysis failed with raw data, so explanatory data 
#' will first be transformed into PCoA axes. 
#' 
#' ### Soil fungi variation partitioning
soil_pcoa <- pcoa_fun(soil, method = "euclidean", corr = "lingoes")
vpdat_its <- list(
    Y = fspe$its %>% filter(region == "BM", field_type == "restored") %>% 
        select(-field_type, -region) %>% 
        data.frame(., row.names = 1),
    X1 = pspe_pcoa_ab$site_vectors,
    X2 = soil_pcoa$site_vectors
)
vpdat_its_zcols <- vpdat_its %>% map(\(df) which(apply(df, 2, sum) == 0))
#+ varpart_its
(vp_its <- varpart(vpdat_its$Y %>% select(-vpdat_its_zcols$Y), vpdat_its$X1, vpdat_its$X2))
#+ varpart_its_plot,fig.align='center'
plot(vp_its, digits = 2, bg = c("tan", "palegreen"))
#' Testing significance of A and B fractions with `rda()`
anova(rda(vpdat_its$Y %>% select(-vpdat_its_zcols$Y) ~ Axis.1 + Axis.2 + Condition(as.matrix(vpdat_its$X2)), data = vpdat_its$X1))
#' Plant data not significant
anova(rda(vpdat_its$Y %>% select(-vpdat_its_zcols$Y) ~ Axis.1 + Axis.2 + Condition(as.matrix(vpdat_its$X1)), data = vpdat_its$X2))
#' Soil data not significant
#' 
#' With other axes accounted for, plant or soil axes explain ~23% of the ITS fungal community variation in Blue Mounds. 
#' Neither is significant. 
#' 
#' ### AMF variation partitioning
vpdat_amf <- list(
    Y = fspe$amf %>% filter(region == "BM", field_type == "restored") %>% 
        select(-field_type, -region) %>% 
        data.frame(., row.names = 1),
    X1 = pspe_pcoa_ab$site_vectors,
    X2 = soil_pcoa$site_vectors
)
vpdat_amf_zcols <- vpdat_amf %>% map(\(df) which(apply(df, 2, sum) == 0))
#+ varpart_amf
(vp_amf <- varpart(vpdat_amf$Y %>% select(-vpdat_amf_zcols$Y), vpdat_amf$X1, vpdat_amf$X2))
#+ varpart_amf_plot,fig.align='center'
plot(vp_amf, digits = 2, bg = c("tan", "palegreen"))
#' Testing significance of A and B fractions with `rda()`
anova(rda(vpdat_amf$Y %>% select(-vpdat_amf_zcols$Y) ~ Axis.1 + Axis.2 + Condition(as.matrix(vpdat_amf$X2)), data = vpdat_amf$X1))
#' Plant data not significant
anova(rda(vpdat_amf$Y %>% select(-vpdat_amf_zcols$Y) ~ Axis.1 + Axis.2 + Condition(as.matrix(vpdat_amf$X1)), data = vpdat_amf$X2))
#' Soil data not significant
#' 
#' With soil axes accounted for, plant axes explain 31% of the ITS fungal community variation in Blue Mounds. 
#' The unique explanation made by soil variables is about half as much. Neither is significant. 
#'
#' ## Response correlations
#' Fungal communities varied with years since restoration, C4 grass, and forbs. How do these predictors affect
#' microbial biomass and function?
#+ func_vars_df
func_vars <-
    sites %>%
    left_join(ptr %>% select(field_name, C4_grass, forb), by = join_by(field_name)) %>%
    left_join(fb %>% select(field_name, fungi, amf), by = join_by(field_name)) %>%
    left_join(wsa %>% select(field_name, wsa), by = join_by(field_name)) %>%
    rename(
        C4_grass_pct = C4_grass,
        forb_pct = forb,
        mass_fungi = fungi,
        mass_amf = amf
    ) %>%
    select(-field_key)
#' Plant traits data are only available in Wisconsin; try these first. 
#+ pairs_wi,fig.align='center'
ggpairs(
    data = func_vars %>% filter(field_type == "restored", region != "FL"),
    columns = 4:9
) +
    theme_bw()
#' In all restored sites, response correlations are still possible without plant traits.
#+ pairs_all,fig.align='center'
ggpairs(
    data = func_vars %>% filter(field_type == "restored") %>% 
        select(-C4_grass_pct, -forb_pct),
    columns = 4:7
) +
    theme_bw()
#' Intercorrelations among responses are obvious and not very useful. Years since restoration correlates 
#' positively with fungal biomass, and with sites from Fermi included, this relationship is fairly strong.
#' Fermi sites, with abundant SOM, are probably influencing this result more than years are, though. 
#' 
#' ** What if we correlate the constrained axes with responses?** This might help a little by 
#' reordering sites based on time since restoration *and* community differences, reducing the leverage of
#' sites from Fermi. Let's arrange the first axes from each dbRDA and join them with response data. 
#+ axis_correlations
axis_corr <- 
    list(
        bm_tr_its = dbrda_bm_tr_its$plot_data$sites,
        bm_ab_its = dbrda_bm_ab_its$plot_data$sites,
        bm_tr_amf = dbrda_bm_tr_amf$plot_data$sites,
        bm_ab_amf = dbrda_bm_ab_amf$plot_data$sites
    ) %>% 
    map( ~ .x %>%
             data.frame() %>% 
             rownames_to_column(var = "field_name") %>% 
             left_join(sites %>% select(field_name, yr_since), by = join_by(field_name)) %>% 
             left_join(fb %>% select(field_name, fungi, amf), by = join_by(field_name)) %>% 
             left_join(wsa %>% select(field_name, wsa), by = join_by(field_name)) %>% 
             rename(mass_fungi = fungi, mass_amf = amf) %>% 
             select(-dbRDA2)
    )
#+ axis_correlation_plots_list
axis_corr_plot <-
    axis_corr %>%
    map( ~ .x %>%
            ggpairs(columns = 2:6) +
            theme_bw())
#' 
#' Pairs panels were investigated (not shown), and nothing new shows up.








