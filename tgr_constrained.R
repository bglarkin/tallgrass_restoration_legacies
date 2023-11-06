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
#' This script will combine and summarize the many single datasets presented so far. 
#' Several multivariate analyses will be considered. Symmetrical analyses will attempt 
#' to discern the relative contributions of plant and soil data on site ordination.
#' Asymmetric constrained analyses will attempt to determine the most significant 
#' contributors of fungal community difference. 
#' 
#' Other data, like microbial biomass and soil water stable aggregates may also be 
#' considered.
#' 
#' All explanatory data are provided as field-based averages. Using them to constrain
#' microbial community data will require alignment with field-based summary data and cannot
#' use subsamples of microbial data.
#' 
#' Note that the plant data are not available for all sites, leading to some
#' complication in the analysis. See further notes below. 
#' 
#' # Packages and libraries
packages_needed = c("tidyverse", "vegan", "conflicted", "knitr", "GGally", "ape", "ade4")
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
#' ## PCoA axes
#' Produce community axes for use in constrained analyses later. 
#+ pcoa_function
pcoa_fun <- function(s, ft=c("restored"), rg, method="bray", binary=FALSE, corr="none") {
    d <- vegdist(
        data.frame(
            s %>% 
                filter(field_type %in% ft, region %in% rg) %>% 
                select(-field_type, -region),
            row.names = 1),
        method = method,
        binary = binary)
    p <- pcoa(d, correction = corr)
    p_vals <- data.frame(p$values) %>% 
        rownames_to_column(var = "Dim") %>% 
        mutate(Dim = as.integer(Dim))
    p_vec <- data.frame(p$vectors)
    if(corr == "none" | ncol(p_vals) == 6) {
        p_ncomp <- with(p_vals, which(Relative_eig < Broken_stick)[1]-1)
    } else {
        p_ncomp <- with(p_vals, which(Rel_corr_eig < Broken_stick)[1]-1)
    }
    ncomp <- if(p_ncomp <= 2) {2} else {p_ncomp}
    # Ordination plot
    scores <- p_vec[, 1:2]
    # Output data
    output <- list(important_components           = p_ncomp,
                   correction_note                = p$note,
                   values                         = p_vals[1:(ncomp+1), ],
                   site_vectors                   = scores)
    return(output)
}
#' 
#' ## Constrained analyses
#' This function performs a db-RDA on microbial data with soil/site covariables and 
#' forward selects on plant communities, agricultural nutrients, and soil carbon.
dbrda_fun <- function(s, pspe_pcoa="none", ft, rg) {
    fspe_bray <- vegdist(
        data.frame(
            s %>% 
                filter(field_type %in% ft, region %in% rg) %>% 
                select(-field_type, -region),
            row.names = 1),
        method = "bray")
    if(is.data.frame(pspe_pcoa) == TRUE) {
        pspe_ax <- pspe_pcoa %>% 
            rownames_to_column("field_name") %>% 
            rename(plant1 = Axis.1, plant2 = Axis.2)
    }
    ptr_norm <- decostand(
        data.frame(
            ptr %>% 
                filter(field_type %in% ft, region %in% rg) %>% 
                select(-field_type, -region),
            row.names = 1),
        "normalize") %>% 
        rownames_to_column("field_name")
    soil_cov_z <- decostand(
        data.frame(
            soil %>% 
                filter(field_type %in% ft, region %in% rg) %>% 
                select(-field_type, -region, -OM, -NO3, -P, -K),
            row.names = 1),
        "standardize")
    soil_cov_sc <- data.frame(scores(rda(soil_cov_z), display = "sites")) %>% 
        rename(soil1 = PC1, soil2 = PC2) %>% 
        rownames_to_column("field_name")
    soil_expl_z <- decostand(
        data.frame(
            soil %>% 
                filter(field_type %in% ft, region %in% rg) %>% 
                select(field_name, OM, NO3, P, K, -field_type, -region),
            row.names = 1),
        "standardize") %>% 
        rownames_to_column("field_name")
    yr_z <- decostand(
        data.frame(
            sites %>% 
                filter(field_type %in% ft, region %in% rg) %>% 
                select(field_name, yr_since),
            row.names = 1),
        "standardize") %>% 
        rownames_to_column("field_name")
    # Create explanatory data frame and covariables matrix
    env <- if(is.data.frame(pspe_pcoa) == FALSE) {
        soil_cov_sc %>% 
            left_join(ptr_norm, by = join_by(field_name)) %>% 
            left_join(soil_expl_z, by = join_by(field_name)) %>% 
            left_join(yr_z, by = join_by(field_name)) %>% 
            column_to_rownames("field_name") %>% 
            drop_na()
    } else {
        soil_cov_sc %>% 
            left_join(pspe_ax, by = join_by(field_name)) %>% 
            left_join(soil_expl_z, by = join_by(field_name)) %>% 
            left_join(yr_z, by = join_by(field_name)) %>% 
            column_to_rownames("field_name") %>% 
            drop_na()
    }
    covars <- as.matrix(env[, 1:2])
    expl <- env[, 3:ncol(env)]
    # Forward select on explanatory data with covariables
    mod_null <- dbrda(fspe_bray ~ 1 + Condition(covars), data = expl, sqrt.dist = TRUE)
    mod_full <- dbrda(fspe_bray ~ . + Condition(covars), data = expl, sqrt.dist = TRUE)
    mod_step <- ordistep(mod_null, 
                         scope = formula(mod_full), 
                         direction = "forward", 
                         permutations = how(nperm = 1999), 
                         trace = FALSE)
    mod_glax <- anova(mod_step, step = 1999)
    mod_inax <- anova(mod_step, by = "axis", step = 1999)
    mod_scor <- scores(mod_step, choices = c(1,2), display = c("bp", "sites"), tidy = FALSE)
    return(list(
        plot_data = mod_scor,
        select_mod = mod_step,
        global_axis_test = mod_glax,
        individual_axis_test = mod_inax
    ))
} 
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
    select(field_name, field_type, region, everything(), -field_key)
#' Plant releve data was available from four Fermi sites, but survey data aren't correctable 
#' between Fermi and other sites. Also, translating counts of species to counts of traits isn't appropriate.
#' Trait data isn't included for the plant presence dataset. 
#' 
#' ### Plant communities
#' Both abundance and presence data are used. 
#' Plant site-species matrices use field names for rows, site metadata will be joined. 
pspe <- list(
    ab = read_csv(paste0(getwd(), "/clean_data/spe_plant_abund.csv"), show_col_types = FALSE),
    pr = read_csv(paste0(getwd(), "/clean_data/spe_plant_presence.csv"), show_col_types = FALSE)
) %>% map(function(x) x %>% 
              rename(field_name = SITE) %>% 
              left_join(sites %>% select(starts_with("field"), region), by = join_by("field_name")) %>% 
              select(field_name, field_type, region, everything(), -field_key))
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
#' ### Fungal biomass
fb <- read_csv(paste0(getwd(), "/clean_data/plfa.csv"), show_col_types = FALSE) %>% 
    left_join(sites %>% select(starts_with("field"), region), by = join_by("field_key", "field_name")) %>% 
    select(field_name, field_type, region, everything(), -field_key)
#' 
#' # Analysis and Results
#' A great number of symmetric and asymmetric comparative and constrained analyses have been attempted
#' with these data. The best and simplest is a partial db-RDA, executed with the custom function
#' `dbrda_fun()`. The strategy is to use most of the soil abiotic data and precipitation as 
#' covariables to remove site-dependent affects on the microbial community. These covariables are numerous,
#' so they are simplified by being normalized and then transformed into two PCA axes. Then, variables 
#' of experimental interest are forward selected. These are the variables that we expect to change 
#' as a result of the restoration, or affect the microbial community directly. 
#' Agricultural nutrients, soil organic matter, and plant community/traits data are in this group. 
#' 
#' As before, the plant data aren't consistent across all regions, so the analysis is adjusted accordingly as described 
#' below. Restricting this analysis to Blue Mounds only is probably the most defensible approach,
#' but with the covariables, other regions can be included too. It's also possible that a precise 
#' approach of the permutation scheme (with `how()` from package `permute`) would properly handle the 
#' design of blocks and replicates used here.  
#' 
#' ## Plant community axes
#' In this analysis, both plant community and traits data can be used. To analyze all regions, community
#' data based on presence/absence is the only possibility. Let's use `pcoa_fun()` to produce 
#' community axes for the abundance and presence/absence plant data. 
#+ plant_pcoa_ab
(pspe_pcoa_ab <- pcoa_fun(pspe$ab, rg = c("BM", "FG", "LP")))
#+ plant_pcoa_pr
(pspe_pcoa_pr <- pcoa_fun(pspe$pr, rg = c("BM", "FG", "LP", "FL"), method = "jaccard", binary = TRUE))
#' 
#' ## Partial db-RDA
#' ### General fungal community (ITS sequence abundance)
#' #### Blue Mounds with plant traits data
#+ dbrda_bm_tr_its
dbrda_fun(s = fspe$its, pspe_pcoa = "none", ft = c("restored"), rg = c("BM"))
#' No explanatory variables were selected
#' 
#' #### Wisconsin sites with plant traits
#+ dbrda_wi_tr_its
(dbrda_wi_tr_its <- dbrda_fun(s = fspe$its, pspe_pcoa = "none", ft = c("restored"), rg = c("BM", "LP", "FG")))
#' Global (p=0.011) and individual (p=0.015, dbRDA1) axis tests were significant. 
#' Years since restoration was selected, and it explains 17% of the variation. Forb and C4 grass are 
#' runners-up but appear highly correlated with years (not shown). 
#' 
#' #### Wisconsin sites with plant community axes
#+ dbrda_wi_ab_its
dbrda_fun(s = fspe$its, pspe_pcoa = pspe_pcoa_ab$site_vectors, ft = c("restored"), rg = c("BM", "LP", "FG"))
#' No explanatory variables were selected
#' 
#' #### All regions with plant community axes
#+ dbrda_all_pr_its
(dbrda_all_pr_its <- dbrda_fun(s = fspe$its, pspe_pcoa = pspe_pcoa_pr$site_vectors, ft = c("restored"), rg = c("BM", "LP", "FG", "FL")))
#' Global and individual axes are significant and strong. Potassium was selected and it explains 13% of the variation here.
#' It seems that this single variable must be driven by Fermi. 
#' 
#' K and plant 1 are runners up
#' 
#' ### AMF community (18S sequence abundance)
#' #### Blue Mounds with plant traits data
#+ dbrda_bm_tr_amf
dbrda_fun(s = fspe$amf, pspe_pcoa = "none", ft = c("restored"), rg = c("BM"))
#' No explanatory variables were selected
#' 
#' #### Wisconsin sites with plant traits data
#+ dbrda_wi_tr_amf
(dbrda_wi_tr_amf <- dbrda_fun(s = fspe$amf, pspe_pcoa = "none", ft = c("restored"), rg = c("BM", "LP", "FG")))
#' Global and single constrained axes are significant in site rank at p<0.05. Forb abundance was the selected
#' explanatory variable, explanaing 22% of the variation in communities. 
#' 
#' Forb and C4 grass are runners up
#' 
#' #### Wisconsin sites with plant community axes
#+ dbrda_wi_ab_amf
dbrda_fun(s = fspe$amf, pspe_pcoa = pspe_pcoa_ab$site_vectors, ft = c("restored"), rg = c("BM", "LP", "FG"))
#' Years
#' 
#' #### All regions with plant community axes
#+ dbrda_all_pr_amf
dbrda_fun(s = fspe$amf, pspe_pcoa = pspe_pcoa_pr$site_vectors, ft = c("restored"), rg = c("BM", "LP", "FG", "FL"))
#' Global and single constrained axes are significant in site rank at p<0.05. As before, Fermi brings 
#' potassium as a significant predictor of microbial communities. 




# Find the fields with high AMF and actinomycetes and check those.


site_sc <- dbrda_wi_tr_its$plot_data$sites
site_bp <- dbrda_wi_tr_its$plot_data$biplot

site_df <- site_sc %>% 
    data.frame() %>% 
    rownames_to_column(var = "field_name") %>% 
    left_join(sites, by = join_by(field_name))
bp_df <- site_bp %>% 
    data.frame() %>% 
    rownames_to_column(var = "envvar") %>% 
    mutate(
        origin = 0,
        m = MDS1 / dbRDA1, 
        d = sqrt(dbRDA1^2 + MDS1^2), 
        dadd = sqrt((max(dbRDA1)-min(MDS1))^2 + (max(MDS1)-min(MDS1))^2)*0.1,
        labx = ((d+dadd)*cos(atan(m)))*(dbRDA1/abs(dbRDA1)), 
        laby = ((d+dadd)*sin(atan(m)))*(dbRDA1/abs(dbRDA1)))
ggplot(site_df, aes(x = dbRDA1, y = MDS1)) +
    geom_text(aes(label = field_name)) +
    geom_segment(data = bp_df, 
                 aes(x = origin, xend = dbRDA1, y = origin, yend = MDS1), 
                 arrow = arrow(length = unit(3, "mm")),
                 color = "blue") +
    geom_text(data = bp_df, 
              aes(x = labx, y = laby, label = envvar),
              color = "blue") +
    theme_bw()
