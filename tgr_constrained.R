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
#' PCoA Function returns basic diagnostics
pcoa_fun <- function(data, ax, corr="none") {
    # ax = number of axes desired in the output
    p <- pcoa(data, correction = corr)
    p_vals <- data.frame(p$values) %>% 
        rownames_to_column(var = "Dim") %>% 
        mutate(Dim = as.integer(Dim))
    p_vec <- data.frame(p$vectors)
    
    if(corr == "none" | ncol(p_vals) == 6) {
        p_bstick <- ggplot(p_vals, aes(x = factor(Dim), y = Relative_eig)) + 
            geom_col(fill = "gray70", color = "gray30") + 
            geom_line(aes(x = Dim, y = Broken_stick), color = "red") +
            geom_point(aes(x = Dim, y = Broken_stick), color = "red") +
            labs(x = "Dimension", title = "PCoA Eigenvalues and Broken Stick Model") +
            theme_bw()
        p_ncomp <- with(p_vals, which(Relative_eig < Broken_stick)[1]-1)
        print(p_bstick)
    } else {
        p_bstick <- ggplot(p_vals, aes(x = factor(Dim), y = Rel_corr_eig)) + 
            geom_col(fill = "gray70", color = "gray30") + 
            geom_line(aes(x = Dim, y = Broken_stick), color = "red") +
            geom_point(aes(x = Dim, y = Broken_stick), color = "red") +
            labs(x = "Dimension", title = "PCoA Eigenvalues and Broken Stick Model") +
            theme_bw()
        p_ncomp <- with(p_vals, which(Rel_corr_eig < Broken_stick)[1]-1)
        print(p_bstick)
    }
    
    output <- list(
        values         = p_vals[1:(ax+1), ], 
        site_vectors   = p_vec[, 1:ax],
        comp_exceed_bs = p_ncomp)
    
    return(output)
}
#' CoIA runs a coinertia analysis and randomization test. 
coia_fun <- function(data1, data2, pco=TRUE) {
    if(isTRUE(pco)) {
        dudi1 <- dudi.pco(data1, scannf = FALSE)
        dudi2 <- dudi.pca(data2, scale = FALSE, center = FALSE, scannf = FALSE)
        coin <- coinertia(dudi1, dudi2, scannf = FALSE, nf = 2)
        print(summary(coin))
        print(procuste.randtest(dudi1$tab, dudi2$tab, 9999))
        
        return(coin)
        
    } else {
        dudi1 <- dudi.pca(data1, scale = FALSE, center = FALSE, scannf = FALSE)
        dudi2 <- dudi.pca(data2, scale = FALSE, center = FALSE, scannf = FALSE)
        coin <- coinertia(dudi1, dudi2, scannf = FALSE, nf = 2)
        
        print(summary(coin))
        print(randtest(coin, 9999))
        return(coin)
    }
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
            levels = c("corn", "restored", "remnant")),
        yr_since = replace(yr_since, which(field_type == "remnant"), "+"),
        yr_since = replace(yr_since, which(field_type == "corn"), "-")) %>%
    select(-lat, -long, -yr_restore, -yr_rank) %>% 
    arrange(field_key)
#' 
#' ## Plant data
#' ### Traits
#' Plant abundance data was only available at 16 sites (none at Fermi).
#' A vector of site names produced to filter species and env tables in later analysis.
plant_func_ab <- read_csv(paste0(getwd(), "/clean_data/plant_trait_abund.csv"), show_col_types = FALSE)
pab_fields <- plant_func_ab$field_name
#' Plant releve data was available from four Fermi sites. This was added to the abundance
#' data set, and then the entire matrix was transformed into presence-absence data. Note 
#' that survey effort are not correctable between Fermi and the other sites in this data
#' table. 
plant_func_pr <- read_csv(paste0(getwd(), "/clean_data/plant_trait_presence.csv"), show_col_types = FALSE)
ppr_fields <- plant_func_pr$field_name
#' 
#' ### Plant communities
#' Plant site-species matrices use field names for rows 
plant_spe_ab <- read_csv(paste0(getwd(), "/clean_data/spe_plant_abund.csv"), show_col_types = FALSE) %>% rename(field_name = SITE)
plant_spe_pr <- read_csv(paste0(getwd(), "/clean_data/spe_plant_presence.csv"), show_col_types = FALSE) %>% rename(field_name = SITE)
#' 
#' ## Environmental data
#' Use precipitation as proxy for soil moisture
ppt = read_csv(paste0(getwd(), "/clean_data/site_precip_normal.csv"), show_col_types = FALSE)
#' Create subsets of soil environmental data to align with plant abundance or presence sites
soil <- read_csv(paste0(getwd(), "/clean_data/soil.csv"), show_col_types = FALSE) %>% 
    left_join(ppt, by = join_by("field_key")) %>% 
    filter(field_key %in% sites$field_key) %>% 
    select(-field_key)
soil_pab <- soil %>% filter(field_name %in% pab_fields)
soil_ppr <- soil %>% filter(field_name %in% ppr_fields)
#' 
#' ## Microbial data
#' ### Fungal communities
#' Fungal species matrices must have field names instead of field keys to align all datasets.
#' Create subsets of fungal species matrices to align with plant abundance or presence sites
#' Sites-species tables contain rarefied sequence abundances. This list includes
#' composition summarized in fields. 
#' CSV files were produced in [process_data.R](process_data.md)
spe <- list(
    its = read_csv(paste0(getwd(), "/clean_data/spe_ITS_rfy.csv"), show_col_types = FALSE),
    amf = read_csv(paste0(getwd(), "/clean_data/spe_18S_rfy.csv"), show_col_types = FALSE)
) %>% map(function(x) x %>% 
              left_join(sites %>% select(field_key, field_name), by = join_by("field_key")) %>% 
              select(field_name, everything(), -field_key))
spe$its_pab <- spe$its %>% filter(field_name %in% pab_fields)
spe$its_ppr <- spe$its %>% filter(field_name %in% ppr_fields)
spe$amf_pab <- spe$amf %>% filter(field_name %in% pab_fields)
spe$amf_ppr <- spe$amf %>% filter(field_name %in% ppr_fields)
#' 
#' ### Species metadata
#' The OTUs and sequence abundances in these files matches the rarefied data in `spe$` above.
#' CSV files were produced in the [microbial diversity script](microbial_diversity.md).
spe_meta <- list(
    its = read_csv(paste0(getwd(), "/clean_data/speTaxa_ITS_rfy.csv"), show_col_types = FALSE),
    amf =  read_csv(paste0(getwd(), "/clean_data/speTaxa_18S_rfy.csv"), show_col_types = FALSE)
)
#' 
#' ### Fungal biomass
fa <- read_csv(paste0(getwd(), "/clean_data/plfa.csv"), show_col_types = FALSE) %>% 
    select(-field_key)
#' 
#' # Analysis and Results
#' This is a stepwise analysis procedure to work through to find the most useful 
#' analyses and visualization products. This can start with general fungi, also any 
#' steps which produce compelling results should proceed with AMF in series. 
#' 
#' Note that the plant presence data set, although it allows the use of more sites, 
#' is more complicated to justify. As a binary matrix, it's fine, but transforming 
#' that into functional group abundance is maybe questionable. 
#' 
#' 1. CoIA
#' This analysis isn't useful and will be removed. 
#'   a. Starting with plant abundance data and the 16 sites where we have this:
#'   b. Symmetrical ordination of plant site-species matrix and soil data will inform 
#'      the relative contributions of each. 
#'   c. Symmetrical ordination of plant abundance in functional groups will inform how 
#'      well the functional groups table aligns with the plant species data.
#'   d. A and B can be repeated with the plant presence-absence data, which expands 
#'      on the number of available sites, to compare results.
#'   e. Finally, comparing the plant abundance and plant presence matrices might help 
#'      justify (or not) their similarity. 
#' 2. Variation Partitioning
#'   a. Asymmetrical ordination but assesses relative contribution and interactions 
#'      among multiple ENV data sets. ENV must have fewer columns than sites. 
#'   b. Use db-RDA to produce the plant species response matrix. Use chord normalization
#'      to handle the soil data.
#'      as comparative ENV sets. 
#'   c. B can be repeated with the plant presence data set.
#' 3. Db-RDA
#'   a. Although this was done initially in variation partitioning, at this point 
#'      I will have more information about which environmental variables to include. 
#'      The key question is whether or not to include plant data, and how, since plant 
#'      data influences the number of available sites. 
#'   b. Step 1 is to forward select on all sites using all available environmental 
#'      variables, which at this point is soil only. 
#'   c. Step 2 is to use the information from 1 and 2 to add an additional constrained 
#'      analysis based on whatever sites are a best choice. 
#'   d. Try envfit or ordi-functions to overlay experimental items like yr_since
#' 4. LDA
#'   a. This might be worth trying to test any selected variables against a priori 
#'      treatment clusters. 
#' X. Biomass 
#'   a. 
#' X. Water stable aggregates
#'   a.
#' 
#' ## CoIA




pabenv <- 
    coia_fun(
        sqrt(vegdist(data.frame(plant_spe_ab, row.names = 1), "bray")),
        decostand(data.frame(soil_pab, row.names = 1), "standardize"),
        pco = TRUE
    )
plot(pabenv)
#' 
#' 
pabfgenv <- 
    coia_fun(
        decostand(data.frame(plant_func_ab, row.names = 1), "normalize"),
        decostand(data.frame(soil_pab, row.names = 1), "standardize"),
        pco = FALSE
    )
plot(pabfgenv)
#' 
#' 
pprenv <- 
    coia_fun(
        sqrt(vegdist(data.frame(plant_spe_pr, row.names = 1), "jaccard")),
        decostand(data.frame(soil_ppr, row.names = 1), "standardize"),
        pco = TRUE
    )
plot(pprenv)
pprfgenv <- 
    coia_fun(
        decostand(data.frame(plant_func_pr, row.names = 1), "normalize"),
        decostand(data.frame(soil_ppr, row.names = 1), "standardize"),
        pco = FALSE
    )
plot(pprfgenv)






#' In blue mounds restored fields, 
#' does SOM increase over time?
#' does WSA increase over time?
#' does biomass increase over time?
#' 
#' What is the effect of controlling for plant functional groups on the above questions?
#' 
#' does the microbial community change over time with soil nutrients controlled?
#' 
#' Are the microbes explained by plants (and field type maybe, resto vs remnant) in the presence of soil nutrients?
#' This could be done with all fields and field types...well in Wisconsin anyway. Could try with the plant presence too.
#' Could remove corn from this because it's too different. 
#' - Filtering soil nutrients into covariates or explanatory. Maybe keep the fertilizer ones as explanatory?
#' 
#' another try.
#' Are microbes explained by plant traits, ppt, and fertilizer nutrients with resto and remnant fields?
#' - Blue Mounds and Wisconsin
#' - ITS and 18S
#' Consider doing this procedure with only restored fields, and then with remnants, and see if they tell the same story
#' (you should redo plant traits without cornfields to see if resto and remnant are different)
#' Then, do axis values correlate with WSA, PLFA/NLFA, yr_since, SOM?
#' Do guilds correlate with axis values?
#' 
#' Could start with variation partitioning:
#' Effect of plant traits and soil nutrients on microbial data, assess contribution of each and together
#' Might have to do this with restored fields only. 
#' Then how do you get an axis to compare with responses? 
#' You'd need to do the forward selection db-RDA with covariables as described above. 
#' 
#' We already have years correlating with microbial communities in restored fields. How does that fact fit in?
#' 
#' 
#' 
#' ## Variation Partitioning
#' 