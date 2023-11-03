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
#' This is a stepwise analysis procedure to work through to find the most useful 
#' analyses and visualization products. This can start with general fungi, also any 
#' steps which produce compelling results should proceed with AMF in series. 
#' 
#' Note that the plant presence data set, although it allows the use of more sites, 
#' is more complicated to justify. As a binary matrix, it's fine, but transforming 
#' that into functional group abundance is maybe questionable. 
#' 
#' 1. Variation Partitioning
#'   a. Asymmetrical ordination but assesses relative contribution and interactions 
#'      among multiple ENV data sets. ENV must have fewer columns than sites. 
#'   b. Use db-RDA to produce the plant species response matrix. Use chord normalization
#'      to handle the soil data.
#'      as comparative ENV sets. 
#'   c. B can be repeated with the plant presence data set.
#' 1. Db-RDA
#'   a. Although this was done initially in variation partitioning, at this point 
#'      I will have more information about which environmental variables to include. 
#'      The key question is whether or not to include plant data, and how, since plant 
#'      data influences the number of available sites. 
#'   b. Step 1 is to forward select on all sites using all available environmental 
#'      variables, which at this point is soil only. 
#'   c. Step 2 is to use the information from 1 and 2 to add an additional constrained 
#'      analysis based on whatever sites are a best choice. 
#'   d. Try envfit or ordi-functions to overlay experimental items like yr_since
#' 1. LDA
#'   a. This might be worth trying to test any selected variables against a priori 
#'      treatment clusters. 
#' 1. Biomass 
#'   a. 
#' 1. Water stable aggregates
#'   a.
#' 







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
#' This analysis will show the relative explanatory power of soil and/or plant trait data on microbial communities. 
#' It's useful here as it will also reveal the more important explanatory variables from either
#' data set through forward selection. 
#' 
#' With so many possible ways to parse the data, it will be easy to get carried away with variation 
#' partitioning. Let's exclude some potentials here. Cornfields won't be included because soil and plant 
#' data are discontinuous with the rest of the fields. Should restored fields be included? Using plant
#' data as explanatory makes sense in restored fields because the plants were seeded as they would be 
#' in a designed experiment. This might be less relevant in remnant fields where plants and soils have 
#' co-filtered through feedbacks. 
#' 
#' Plant trait data restricts this to only those fields which have plant abundance because counts 
#' traits isn't appropriate to use. If the larger 
#' number of fields is important to look at later, variation partitioning could be run with the plant
#' presence community data transformed into axes instead. Finally, it must be decided whether to include
#' all the plant abundance sites or just use Blue Mounds because it's closest to a chronosequence. 
#' Since years since restoration doesn't fit well in variation partitioning as a variable, I don't think
#' we need to stick with the best chronosequence, though. But, if we end up doing that later, it would
#' be good to know what the outcome of variation partitioning was.  
#' 
#' Maybe it's best to create a list of possible analyses and rank them from most to least important.
#' 
#' - Plant abundance regions, restored sites only, 18S and ITS
#'    - If warranted, Blue Mounds, restored sites only, 18S and ITS
#' - Plant abundance regions, restored and remnant sites, 18S and ITS
#' - Plant presence regions, restored and remnant sites, 18S and ITS
#'    - Use significant axes from PCoA for plant community data. 


# LETS MAKE THIS INTO A FUNCTION! Possible, or is each analysis too different?
# First, just fiddle around with this and keep a list of what works and what doesn't. 



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

pspe_pcoa_ab <- pcoa_fun(pspe$ab, rg = c("BM", "FG", "LP"))
pspe_pcoa_pr <- pcoa_fun(pspe$pr, rg = c("BM", "FG", "LP", "FL"), method = "jaccard", binary = TRUE)


s = fspe$its
pspe_pcoa = "none"
ft = c("restored")
rg = c("BM")


#' 
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
    # Create explanatory data frame and covariables matrix
    env <- if(is.data.frame(pspe_pcoa) == FALSE) {
        soil_cov_sc %>% 
            left_join(ptr_norm, by = join_by(field_name)) %>% 
            left_join(soil_expl_z, by = join_by(field_name)) %>% 
            left_join(sites %>% select(field_name, yr_since), by = join_by(field_name)) %>% 
            column_to_rownames("field_name") %>% 
            drop_na()
    } else {
        soil_cov_sc %>% 
            left_join(pspe_ax, by = join_by(field_name)) %>% 
            left_join(soil_expl_z, by = join_by(field_name)) %>% 
            left_join(sites %>% select(field_name, yr_since), by = join_by(field_name)) %>% 
            column_to_rownames("field_name") %>% 
            drop_na()
    }
    covars <- as.matrix(env[, 1:2])
    expl <- env[, 3:ncol(env)]
    # Forward select on explanatory data with covariables
    mod_null <- dbrda(fspe_bray ~ 1 + Condition(covars), data = expl, sqrt.dist = TRUE)
    mod_full <- dbrda(fspe_bray ~ . + Condition(covars), data = expl, sqrt.dist = TRUE)
    mod_step <- ordistep(mod_null, scope = formula(mod_full), direction = "forward", permutations = 1999, trace = FALSE)
    return(list(
        select_mod = mod_step,
        covars,
        expl
    ))
}


#' dbrda
#' four approaches with restored only fields, with its and amf
#' plant abundance sites with soil and plant functional trait vars
#' plant abundance sites with soil and plant community axes (no axes significant)
#' Blue Mounds sites with soil and plant functional trait vars
#' plant presence sites with soil and plant community axes (no axes significant)

dbrda_fun(s = fspe$its, pspe_pcoa = "none", ft = c("restored"), rg = c("BM")) # none
dbrda_fun(s = fspe$its, pspe_pcoa = "none", ft = c("restored"), rg = c("BM", "LP", "FG")) # forb, c4grass depending on model run
dbrda_fun(s = fspe$its, pspe_pcoa = pspe_pcoa_ab$site_vectors, ft = c("restored"), rg = c("BM", "LP", "FG")) # none
dbrda_fun(s = fspe$its, pspe_pcoa = pspe_pcoa_pr$site_vectors, ft = c("restored"), rg = c("BM", "LP", "FG", "FL")) # plant2, K depending on run





dbrda_fun(s = fspe$amf, pspe_pcoa = "none", ft = c("restored"), rg = c("BM")) # none
dbrda_fun(s = fspe$amf, pspe_pcoa = "none", ft = c("restored"), rg = c("BM", "LP", "FG")) # forb
dbrda_fun(s = fspe$amf, pspe_pcoa = pspe_pcoa_ab$site_vectors, ft = c("restored"), rg = c("BM", "LP", "FG")) # none
dbrda_fun(s = fspe$amf, pspe_pcoa = pspe_pcoa_pr$site_vectors, ft = c("restored"), rg = c("BM", "LP", "FG", "FL")) # none, K depending on run


# Find the fields with high AMF and actinomycetes and check those.
