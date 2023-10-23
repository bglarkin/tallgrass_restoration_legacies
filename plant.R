#' ---
#' title: "Plant data: communities and traits"
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
#' Plant data comprises two separate data sets. Mike Healy did quadrat surveys at 
#' all sites except Fermi, recording plant abundance. At Fermi, we only 
#' have relevé data with presence/absence. These data were provided by Mike Miller. 
#' Wisconsin sites were surveyed in Aug-Sept 2016. Fermi sites were surveyed in summer 2017
#' 
#' Plant metadata includes taxonomy and life history traits, and should cover
#' both plant data sets. With the abundance-data sites, trait data are reported in percent cover. 
#' In the presence-data sites, traits are in counts of species with that trait per field. 
#' 
#' This script produces basic visualization and diagnostic views of the plant data. Two traits matrices
#' are output, one for sites with abundance data (16 sites), and one for sites with presence data only (20 sites).  
#' 
#' Plant data may be used in matrix form for ordination, etc. Abundance and presence/absence
#' matrices are available. The data may also be transposed and summarized as abundance
#' in taxonomic or life history classes. This makes sense for the abundance data, but
#' maybe makes less sense with the presence/absence data because the summary would
#' really be richness, not abundance. 
#' 
#' Plants data were obtained in the [TRY Plant Trait Database](https://www.try-db.org/TryWeb/Home.php) 
#' ([Kattge et al. 2010](https://besjournals.onlinelibrary.wiley.com/doi/10.1111/j.2041-210X.2010.00067.x), 
#' [Kattge et al. 2011](https://onlinelibrary.wiley.com/doi/10.1111/j.1365-2486.2011.02451.x)) in 2016.
#' 
#' # Packages and libraries
packages_needed = c("GGally", "tidyverse", "vegan", "colorspace")
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
#' # Functions
#' Cleanplot PCA produces informative visualizations of PCA ordinations 
#' [(Borcard et al. 2018)](http://link.springer.com/10.1007/978-3-319-71404-2)
#+ cleanplot_pca
source(paste0(getwd(), "/supporting_files/cleanplot_pca.txt"))
#' 
#' # Data
#' Plant community data includes:
#' 
#' - Metadata, taxonomy and traits
#' - Abundance data, surveyed in 2016, sites limited to Wisconsin only
#' - Presence data, from Fermi in 2015, all other sites converted to presence data, does not 
#' include Fermi switchgrass or corn fields.
#+ plant_data_list
plant <- list(
    meta = read_csv(paste0(getwd(), "/clean_data/spe_plant_meta.csv"), show_col_types = FALSE) %>% 
        rename_with(tolower),
    ab   = read_csv(paste0(getwd(), "/clean_data/spe_plant_abund.csv"), show_col_types = FALSE) %>% 
        rename(field_name = SITE),
    pr   = read_csv(paste0(getwd(), "/clean_data/spe_plant_presence.csv"), show_col_types = FALSE) %>% 
        rename(field_name = SITE)
)
#' Metadata from sites, as in previous
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
#+ sites_no_cornfields
sites_noc <- sites %>% 
    filter(field_type != "corn")
#' 
#' ## Data wrangling: traits data
#' ### Abundance data (16 sites)
#' The alignment of plant codes among files must be confirmed. Show any mismatched codes.
which(!(colnames(plant$ab[, -1]) %in% plant$meta$code))
which(!(colnames(plant$pr[, -1]) %in% plant$meta$code))
#' The codes are aligned. Factor levels in metadata must be cleaned up and simplified. 
recode_group <- c(`annual/biennial` = "biennial", `annual/perennial` = "perennial", `C4-grass` = "C4_grass", `C3-grass` = "C3_grass")
filter_group <- c("bare", "litter", "unknown", "fern", "moss", "sedge", "Unknown", "shrub", "tree", "C3/C4", "rush")
#' Plant traits are summarized among the set of sites with plant abundance data
#+ p_ab_trait
p_ab_trait <-     
    plant$ab %>% 
    pivot_longer(-field_name, names_to = "code", values_to = "pct_cvr") %>% 
    filter(pct_cvr > 0) %>% 
    left_join(plant$meta %>% select(-(genus:name), -lifeform, -family), by = "code") %>% 
    pivot_longer(lifehist:category, names_to = "variable", values_to = "group", values_drop_na = TRUE) %>% 
    mutate(group = recode(group, !!!recode_group)) %>% 
    group_by(field_name, variable, group) %>% 
    summarize(pct_cvr = sum(pct_cvr), .groups = "drop") %>% 
    filter(!(group %in% filter_group)) %>% 
    select(-variable) %>% 
    arrange(group) %>%
    pivot_wider(names_from = group, values_from = pct_cvr, values_fill = 0) %>% 
    select(field_name, annual, biennial, perennial, native, nonnative, C3_grass, C4_grass, forb, legume, shrubTree)
#' 
#' ### Presence data (20 sites)
#+ p_pr_trait
p_pr_trait <- 
    plant$pr %>% 
    pivot_longer(-field_name, names_to = "code", values_to = "count") %>% 
    filter(count > 0) %>% 
    left_join(plant$meta %>% select(-(genus:name), -lifeform, -family), by = "code") %>% 
    pivot_longer(lifehist:category, names_to = "variable", values_to = "group", values_drop_na = TRUE) %>% 
    mutate(group = recode(group, !!!recode_group)) %>% 
    group_by(field_name, variable, group) %>% 
    summarize(count = sum(count), .groups = "drop") %>% 
    filter(!(group %in% filter_group)) %>% 
    select(-variable) %>% 
    arrange(group) %>%
    pivot_wider(names_from = group, values_from = count, values_fill = 0) %>% 
    select(field_name, annual, biennial, perennial, native, nonnative, C3_grass, C4_grass, forb, legume, shrubTree)
#' 
#' # Results
#' ## Visualization and output of abundance data
#' Run a PCA on chord-transformed traits data from sites with abundance data, perform 
#' typical basic diagnostics. This should be done without corn fields because they exert 
#' too strong a difference on everything else. 
p_ab_trait_ch <- decostand(data.frame(p_ab_trait %>% filter(field_name %in% sites_noc$field_name), row.names = 1), "normalize")
p_ab_trait_pca <- rda(p_ab_trait_ch)
p_ab_trait_pca %>% summary(., display = NULL)
#+ p_ab_trait_pca_scree,fig.align='center'
screeplot(p_ab_trait_pca, bstick = TRUE)
#+ p_ab_trait_pca_cleanplot,fig.align='center'
cleanplot.pca(p_ab_trait_pca)
#' Axis 1 & 2 explain 77% of the variation, and both eigenvalues exceed the broken stick model.
#' Traits forb, perennial, annual, nonnative, and C4 grass exceed the unit circle, suggesting a strong correlation
#' with site differences. Traits appear collinear, explore which ones produce high VIF. 
sort(diag(solve(cor(data.frame(p_ab_trait, row.names = 1)))), decreasing = TRUE) 
#' Many are very high. Traits levels need not be mutually exclusive, but it appears here that they are. 
#' Native opposes nonnative abundance, perennial opposes annual, and forb opposes C4 grass. The last case
#' may be driven by Karla Ott's field, but in other cases. What happens to VIF when one of these opposing factor levels 
#' are removed?
sort(diag(solve(cor(data.frame(p_ab_trait, row.names = 1) %>% select(-annual, -nonnative, -C4_grass)))), decreasing = TRUE)
#' It's clear that some are correlated and there is good reason to remove them. It's useful here because 
#' even after these are removed, they suggest their opposite factor level. This will help inform forward selection later. 
#' Export the traits matrix for sites with abundance data:
#+ p_ab_trait_export
write_csv(p_ab_trait, paste0(getwd(), "/clean_data/plant_trait_abund.csv"))
#' 
#' ## Visualization and output of presence data
#' Run a PCA on chord-transformed traits data from sites with abundance data, perform 
#' typical basic diagnostics.
p_pr_trait_ch <- decostand(data.frame(p_pr_trait %>% filter(field_name %in% sites_noc$field_name), row.names = 1), "normalize")
p_pr_trait_pca <- rda(p_pr_trait_ch)
p_pr_trait_pca %>% summary(., display = NULL)
#+ p_pr_trait_pca_scree,fig.align='center'
screeplot(p_pr_trait_pca, bstick = TRUE)
#+ p_pr_trait_pca_cleanplot,fig.align='center'
cleanplot.pca(p_pr_trait_pca)
#' Axis 1 & 2 explain 79% of the variation. The axis 1 eigenvalue is the only one which exceeds a broken stick model.
#' Traits forb, biennial, perennial, native, and nonnative exceed the unit circle, suggesting a strong correlation
#' with site differences. Traits appear collinear, explore which ones produce high VIF. 
sort(diag(solve(cor(data.frame(p_pr_trait, row.names = 1)))), decreasing = TRUE) 
#' Many are very high, and lie in opposition by factor levels as in the abundance data, but possibly 
#' a little less strong in terms of pure linear correlation. Biennial probably describes forbs, and can 
#' easily be discarded. It also looks like most perennials are also native? 
sort(diag(solve(cor(data.frame(p_pr_trait, row.names = 1) %>% select(-nonnative, -biennial, -native)))), decreasing = TRUE)
#' It's clear that some are correlated and there is good reason to remove them. It's useful here because 
#' even after these are removed, they suggest their opposite factor level. This will help inform forward selection later. 
#' Export the traits matrix for sites with abundance data:
#+ p_pr_trait_export
write_csv(p_ab_trait, paste0(getwd(), "/clean_data/plant_trait_abund.csv"))
#' 