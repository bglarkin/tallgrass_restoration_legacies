#' ---
#' title: "Microbial data: fatty acids (biomass)"
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
#' Microbial data include site-species tables derived from high-throughput sequencing 
#' and PLFA/NLFA data which Ylva did. 
#' 
#' This presents basic visualizations of microbial biomass inferred with PLFA/NLFA 
#' quantification.
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
#+ plfa
fa <- read_csv(paste0(getwd(), "/clean_data/plfa.csv"), show_col_types = FALSE)
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
#+ summary_data
fa_meta <- 
    sites %>% 
    left_join(fa, by = c("field_key", "field_name"))
fa_grp <- 
    fa_meta %>% 
    select(-starts_with("fa_")) %>% 
    pivot_longer(cols = gram_pos:nlfa_plfa_ratio, names_to = "group", values_to = "qty")
#' 
#' # Results
#' ## Biomass in field types and regions
#' Let's first visualize the data across regions and field types
#+ wsa_visual_fig,fig.align='center'
ggplot(
    fa_grp %>% 
        group_by(field_type, region, group) %>% 
        summarize(qty_avg = mean(qty), .groups = "drop") %>% 
        filter(!(group %in% c("pseudo_amf", "nlfa_plfa_ratio"))),
    aes(x = field_type, y = qty_avg, group = region)) +
    facet_wrap(vars(group), scales = "free_y") +
    geom_point(data = fa_grp %>% filter(!(group %in% c("pseudo_amf", "nlfa_plfa_ratio"))), 
                                        aes(x = field_type, y = qty, fill = region), shape = 21) +
    geom_line(linetype = "dotted") +
    geom_point(aes(fill = region), shape = 21, size = 4) +
    scale_fill_discrete_qualitative(name = "Region", palette = "Dark2") +
    labs(x = "", y = "Biomass", 
         caption = "Large points show regional means. Small points show values from individual fields.") +
    theme_bw()
#' We see a variety of patterns across field types. Most often, biomass is highest in restored fields, with 
#' notable exceptions for actionmycetes and Fermilab. In Blue Mounds, the pattern is consistent across field types,
#' and the magnitude of difference isn't large. 
#' 
#' ## Ordination with PCA
fa_z <- decostand(data.frame(fa_meta %>% select(field_name, starts_with("fa")), row.names = 1), "standardize")
fa_pca <- rda(fa_z)
fa_pca %>% summary(., display = NULL)
#' Axes 1 and 2 explain 85% of the variation in sites. 
#+ fa_screeplot_fig,fig.align='center'
screeplot(fa_pca, bstick = TRUE)
#' Only the first eigenvalue exceeds the broken stick model, which is unsurprising because Axis 1 explains
#' 77% of the total variation. This will make these data very easy to use if it holds up. 
#+ fa_cleanplot_fig,fig.align='center'
cleanplot.pca(fa_pca)
#' This shows that all the biomass indices increase along Axis 1. Cornfields appear on the left with low biomass,
#' but as biomass increases, the signal is based on regions, with Blue Mounds, then Faville Grove, and then
#' Fermi fields appearing in order. Let's look at the summary fields to see if there's anything interesting 
#' left. 
fa_sum_z <- decostand(data.frame(fa_meta %>% 
                                 select(field_name, gram_pos, gram_neg, bacteria, fungi, actinomycetes, amf), 
                             row.names = 1), 
                  "standardize")
fa_sum_pca <- rda(fa_sum_z)
fa_sum_pca %>% summary(., display = NULL)
#' Axes 1 and 2 explain 96% of the variation in sites. 
#+ fa_sum_screeplot_fig,fig.align='center'
screeplot(fa_sum_pca, bstick = TRUE)
#' Only the first eigenvalue exceeds the broken stick model, which is unsurprising because Axis 1 explains
#' 81% of the total variation. This will make these data very easy to use if it holds up. 
#+ fa_sum_cleanplot_fig,fig.align='center'
cleanplot.pca(fa_sum_pca)
#' AMF have a strong correlation on a weak axis. It appears that some of the restored fields have the 
#' most AMF, which is what we saw in the visualization earlier. In this PCA, columns were standardized
#' to highlight effect sizes and de-emphasize magnitudes. This could be argued because gram_pos and gram_neg
#' are derivatives of bacteria, for example. This was run once with raw data (not shown), and the patterns
#' were the same. The second axis had more importance, showing that the signal with AMF is probably 
#' worth holding on to. 
#' 
#' ## Nutrients and microbial biomass
#' See [Soil abiotic properties](soil_properties.md) for a basic comparison and correlation. 
#' 
#' ## Biomass over field age
#' We can compare biomass with restoration age in Blue Mounds fields.
fa_meta %>% 
    filter(region == "BM", field_type == "restored") %>% 
    mutate(yr_since = as.numeric(yr_since)) %>% 
    select(field_name, yr_since, gram_pos, gram_neg, bacteria, fungi, actinomycetes, amf) %>% 
    ggpairs(columns = 2:ncol(.)) + theme_bw()
#' AMF decline with moderate strength as fields age. 
