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
#' This presents basic visualizations of microbial biomass inferred with PLFA/NLFA 
#' quantification done by YL.
#' 
#' Biomass may also correlate with water stable aggregation, so we'll look at that too. 
#'
#' **Note:** Only fatty acid 18.2 is used for fungi because 18.2w9 is also found
#' in gram-negative bacteria. 
#' 
#' # Packages and libraries
packages_needed = c("GGally", "tidyverse", "vegan", "colorspace", "ggbeeswarm")
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
    rename(fungi_18.2 = fa_18.2) %>% 
    select(-starts_with("fa_"), -fungi) %>% 
    pivot_longer(cols = fungi_18.2:nlfa_plfa_ratio, names_to = "group", values_to = "qty")
#+ wsa_data
# Remove rows from old field sites (26 and 27)
wsa <- read_csv(paste0(getwd(), "/clean_data/wsa.csv"), show_col_types = FALSE)[-c(26:27), ] %>% 
    left_join(sites, by = "field_key")
#' 
#' # Results
#' ## Biomass in field types and regions
#' Let's first visualize the data across regions and field types
#+ fa_visual_fig,fig.align='center'
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
#' Let's make a figure that can work in the summary.
#+ fa_boxplot,fig.align='center',fig.width=6,fig.height=3.5
ggplot(fa_grp %>% filter(group %in% c("amf", "fungi")), aes(x = field_type, y = qty)) +
    facet_wrap(vars(group), scales = "free_y") +
    geom_boxplot(fill = "gray90", varwidth = FALSE, outlier.shape = NA) +
    geom_beeswarm(aes(shape = region, fill = field_type), size = 2, dodge.width = 0.3) +
    labs(y = expression(Biomass~(nmol%*%g[soil]^-1))) +
    scale_fill_discrete_qualitative(name = "Field Type", palette = "Harmonic") +
    scale_shape_manual(name = "Region", values = c(21, 22, 23, 24)) +
    theme_bw() +
    theme(axis.title.x = element_blank()) +
    guides(fill = guide_legend(override.aes = list(shape = 21)))
#' 
#' And one that shows the correlation between AMF and time. Note: there is no relationship between fungi and time.
#+ fa_amf_reg
fa_grp %>% 
    filter(group == "amf", field_type == "restored", region == "BM") %>% 
    lm(yr_since ~ qty, data = .) %>% 
    summary()
#+ fa_amf_yrs,fig.align='center',fig.width=3,fig.height=3.5
fa_grp %>% 
    filter(group == "amf", field_type == "restored", region == "BM") %>% 
    ggplot(aes(x = as.numeric(yr_since), y = qty)) +
    geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.4) +
    geom_point(fill = "#5CBD92", shape = 21, size = 2.5) +
    labs(x = "Years since restoration", y = expression(Biomass~(nmol%*%g[soil]^-1))) +
    theme_bw()
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
#+ fa_yr_since_fig,fig.align='center'
fa_meta %>% 
    filter(region == "BM", field_type == "restored") %>% 
    mutate(yr_since = as.numeric(yr_since)) %>% 
    select(field_name, yr_since, gram_pos, gram_neg, bacteria, fungi, actinomycetes, amf) %>% 
    ggpairs(columns = 2:ncol(.)) + theme_bw()
#' AMF decline with moderate strength as fields age. 
#' 
#' ## Biomass and WSA
#+ biomass_wsa_wisconsin_fig,fig.align='center'
fa_grp %>% 
    filter(group %in% c("fungi_18.2", "amf"), 
           region != "FL") %>% 
    pivot_wider(names_from = "group", values_from = "qty", names_prefix = "mass_") %>% 
    left_join(wsa %>% select(field_key, wsa), by = join_by(field_key)) %>% 
    mutate(field_type = factor(field_type, ordered = FALSE),
           yr_since = as.integer(yr_since)) %>% 
    ggpairs(columns = 5:8, ggplot2::aes(color = field_type, shape = region))
#' Whatever is going on between fungal biomass and wsa isn't related to age, it might be related to 
#' something in the plant community or another site variable, but it's not strong and doesn't 
#' immediately appear to be related to any of the primary questions here. 
#+ biomass_wsa_bm_fig,fig.align='center'
fa_grp %>% 
    filter(group %in% c("fungi_18.2", "amf"), 
           region == "BM", 
           field_type == "restored") %>% 
    pivot_wider(names_from = "group", values_from = "qty", names_prefix = "mass_") %>% 
    left_join(wsa %>% select(field_key, wsa), by = join_by(field_key)) %>% 
    mutate(field_type = factor(field_type, ordered = FALSE),
           yr_since = as.integer(yr_since)) %>% 
    ggpairs(columns = 5:8)
#' Again, no relationship between biomass and wsa. 
