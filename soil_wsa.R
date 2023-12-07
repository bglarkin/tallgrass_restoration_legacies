#' ---
#' title: "Soil: water stable aggregates"
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
#' Water stable aggregates were determined. WSA is a functional attribute of soil
#' created by fungi. WSA could be added to the soil chemical analysis data table. 
#' WSA isn't precisely a soil chemical feature (these are based more on parent material
#' and weathering, while WSA are based on fungal activity), but they are related. 
#' Whether to join them or not probably depends on the test desired. The first thing
#' to do there is just see if WSAs vary in any meaningful way that could help 
#' other analyses. 
#' 
#' I don't know who obtained these data or how it was done. Need to follow up.  
#' 
#' # Packages and libraries
packages_needed = c("GGally", "rsq", "lme4", "multcomp", "tidyverse", "ggbeeswarm", "knitr", "conflicted", "colorspace")
packages_installed = packages_needed %in% rownames(installed.packages())
#+ packages,message=FALSE
if (any(!packages_installed)) {
    install.packages(packages_needed[!packages_installed])
}
#+ libraries,message=FALSE
for (i in 1:length(packages_needed)) {
    library(packages_needed[i], character.only = T)
}
#+ conflicts
conflicts_prefer(dplyr::select)
conflicts_prefer(dplyr::filter)
#' 
#' # Data
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
#+ water_stable_aggregates
# Remove rows from old field sites (26 and 27)
wsa <- read_csv(paste0(getwd(), "/clean_data/wsa.csv"), show_col_types = FALSE)[-c(26:27), ] %>% 
    left_join(sites, by = "field_key")
#' 
#' # Results
#' ## WSA in field types and regions
#' Let's first visualize the data across regions and field types
#+ wsa_visual_fig,fig.align='center'
ggplot(wsa %>% group_by(field_type, region) %>% summarize(wsa_avg = mean(wsa), .groups = "drop"),
       aes(x = field_type, y = wsa_avg, group = region)) +
    geom_point(data = wsa, aes(x = field_type, y = wsa, fill = region), shape = 21) +
    geom_line(linetype = "dotted") +
    geom_point(aes(fill = region), shape = 21, size = 5) +
    scale_fill_discrete_qualitative(name = "Region", palette = "Dark2") +
    labs(x = "", y = "Water stable aggregates (%)", 
         caption = "Large points show regional means. Small points show values from individual fields.") +
    theme_bw()
#' A minor interaction appears, with WSA in restored fields being higher in all regions except Lake Petite. 
#' Let's try a mixed model with region as a random effect and test the difference in WSAs across field types.
wsa_mod <- lmer(wsa ~ field_type + (1 | region), data = wsa)
wsa_mod_null <- lmer(wsa ~ 1 + (1 | region), data = wsa)
wsa_mod_tuk <- glht(wsa_mod, linfct = mcp(field_type = "Tukey"), test = adjusted("holm"))
#' Model summaries: display results of fitted model, null model, and the likelihood ratio test of the 
#' null vs. tested model to assess significance of field type. 
summary(wsa_mod)
#+ lik_ratio_test,message=FALSE
print(anova(wsa_mod, wsa_mod_null, REML = FALSE))
summary(wsa_mod_tuk)
cld(wsa_mod_tuk)
#+ wsa_effect_table
wsa %>% 
    group_by(field_type) %>% 
    summarize(wsa_avg = round(mean(wsa), 1), .groups = "drop") %>% 
    mutate(sig = c("a", "b", "a")) %>% 
    kable(format = "pandoc")
#' Percent water stable aggregates differed by field type based on a likelihood ratio 
#' test ($\chi^2(2)=12.23,~p<0.05$). In restored fields, water stable aggregates were 14.2% higher 
#' than in corn fields and 13% higher than in remnant fields based on Tukey's post hoc test ($p<0.05$). 
#' Let's view a boxplot of the result. 
sig_letters <- data.frame(
    lab = c("a", "b", "a"),
    xpos = c(1,2,3),
    ypos = rep(85,3)
)
#+ wsa_regions_fieldtypes_fig,fig.align='center',fig.width=4.5,fig.height=3.5
ggplot(wsa, 
       aes(x = field_type, y = wsa)) +
    geom_boxplot(fill = "gray90", varwidth = FALSE, outlier.shape = NA) +
    geom_beeswarm(aes(fill = region), shape = 21, size = 2, dodge.width = 0.2) +
    geom_text(data = sig_letters, aes(x = xpos, y = ypos, label = lab)) +
    labs(x = "", y = "Water stable aggregates (%)") +
    scale_fill_discrete_qualitative(name = "Region", palette = "Dark2") +
    theme_bw()
#' In the figure above, Linear mixed models with region as random effect. 
#' Letters show differences based on Tukey's post hoc with Holm correction at p<0.05
#' ## WSA over time in restored fields
#' Percent WSA varies greatly across restored fields. Can some of this variation be explained by 
#' how long it's been since the field was restored? This comparison can only be justified in the Blue 
#' Mounds area. 
wsa_bm_resto <- 
    wsa %>% 
    filter(field_type == "restored", region == "BM") %>% 
    mutate(yr_since = as.numeric(yr_since))
wsa_bm_resto_lm <- 
    lm(wsa ~ yr_since, data = wsa_bm_resto)
summary(wsa_bm_resto_lm)
#' WSA doesn't change based on time in the Blue Mounds restored fields. 