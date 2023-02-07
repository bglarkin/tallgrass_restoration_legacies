#' ---
#' title: "Microbial data: overview of data, diversity statistics"
#' author: "Beau Larkin\n"
#' date: "Last updated: `r format(Sys.time(), '%d %B, %Y')`"
#' output:
#'   github_document:
#'     toc: true
#'     toc_depth: 3
#'     df_print: paged
#' ---
#'
#' # Description
#' Microbial data include site-species tables derived from high-throughput sequencing 
#' and PLFA/NLFA extractions and measurement. Lipid workflow was completed by Ylva Lekberg. 
#' 
#' The overview here presents basic statistics and visualizations of diversity in the  
#' microbial species data. 
#' 
#' - Diversity and evenness of microbial communities
#' - Interpretation of differences in diversity among regions and field types, and over years.
#' 
#' # Packages and libraries
packages_needed = c(
    "rsq",
    "lme4",
    "multcomp",
    "tidyverse",
    "vegan",
    "ggbeeswarm",
    "knitr",
    "conflicted",
    "colorspace"
)
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
#' # Data
#' ## Sites-species tables
#' CSV files were produced in `process_data.R`
its_otu_all <- read_csv(paste0(getwd(), "/clean_data/spe_ITS_otu_siteSpeMatrix_allReps.csv"), 
                        show_col_types = FALSE)
its_sv_all  <- read_csv(paste0(getwd(), "/clean_data/spe_ITS_sv_siteSpeMatrix_allReps.csv"), 
                        show_col_types = FALSE)
amf_otu_all <- read_csv(paste0(getwd(), "/clean_data/spe_18S_otu_siteSpeMatrix_allReps.csv"), 
                        show_col_types = FALSE)
amf_sv_all  <- read_csv(paste0(getwd(), "/clean_data/spe_18S_sv_siteSpeMatrix_allReps.csv"), 
                        show_col_types = FALSE)
#' 
#' ## Average sequence abundances in each field
#' We examine diversity at the field level, so diversity obtained at samples should be averaged 
#' for each field. We collected ten samples from each field, but processing failed for some samples
#' at one step or another in the pipeline. Samples must be randomly resampled to the smallest 
#' number obtained in a series to produce comparable diversity metrics. The following function 
#' will resample the site-species data to the correct number of samples. 
#+ resample_fields_function 
resample_fields <- function(data, min, cluster_type) {
    set.seed(482)
    data %>% 
        group_by(site_key) %>% 
        slice_sample(n = min) %>%
        summarize(across(starts_with(cluster_type), list(avg = mean)))
}
#' 
#' The minimum number of samples in a field for each gene is:
#' 
#' - ITS = 8 samples
#' - 18S = 7 samples
#' 
#' With this, we can run the function for each dataset:
its_otu_avg <- resample_fields(its_otu_all, 8, "otu")
its_sv_avg  <- resample_fields(its_sv_all,  8, "sv")
amf_otu_avg <- resample_fields(amf_otu_all, 7, "otu")
amf_sv_avg  <- resample_fields(amf_sv_all,  7, "sv")
#' 
#' ## Site metadata and design
#' Set remnants to 50 years old as a placeholder. This number will not be used in 
#' a quantitative sense, for example in models. 
#' Oldfields are filtered out because they could not be replicated in regions. 
rem_age <- 50
sites   <- read_csv(paste0(getwd(), "/clean_data/site.csv"), show_col_types = FALSE) %>% 
    mutate(site_type = factor(site_type, ordered = TRUE, levels = c("corn", "restored", "remnant")),
           yr_since = replace(yr_since, which(site_type == "remnant"), rem_age)) %>% 
    filter(site_type != "oldfield") %>% 
    rename(field_type = site_type)
#' 
#' # Analysis and Results
#' ## Diversity
#' Microbial diversity is considered for each of four datasets: OTU or SV clustering for 18S or ITS gene 
#' sequencing. For each set, Hill's numbers are produced ([Hill 1973](http://doi.wiley.com/10.2307/1934352), 
#' [Borcard and Legendere 2018, p. 373](http://link.springer.com/10.1007/978-3-319-71404-2)) and plotted,
#' with means differences tested using mixed-effects linear models in `lmer` ([Bates et al. 2015](https://doi.org/10.18637/jss.v067.i01)).
#' Correlations are then produced to visualize change in diversity trends over time, 
#' with similar mixed-effects tests performed. 
#' 
#' Hill's numbers, brief description:
#' 
#' - $N_{0}$  = species richness
#' - $N_{1}$  = Shannon's diversity ($e^H$; excludes rarest species, considers the number of "functional" species)
#' - $N_{2}$  = Simpson's diversity ($1 / \lambda$; number of "codominant" species)
#' - $E_{10}$ = Shannon's evenness (Hill's ratio $N_{1} / N_{0}$)
#' - $E_{20}$ = Simpson's evenness (Hill's ratio $N_{2} / N_{0}$)
#' 
#' ### Functions and variables
#' The following functions are used to streamline code and reduce errors:
#' 
#' #### Calculate Hill's series on a samples-species matrix
#+ calc_diversity_function
calc_diversity <- function(spe) {
    spe_mat <- data.frame(spe, row.names = 1)
    
    N0  <- apply(spe_mat > 0, MARGIN = 1, FUN = sum)
    N1  <- exp(diversity(spe_mat))
    N2  <- diversity(spe_mat, "inv")
    E10 <- N1 / N0
    E20 <- N2 / N0
    
    return(
        data.frame(N0, N1, N2, E10, E20) %>%
            rownames_to_column(var = "site_key") %>%
            mutate(site_key = as.integer(site_key)) %>%
            left_join(
                sites %>% select(starts_with("site"), field_type, region, yr_rank, yr_since),
                by = "site_key"
            ) %>%
            pivot_longer(
                cols = N0:E20,
                names_to = "hill_index",
                values_to = "value"
            ) %>%
            mutate(hill_index = factor(
                hill_index,
                ordered = TRUE,
                levels = c("N0", "N1", "N2", "E10", "E20")
            ))
    )
}
#' 
#' #### Test diversity measures across site types with mixed model
#+ test_diversity_function
test_diversity <- function(data) {
    hills <- levels(data$hill_index)
    for(i in 1:length(hills)) {
        print(hills[i])
        mod_data <- data %>% 
            filter(hill_index == hills[i]) %>% 
            mutate(field_type = factor(field_type, ordered = FALSE))
        mmod <- lmer(value ~ field_type + (1 | region), data = mod_data, REML = FALSE)
        mmod_null <- lmer(value ~ 1 + (1 | region), data = mod_data, REML = FALSE)
        print(anova(mmod, mmod_null))
        mod_tuk <- glht(mmod, linfct = mcp(field_type = "Tukey"), test = adjusted("holm"))
        print(mod_tuk)
        print(cld(mod_tuk))
    }
}
#' #### Change in diversity over time 
#' Do Hill's numbers correlate with years since restoration?
#' This is only appropriate to attempt in the Blue Mounds region, and even there, it will be difficult
#' to justify that the area meets the criteria for a chronosequence. 
test_age <- function(data, caption=NULL) {
    temp_df <-
        data %>%
        filter(field_type == "restored", region == "BM") %>%
        pivot_wider(names_from = hill_index, values_from = value) %>%
        select(-starts_with("site"),-field_type,-region,-yr_rank)
    lapply(temp_df %>% select(-yr_since), function(z) {
        test <-
            cor.test(temp_df$yr_since,
                     z,
                     alternative = "two.sided",
                     method = "pearson")
        return(data.frame(
            cor = round(test$estimate, 2),
            R2 = round(test$estimate^2, 2),
            pval = round(test$p.value, 3)
        ))
    }) %>%
        bind_rows(.id = "hill_num") %>%
        remove_rownames() %>%
        mutate(sig = case_when(pval <= 0.05 ~ "*", TRUE ~ "")) %>%
        kable(format = "pandoc", caption = caption)
}
#' 
#' #### Calculate diversity of all samples-species matrices
#' Create list of matrices and process with `calc_diversity()`. Naming list objects as
#' their desired output names will enhance understanding later. 
#+ spe_avg_temp_list
spe_list <- list(
    its_otu = its_otu_avg,
    its_sv = its_sv_avg,
    amf_otu = amf_otu_avg,
    amf_sv = amf_sv_avg
)
#+ diversity_calculations
div <- lapply(spe_list, calc_diversity)
#' 
#' ### Fungi (ITS gene) in OTU clusters, averaged to 8 samples per field.
#' Run the linear model and test differences among field types for diversity.
#+ test_div_its_otu
test_diversity(div$its_otu)
#' Key model results depend on which sites were sampled using `resample_fields()` above. 
#' Changing the value of `set.seed()` in that function may alter the results of this model. 
#' 
#' - $N_{0}$: field type is significant by likelihood ratio test at p<0.001 with region as
#' a random effect. Species richness in corn fields was less than restored or remnants, which 
#' didn't differ at p=0.05. 
#' - $N_{1}$: model fit is questionable due to singular fit, but field type is significant
#' by likelihood ratio test at p<0.01 with region as a random effect. Species richness in corn fields was less 
#' than restored or remnants, which didn't differ at p=0.05. 
#' - $N_{2}$, $E_{10}$, and $E_{20}$: model fits for both null and full models were singular and NS at p<0.05.
#' 
#' Figure labels are generated and the diversity data are plotted below. An interaction plot follows, 
#' and is useful to consider what the model can and cannot say about differences in regions and field types. 
#+ div_its_otu_labs
labs_its_otu <- data.frame(
    hill_index = factor(c(rep("N0", 3), rep("N1", 3)), ordered = TRUE, levels = c("N0", "N1", "N2", "E10", "E20")),
    lab = c("a", "b", "b", "a", "b", "b"),
    xpos = rep(c(1,2,3), 2),
    ypos = rep(c(620, 170), each = 3)
)
#+ plot_div_its_otu,fig.width=9,fig.height=7
ggplot(div$its_otu, aes(x = field_type, y = value)) +
    facet_wrap(vars(hill_index), scales = "free_y") +
    geom_boxplot(varwidth = TRUE, fill = "gray90", outlier.shape = NA) +
    geom_beeswarm(aes(fill = region), shape = 21, size = 2, dodge.width = 0.2) +
    geom_text(data = labs_its_otu, aes(x = xpos, y = ypos, label = lab)) +
    labs(x = "", y = "Index value", title = "TGP microbial diversity (Hill's), ITS, 97% OTU",
         caption = "N0-richness, N1-e^Shannon, N2-Simpson, E10=N1/N0, E20=N2/N0, width=n") +
    scale_fill_discrete_qualitative(palette = "Dark3") +
    theme_bw()
#' Richness and evenness parameters increase from corn, to restored, to remnant fields, and some 
#' support exists for this pattern to occur across regions. 
#+ plot_div_its_otu_interaction,fig.width=9,fig.height=7
ggplot(
    div$its_otu %>% 
        group_by(field_type, region, hill_index) %>% 
        summarize(avg_value = mean(value), .groups = "drop"),
    aes(x = field_type, y = avg_value, group = region)) +
    facet_wrap(vars(hill_index), scales = "free_y") +
    geom_line(aes(linetype = region)) +
    geom_point(aes(fill = region), size = 2, shape = 21) +
    labs(x = "", y = "Average value", title = "Interaction plot of Hill's numbers, ITS, 97% OTU") +
    scale_fill_discrete_qualitative(palette = "Dark3") +
    theme_bw()
#' 
#' Key observations:
#' 
#' - The restored field at LP contains very high diversity, co-dominance, and evenness of fungi.
#' - The restored field at FG contains low diversity, co-dominance, and evenness. 
#' - Interactions are less an issue with $N_{0}$ and $N_{1}$
#' 
#' Next, trends in diversity are correlated with years since restoration, with 0 used for corn fields 
#' and 50 used for remnants. Statistical testing of this relationship is not valid because the ages for 
#' corn and remnant aren't justified, and the fields aren't justifiable as a chronosequence. 
#+ plot_yrs_since_resto,fig.width=9,fig.height=7
ggplot(div$its_otu, aes(x = yr_since, y = value)) +
facet_wrap(vars(hill_index), scales = "free_y") +
    geom_point(aes(fill = region, shape = field_type), size = 2) +
    labs(x = "Years since restoration", y = "index value", title = "Change in TGP microbial diversity (Hill's), ITS, 97% OTU",
         caption = "N0-richness, N1-e^Shannon, N2-Simpson, E10=N1/N0, E20=N2/N0") +
    scale_shape_manual(name = "field type", values = c(21:23)) +
    scale_fill_discrete_qualitative(name = "region", palette = "Dark3") +
    guides(fill = guide_legend(override.aes = list(shape = 21)),
           shape = guide_legend(override.aes = list(fill = NA))) +
    theme_bw()
#' 
#' Possibly, it's justified to correlate restoration age with diversity at Blue Mounds only, 
#' and with restored fields only. A Pearson's correlation is used:
#+ test_age_its_otu
test_age(div$its_otu, caption = "Correlation between Hill's numbers and field age in the Blue Mounds region")
#' 
#' Hill's $N_{1}$ decreases with age since restoration in the Blue Mounds area ($R^2$=0.60, p<0.05). 
#' This is odd and points to a confounding effect driven by difference in restoration strategy over time.
#' It's also possible that site differences (soils, etc.) also confound this relationship. It's possible
#' that we cannot attempt to present this as a time-based result at all. 
#' 
#' In any case, let's take a look at Shannon's diversity over time in Blue Mounds's restored fields.
#+ bm_test_age,fig.width=7,fig.height=6,fig.align="center"
div$its_otu %>% 
    filter(region == "BM", field_type == "restored", hill_index == "N1") %>% 
    ggplot(aes(x = yr_since, y = value)) +
    geom_text(aes(label = site_name)) +
    labs(x = "Years since restoration", y = expression("Shannon's diversity"~(N[1]))) +
    theme_classic()
#' Karla Ott's field was almost exclusively dominated by big bluestem, possibly leading to
#' a simpler microbial community. Right now, my interpretation is that restoration strategies changed 
#' over time and although restored plant communities persisted, microbial communities simplified over time.
#' Immediately after restoration, microbial diversity increased rapidly and was not sustained 
#' because the soil properties ultimately didn't change very much. 
#' 
#' Site factors (soil type) are hard to tease out, but in later analyses we will try using measured 
#' soil chemical properties. 
    


div





#' ### Fungi (ITS gene) in sequence-variant (SV) clusters, averaged to 8 samples per field.
#' Run the linear model and test differences among field types for diversity.
#+ test_div_its_otu
test_diversity(div$its_otu)
#' Key model results depend on which sites were sampled using `resample_fields()` above. 
#' Changing the value of `set.seed()` in that function may alter the results of this model. 
#' 
#' - $N_{0}$: field type is significant by likelihood ratio test at p<0.001 with region as
#' a random effect. Species richness in corn fields was less than restored or remnants, which 
#' didn't differ at p=0.05. 
#' - $N_{1}$: model fit is questionable due to singular fit, but field type is significant
#' by likelihood ratio test at p<0.01 with region as a random effect. Species richness in corn fields was less 
#' than restored or remnants, which didn't differ at p=0.05. 
#' - $N_{2}$, $E_{10}$, and $E_{20}$: model fits for both null and full models were singular and NS at p<0.05.
#' 
#' Figure labels are generated and the diversity data are plotted below. An interaction plot follows, 
#' and is useful to consider what the model can and cannot say about differences in regions and field types. 
#+ div_its_otu_labs
labs_its_otu <- data.frame(
    hill_index = factor(c(rep("N0", 3), rep("N1", 3)), ordered = TRUE, levels = c("N0", "N1", "N2", "E10", "E20")),
    lab = c("a", "b", "b", "a", "b", "b"),
    xpos = rep(c(1,2,3), 2),
    ypos = rep(c(620, 170), each = 3)
)
#+ plot_div_its_otu,fig.width=9,fig.height=7
ggplot(div$its_otu, aes(x = field_type, y = value)) +
    facet_wrap(vars(hill_index), scales = "free_y") +
    geom_boxplot(varwidth = TRUE, fill = "gray90", outlier.shape = NA) +
    geom_beeswarm(aes(fill = region), shape = 21, size = 2, dodge.width = 0.2) +
    geom_text(data = labs_its_otu, aes(x = xpos, y = ypos, label = lab)) +
    labs(x = "", y = "Index value", title = "TGP microbial diversity (Hill's), ITS, 97% OTU",
         caption = "N0-richness, N1-e^Shannon, N2-Simpson, E10=N1/N0, E20=N2/N0, width=n") +
    scale_fill_discrete_qualitative(palette = "Dark3") +
    theme_bw()
#' Richness and evenness parameters increase from corn, to restored, to remnant fields, and some 
#' support exists for this pattern to occur across regions. 
#+ plot_div_its_otu_interaction,fig.width=9,fig.height=7
ggplot(
    div$its_otu %>% 
        group_by(field_type, region, hill_index) %>% 
        summarize(avg_value = mean(value), .groups = "drop"),
    aes(x = field_type, y = avg_value, group = region)) +
    facet_wrap(vars(hill_index), scales = "free_y") +
    geom_line(aes(linetype = region)) +
    geom_point(aes(fill = region), size = 2, shape = 21) +
    labs(x = "", y = "Average value", title = "Interaction plot of Hill's numbers, ITS, 97% OTU") +
    scale_fill_discrete_qualitative(palette = "Dark3") +
    theme_bw()
#' 
#' Key observations:
#' 
#' - The restored field at LP contains very high diversity, co-dominance, and evenness of fungi.
#' - The restored field at FG contains low diversity, co-dominance, and evenness. 
#' - Interactions are less an issue with $N_{0}$ and $N_{1}$
#' 
#' Next, trends in diversity are correlated with years since restoration, with 0 used for corn fields 
#' and 50 used for remnants. Statistical testing of this relationship is not valid because the ages for 
#' corn and remnant aren't justified, and the fields aren't justifiable as a chronosequence. 
#+ plot_yrs_since_resto,fig.width=9,fig.height=7
ggplot(div$its_otu, aes(x = yr_since, y = value)) +
    facet_wrap(vars(hill_index), scales = "free_y") +
    geom_point(aes(fill = region, shape = field_type), size = 2) +
    labs(x = "Years since restoration", y = "index value", title = "Change in TGP microbial diversity (Hill's), ITS, 97% OTU",
         caption = "N0-richness, N1-e^Shannon, N2-Simpson, E10=N1/N0, E20=N2/N0") +
    scale_shape_manual(name = "field type", values = c(21:23)) +
    scale_fill_discrete_qualitative(name = "region", palette = "Dark3") +
    guides(fill = guide_legend(override.aes = list(shape = 21)),
           shape = guide_legend(override.aes = list(fill = NA))) +
    theme_bw()
#' 
#' Possibly, it's justified to correlate restoration age with diversity at Blue Mounds only, 
#' and with restored fields only. A Pearson's correlation is used:
#+ test_age_its_otu
test_age(div$its_otu, caption = "Correlation between Hill's numbers and field age in the Blue Mounds region")
#' 
#' Hill's $N_{1}$ decreases with age since restoration in the Blue Mounds area ($R^2$=0.60, p<0.05). 
#' This is odd and points to a confounding effect driven by difference in restoration strategy over time.
#' It's also possible that site differences (soils, etc.) also confound this relationship. It's possible
#' that we cannot attempt to present this as a time-based result at all. 
#' 
#' In any case, let's take a look at Shannon's diversity over time in Blue Mounds's restored fields.
#+ bm_test_age,fig.width=7,fig.height=6,fig.align="center"
div$its_otu %>% 
    filter(region == "BM", field_type == "restored", hill_index == "N1") %>% 
    ggplot(aes(x = yr_since, y = value)) +
    geom_text(aes(label = site_name)) +
    labs(x = "Years since restoration", y = expression("Shannon's diversity"~(N[1]))) +
    theme_classic()