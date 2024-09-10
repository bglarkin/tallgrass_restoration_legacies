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
#' Microbial data analyzed here include site-species tables derived from high-throughput sequencing of 
#' ITS and 18S genes and clustering into 97% similar OTUs.
#' This report presents basic statistics and visualizations of species richness, Shannon's 
#' diversity/evenness, and Simpson's diversity/evenness in the microbial species data across field types. 
#' 
#' - Diversity and evenness of microbial communities
#' - Interpretation of differences in diversity among regions and field types, and over years
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
#' List *spe* holds average sequence abundances for the top 6 samples per field. 
#' CSV files were produced in `process_data.R`
spe <- list(
    its_raw = read_csv(paste0(getwd(), "/clean_data/spe_ITS_raw.csv"), 
                        show_col_types = FALSE),
    its_rfy = read_csv(paste0(getwd(), "/clean_data/spe_ITS_rfy.csv"), 
                       show_col_types = FALSE),
    amf_raw = read_csv(paste0(getwd(), "/clean_data/spe_18S_raw.csv"), 
                        show_col_types = FALSE),
    amf_rfy = read_csv(paste0(getwd(), "/clean_data/spe_18S_rfy.csv"), 
                       show_col_types = FALSE)
)
#' 
#' ## Site metadata and design
sites <- read_csv(paste0(getwd(), "/clean_data/sites.csv"), show_col_types = FALSE) %>% 
    mutate(field_type = factor(field_type, ordered = TRUE, levels = c("corn", "restored", "remnant"))) %>% 
    select(-lat, -long, -yr_restore, -yr_rank)
#' 
#' # Functions
#' The following functions are used to streamline code and reduce errors:
#' 
#' ## Calculate Hill's series on a samples-species matrix
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
            rownames_to_column(var = "field_key") %>%
            mutate(field_key = as.integer(field_key)) %>%
            left_join(sites, by = join_by(field_key)) %>%
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
#' ## Test diversity measures across site types with mixed model
#+ test_diversity_function
test_diversity <- function(data) {
    hills <- levels(data$hill_index)
    for(i in 1:length(hills)) {
        cat("---------------------------------\n")
        print(hills[i])
        cat("---------------------------------\n\n")
        mod_data <- data %>% 
            filter(hill_index == hills[i]) %>% 
            mutate(field_type = factor(field_type, ordered = FALSE))
        mmod <- lmer(value ~ field_type + (1 | region), data = mod_data, REML = FALSE)
        print(mmod)
        cat("\n---------------------------------\n\n")
        mmod_null <- lmer(value ~ 1 + (1 | region), data = mod_data, REML = FALSE)
        print(mmod_null)
        cat("\n---------------------------------\n\n")
        print(anova(mmod, mmod_null))
        mod_tuk <- glht(mmod, linfct = mcp(field_type = "Tukey"), test = adjusted("holm"))
        print(mod_tuk)
        print(cld(mod_tuk))
        cat("\n\n\n")
    }
}
#' ## Change in diversity over time 
#' Do Hill's numbers correlate with years since restoration?
#' This may only be appropriate to attempt in the Blue Mounds region, and even there, it may be difficult
#' to justify that the area meets the criteria for a chronosequence. 
test_age <- function(data, caption=NULL) {
    temp_df <-
        data %>%
        filter(field_type == "restored", region == "BM") %>%
        pivot_wider(names_from = hill_index, values_from = value) %>%
        select(-starts_with("field"),-region)
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
#' # Analysis and Results
#' 
#' ## OTUs pre/post corrections
#' 
#' Datasets were corrected for survey effort (min samples per field) and sequencing depth 
#' (min sequences per field). What was the effect of these actions on the number of OTUs recovered?
#' After rarefying, zero-abundance OTUs were removed.  
#' Few were lost due to rarefying, as we can see by counting columns (less column 1 because it 
#' has field site keys):
Map(function(x) ncol(x)-1, spe) 
#' It appears that little will be lost in terms of richness or diversity by rarefying.
#' 
#' ## Microbial diversity
#' 
#' Microbial diversity is considered for 18S or ITS gene 
#' datasets. For each set, Hill's numbers are produced ([Hill 1973](http://doi.wiley.com/10.2307/1934352), 
#' [Borcard and Legendere 2018, p. 373](http://link.springer.com/10.1007/978-3-319-71404-2)) and plotted,
#' with means differences tested using mixed-effects linear models in `lmer` ([Bates et al. 2015](https://doi.org/10.18637/jss.v067.i01)).
#' Correlations are then produced to visualize change in diversity trends over time, 
#' with similar mixed-effects tests performed. 
#' 
#' The statistical tests are on shaky ground due to imbalance, but are presented here as an 
#' attempt to at least explore some differences and think more later about how we could 
#' present them in a more valid way. 
#' 
#' Hill's numbers, brief description:
#' 
#' - $N_{0}$  = species richness
#' - $N_{1}$  = Shannon's diversity ($e^H$; excludes rarest species, considers the number of "functional" species)
#' - $N_{2}$  = Simpson's diversity ($1 / \lambda$; number of "codominant" species)
#' - $E_{10}$ = Shannon's evenness (Hill's ratio $N_{1} / N_{0}$)
#' - $E_{20}$ = Simpson's evenness (Hill's ratio $N_{2} / N_{0}$)
#' 
#+ diversity_calculations
# Rarefied tables only
div <- Map(calc_diversity, spe[c(2,4)])
#' 
#' ### Fungi (ITS gene)
#' #### Diversity across field types
#' Run the linear model and test differences among field types for diversity.
#+ test_div_its_otu
test_diversity(div$its_rfy)
#' 
#' #### Result: ITS diversity
#' 
#' - $N_{0}$: field type is significant by likelihood ratio test at p<0.001 with region as
#' a random effect. Species richness in corn fields was less than restored or remnants, which 
#' didn't differ at p=0.05. 
#' - $N_{1}$: field type is significant
#' by likelihood ratio test at p<0.05 with region as a random effect. Shannon's diversity in corn fields was less 
#' than restored or remnants, which didn't differ at p=0.05. 
#' - $N_{2}$, $E_{10}$, and $E_{20}$: model fits for both null and full models were singular and NS at p<0.05.
#' 
#' Figure labels are generated and the diversity data are plotted below. An interaction plot follows, 
#' and is useful to consider what the model can and cannot say about differences in regions and field types. 
#+ div_its_otu_labs
labs_its <- data.frame(
    hill_index = factor(c(rep("N0", 3), rep("N1", 3)), 
                        ordered = TRUE, 
                        levels = c("N0", "N1", "N2", "E10", "E20")),
    lab = c("a", "b", "b", "a", "b", "b"),
    xpos = rep(c(1,2,3), 2),
    ypos = rep(c(590, 170), each = 3)
)
#+ plot_div_its_rfy,fig.width=9,fig.height=7,fig.align='center'
ggplot(div$its_rfy, aes(x = field_type, y = value)) +
    facet_wrap(vars(hill_index), scales = "free_y") +
    geom_boxplot(varwidth = TRUE, fill = "gray90", outlier.shape = NA) +
    geom_beeswarm(aes(fill = region), shape = 21, size = 2, dodge.width = 0.2) +
    geom_label(data = labs_its, aes(x = xpos, y = ypos, label = lab), label.size = NA) +
    labs(x = "", y = "Index value", title = "TGP microbial diversity (Hill's), ITS, 97% OTU",
         caption = "N0-richness, N1-e^Shannon, N2-Simpson, E10=N1/N0, E20=N2/N0, width=n,\nletters indicate significant differences at p<0.05") +
    scale_fill_discrete_qualitative(palette = "Dark3") +
    theme_bw()
#' Let's also make a figure to show just the diversity measures, and with better labels. 
hill_labs <- c("Richness", "Shannon's", "Simpson's")
names(hill_labs) <- c("N0", "N1", "N2")
#+ div_its_box,fig.width=7,fig.height=3
div$its_rfy %>% 
    filter(hill_index %in% c("N0", "N1", "N2")) %>% 
    ggplot(aes(x = field_type, y = value)) +
    facet_wrap(vars(hill_index), scales = "free_y", labeller = labeller(hill_index = hill_labs)) +
    geom_boxplot(varwidth = TRUE, fill = "gray90", outlier.shape = NA) +
    geom_beeswarm(aes(shape = region, fill = field_type), size = 2, dodge.width = 0.3) +
    geom_label(data = labs_its, aes(x = xpos, y = ypos, label = lab), label.size = NA) +
    labs(y = "Species equivalent (n)") +
    scale_fill_discrete_qualitative(name = "Field Type", palette = "Harmonic") +
    scale_shape_manual(name = "Region", values = c(21, 22, 23, 24)) +
    theme_bw() +
    theme(axis.title.x = element_blank()) +
    guides(fill = guide_legend(override.aes = list(shape = 21)))

#' Richness and evenness parameters increase from corn, to restored, to remnant fields, and some 
#' support exists for this pattern to occur across regions. 
#+ plot_div_its_otu_interaction,fig.width=9,fig.height=7,fig.align='center'
ggplot(
    div$its_rfy %>% 
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
#' #### Key observations:
#' 
#' - The restored field at LP contains very high diversity, co-dominance, and evenness of fungi.
#' - The restored field at FG contains low diversity, co-dominance, and evenness. 
#' - Interactions are less an issue with $N_{0}$ and $N_{1}$
#' 
#' #### Diversity over time
#' Next, trends in diversity are correlated with years since restoration. This can only be attempted
#' with Fermi and Blue Mounds sites; elsewhere, blocks cannot be statistically accounted for because
#' treatments aren't replicated within them. 
#+ plot_yrs_since_resto_FLBM,fig.width=9,fig.height=7,fig.align='center'
div$its_rfy %>% 
    filter(field_type == "restored", 
           hill_index %in% c("N0", "N1", "N2"),
           region %in% c("BM", "FL")) %>% 
    ggplot(aes(x = yr_since, y = value)) +
    facet_grid(cols = vars(region), rows = vars(hill_index), scales = "free_y") +
    geom_point(shape = 21, fill = "gray50", size = 2) +
    labs(x = "Years since restoration", y = "index value", title = "TGP microbial diversity (Hill's) over time, ITS, 97% OTU",
         caption = "N0-richness, N1-e^Shannon, N2-Simpson") +
    theme_bw()
#' If the relationship exists, it is in Blue Mounds only. 
#' 
#' It's probably justified to correlate diversity with field age in Blue Mounds's restored fields.
#' A Pearson's correlation is used:
#+ test_age_its_otu
test_age(div$its_rfy, 
         caption = "Correlation between Hill's numbers and field age in the Blue Mounds region: ITS, 97% OTU")
#' 
#' Hill's $N_{1}$ decreases with age since restoration in the Blue Mounds area ($R$=-0.84, p<0.05). 
#' It's driven primarily by Karla Ott's older restoration which differs in plant community. 
#' Removing this field, the correlation is no longer significant, and it not even close.   
#' This is odd and points to a confounding effect driven by difference in restoration strategy over time.
#' Karla Ott's field is heavily dominated by a C4 grass and has relatively poor plant species diversity.
#' It's possible that other site differences (soils, etc.) also confound this relationship. It's possible
#' that we cannot attempt to present this as a time-based result at all. Or, maybe the number of
#' functionally dominant species slowly declines over time due to lack of disturbance and substrate 
#' diversity.
#' 
#' That diversity metrics aren't changing over time, or are possibly declining over time
#' after restoration is concerning and worth mentioning. 
#' 
#' In any case, let's take a look at Shannon's diversity over time in Blue Mounds' restored fields.
#+ bm_test_age,message=FALSE,fig.width=7,fig.height=6,fig.align='center'
div$its_rfy %>% 
    filter(region == "BM", field_type == "restored", hill_index == "N1") %>% 
    ggplot(aes(x = yr_since, y = value)) +
    geom_smooth(method = "lm", se = TRUE, linewidth = 0.8, fill = "gray70") +
    geom_label(aes(label = field_name)) +
    labs(x = "Years since restoration", y = expression("Shannon's diversity"~(N[1]))) +
    theme_classic()
#' Karla Ott's field was almost exclusively dominated by big bluestem, possibly leading to
#' a simpler microbial community. That field has too much leverage on this plot.  
#' Right now, my interpretation is that restoration strategies changed 
#' over time and although restored plant communities persisted, microbial communities simplified over time.
#' Immediately after restoration, microbial diversity increased rapidly and was not sustained 
#' because the soil properties ultimately didn't change very much. 
#' 
#' Site factors (soil type) are hard to tease out, but in later analyses we will try using measured 
#' soil chemical properties. 
#'  
#'  
#' ### AMF (18S gene)
#' Run the linear model and test differences among field types for diversity.
#+ test_div_amf_otu
test_diversity(div$amf_rfy)
#' 
#' #### Result: AMF diversity
#' Despite apparent trends across field types, variances are large and interactions apparent. All
#' model fits are questionable due to [singularity](https://rdrr.io/cran/lme4/man/isSingular.html). 
#' The following results and plots are provisional and included here for consideration only.  
#' 
#' - $N_{0}$: Species richness differences by field type are significant by likelihood ratio test at p<0.005 with region as
#' a random effect. Species richness in corn fields was less than restored, and neither differed
#' from remnants at p=0.05. 
#' - $N_{1}$: field type is significant by likelihood ratio test at p<0.001 with region as a random effect. 
#' Shannon's diversity in corn fields was less than restored or remnants, which didn't differ at p=0.05. 
#' - $N_{2}$: field type is significant by likelihood ratio test at p<0.001 with region as a random effect.
#' Restored and remnant fields host a larger group of co-dominant AMF than are found in cornfields 
#' (p=0.05).
#' - $E_{10}$: field type is significant by likelihood ratio test at p<0.05 with region as a random effect. 
#' Weak support for higher evenness of functionally abundant species in remnant fields was found 
#' (p<0.05).
#' - $E_{20}$: Similar trend as $E_{10}$ but NS (p>0.05). 
#' 
#' Figure labels are generated and the diversity data are plotted below. An interaction plot follows, 
#' and is useful to consider what the model can and cannot say about differences in regions and field types. 
#+ div_amf_otu_labs
labs_amf <- data.frame(
    hill_index = factor(c(rep("N0", 3), rep("N1", 3), rep("N2", 3), rep("E10", 3)), 
                        ordered = TRUE, 
                        levels = c("N0", "N1", "N2", "E10", "E20")),
    lab = c("a", "b", "b", "a", "b", "b", "a", "b", "b", "a", "a", "b"),
    xpos = rep(c(1,2,3), 4),
    ypos = rep(c(66, 35, 26, 0.59), each = 3)
)
#+ plot_div_amf_otu,fig.width=9,fig.height=7,fig.align='center'
ggplot(div$amf_rfy, aes(x = field_type, y = value)) +
    facet_wrap(vars(hill_index), scales = "free_y") +
    geom_boxplot(varwidth = TRUE, fill = "gray90", outlier.shape = NA) +
    geom_beeswarm(aes(fill = region), shape = 21, size = 2, dodge.width = 0.2) +
    geom_label(data = labs_amf, aes(x = xpos, y = ypos, label = lab), label.size = NA) +
    labs(x = "", y = "Index value", title = "TGP microbial diversity (Hill's), AMF (18S), 97% OTU",
         caption = "N0-richness, N1-e^Shannon, N2-Simpson, E10=N1/N0, E20=N2/N0, width=n,\nletters indicate significant differences at p<0.05") +
    scale_fill_discrete_qualitative(palette = "Dark3") +
    theme_bw()
#' Let's also make a figure to show just the diversity measures, and with better labels. 
#+ div_amf_box,fig.width=7,fig.height=3
div$amf_rfy %>% 
    filter(hill_index %in% c("N0", "N1", "N2")) %>% 
    ggplot(aes(x = field_type, y = value)) +
    facet_wrap(vars(hill_index), scales = "free_y", labeller = labeller(hill_index = hill_labs)) +
    geom_boxplot(varwidth = TRUE, fill = "gray90", outlier.shape = NA) +
    geom_beeswarm(aes(shape = region, fill = field_type), size = 2, dodge.width = 0.3) +
    geom_label(data = labs_amf %>% filter(hill_index != "E10"), aes(x = xpos, y = ypos, label = lab), label.size = NA) +
    labs(y = "Species equivalent (n)") +
    scale_fill_discrete_qualitative(name = "Field Type", palette = "Harmonic") +
    scale_shape_manual(name = "Region", values = c(21, 22, 23, 24)) +
    theme_bw() +
    theme(axis.title.x = element_blank()) +
    guides(fill = guide_legend(override.aes = list(shape = 21)))
#' Richness increases from corn, to restored and remnant fields, and some 
#' support exists for this pattern to occur across regions. The trend is weakest with $N_{0}$, suggesting
#' that both restored and remnant soils contain more functionally abundant and co-dominant species than are found in cornfields, 
#' but some cornfields have "long tails" of rare species. Wide variances stifle inferences. The 
#' trend detected in evenness suggests that a few weedy AMF species dominate cornfields but most restored 
#' fields host more balanced communities that are more similar to remnants. 
#+ plot_div_amf_otu_interaction,fig.width=9,fig.height=7,fig.align='center'
ggplot(
    div$amf_rfy %>% 
        group_by(field_type, region, hill_index) %>% 
        summarize(avg_value = mean(value), .groups = "drop"),
    aes(x = field_type, y = avg_value, group = region)) +
    facet_wrap(vars(hill_index), scales = "free_y") +
    geom_line(aes(linetype = region)) +
    geom_point(aes(fill = region), size = 2, shape = 21) +
    labs(x = "", y = "Average value", title = "Interaction plot of Hill's numbers, AMF (18S), 97% OTU") +
    scale_fill_discrete_qualitative(palette = "Dark3") +
    theme_bw()
#' 
#' #### Key observations:
#' 
#' - Cornfields differ from restored and remnant fields, with lower richness, fewer dominant species,
#' and greater dominance of those species. 
#' - Differences between restored and remnant fields change directions based on region, with FL
#' and LP matching the hypothesized pattern but BM and FG reversing it. 
#' - Particular species may be strong interactors here. 
#' 
#' #### Diversity over time at Blue Mounds (AMF)
#' It's probably justified to correlate diversity with field age in Blue Mounds's restored fields.
#' A Pearson's correlation is used:
#+ test_age_amf_otu
test_age(div$amf_rfy, caption = "Correlation between Hill's numbers and field age in the Blue Mounds region: AMF (18S), 97% OTU")
#' 
#' The relationships are too weak to examine further. 
#' 

