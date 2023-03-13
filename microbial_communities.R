#' ---
#' title: "Microbial data: community differences"
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
#' Microbial data include site-species tables derived from high-throughput sequencing and 
#' clustering in QIIME by Lorinda Bullington and PLFA/NLFA data which Ylva Lekberg did. 
#' 
#' This presents basic visualizations of community differences among sites/regions
#' based on ITS data. 
#' 
#' One goal here is to see whether choosing OTU or SV clusters presents qualitatively different
#' outcomes in ordinations. We will choose one (OTUs) if they start to look similar, as they 
#' have so far.
#' 
#' Species distance matrices are resampled to the minimum number which successfully amplified per
#' field. This was done to equalize sampling effort. This procedure can easily be undone in 
#' the [process_data script](process_data.md)
#' 
#' # Packages and libraries
packages_needed = c("tidyverse", "vegan", "colorspace", "ape", "knitr")
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
#' # Data
#' ## Sites-species tables
#' Samples-species tables with rarefied sequence abundances. 
#' CSV files were produced in [process_data.R](process_data.md)
spe <- list(
    its = read_csv(
        paste0(getwd(), "/clean_data/spe_ITS_rfy.csv"),
        show_col_types = FALSE
    ),
    amf = read_csv(
        paste0(getwd(), "/clean_data/spe_18S_rfy.csv"),
        show_col_types = FALSE
    )
)
#' ## Species metadata
#' Needed to make inset figures showing most important categories of species. The OTUs
#' and sequence abundances in these files matches the rarefied data in `spe$` above.
#' CSV files were produced in the [microbial diversity script](microbial_diversity.md).
spe_meta <- list(
    its =
        read_csv(
            paste0(getwd(), "/clean_data/speTaxa_ITS_rfy.csv"),
            show_col_types = FALSE
        ),
    amf = 
        read_csv(
            paste0(getwd(), "/clean_data/speTaxa_18S_rfy.csv"),
            show_col_types = FALSE
        )
)
#' 
#' ## Site metadata
#' Needed for figure interpretation and permanova designs.
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
#' ## Distance tables
#' Creating distance objects from the samples-species tables is done with the typical 
#' process of `vegdist()` in vegan.
#' Bray-Curtis or Ruzicka distance are both appropriate methods for these data, but Bray-Curtis has 
#' produced axes with better explanatory power (Ruzicka is used with method="jaccard")
#' With the 18S data, we can take advantage of phylogenetic relationships in a UNIFRAC distance
#' matrix. The UNIFRAC distance was produced in QIIME II and needs some wrangling to 
#' conform to the standards of a distance object in R. The following list contains vegdist-produced
#' distance objects for ITS and 18S, and it includes UNIFRAC distance for 18S. 
#' 
distab <- list(
    its = vegdist(data.frame(spe$its, row.names = 1), method = "bray"),
    amf_bray = vegdist(data.frame(spe$amf, row.names = 1), method = "bray"),
    amf_uni = sites %>%
        select(field_name, field_key) %>%
        left_join(read_delim(
            paste0(getwd(), "/otu_tables/18S/18S_weighted_Unifrac.tsv"),
            show_col_types = FALSE
        ),
        by = join_by(field_name)) %>%
        select(field_key, everything(),-field_name) %>%
        data.frame(row.names = 1) %>%
        as.dist()
) 
#' 
#' ## Functions
#' A function handles the Principal Components Analysis (PCoA) diagnostics, with outputs and figures 
#' saved to a list for later use. 
#+ pcoa_function
pcoa_fun <- function(d, env=sites, corr="none", df_name, nperm=1999) {
    # Multivariate analysis
    p <- pcoa(d, correction = corr)
    p_vals <- data.frame(p$values) %>% 
        rownames_to_column(var = "Dim") %>% 
        mutate(Dim = as.integer(Dim))
    p_vec <- data.frame(p$vectors)
    # Permutation tests (PERMANOVA)
    p_permtest <-
        with(env,
             adonis2(
                 d ~ field_type,
                 data = env,
                 permutations = nperm,
                 strata = region
             ))
    # Diagnostic plots
    if(corr == "none" | ncol(p_vals) == 6) {
        p_bstick <- ggplot(p_vals, aes(x = factor(Dim), y = Relative_eig)) + 
            geom_col(fill = "gray70", color = "gray30") + 
            geom_line(aes(x = Dim, y = Broken_stick), color = "red") +
            geom_point(aes(x = Dim, y = Broken_stick), color = "red") +
            labs(x = "Dimension", 
                 title = paste0("PCoA Eigenvalues and Broken Stick Model (", df_name, ")")) +
            theme_bw()
        p_ncomp <- with(p_vals, which(Relative_eig < Broken_stick)[1]-1)
        eig <- round(p_vals$Relative_eig[1:2] * 100, 1)
    } else {
        p_bstick <- ggplot(p_vals, aes(x = factor(Dim), y = Rel_corr_eig)) + 
            geom_col(fill = "gray70", color = "gray30") + 
            geom_line(aes(x = Dim, y = Broken_stick), color = "red") +
            geom_point(aes(x = Dim, y = Broken_stick), color = "red") +
            labs(x = "Dimension", 
                 title = paste0("PCoA Eigenvalues and Broken Stick Model (", df_name, ")")) +
            theme_bw()
        p_ncomp <- with(p_vals, which(Rel_corr_eig < Broken_stick)[1]-1)
        eig <- round(p_vals$Rel_corr_eig[1:2] * 100, 1)
    }
    ncomp <- if(p_ncomp <= 2) {2} else {p_ncomp}
    # Ordination plot
    scores <-
        p_vec[, 1:ncomp] %>%
        rownames_to_column(var = "field_key") %>%
        mutate(field_key = as.integer(field_key)) %>%
        left_join(sites, by = "field_key") %>% 
        select(-field_name)
    # Output data
    output <- list(dataset = df_name,
                   components_exceed_broken_stick = p_ncomp,
                   correction_note = p$note,
                   values = p_vals[1:(ncomp+1), ], 
                   eigenvalues = eig,
                   site_vectors = scores,
                   broken_stick_plot = p_bstick,
                   permanova = p_permtest)
    return(output)
}
#' # Results
#' ## Ordinations
#' Bray-Curtis or Ruzicka distance are both appropriate, but Bray-Curtis has 
#' produced axes with better explanatory power (Ruzicka is used with method="jaccard")
#' 
#' In trial runs, no negative eigenvalues were observed (not shown). No 
#' correction is needed for these ordinations.
#' 
#' ### PCoA with ITS gene, OTU clusters
#+ pcoa_its_otu,fig.align='center'
(pcoa_its <- pcoa_fun(distab$its, df_name = "ITS gene, 97% OTU"))
#' 
#' Axis 1 explains `r pcoa_its$eigenvalues[1]`% of the variation and is the only eigenvalue that exceeds a 
#' broken stick model. The most substantial variation here will be on the first axis,
#' although axis 2 explains `r pcoa_its$eigenvalues[2]`% of the variation and was very close to the broken
#' stick value. Testing the design factor *field_type* (with *region* treated as a block 
#' using the `strata` argument of `adonis2`) revealed a significant
#' clustering $(R^2=`r round(pcoa_its$permanova$R2[1], 2)`, p=`r pcoa_its$permanova$Pr[1]`)$.  
#' 
#' Let's view a plot with abundances of community subgroups inset.
pcoa_its$ord <-
    ggplot(pcoa_its$site_vectors, aes(x = Axis.1, y = Axis.2)) +
    geom_point(aes(fill = field_type, shape = region), size = 10) +
    geom_text(aes(label = yr_since)) +
    scale_fill_discrete_qualitative(palette = "harmonic") +
    scale_shape_manual(values = c(21, 22, 23, 24)) +
    labs(
        x = paste0("Axis 1 (", pcoa_its$eig[1], "%)"),
        y = paste0("Axis 2 (", pcoa_its$eig[2], "%)"),
        title = paste0(
            "PCoA Ordination of field-averaged species data (",
            pcoa_its$dataset,
            ")"
        ),
        caption = "Text indicates years since restoration, with corn (-) and remnants (+) never restored."
    ) +
    lims(y = c(-0.5,0.32)) +
    theme_bw() +
    guides(fill = guide_legend(override.aes = list(shape = 21)))
pcoa_its$inset <-
    spe_meta$its %>%
    filter(primary_lifestyle %in% c("soil_saprotroph", "wood_saprotroph", "plant_pathogen")) %>%
    mutate(field_type = factor(
        field_type,
        ordered = TRUE,
        levels = c("corn", "restored", "remnant")
    )) %>%
    group_by(primary_lifestyle, field_type) %>%
    summarize(avg_seq_abund = mean(seq_abund), .groups = "drop") %>%
    ggplot(aes(x = primary_lifestyle, y = avg_seq_abund)) +
    geom_col(aes(fill = field_type)) +
    labs(x = "",
         y = "Seq. abund. (avg)") +
    scale_fill_discrete_qualitative(palette = "Harmonic") +
    scale_x_discrete(label = c("soil sapr", "wood sapr", "plnt path")) +
    # coord_flip() +
    theme_classic() +
    theme(legend.position = "none")
#+ its_guilds_fig,fig.align='center'
pcoa_its$ord +
    annotation_custom(
        ggplotGrob(pcoa_its$inset + theme(
            plot.background = element_rect(colour = "black", fill = "gray90")
        )),
        xmin = -0.38,
        xmax = 0,
        ymin = -0.53,
        ymax = -0.18
    )
#' 
#' Community trajectories revealed in the ordination clearly depend on both region and field type.
#' Faville Grove shows a linear progression
#' from corn to remnant and Lake Petite does as well, although with few sites and only 
#' single restoration ages these are weak supports. With Blue Mounds sites, the general 
#' progression along Axis 1 is to increase in age from left to right, but the remnant 
#' doesn't seem representative because it clusters far from everything else and associates
#' most strongly with the neighboring restored field (both on Merel Black's property). Restored fields 
#' at Fermi separate well away from cornfields, but less age structure is found. Instead, the 
#' old restorations in the ring most resemble the Railroad Remnant (which is in a different soil...), 
#' the switchgrass restored fields take a potentially novel path toward distant remnants. 
#' 
#' On axis 1, four clusters are apparent in at least two partitioning schemes. 
#' It will be interesting to see if we can pull those apart with explanatory variables. 
#' 
#' Restoration age will be explored in-depth with the subset of restoration fields. 
#' 
#' Here we can also begin considering what an inset plot to display metadata might look like. 
#' Let's plot and test the relationship between age and community axis scores with restored fields 
#' only.
#+ its_resto_scores_data
its_resto_scores <-
    pcoa_its$site_vectors %>%
    filter(field_type == "restored") %>%
    mutate(yr_since = as.numeric(yr_since))
#+ its_resto_scores_lm
summary(lm(Axis.1 ~ yr_since,
           data = its_resto_scores))
#+ its_resto_scores_fig,fig.align='center',message=FALSE
its_resto_scores %>%
    pivot_longer(Axis.1:Axis.2, names_to = "axis", values_to = "score") %>%
    ggplot(aes(x = yr_since, y = score)) +
    facet_wrap(vars(axis), scales = "free") +
    geom_smooth(aes(linetype = axis), method = "lm", se = FALSE, linewidth = 0.5) +
    geom_point(aes(shape = region), fill = "grey", size = 2) +
    labs(x = "Years since restoration",
         y = "PCoA axis score",
         title = "Correlations, axis scores and years since restoration (ITS, 97% OTU)",
         caption = "Blue lines show linear model fit; solid line is significant at p<0.05") +
    scale_shape_manual(values = c(21, 22, 23, 24)) +
    scale_linetype_manual(values = c('solid', 'dashed'), guide = "none") +
    theme_bw()
#' 
#' Indeed, Axis 1 does correlate well with age $(R^2_{Adj}=0.51, p<0.005)$.
#' 
#' It's probably better to do with with a new ordination of just restoration sites
#' and a constrained ordination with years and other environmental variables.
#' 
#' ### PCoA with 18S gene, uninformed distance
#' "Uninformed" refers to Bray-Curtis or another distance informational metric available in `vegdist()`.
#' It is as opposed to UNIFRAC distance, which is informed by phylogeny.
#+ pcoa_amf_otu,fig.align='center'
(pcoa_amf_bray <- pcoa_fun(distab$amf_bray, df_name = "18S gene, 97% OTU, Bray-Curtis distance"))
#' 
#' Four axes are significant by a broken stick model, between them explaining
#' `r round(sum(pcoa_amf_bray$values$Relative_eig[1:4])*100, 1)`% of the 
#' variation in AMF among fields. It may be worthwhile to examine structure on Axes 3 and 4 
#' sometime. The most substantial variation here is on the first axis (`r pcoa_amf_bray$eigenvalues[1]`%) with Axis 2
#' explaining `r pcoa_amf_bray$eigenvalues[2]`% of the variation in AMF abundances. 
#' Testing the design factor *field_type* (with *region* treated as a block 
#' using the `strata` argument of `adonis2`) revealed a significant
#' clustering $(R^2=`r round(pcoa_amf_bray$permanova$R2[1], 2)`, p=`r round(pcoa_amf_bray$permanova$Pr[1], 3)`)$. 
#' 
#' Let's view a plot with abundances of community subgroups inset. 
pcoa_amf_bray$ord <- 
    ggplot(pcoa_amf_bray$site_vectors, aes(x = Axis.1, y = Axis.2)) +
    geom_point(aes(fill = field_type, shape = region), size = 10) +
    geom_text(aes(label = yr_since)) +
    scale_fill_discrete_qualitative(palette = "harmonic") +
    scale_shape_manual(values = c(21, 22, 23, 24)) +
    labs(x = paste0("Axis 1 (", pcoa_amf_bray$eig[1], "%)"), 
         y = paste0("Axis 2 (", pcoa_amf_bray$eig[2], "%)"), 
         title = paste0("PCoA Ordination of field-averaged species data (", pcoa_amf_bray$dataset, ")"),
         caption = "Text indicates years since restoration, with corn (-) and remnants (+) never restored.") +
    lims(x = c(-0.6,0.35)) +
    theme_bw() +
    guides(fill = guide_legend(override.aes = list(shape = 21)))
pcoa_amf_bray$inset <- 
    spe_meta$amf %>% 
    filter(family %in% c("Claroideoglomeraceae", "Paraglomeraceae", "Diversisporaceae", "Gigasporaceae")) %>% 
    mutate(field_type = factor(field_type, ordered = TRUE, levels = c("corn", "restored", "remnant"))) %>% 
    group_by(family, field_type) %>% 
    summarize(avg_seq_abund = mean(seq_abund), .groups = "drop") %>% 
    ggplot(aes(x = family, y = avg_seq_abund)) +
    geom_col(aes(fill = field_type)) +
    labs(x = "",
         y = "Seq. abundance (avg)") +
    scale_fill_discrete_qualitative(palette = "Harmonic") +
    scale_x_discrete(label = function(x) abbreviate(x, minlength = 6)) +
    # coord_flip() +
    theme_classic() +
    theme(legend.position = "none")
#+ amf_families_fig,fig.align='center'
pcoa_amf_bray$ord +
    annotation_custom(
        ggplotGrob(pcoa_amf_bray$inset + theme(
            plot.background = element_rect(colour = "black", fill = "gray90")
        )),
        xmin = -0.63,
        xmax = -0.2,
        ymin = -0.32,
        ymax = -0.04
    )
#'
#' Community trajectories revealed in the ordination correlate with field type. 
#' Corn fields stand well apart with AMF communities,
#' with restored and remnant fields clustering closer than we had seen with ITS-identified fungi. 
#' Restoration age along Axis 1 follows a near-linear progression in Blue Mounds fields; with Fermi,
#' we see a weaker age progression and instead a strong separation between "ring fields" and switchgrass
#' plots as before. Restored fields' fidelity to remnants seems stronger with AMF than we 
#' had seen with general fungi. 
#' 
#' What's becoming apparent here is that Axis 1 separates strongly on *field_type* and years 
#' since restoration, and Axis 2 further separates on years since restoration. A consistent signal
#' of region isn't obvious. 
#' 
#' Let's test the relationship between age and community axis scores with restored fields 
#' only.
#+ amf_yrs_scores_data
amf_resto_scores <-
    pcoa_amf_bray$site_vectors %>%
    filter(field_type == "restored") %>%
    mutate(yr_since = as.numeric(yr_since))
#+ amf_yrs_scores_lm_1
summary(lm(Axis.1 ~ yr_since,
           data = amf_resto_scores))
#+ amf_yrs_scores_lm_2
summary(lm(Axis.2 ~ yr_since,
           data = amf_resto_scores))
#+ amf_yrs_scores_fig,fig.align='center',message=FALSE
amf_resto_scores %>%
    pivot_longer(Axis.1:Axis.2, names_to = "axis", values_to = "score") %>%
    ggplot(aes(x = yr_since, y = score)) +
    facet_wrap(vars(axis), scales = "free") +
    geom_smooth(method = "lm", se = FALSE, linewidth = 0.5) +
    geom_point(aes(shape = region), fill = "grey", size = 2) +
    labs(x = "Years since restoration",
         y = "PCoA axis score",
         title = "Correlations, axis scores and years since restoration (18S gene, 97% OTU, Bray-Curtis distance)",
         caption = "Blue lines show linear model fit; solid line is significant at p<0.05") +
    scale_shape_manual(values = c(21, 22, 23, 24)) +
    theme_bw()
#' 
#' Both axes correlate significantly and strongly with years since restoration. 
#' Axis 2 shows a stronger relationship $(R^2_{Adj}=0.61, p<0.001)$, and Axis 1 
#' shows a moderately strong relationship $(R^2_{Adj}=0.46, p<0.005)$ 
#' 
#' 
#' ### PCoA with 18S gene, UNIFRAC distance
#' 
#+ pcoa_amf_uni,fig.align='center'
(pcoa_amf_uni <- pcoa_fun(distab$amf_uni, df_name = "18S gene, 97% OTU, UNIFRAC distance", corr = "lingoes"))
#' 
#' Three axes are significant by a broken stick model, between them explaining
#' `r round(sum(pcoa_amf_uni$values$Rel_corr_eig[1:3])*100, 1)`% of the 
#' variation in AMF among fields. The most substantial variation here is on the first axis 
#' (`r pcoa_amf_uni$eigenvalues[1]`%) with Axis 2
#' explaining `r pcoa_amf_uni$eigenvalues[2]`% of the variation in AMF abundances. 
#' Testing the design factor *field_type* (with *region* treated as a block 
#' using the `strata` argument of `adonis2`) revealed a significant
#' clustering $(R^2=`r round(pcoa_amf_uni$permanova$R2[1], 2)`, p=`r round(pcoa_amf_uni$permanova$Pr[1], 3)`)$. 
#' 
#' Let's view a plot with abundances of community subgroups inset. 
pcoa_amf_uni$ord <- 
    ggplot(pcoa_amf_uni$site_vectors, aes(x = Axis.1, y = Axis.2)) +
    geom_point(aes(fill = field_type, shape = region), size = 10) +
    geom_text(aes(label = yr_since)) +
    scale_fill_discrete_qualitative(palette = "harmonic") +
    scale_shape_manual(values = c(21, 22, 23, 24)) +
    labs(x = paste0("Axis 1 (", pcoa_amf_uni$eig[1], "%)"), 
         y = paste0("Axis 2 (", pcoa_amf_uni$eig[2], "%)"), 
         title = paste0("PCoA Ordination of field-averaged species data (", pcoa_amf_uni$dataset, ")"),
         caption = "Text indicates years since restoration, with corn (-) and remnants (+) never restored.") +
    # lims(x = c(-0.6,0.35)) +
    theme_bw() +
    guides(fill = guide_legend(override.aes = list(shape = 21)))
#+ amf_uni_families_fig,fig.align='center'
# pcoa_amf_bray$inset reused here because it doesn't change
pcoa_amf_uni$ord +
    annotation_custom(
        ggplotGrob(pcoa_amf_bray$inset + theme(
            plot.background = element_rect(colour = "black", fill = "gray90")
        )),
        xmin = 0.07,
        xmax = 0.19,
        ymin = 0.015,
        ymax = 0.09
    )
#'
#' Community trajectories revealed in the ordination separate cornfields from everything else. 
#' Using UNIFRAC distance has really dissolved most of what was apparent with the Bray-Curtis distance.  
#' Corn fields stand well apart with AMF communities, but no signal appears for other field types
#' or for years since restoration. I guess what this shows is that for AMF, restored fields 
#' almost immediately resemble remnants (but there must be some outlier taxa in Eric Rahnheim's place). 
#' 
#' What's becoming apparent here is that Axis 1 separates strongly on *field_type* and years 
#' since restoration, and Axis 2 further separates on years since restoration. A consistent signal
#' of region isn't obvious. 
#' 
#' Let's test the relationship between age and community axis scores with restored fields 
#' only. I don't expect much. 
#+ amf_uni_yrs_scores_data
amf_uni_resto_scores <-
    pcoa_amf_uni$site_vectors %>%
    filter(field_type == "restored") %>%
    mutate(yr_since = as.numeric(yr_since))
#+ amf_uni_yrs_scores_lm_1
summary(lm(
    Axis.1 ~ yr_since,
    data = amf_uni_resto_scores
))
#+ amf_uni_yrs_scores_lm_2
summary(lm(
    Axis.2 ~ yr_since,
    data = amf_uni_resto_scores
))
#+ amf_uni_yrs_scores_fig,fig.align='center',message=FALSE
amf_uni_resto_scores %>%
    pivot_longer(Axis.1:Axis.2, names_to = "axis", values_to = "score") %>%
    ggplot(aes(x = yr_since, y = score)) +
    facet_wrap(vars(axis), scales = "free") +
    geom_smooth(method = "lm", se = FALSE, linewidth = 0.5) +
    geom_point(aes(shape = region), fill = "grey", size = 2) +
    labs(x = "Years since restoration",
         y = "PCoA axis score",
         title = "Correlations, axis scores and years since restoration (18S gene, 97% OTU, UNIFRAC distance)",
         caption = "Blue lines show linear model fit; solid line is significant at p<0.05") +
    scale_shape_manual(values = c(21, 22, 23, 24)) +
    theme_bw()
#' 
#' Both axes correlate significantly but with less than moderate strength with years since restoration. 
#' Axis 2 again shows a stronger relationship $(R^2_{Adj}=0.31, p<0.05)$, and Axis 1 
#' is close with $(R^2_{Adj}=0.30, p<0.05)$ 
#' 

#' Try fitting years since restoration with envfit...
