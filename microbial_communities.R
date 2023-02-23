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
#'     fig_height: 6
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
packages_needed = c("tidyverse", "vegan", "colorspace", "ape", "knitr", "gridExtra")
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
#' CSV files were produced in [process_data.R](process_data.md)
spe <- list(
    its_otu = read_csv(
        paste0(getwd(), "/clean_data/spe_ITS_otu_siteSpeMatrix_avg.csv"),
        show_col_types = FALSE
    ),
    its_sv  = read_csv(
        paste0(getwd(), "/clean_data/spe_ITS_sv_siteSpeMatrix_avg.csv"),
        show_col_types = FALSE
    ),
    amf_otu = read_csv(
        paste0(getwd(), "/clean_data/spe_18S_otu_siteSpeMatrix_avg.csv"),
        show_col_types = FALSE
    ),
    amf_sv  = read_csv(
        paste0(getwd(), "/clean_data/spe_18S_sv_siteSpeMatrix_avg.csv"),
        show_col_types = FALSE
    )
)
#' ## Species metadata
#' Needed to make inset figures showing most important categories of species. 
#' CSV files were produced in the [microbial diversity script](microbial_diversity.md).
spe_meta <- list(
    its_otu =
        read_csv(
            paste0(getwd(), "/clean_data/speGuild_ITS_otu.csv"),
            show_col_types = FALSE
        ),
    its_sv  =
        read_csv(
            paste0(getwd(), "/clean_data/speGuild_ITS_sv.csv"),
            show_col_types = FALSE
        ),
    amf_otu = 
        read_csv(
            paste0(getwd(), "/clean_data/speTaxa_18S_otu.csv"),
            show_col_types = FALSE
        ),
    amf_sv = 
        read_csv(
            paste0(getwd(), "/clean_data/speTaxa_18S_sv.csv"),
            show_col_types = FALSE
        )
)
#' 
#' ## Site metadata
#' Needed for figure interpretation and permanova designs.
sites <- read_csv(paste0(getwd(), "/clean_data/site.csv"), show_col_types = FALSE) %>% 
    mutate(field_type = factor(site_type, ordered = TRUE, levels = c("corn", "restored", "remnant")),
           yr_since = replace(yr_since, which(site_type == "remnant"), "+"),
           yr_since = replace(yr_since, which(site_type == "corn"), "-")) %>% 
    filter(site_type != "oldfield") %>% 
    select(-lat, -long, -yr_restore, -site_type)
#' 
#' ## Functions
#' A function handles the Principal Components Analysis (PCoA) diagnostics, with outputs and figures 
#' saved to a list for later use. 
#' 
#' USE THE STRATA ARGUMENT IN ADONIS2
#' 
pcoa_fun <- function(data, env=sites, corr="none", d_method="bray", df_name, nperm=1999) {
    # Multivariate analysis
    d <- vegdist(data.frame(data, row.names = 1), d_method)
    p <- pcoa(d, correction = corr)
    p_vals <- data.frame(p$values) %>% 
        rownames_to_column(var = "Dim") %>% 
        mutate(Dim = as.integer(Dim))
    p_vec <- data.frame(p$vectors)
    # Permutation tests (PERMANOVA)
    glo_perm <-
        adonis2(
            data.frame(data, row.names = 1) ~ field_type + region,
            data = env,
            by = NULL,
            permutations = nperm,
            method = d_method
        )
    mar_perm <-
        adonis2(
            data.frame(data, row.names = 1) ~ field_type + region,
            data = env,
            by = "margin",
            permutations = nperm,
            method = d_method
        )
    int_perm <-
        adonis2(
            data.frame(data, row.names = 1) ~ field_type * region,
            data = env,
            by = "margin",
            permutations = nperm,
            method = d_method
        )
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
    } else {
        p_bstick <- ggplot(p_vals, aes(x = factor(Dim), y = Rel_corr_eig)) + 
            geom_col(fill = "gray70", color = "gray30") + 
            geom_line(aes(x = Dim, y = Broken_stick), color = "red") +
            geom_point(aes(x = Dim, y = Broken_stick), color = "red") +
            labs(x = "Dimension", 
                 title = paste0("PCoA Eigenvalues and Broken Stick Model (", df_name, ")")) +
            theme_bw()
        p_ncomp <- with(p_vals, which(Rel_corr_eig < Broken_stick)[1]-1)
    }
    ncomp <- if(p_ncomp <= 2) {2} else {p_ncomp}
    # Ordination plot
    scores <-
        p_vec[, 1:ncomp] %>%
        rownames_to_column(var = "site_key") %>%
        mutate(site_key = as.integer(site_key)) %>%
        left_join(sites, by = "site_key") %>% 
        select(-site_name, -yr_rank)
    eig <- round(p_vals$Relative_eig[1:2] * 100, 1)
    ord <- ggplot(scores, aes(x = Axis.1, y = Axis.2)) +
        geom_point(aes(fill = field_type, shape = region), size = 10) +
        geom_text(aes(label = yr_since)) +
        scale_fill_discrete_qualitative(palette = "harmonic") +
        scale_shape_manual(values = c(21, 22, 23, 24)) +
        labs(x = paste0("Axis 1 (", eig[1], "%)"), 
             y = paste0("Axis 2 (", eig[2], "%)"), 
             title = paste0("PCoA Ordination of field-averaged species data (", df_name, ")"),
             caption = "Text indicates years since restoration, with corn (-) and remnants (+) never restored.") +
        theme_bw() +
        guides(fill = guide_legend(override.aes = list(shape = 21)))
    # Output data
    output <- list(components_exceed_broken_stick = p_ncomp,
                   correction_note = p$note,
                   values = p_vals[1:(ncomp+1), ], 
                   site_vectors = scores,
                   broken_stick_plot = p_bstick,
                   global_permanova = glo_perm,
                   margin_terms_permanova = mar_perm,
                   interaction_terms_permanova = int_perm,
                   ordination_plot = ord)
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
(pcoa_its_otu <- pcoa_fun(spe$its_otu, df_name = "ITS gene, 97% OTU"))
#' 
#' Axis 1 explains 19% of the variation and is the only eigenvalue that exceeds a 
#' broken stick model. The most substantial variation here will be on the first axis,
#' although axis 2 explains 11% of the variation and was very close to the broken
#' stick value. Testing the design factors *region* and *field_type* revealed a significant
#' global ordination ($R^2=0.36, p<0.001$), with independent factors both significant 
#' (*region:* $R^2=0.18, p<0.001$ and *field_type:* $R^2=0.17, p<0.001$). An interaction between
#' these factors was not supported. 
#' 
#' Community trajectories revealed in the ordination depend on both region and field type, so 
#' it's not surprising that both were significant. Faville Grove shows a linear progression
#' from corn to remnant and Lake Petite does as well, although with few sites and only 
#' single restoration ages these are weak supports. With Blue Mounds sites, the general 
#' progression along Axis 1 is to increase in age from left to right, but the remnant 
#' doesn't seem representative because it clusters far from everything else and associates
#' most strongly with the neighboring restored field (both on Merel Black's property). Restored fields 
#' at Fermi separate well away from cornfields, but less age structure is found. Instead, the 
#' old restorations in the ring most resemble the Railroad Remnant (which is in a different soil...), 
#' the switchgrass restored fields take a potentially novel path toward distant remnants. 
#' 
#' Here we can also begin considering what an inset plot to display metadata might look like. 
#' There isn't much room for it...maybe it goes along side. 
#+ its_guilds_fig,fig.align='center',fig.width=10,fig.height=5
grid.arrange(
    spe_meta$its_otu %>% 
        filter(guild %in% c("Undefined Saprotroph", "Wood Saprotroph", "Plant Pathogen", "Ectomycorrhizal", "Endophyte")) %>% 
        mutate(field_type = factor(field_type, ordered = TRUE, levels = c("corn", "restored", "remnant"))) %>% 
        group_by(guild, field_type) %>% 
        summarize(avg_seq_abund = mean(seq_abund), .groups = "drop") %>% 
        ggplot(aes(x = guild, y = avg_seq_abund)) +
        geom_col(aes(fill = field_type)) +
        labs(x = "",
             y = "Average sequence abundance",
             title = "Guild abundances") +
        scale_fill_discrete_qualitative(palette = "Harmonic") +
        coord_flip() +
        theme_classic() +
        theme(legend.position = "none"),
    pcoa_its_otu$ordination_plot,
    ncol = 2)
#' This is messy but a place to start thinking about Fig 1. 
#' 
#' Let's plot and test the relationship between age and community axis scores with restored fields 
#' only.
#+ its_yrs_scores_data
its_yrs_scores <-
    pcoa_its_otu$site_vectors %>%
    filter(field_type == "restored") %>%
    mutate(yr_since = as.numeric(yr_since))
#+ its_yrs_scores_fig,fig.align='center',message=FALSE
its_yrs_scores %>%
    pivot_longer(Axis.1:Axis.2, names_to = "axis", values_to = "score") %>%
    ggplot(aes(x = score, y = yr_since)) +
    facet_wrap(vars(axis), scales = "free") +
    geom_smooth(aes(linetype = axis), method = "lm", se = FALSE, linewidth = 0.5) +
    geom_point(aes(shape = region), fill = "grey", size = 2) +
    labs(x = "PCoA axis score",
         y = "Years since restoration",
         title = "Correlations, axis scores and years since restoration (ITS, 97% OTU)",
         caption = "Blue lines show linear model fit; solid line is significant at p<0.05") +
    scale_shape_manual(values = c(21, 22, 23, 24)) +
    scale_linetype_manual(values = c('solid', 'dashed'), guide = "none") +
    theme_bw()
#+ its_yrs_scores_lm
summary(lm(
    yr_since ~ Axis.1,
    data = its_yrs_scores
))
#' 
#' Indeed, Axis 1 does correlate well with age ($R^2_{Adj}=0.51, p<0.005$). 
#' 
#' ### PCoA with ITS gene, SV clusters
#+ pcoa_its_sv,fig.align='center'
(pcoa_its_sv <- pcoa_fun(spe$its_sv, df_name = "ITS gene, 100% SV"))
#' 
#' No axes are justified by the broken stick model. Otherwise, the results match what 
#' was found with OTU-based data with weaker support. 
#' 
#' ### PCoA with 18S gene, OTU clusters
#+ pcoa_amf_otu,fig.align='center'
(pcoa_amf_otu <- pcoa_fun(spe$amf_otu, df_name = "18S gene, 97% OTU"))
#' 
#' Four axes are significant by a broken stick model, between them explaining 64% of the 
#' variation in AMF among fields. It may be worthwhile to examine structure on Axes 3 and 4 
#' sometime. The most substantial variation here is on the first axis (28%) with Axis 2
#' explaining 18% of the variation in AMF abundances. 
#' Testing the design factors *region* and *field_type* revealed a significant
#' global ordination ($R^2=0.42, p<0.001$), with independent factors both significant 
#' (*region:* $R^2=0.16, p<0.05$ and *field_type:* $R^2=0.25, p<0.001$). An interaction between
#' these factors was not supported. 
#' 
#' Community trajectories revealed in the ordination depend on both region and field type, so 
#' it's not surprising that both were significant. Corn fields stand well apart with AMF communities,
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
#' Here we can also begin considering what an inset plot to display metadata might look like. 
#' There isn't much room for it...maybe it goes along side. 
#+ amf_guilds_fig,fig.align='center',fig.width=10,fig.height=5
grid.arrange(
spe_meta$amf_otu %>% 
    mutate(field_type = factor(field_type, ordered = TRUE, levels = c("corn", "restored", "remnant"))) %>% 
    group_by(family, field_type) %>% 
    summarize(avg_seq_abund = mean(seq_abund), .groups = "drop") %>% 
    ggplot(aes(x = family, y = avg_seq_abund)) +
    geom_col(aes(fill = field_type)) +
    labs(x = "",
         y = "Average sequence abundance",
         title = "Family abundances") +
    scale_fill_discrete_qualitative(palette = "Harmonic") +
    coord_flip() +
    theme_classic() +
    theme(legend.position = "none"),
pcoa_amf_otu$ordination_plot,
ncol = 2)
#' Again, messy but useful.
#' 
#' Let's test the relationship between age and community axis scores with restored fields 
#' only.
#+ amf_yrs_scores_data
amf_yrs_scores <-
    pcoa_amf_otu$site_vectors %>%
    filter(field_type == "restored") %>%
    mutate(yr_since = as.numeric(yr_since))
#+ amf_yrs_scores_fig,fig.align='center',message=FALSE
amf_yrs_scores %>%
    pivot_longer(Axis.1:Axis.2, names_to = "axis", values_to = "score") %>%
    ggplot(aes(x = score, y = yr_since)) +
    facet_wrap(vars(axis), scales = "free") +
    geom_smooth(method = "lm", se = FALSE, linewidth = 0.5) +
    geom_point(aes(shape = region), fill = "grey", size = 2) +
    labs(x = "PCoA axis score",
         y = "Years since restoration",
         title = "Correlations, axis scores and years since restoration (AMF, 97% OTU)",
         caption = "Blue lines show linear model fit; solid lines are significant at p<0.05") +
    scale_shape_manual(values = c(21, 22, 23, 24)) +
    theme_bw()
#+ amf_yrs_scores_lm_1
summary(lm(
    yr_since ~ Axis.1,
    data = amf_yrs_scores
))
#+ amf_yrs_scores_lm_2
summary(lm(
    yr_since ~ Axis.2,
    data = amf_yrs_scores
))
#' 
#' Both axes correlate significantly and strongly with years since restoration. 
#' Axis 2 shows a stronger relationship ($R^2_{Adj}=0.68, p<0.001), and Axis 1 
#' shows a moderate relationship ($R^2_{Adj}=0.38, p<0.01) 
#' 
#' ### PCoA with 18S gene, SV clusters
#+ pcoa_amf_sv,fig.align='center'
(pcoa_amf_sv <- pcoa_fun(spe$amf_sv, df_name = "18S gene, 100% SV"))
#' 
#' These results align strongly with those obtained from OTU clusters, with just slightly weaker
#' support. Suggest that we continue with only OTUs from here on out. 

