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
#' based on ITS and 18S data. 
#' 
#' During data processing, not all samples were retained. Some had failed to amplify and others
#' had very few sequences, leading to the potential for a loss of information during rarefication. 
#' With the loss of some samples, all fields were resampled to the same lower number of samples. 
#' This was done to equalize sampling effort (from a statistical perspective). 
#' This procedure can easily be undone in the [process_data script](process_data.md)
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
#' Sites-species tables with rarefied sequence abundances. This list includes
#' composition summarized by fields or unsummarized (all samples). 
#' CSV files were produced in [process_data.R](process_data.md)
spe <- list(
    its = read_csv(
        paste0(getwd(), "/clean_data/spe_ITS_rfy.csv"),
        show_col_types = FALSE
    ),
    its_samps = read_csv(
        paste0(getwd(), "/clean_data/spe_ITS_rfy_samples.csv"),
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
#' Needed for figure interpretation and permanova designs. The subset of restored fields in Blue Mounds 
#' only will also be used and is parsed here.
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
sites_resto_bm <- 
    sites %>% 
    filter(field_type == "restored",
           region == "BM") %>% 
    select(-field_name, -region) %>% 
    mutate(yr_since = as.numeric(yr_since))   
#' 
#' ## Distance tables
#' Creating distance objects from the samples-species tables is done with the typical 
#' process of `vegdist()` in vegan.
#' Bray-Curtis or Ruzicka (used with method="jaccard") distance are both appropriate 
#' methods for these data, but Bray-Curtis has produced axes with better explanatory power. 
#' With the 18S data, we can take advantage of phylogenetic relationships in a UNIFRAC distance
#' matrix. The UNIFRAC distance was produced in QIIME II and needs some wrangling to 
#' conform to the standards of a distance object in R. The following list contains vegdist-produced
#' distance objects for ITS and 18S, and it includes UNIFRAC distance for 18S. 
#' 
#' **List of objects in `distab`**
#' 
#' - its: the rarefied data, summed from 8 samples in each field
#' - its_samps: rarefied data from 8 samples per field, all fields retained
#' - its_resto_bm: rarefied data, summed from 8 samples in each field, filtered to include Blue Mounds region only
#' - its_resto_samps_bm: rarefied data from 8 samples in each field, not summed, filtered to include Blue Mounds region only
#' - amf_bray: rarefied data, summed from 7 samples from each field, bray-curtis distance
#' - amf_uni: rarefied data, summed from 7 samples from each field, UNIFRAC distance
#' 
#+ distab_list
distab <- list(
    its       = vegdist(data.frame(spe$its, row.names = 1), method = "bray"),
    its_samps = vegdist(
        data.frame(
            spe$its_samps %>% 
                mutate(field_sample = paste(field_key, sample, sep = "_")) %>% 
                column_to_rownames(var = "field_sample") %>% 
                select(-field_key, -sample)
        )
    ),
    its_resto_bm = vegdist(
        data.frame(
            spe$its %>% 
                filter(field_key %in% sites_resto_bm$field_key), 
            row.names = 1
            ), 
        method = "bray"),
    its_resto_samps_bm = vegdist(
        data.frame(
            spe$its_samps %>% 
                filter(field_key %in% sites_resto_bm$field_key) %>% 
                mutate(field_sample = paste(field_key, sample, sep = "_")) %>% 
                column_to_rownames(var = "field_sample") %>% 
                select(-field_key, -sample)
            ),
        method = "bray"),
    amf_bray  = vegdist(data.frame(spe$amf, row.names = 1), method = "bray"),
    amf_uni   = sites %>%
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
#' Functions handle the Principal Components Analysis (PCoA) diagnostics, with outputs and figures 
#' saved to a list for later use. 
#' 
#' - `pcoa_fun()` is used with data where samples have been summed in fields. 
#' - `pcoa_samps_fun()` is used with rarefied subsample data from all fields.
#' - `pcoa_samps_bm_fun()` is used for the subsample data from Blue Mounds restored fields. The variable 
#' **yr_since** is continuous with this dataset and is tested with `envfit()`.
#' 
#+ pcoa_function
pcoa_fun <- function(d, env=sites, corr="none", df_name, nperm=1999) {
    set.seed <- 397
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
    output <- list(dataset                        = df_name,
                   components_exceed_broken_stick = p_ncomp,
                   correction_note                = p$note,
                   values                         = p_vals[1:(ncomp+1), ], 
                   eigenvalues                    = eig,
                   site_vectors                   = scores,
                   broken_stick_plot              = p_bstick,
                   permanova                      = p_permtest)
    return(output)
}
#' 
#+ pcoa_samps_function
pcoa_samps_fun <- function(s, d, env=sites, corr="none", df_name, nperm=1999) {
    set.seed <- 438
    # Multivariate analysis
    p <- pcoa(d, correction = corr)
    p_vals <- data.frame(p$values) %>% 
        rownames_to_column(var = "Dim") %>% 
        mutate(Dim = as.integer(Dim))
    p_vec <- data.frame(p$vectors)
    # Wrangle site data
    env_w <- env %>% 
        left_join(s %>% select(field_key, sample), by = join_by(field_key)) %>% 
        mutate(field_sample = paste(field_key, sample, sep = "_")) %>% 
        column_to_rownames(var = "field_sample")
    # Permutation tests (PERMANOVA)
    # Fields as replicate strata with subsamples
    # Regions as blocks
    h <- with(env_w, 
              how(within = Within(type="free"), 
                  plots  = Plots(strata=field_key, type="free"),
                  blocks = region,
                  nperm  = nperm))
    p_permtest <- adonis2(
        d ~ field_type,
        data = env_w,
        permutations = h)
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
        rownames_to_column(var = "field_sample") %>%
        separate_wider_delim(field_sample, delim = "_", names = c("field_key", "sample_key"), cols_remove = TRUE) %>% 
        mutate(field_key = as.integer(field_key)) %>%
        left_join(env, by = join_by(field_key)) %>% 
        select(-field_type)
    # Output data
    output <- list(dataset                        = df_name,
                   components_exceed_broken_stick = p_ncomp,
                   correction_note                = p$note,
                   values                         = p_vals[1:(ncomp+1), ], 
                   eigenvalues                    = eig,
                   site_vectors                   = scores,
                   broken_stick_plot              = p_bstick,
                   permanova                      = p_permtest)
    return(output)
}
#' 
#+ pcoa_samps_bm_function
pcoa_samps_bm_fun <- function(s, d, env=sites_resto_bm, corr="none", df_name, nperm=1999) {
    set.seed <- 845
    # Multivariate analysis
    p <- pcoa(d, correction = corr)
    p_vals <- data.frame(p$values) %>% 
        rownames_to_column(var = "Dim") %>% 
        mutate(Dim = as.integer(Dim))
    p_vec <- data.frame(p$vectors)
    # Wrangle site data
    env_w <- env %>% 
        left_join(s %>% select(field_key, sample), by = join_by(field_key)) %>% 
        mutate(field_sample = paste(field_key, sample, sep = "_")) %>% 
        column_to_rownames(var = "field_sample")
    # Permutation test (PERMANOVA)
    p_permtest <- adonis2(d ~ field_key, data = env_w, permutations = nperm)
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
    # Permutation test (ENVFIT)
    # Fields as strata
    h = with(env_w, 
             how(within = Within(type="none"), 
                 plots = Plots(strata=field_key, type="free"), 
                 nperm = nperm))
    fit <- envfit(p_vec ~ yr_since, env_w, permutations = h, choices = c(1:ncomp))
    fit_sc <- scores(fit, c("vectors"))
    # Ordination plotting data
    scores <-
        p_vec[, 1:ncomp] %>%
        rownames_to_column(var = "field_sample") %>%
        separate_wider_delim(field_sample, delim = "_", names = c("field_key", "sample_key"), cols_remove = TRUE) %>% 
        mutate(field_key = as.integer(field_key)) %>%
        left_join(env, by = join_by(field_key)) %>% 
        select(-field_type)
    # Output data
    output <- list(dataset                        = df_name,
                   components_exceed_broken_stick = p_ncomp,
                   correction_note                = p$note,
                   values                         = p_vals[1:(ncomp+1), ], 
                   eigenvalues                    = eig,
                   site_vectors                   = scores,
                   broken_stick_plot              = p_bstick,
                   permanova                      = p_permtest,
                   vector_fit                     = fit,
                   vector_fit_scores              = fit_sc)
    return(output)
}
#' 
#' # Results
#' ## Ordinations
#' Bray-Curtis or Ruzicka distance are both appropriate, but Bray-Curtis has 
#' produced axes with better explanatory power.
#' 
#' ### PCoA with ITS gene, OTU clusters
#' In trial runs, no negative eigenvalues were observed (not shown). No 
#' correction is needed for these ordinations.
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
    lims(y = c(-0.35,0.44)) +
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
        xmax = -0.05,
        ymin = 0.17,
        ymax = 0.46
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
#' Indeed, Axis 1 does correlate well with age $(R^2_{Adj}=0.51, p<0.005)$. But it isn't appropriate
#' to use these scores for the correlation because they were created with the corn and remnant fields
#' in the ordination as well. 
#' 
#' The most appropriate way to look at communities vs. field age is with the Blue Mounds restored
#' fields. The function `pcoa_its_samps_bm()` will take care of this. Field age will be fitted 
#' to the ordination and tested using `envfit()`.  
#' 
#' ### PCoA with ITS gene, OTU clusters, Blue Mounds restored fields, all subsamples
#' In trial runs, no negative eigenvalues were observed (not shown). No 
#' correction is needed for these ordinations.
#+ pcoa_its_samps_bm,fig.align='center'
(pcoa_its_samps_bm <- pcoa_samps_bm_fun(spe$its_samps, 
                                        distab$its_resto_samps_bm, 
                                        sites_resto_bm, 
                                        df_name="BM restored, ITS gene, 97% OTU"))
#' 
#' Axis 1 explains `r pcoa_its_samps_bm$eigenvalues[1]`% and axis 2 
#' explains `r pcoa_its_samps_bm$eigenvalues[2]`% of the variation in the community data. Both axes are important
#' based on the broken stick model. The relatively low percent variation explained is partly due to the 
#' high number of dimensions used when all samples from fields are included. 
#' The fidelity of samples to fields was significant based on a permutation test
#' $(R^2=`r round(pcoa_its_samps_bm$permanova$R2[1], 2)`, p=`r pcoa_its_samps_bm$permanova$Pr[1]`)$. 
#' In this case, the partial $R^2$ shows the proportion of sum of squares from the total. It is a low number
#' here because so much unexplained variation exists, resulting in a high sum of squares that is outside 
#' the assignment of subsamples to fields.
#' 
#' Years since restoration has a moderately strong correlation with communities and was significant 
#' with a permutation test where samples were constrained within
#' fields to account for lack of independence $(R^2=`r round(pcoa_its_samps_bm$vector_fit$vectors$r, 2)`)$.
#' 
#' Let's view an ordination plot with hulls around subsamples and a fitted vector for field age overlaid.  
#+ its_samps_bm_plotdata
centroid_its_bm <- aggregate(cbind(Axis.1, Axis.2) ~ field_key, data = pcoa_its_samps_bm$site_vectors, mean) %>% 
    left_join(sites %>% select(field_key, yr_since), by = join_by(field_key))
hull_its_bm <- pcoa_its_samps_bm$site_vectors %>% 
    group_by(field_key) %>% 
    slice(chull(Axis.1, Axis.2))
#+ its_samps_bm_fig,fig.align='center',message=FALSE
ggplot(pcoa_its_samps_bm$site_vectors, aes(x = Axis.1, y = Axis.2)) +
    geom_point(fill = "#5CBD92", shape = 21) +
    geom_polygon(data = hull_its_bm, aes(group = as.character(field_key)), fill = "#5CBD92", alpha = 0.3) +
    geom_point(data = centroid_its_bm, fill = "#5CBD92", size = 8, shape = 21) +
    geom_text(data = centroid_its_bm, aes(label = yr_since)) +
    geom_segment(aes(x = 0, 
                     y = 0, 
                     xend = pcoa_its_samps_bm$vector_fit_scores[1] * 0.4, 
                     yend = pcoa_its_samps_bm$vector_fit_scores[2] * 0.4),
                 color = "blue", 
                 arrow = arrow(length = unit(3, "mm"))) +
    labs(
        x = paste0("Axis 1 (", pcoa_its_samps_bm$eigenvalues[1], "%)"),
        y = paste0("Axis 2 (", pcoa_its_samps_bm$eigenvalues[2], "%)"),
        title = paste0(
            "PCoA Ordination (",
            pcoa_its_samps_bm$dataset,
            ")"
        ),
        caption = "Text indicates years since restoration.\nYears since restoration significant at p<0.05."
    ) +
    theme_bw() +
    theme(legend.position = "none")
#' 
#' ### PCoA with ITS gene, OTU clusters, all fields and regions, all subsamples
#' This leverages the information from all subsamples. Modifications to `how()` from 
#' package [permute](https://cran.r-project.org/package=permute) allow for the more complex design. 
#' 
#' Negative eigenvalues were produced in trial runs (not shown). A Lingoes correction was applied.
#+ pcoa_its_samps,fig.align='center'
(pcoa_its_samps <- pcoa_samps_fun(spe$its_samps, 
                                  distab$its_samps, 
                                  corr="lingoes", 
                                  df_name = "ITS gene, 97% OTU"))
#' 
#' Axis 1 explains `r pcoa_its_samps$eigenvalues[1]`% and axis 2 
#' explains `r pcoa_its_samps$eigenvalues[2]`% of the variation in the community data. Both axes are important
#' based on the broken stick model, in fact, the broken stick model shows that `r pcoa_its_samps$components_exceed_broken_stick`
#' axes are important in explaining variation with this dataset. 
#' The relatively low percent variation explained on axes 1 and 2 is partly due to the 
#' high number of dimensions used when all samples from fields are included. 
#' The fidelity of samples to fields was strong based on a permutation test when restricting permutations to
#' fields (=plots in `how()`) within regions (=blocks in `how()`) 
#' $(R^2=`r round(pcoa_its_samps$permanova$R2[1], 2)`, p=`r pcoa_its_samps$permanova$Pr[1]`)$. 
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


