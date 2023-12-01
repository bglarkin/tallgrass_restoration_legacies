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
#' Plant data comprises two data sets. Mike Healy did quadrat surveys at 
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
packages_needed = c("GGally", "tidyverse", "vegan", "colorspace", "ape")
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
#' ## Cleanplot PCA
#' Cleanplot PCA produces informative visualizations of PCA ordinations 
#' [(Borcard et al. 2018)](http://link.springer.com/10.1007/978-3-319-71404-2)
#+ cleanplot_pca
source(paste0(getwd(), "/supporting_files/cleanplot_pca.txt"))
#' 
#' ## PCoA
#+ pcoa_function
pcoa_fun <- function(s, d, env=sites, corr="none", df_name, nperm=1999) {
    set.seed <- 397
    # Multivariate analysis
    p <- pcoa(d, correction = corr)
    p_vals <- data.frame(p$values) %>% 
        rownames_to_column(var = "Dim") %>% 
        mutate(Dim = as.integer(Dim))
    p_vec <- data.frame(p$vectors)
    # Wrangle site data
    env_w <- env %>% filter(field_name %in% s$field_name)
    # Permutation tests (PERMANOVA)
    h <- with(env_w, 
              how(within = Within(type="none"), 
                  plots  = Plots(strata=field_name, type="free"),
                  blocks = region,
                  nperm  = nperm))
    p_permtest <- adonis2(d ~ field_type, data = env_w, permutations = h)
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
        rownames_to_column(var = "field_name") %>%
        left_join(sites, by = "field_name")
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
pcoa_samps_fun <- function(s, m, bi=FALSE, env, corr="none", df_name, nperm=1999) {
    set.seed <- 438
    # Create distance object, check for zero columns
    s_df <- s %>% 
        mutate(field_sample = paste(field_name, sample, sep = "_")) %>% 
        select(field_sample, everything(), -field_name, -sample, -region, -field_type) %>% 
        data.frame(row.names = 1)
    s_nz <- s_df[, which(apply(s_df, 2, sum) > 0)]
    d <- vegdist(s_nz, method = m, binary = bi)
    # Multivariate analysis
    p <- pcoa(d, correction = corr)
    p_vals <- data.frame(p$values) %>% 
        rownames_to_column(var = "Dim") %>% 
        mutate(Dim = as.integer(Dim))
    p_vec <- data.frame(p$vectors)
    # Wrangle site data
    env_w <- data.frame(field_sample = rownames(s_nz)) %>% 
        separate_wider_delim(field_sample, delim = "_", names = c("field_name", "sample"), cols_remove = FALSE) %>% 
        left_join(env, by = join_by(field_name)) %>% 
        column_to_rownames(var = "field_sample")
    # Permutation tests (PERMANOVA)
    # Fields as replicate strata with subsamples
    # Regions as blocks
    h <- with(env_w, 
              how(within = Within(type="free"), 
                  blocks = region,
                  nperm  = nperm))
    p_permtest <- adonis2(
        d ~ field_name,
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
    fit <- NULL
    fit_sc <- NULL
    if(m == "bray") {
        # Permutation test (ENVFIT)
        # Fields as strata
        h = with(env_w, 
                 how(within = Within(type="none"), 
                     plots = Plots(strata=field_key, type="free"), 
                     nperm = nperm))
        fit <- envfit(p_vec ~ as.numeric(yr_since), env_w, permutations = h, choices = c(1:ncomp))
        fit_sc <- scores(fit, c("vectors"))
    }
    # Ordination plot
    scores <-
        p_vec[, 1:ncomp] %>%
        rownames_to_column(var = "field_sample") %>%
        separate_wider_delim(field_sample, delim = "_", names = c("field_name", "sample"), cols_remove = TRUE) %>% 
        left_join(env, by = join_by(field_name))
    # Output data
    output <- list(dataset           = df_name,
                   important_comp    = p_ncomp,
                   correction        = p$note,
                   values            = p_vals[1:(ncomp+1), ], 
                   eigenvalues       = eig,
                   site_vectors      = scores,
                   broken_stick_plot = p_bstick,
                   permanova         = p_permtest,
                   vector_fit        = fit,
                   vector_fit_scores = fit_sc)
    return(output)
}
#' 
#' # Data
#' ## Metadata from sites, as in previous
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
#' ## Plant communities, subsample data
#' Raw field data can be used to compare among replicates (fields). Long-form data must be wrangled into 
#' matrices and then converted to distance objects. 
plant_field_data <- read_csv(paste0(getwd(), "/plant_field_data/plant_cover.csv"), show_col_types = FALSE) %>% 
    select(-SITE, -TYPE, -FIELDNUM) %>% 
    rename(field_name = FIELDID, sample = PLOT, field_sample = PLOTID, code = PLANTCODE, cover_pct = PCTCVR, presence = PRESENCE) %>% 
    filter(field_name %in% sites$field_name) %>% 
    left_join(sites %>% select(field_name, region, field_type, field_key), by = join_by(field_name)) %>% 
    arrange(field_key, sample, code) %>% 
    select(region, field_type, field_name, sample, field_sample, code, cover_pct, presence)
#+ plant_abiotic_csv
plant_abiotic <- plant_field_data %>% 
    filter(code %in% c("BARESOIL", "LITTER")) %>% 
    write_csv(paste0(getwd(), "/clean_data/plant_abiotic.csv"))
#+ plant_abundance_samples_csv
plant_abund_samples <- plant_field_data %>% 
    filter(!(code %in% c("BARESOIL", "LITTER")),
           region != "FL") %>% 
    select(region, field_type, field_name, sample, code, cover_pct) %>% 
    arrange(code) %>% 
    pivot_wider(names_from = code, values_from = cover_pct, values_fill = 0) %>% 
    arrange(field_name, sample) %>% 
    write_csv(paste0(getwd(), "/clean_data/plant_abund_samples.csv"))
#' 
#' ## Plant communities data
#' 
#' - Metadata, taxonomy and traits
#' - Abundance data, surveyed in 2016, sites limited to Wisconsin only
#' - Abundance data at the level of field samples (to compare among fields)
#' - Presence data, from Fermi in 2015, all other sites converted to presence data, does not 
#' include Fermi switchgrass or corn fields.
#+ plant_data_list
plant <- list(
    meta    = read_csv(paste0(getwd(), "/clean_data/spe_plant_meta.csv"), show_col_types = FALSE) %>% 
        rename_with(tolower),
    ab      = read_csv(paste0(getwd(), "/clean_data/spe_plant_abund.csv"), show_col_types = FALSE) %>% 
        rename(field_name = SITE) %>% select(-BARESOIL, -LITTER),
    ab_samp = read_csv(paste0(getwd(), "/clean_data/plant_abund_samples.csv"), show_col_types = FALSE),
    pr      = read_csv(paste0(getwd(), "/clean_data/spe_plant_presence.csv"), show_col_types = FALSE) %>% 
        rename(field_name = SITE) %>% select(-BARESOIL, -LITTER)
)
#' 
#' ## Distance matrices
#' Distance matrices are needed for ordinations of the plant data. Bray-Curtis distance is used 
#' for abundance data, the Jaccard similarity is used for binary data.
#+ distab_list
distab = list(
    p_ab = vegdist(data.frame(plant$ab, row.names = 1), method = "bray"),
    p_pr = vegdist(data.frame(plant$pr, row.names = 1), method = "jac", binary = TRUE)
)
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
#' These data are problematic because they are essentially counts of traits, and these counts
#' aren't related to cover. 
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
#' ## Traits: sites with abundance data
#' Run a PCA on chord-transformed traits data from sites with abundance data, perform 
#' typical basic diagnostics. This should be done without corn fields because they exert 
#' too strong a difference on everything else. 
p_ab_trait_ch <- decostand(
    data.frame(
        p_ab_trait %>% 
            filter(field_name %in% sites_noc$field_name), 
        row.names = 1), "normalize")
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
#' ## Traits: sites with presence data
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
write_csv(p_pr_trait, paste0(getwd(), "/clean_data/plant_trait_presence.csv"))
#' 
#' ## Plant Communities
#' In restored fields, plant communities don't reflect natural community assembly. Still, it's useful to 
#' examine an ordination of sites to develop an understanding of how they differ. Plant traits data are probably 
#' more useful to looking at development of plant communities over time after restoration. Traits may be filtered
#' more than species, and species' occurrence may not be uniform across sites, though a species' realized niche may 
#' include sites where it is not found. 
#' ### Sites with abundance data
#' An ordiation is run on plant abundance data using `pcoa_fun()`.
#+ pcoa_ab
(pcoa_ab <- pcoa_fun(plant$ab, distab$p_ab, corr="lingoes", df_name = "plant abundance data, 16 sites"))
#' Axis 1 explains `r pcoa_ab$eigenvalues[1]`% of the variation and is the only eigenvalue that exceeds a 
#' broken stick model. The most substantial variation here will be on the first axis.
#' Axis 2 explains `r pcoa_ab$eigenvalues[2]`% of the variation and was not very close to the broken
#' stick value. Testing the design factor *field_type* (with *region* treated as a block 
#' using arguments to `how()` revealed a significant
#' clustering $(R^2=`r round(pcoa_ab$permanova$R2[1], 2)`,~p=`r pcoa_ab$permanova$Pr[1]`)$. 
#' Let's view a plot of these results. 
#+ pcoa_ab_fig,fig.align='center'
ggplot(pcoa_ab$site_vectors, aes(x = Axis.1, y = Axis.2)) +
    geom_point(aes(fill = field_type, shape = region), size = 10) +
    geom_text(aes(label = yr_since)) +
    scale_fill_discrete_qualitative(name = "Field Type", palette = "harmonic") +
    scale_shape_manual(name = "Region", values = c(21, 22, 23, 24)) +
    labs(
        x = paste0("Axis 1 (", pcoa_ab$eig[1], "%)"),
        y = paste0("Axis 2 (", pcoa_ab$eig[2], "%)"),
        title = paste0(
            "PCoA Ordination of field-averaged species data (",
            pcoa_ab$dataset,
            ")"
        ),
        caption = "Text indicates years since restoration, with corn (-) and remnants (+) never restored."
    ) +
    theme_bw() +
    guides(fill = guide_legend(override.aes = list(shape = 21)))
#' 
#' ### Sites with presence data
#' An ordiation is run on plant presence data using `pcoa_fun()`. The dataset includes 20 sites. This 
#' analysis isn't appropriate because the blocks are unbalanced (no cornfield data from Fermi), but it 
#' still shows differences with plant data. 
#+ pcoa_pr
(pcoa_pr <- pcoa_fun(plant$pr, distab$p_pr, corr="none", df_name = "plant presence data, 20 sites"))
#' Axis 1 explains `r pcoa_pr$eigenvalues[1]`% of the variation and axis 2 explains `r pcoa_ab$eigenvalues[2]`% 
#' of the variation. These two eigenvalues exceed the broken stick value. 
#' stick value. Testing the design factor *field_type* (with *region* treated as a block 
#' using arguments to `how()` revealed a significant
#' clustering $(R^2=`r round(pcoa_pr$permanova$R2[1], 2)`,~p=`r pcoa_pr$permanova$Pr[1]`)$. 
#' Let's view a plot of these results. 
#+ pcoa_pr_fig,fig.align='center'
ggplot(pcoa_pr$site_vectors, aes(x = Axis.1, y = Axis.2)) +
    geom_point(aes(fill = field_type, shape = region), size = 10) +
    geom_text(aes(label = yr_since)) +
    scale_fill_discrete_qualitative(name = "Field Type", palette = "harmonic") +
    scale_shape_manual(name = "Region", values = c(21, 22, 23, 24)) +
    labs(
        x = paste0("Axis 1 (", pcoa_pr$eig[1], "%)"),
        y = paste0("Axis 2 (", pcoa_pr$eig[2], "%)"),
        title = paste0(
            "PCoA Ordination of field-averaged species data (",
            pcoa_pr$dataset,
            ")"
        ),
        caption = "Text indicates years since restoration, with corn (-) and remnants (+) never restored."
    ) +
    theme_bw() +
    guides(fill = guide_legend(override.aes = list(shape = 21)))
#' The regional signal is most obvious here. 
#' 
#' ### Restoration sites
#' #### Wisconsin regions
#' Plant communities in restoration fields are artificial. It's important to know how much
#' they differ among fields and if differences are related to years since restoration. We will
#' look at the Wisconsin regions first and permute within regions. 
#+ pcoa_ab_samps_wi,message=FALSE,warnings=FALSE
(pcoa_ab_samps_wi <- 
    pcoa_samps_fun(s = plant$ab_samp %>% filter(field_type == "restored", region != "FL"),
               m = "bray",
               env = sites %>% filter(field_type == "restored", region != "FL"),
               corr = "lingoes",
               df_name = "Plant abundance data, subsamples, Wisconsin regions",
               nperm = 1999))
#' Let's view an ordination plot with hulls around subsamples.  
#+ pcoa_ab_samps_wi_plotdata
centroid_ab_samps_wi <- aggregate(cbind(Axis.1, Axis.2) ~ field_name, data = pcoa_ab_samps_wi$site_vectors, mean) %>% 
    left_join(sites, by = join_by(field_name))
hull_ab_samps_wi <- pcoa_ab_samps_wi$site_vectors %>% 
    group_by(field_name) %>% 
    slice(chull(Axis.1, Axis.2))
#+ pcoa_ab_samps_wi_fig,fig.align='center',message=FALSE
ggplot(pcoa_ab_samps_wi$site_vectors, aes(x = Axis.1, y = Axis.2)) +
    geom_point(aes(fill = region), shape = 21) +
    geom_polygon(data = hull_ab_samps_wi, aes(group = field_name, fill = region), alpha = 0.3) +
    geom_point(data = centroid_ab_samps_wi, aes(fill = region), shape = 21, size = 8) +
    geom_text(data = centroid_ab_samps_wi, aes(label = yr_since)) +
    geom_segment(aes(x = 0, 
                     y = 0, 
                     xend = pcoa_ab_samps_wi$vector_fit_scores[1] * 0.65, 
                     yend = pcoa_ab_samps_wi$vector_fit_scores[2] * 0.65),
                 color = "blue", 
                 arrow = arrow(length = unit(3, "mm"))) +
    labs(
        x = paste0("Axis 1 (", pcoa_ab_samps_wi$eigenvalues[1], "%)"),
        y = paste0("Axis 2 (", pcoa_ab_samps_wi$eigenvalues[2], "%)"),
        title = paste0(
            "PCoA Ordination (",
            pcoa_ab_samps_wi$dataset,
            ")"
        ),
        caption = "Text indicates years since restoration\nYears since restoration significant at p<0.05"
    ) +
    scale_fill_discrete_qualitative(name = "Region", palette = "Dynamic") +
    theme_bw() +
    guides(fill = guide_legend(override.aes = list(shape = 21)))
#' 
#' The permutation test reveals that subsamples cluster to fields based on plant communities, and 
#' that years since restoration is significantly related community difference. This shows that 
#' plant communities and time since restoration are potentially confounded as explanatory variables
#' of soil microbial communities.  
#' 
#' #### Blue Mounds fields
#+ pcoa_ab_samps_bm,message=FALSE,warnings=FALSE
(pcoa_ab_samps_bm <- 
    pcoa_samps_fun(s = plant$ab_samp %>% filter(field_type == "restored", region == "BM"),
               m = "bray",
               env = sites %>% filter(field_type == "restored", region == "BM"),
               corr = "lingoes",
               df_name = "Plant abundance data, subsamples, Blue Mounds",
               nperm = 1999))
#' Let's view an ordination plot with hulls around subsamples.  
#+ pcoa_ab_samps_bm_plotdata
centroid_ab_samps_bm <- aggregate(cbind(Axis.1, Axis.2) ~ field_name, data = pcoa_ab_samps_bm$site_vectors, mean) %>% 
    left_join(sites, by = join_by(field_name))
hull_ab_samps_bm <- pcoa_ab_samps_bm$site_vectors %>% 
    group_by(field_name) %>% 
    slice(chull(Axis.1, Axis.2))
#+ pcoa_ab_samps_bm_fig,fig.align='center',message=FALSE
ggplot(pcoa_ab_samps_bm$site_vectors, aes(x = Axis.1, y = Axis.2)) +
    geom_point(aes(fill = region), shape = 21) +
    geom_polygon(data = hull_ab_samps_bm, aes(group = field_name, fill = region), alpha = 0.3) +
    geom_point(data = centroid_ab_samps_bm, aes(fill = region), shape = 21, size = 8) +
    geom_text(data = centroid_ab_samps_bm, aes(label = yr_since)) +
    geom_segment(aes(x = 0, 
                     y = 0, 
                     xend = pcoa_ab_samps_bm$vector_fit_scores[1] * 0.65, 
                     yend = pcoa_ab_samps_bm$vector_fit_scores[2] * 0.65),
                 color = "blue", 
                 arrow = arrow(length = unit(3, "mm"))) +
    labs(
        x = paste0("Axis 1 (", pcoa_ab_samps_bm$eigenvalues[1], "%)"),
        y = paste0("Axis 2 (", pcoa_ab_samps_bm$eigenvalues[2], "%)"),
        title = paste0(
            "PCoA Ordination (",
            pcoa_ab_samps_bm$dataset,
            ")"
        ),
        caption = "Text indicates years since restoration\nYears since restoration significant at p<0.05"
    ) +
    scale_fill_discrete_qualitative(name = "Field Type", palette = "Dynamic") +
    theme_bw() +
    guides(fill = guide_legend(override.aes = list(shape = 21)))
#' 
#' The Blue Mounds area is our most defensible chronosequence. But even here, 
#' the permutation test reveals that subsamples cluster to fields based on plant communities, and 
#' that years since restoration is significantly related community difference. This shows that 
#' plant communities and time since restoration are potentially confounded as explanatory variables
#' of soil microbial communities.  
#' 
#' ## Plant traits
#' The inverse relationship of forbs and C4 grasses over time is a common feature of tallgrass 
#' prairie restoration. In this case, some of that could be happening naturally, but we also know
#' that the oldest field was planted heavily with bluestem. 
#' 
#' We saw above that plant community change correlated with years since restoration, and it's a 
#' strong relationship with the shift from forbs to C4 grasses over time, for restored sites in Wisconsin:
#+ trait_time_cor_fig_wi,fig.align='center'
p_ab_trait %>% 
    left_join(sites, by = join_by(field_name)) %>% 
    filter(region != "FL", field_type == "restored") %>% 
    mutate(yr_since = as.numeric(yr_since)) %>% 
    select(yr_since, C4_grass, forb) %>% 
    ggpairs(title = "Trait correlations with time in Wisconsin restored fields") +
    theme_bw()
#' 
#' And it's even stronger in just the Blue Mounds:
#+ trait_time_cor_fig_bm,fig.align='center'
p_ab_trait %>% 
    left_join(sites, by = join_by(field_name)) %>% 
    filter(region == "BM", field_type == "restored") %>% 
    mutate(yr_since = as.numeric(yr_since)) %>% 
    select(yr_since, C4_grass, forb) %>% 
    ggpairs(title = "Trait correlations with time in Blue Mounds restored fields") +
    theme_bw()