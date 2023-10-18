Microbial data: community differences
================
Beau Larkin

Last updated: 18 October, 2023

- [Description](#description)
- [Packages and libraries](#packages-and-libraries)
- [Data](#data)
  - [Sites-species tables](#sites-species-tables)
  - [Species metadata](#species-metadata)
  - [Site metadata](#site-metadata)
  - [Distance tables](#distance-tables)
  - [Functions](#functions)
- [Results](#results)
  - [Ordinations](#ordinations)
    - [PCoA with ITS gene, OTU
      clusters](#pcoa-with-its-gene-otu-clusters)
    - [PCoA with ITS gene, OTU clusters, Blue Mounds restored fields,
      all
      subsamples](#pcoa-with-its-gene-otu-clusters-blue-mounds-restored-fields-all-subsamples)
    - [PCoA with ITS gene, OTU clusters, all fields and regions, all
      subsamples](#pcoa-with-its-gene-otu-clusters-all-fields-and-regions-all-subsamples)
    - [PCoA with 18S gene, uninformed
      distance](#pcoa-with-18s-gene-uninformed-distance)
    - [PCoA with 18S gene, UNIFRAC
      distance](#pcoa-with-18s-gene-unifrac-distance)

# Description

Microbial data include site-species tables derived from high-throughput
sequencing and clustering in QIIME by Lorinda Bullington and PLFA/NLFA
data which Ylva Lekberg did.

This presents basic visualizations of community differences among
sites/regions based on ITS and 18S data.

During data processing, not all samples were retained. Some had failed
to amplify and others had very few sequences, leading to the potential
for a loss of information during rarefication. With the loss of some
samples, all fields were resampled to the same lower number of samples.
This was done to equalize sampling effort (from a statistical
perspective). This procedure can easily be undone in the [process_data
script](process_data.md)

# Packages and libraries

``` r
packages_needed = c("tidyverse", "vegan", "colorspace", "ape", "knitr")
packages_installed = packages_needed %in% rownames(installed.packages())
```

``` r
if (any(!packages_installed)) {
    install.packages(packages_needed[!packages_installed])
}
```

``` r
for (i in 1:length(packages_needed)) {
    library(packages_needed[i], character.only = T)
}
```

# Data

## Sites-species tables

Sites-species tables with rarefied sequence abundances. This list
includes composition summarized by fields or unsummarized (all samples).
CSV files were produced in [process_data.R](process_data.md)

``` r
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
```

## Species metadata

Needed to make inset figures showing most important categories of
species. The OTUs and sequence abundances in these files matches the
rarefied data in `spe$` above. CSV files were produced in the [microbial
diversity script](microbial_diversity.md).

``` r
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
```

## Site metadata

Needed for figure interpretation and permanova designs. The subset of
restored fields in Blue Mounds only will also be used and is parsed
here.

``` r
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
```

## Distance tables

Creating distance objects from the samples-species tables is done with
the typical process of `vegdist()` in vegan. Bray-Curtis or Ruzicka
(used with method=“jaccard”) distance are both appropriate methods for
these data, but Bray-Curtis has produced axes with better explanatory
power. With the 18S data, we can take advantage of phylogenetic
relationships in a UNIFRAC distance matrix. The UNIFRAC distance was
produced in QIIME II and needs some wrangling to conform to the
standards of a distance object in R. The following list contains
vegdist-produced distance objects for ITS and 18S, and it includes
UNIFRAC distance for 18S.

**List of objects in `distab`**

- its: the rarefied data, summed from 8 samples in each field
- its_samps: rarefied data from 8 samples per field, all fields retained
- its_resto_bm: rarefied data, summed from 8 samples in each field,
  filtered to include Blue Mounds region only
- its_resto_samps_bm: rarefied data from 8 samples in each field, not
  summed, filtered to include Blue Mounds region only
- amf_bray: rarefied data, summed from 7 samples from each field,
  bray-curtis distance
- amf_uni: rarefied data, summed from 7 samples from each field, UNIFRAC
  distance

``` r
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
```

## Functions

Functions handle the Principal Components Analysis (PCoA) diagnostics,
with outputs and figures saved to a list for later use.

- `pcoa_fun()` is used with data where samples have been summed in
  fields.
- `pcoa_samps_fun()` is used with rarefied subsample data from all
  fields.
- `pcoa_samps_bm_fun()` is used for the subsample data from Blue Mounds
  restored fields. The variable **yr_since** is continuous with this
  dataset and is tested with `envfit()`.

``` r
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
```

``` r
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
```

``` r
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
```

# Results

## Ordinations

Bray-Curtis or Ruzicka distance are both appropriate, but Bray-Curtis
has produced axes with better explanatory power.

### PCoA with ITS gene, OTU clusters

In trial runs, no negative eigenvalues were observed (not shown). No
correction is needed for these ordinations.

``` r
(pcoa_its <- pcoa_fun(distab$its, df_name = "ITS gene, 97% OTU"))
```

    ## $dataset
    ## [1] "ITS gene, 97% OTU"
    ## 
    ## $components_exceed_broken_stick
    ## [1] 1
    ## 
    ## $correction_note
    ## [1] "There were no negative eigenvalues. No correction was applied"
    ## 
    ## $values
    ##   Dim Eigenvalues Relative_eig Broken_stick Cumul_eig Cumul_br_stick
    ## 1   1   1.2815624   0.18674310   0.15733159 0.1867431      0.1573316
    ## 2   2   0.7882559   0.11486085   0.11566492 0.3016039      0.2729965
    ## 3   3   0.5893471   0.08587682   0.09483159 0.3874808      0.3678281
    ## 
    ## $eigenvalues
    ## [1] 18.7 11.5
    ## 
    ## $site_vectors
    ##    field_key      Axis.1       Axis.2 region field_type yr_since
    ## 1          1  0.22530898 -0.053082773     BM   restored       16
    ## 2          2 -0.10355568  0.007898947     BM   restored        3
    ## 3          3 -0.33325192  0.024030077     FG       corn        -
    ## 4          4  0.10206801 -0.300619212     FG    remnant        +
    ## 5          5 -0.05803485 -0.286262691     FG   restored       15
    ## 6          6 -0.29490889  0.127468201     FL       corn        -
    ## 7          7 -0.32755915  0.024346959     FL       corn        -
    ## 8          8  0.09500886 -0.221442168     FL    remnant        +
    ## 9          9  0.24859534 -0.194513948     FL   restored       40
    ## 10        10  0.32311444 -0.082553592     FL   restored       36
    ## 11        11  0.21983056 -0.168618940     FL   restored       35
    ## 12        12  0.21443506  0.239228368     FL   restored       10
    ## 13        13  0.18556349  0.116510972     FL   restored       10
    ## 14        14  0.19523115  0.253841733     FL   restored       10
    ## 15        15  0.24572587 -0.052382971     BM   restored       28
    ## 16        16 -0.37667427  0.119738544     LP       corn        -
    ## 17        17  0.07394546  0.183054802     LP    remnant        +
    ## 18        18 -0.24265416  0.065593993     LP   restored        4
    ## 19        19 -0.14439863  0.090693568     LP   restored        4
    ## 20        20  0.24287879  0.312070983     BM    remnant        +
    ## 21        21  0.19143712  0.279206947     BM   restored       18
    ## 22        22 -0.10909894 -0.215513221     BM   restored        7
    ## 23        23 -0.22079818 -0.087178948     BM   restored        2
    ## 24        24 -0.33904077  0.030270315     BM       corn        -
    ## 25        25 -0.01316768 -0.211785946     BM   restored       11
    ## 
    ## $broken_stick_plot

<img src="microbial_communities_files/figure-gfm/pcoa_its_otu-1.png" style="display: block; margin: auto;" />

    ## 
    ## $permanova
    ## Permutation test for adonis under reduced model
    ## Terms added sequentially (first to last)
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1999
    ## 
    ## adonis2(formula = d ~ field_type, data = env, permutations = nperm, strata = region)
    ##            Df SumOfSqs      R2      F Pr(>F)    
    ## field_type  2   1.2189 0.17762 2.3758  5e-04 ***
    ## Residual   22   5.6438 0.82238                  
    ## Total      24   6.8627 1.00000                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Axis 1 explains 18.7% of the variation and is the only eigenvalue that
exceeds a broken stick model. The most substantial variation here will
be on the first axis, although axis 2 explains 11.5% of the variation
and was very close to the broken stick value. Testing the design factor
*field_type* (with *region* treated as a block using the `strata`
argument of `adonis2`) revealed a significant clustering
$(R^2=0.18, p=5\times 10^{-4})$.

Let’s view a plot with abundances of community subgroups inset.

``` r
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
```

``` r
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
```

<img src="microbial_communities_files/figure-gfm/its_guilds_fig-1.png" style="display: block; margin: auto;" />

Community trajectories revealed in the ordination clearly depend on both
region and field type. Faville Grove shows a linear progression from
corn to remnant and Lake Petite does as well, although with few sites
and only single restoration ages these are weak supports. With Blue
Mounds sites, the general progression along Axis 1 is to increase in age
from left to right, but the remnant doesn’t seem representative because
it clusters far from everything else and associates most strongly with
the neighboring restored field (both on Merel Black’s property).
Restored fields at Fermi separate well away from cornfields, but less
age structure is found. Instead, the old restorations in the ring most
resemble the Railroad Remnant (which is in a different soil…), the
switchgrass restored fields take a potentially novel path toward distant
remnants.

On axis 1, four clusters are apparent in at least two partitioning
schemes. It will be interesting to see if we can pull those apart with
explanatory variables.

Restoration age will be explored in-depth with the subset of restoration
fields.

Here we can also begin considering what an inset plot to display
metadata might look like. Let’s plot and test the relationship between
age and community axis scores with restored fields only.

``` r
its_resto_scores <-
    pcoa_its$site_vectors %>%
    filter(field_type == "restored") %>%
    mutate(yr_since = as.numeric(yr_since))
```

``` r
summary(lm(Axis.1 ~ yr_since,
           data = its_resto_scores))
```

    ## 
    ## Call:
    ## lm(formula = Axis.1 ~ yr_since, data = its_resto_scores)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.18186 -0.08899 -0.03211  0.10525  0.20614 
    ## 
    ## Coefficients:
    ##              Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept) -0.106858   0.053293  -2.005  0.06468 .  
    ## yr_since     0.011515   0.002702   4.262  0.00079 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.131 on 14 degrees of freedom
    ## Multiple R-squared:  0.5647, Adjusted R-squared:  0.5336 
    ## F-statistic: 18.16 on 1 and 14 DF,  p-value: 0.0007897

``` r
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
```

<img src="microbial_communities_files/figure-gfm/its_resto_scores_fig-1.png" style="display: block; margin: auto;" />

Indeed, Axis 1 does correlate well with age $(R^2_{Adj}=0.51, p<0.005)$.
But it isn’t appropriate to use these scores for the correlation because
they were created with the corn and remnant fields in the ordination as
well.

The most appropriate way to look at communities vs. field age is with
the Blue Mounds restored fields. The function `pcoa_its_samps_bm()` will
take care of this. Field age will be fitted to the ordination and tested
using `envfit()`.

### PCoA with ITS gene, OTU clusters, Blue Mounds restored fields, all subsamples

In trial runs, no negative eigenvalues were observed (not shown). No
correction is needed for these ordinations.

``` r
(pcoa_its_samps_bm <- pcoa_samps_bm_fun(spe$its_samps, 
                                        distab$its_resto_samps_bm, 
                                        sites_resto_bm, 
                                        df_name="BM restored, ITS gene, 97% OTU"))
```

    ## Set of permutations < 'minperm'. Generating entire set.
    ## Set of permutations < 'minperm'. Generating entire set.

    ## $dataset
    ## [1] "BM restored, ITS gene, 97% OTU"
    ## 
    ## $components_exceed_broken_stick
    ## [1] 2
    ## 
    ## $correction_note
    ## [1] "There were no negative eigenvalues. No correction was applied"
    ## 
    ## $values
    ##   Dim Eigenvalues Relative_eig Broken_stick Cumul_eig Cumul_br_stick
    ## 1   1    2.133078   0.11483922   0.08352022 0.1148392     0.08352022
    ## 2   2    1.427775   0.07686758   0.06533840 0.1917068     0.14885863
    ## 3   3    1.036148   0.05578344   0.05624749 0.2474902     0.20510612
    ## 
    ## $eigenvalues
    ## [1] 11.5  7.7
    ## 
    ## $site_vectors
    ## # A tibble: 56 × 5
    ##    field_key sample_key  Axis.1  Axis.2 yr_since
    ##        <dbl> <chr>        <dbl>   <dbl>    <dbl>
    ##  1         1 1          -0.194  -0.194        16
    ##  2         1 2          -0.219  -0.0553       16
    ##  3         1 4          -0.183  -0.0668       16
    ##  4         1 5          -0.280  -0.137        16
    ##  5         1 6          -0.331  -0.159        16
    ##  6         1 7          -0.103  -0.0802       16
    ##  7         1 9          -0.208  -0.153        16
    ##  8         1 10         -0.242  -0.0266       16
    ##  9         2 1           0.242   0.216         3
    ## 10         2 2           0.0928  0.0872        3
    ## # ℹ 46 more rows
    ## 
    ## $broken_stick_plot

<img src="microbial_communities_files/figure-gfm/pcoa_its_samps_bm-1.png" style="display: block; margin: auto;" />

    ## 
    ## $permanova
    ## Permutation test for adonis under reduced model
    ## Terms added sequentially (first to last)
    ## Permutation: free
    ## Number of permutations: 1999
    ## 
    ## adonis2(formula = d ~ field_key, data = env_w, permutations = nperm)
    ##           Df SumOfSqs      R2     F Pr(>F)    
    ## field_key  1   0.7907 0.04257 2.401  5e-04 ***
    ## Residual  54  17.7837 0.95743                 
    ## Total     55  18.5745 1.00000                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## $vector_fit
    ## 
    ## ***VECTORS
    ## 
    ##             Axis.1    Axis.2     r2 Pr(>r)  
    ## yr_since -0.996630  0.081979 0.7285 0.0205 *
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Plots: field_key, plot permutation: free
    ## Permutation: none
    ## Number of permutations: 5039
    ## 
    ## 
    ## 
    ## $vector_fit_scores
    ##             Axis.1     Axis.2
    ## yr_since -0.850678 0.06997293

Axis 1 explains 11.5% and axis 2 explains 7.7% of the variation in the
community data. Both axes are important based on the broken stick model.
The relatively low percent variation explained is partly due to the high
number of dimensions used when all samples from fields are included. The
fidelity of samples to fields was significant based on a permutation
test $(R^2=0.04, p=5\times 10^{-4})$. In this case, the partial $R^2$
shows the proportion of sum of squares from the total. It is a low
number here because so much unexplained variation exists, resulting in a
high sum of squares that is outside the assignment of subsamples to
fields.

Years since restoration has a moderately strong correlation with
communities and was significant with a permutation test where samples
were constrained within fields to account for lack of independence
$(R^2=0.73)$.

Let’s view an ordination plot with hulls around subsamples and a fitted
vector for field age overlaid.

``` r
centroid_its_bm <- aggregate(cbind(Axis.1, Axis.2) ~ field_key, data = pcoa_its_samps_bm$site_vectors, mean) %>% 
    left_join(sites %>% select(field_key, yr_since), by = join_by(field_key))
hull_its_bm <- pcoa_its_samps_bm$site_vectors %>% 
    group_by(field_key) %>% 
    slice(chull(Axis.1, Axis.2))
```

``` r
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
```

<img src="microbial_communities_files/figure-gfm/its_samps_bm_fig-1.png" style="display: block; margin: auto;" />

### PCoA with ITS gene, OTU clusters, all fields and regions, all subsamples

This leverages the information from all subsamples. Modifications to
`how()` from package
[permute](https://cran.r-project.org/package=permute) allow for the more
complex design.

Negative eigenvalues were produced in trial runs (not shown). A Lingoes
correction was applied.

``` r
(pcoa_its_samps <- pcoa_samps_fun(spe$its_samps, 
                                  distab$its_samps, 
                                  corr="lingoes", 
                                  df_name = "ITS gene, 97% OTU"))
```

    ## $dataset
    ## [1] "ITS gene, 97% OTU"
    ## 
    ## $components_exceed_broken_stick
    ## [1] 10
    ## 
    ## $correction_note
    ## [1] "Lingoes correction applied to negative eigenvalues: D' = -0.5*D^2 - 0.0573041578256672 , except diagonal elements"
    ## 
    ## $values
    ##    Dim Eigenvalues Corr_eig Rel_corr_eig Broken_stick Cum_corr_eig Cum_br_stick
    ## 1    1    6.688218 6.745522   0.08173811   0.02963639   0.08173811   0.02963639
    ## 2    2    4.136197 4.193501   0.05081428   0.02458589   0.13255240   0.05422228
    ## 3    3    3.315460 3.372764   0.04086909   0.02206064   0.17342149   0.07628292
    ## 4    4    2.632045 2.689349   0.03258789   0.02037713   0.20600938   0.09666005
    ## 5    5    2.132678 2.189982   0.02653686   0.01911451   0.23254624   0.11577456
    ## 6    6    2.025277 2.082582   0.02523545   0.01810441   0.25778169   0.13387896
    ## 7    7    1.691608 1.748912   0.02119225   0.01726266   0.27897395   0.15114162
    ## 8    8    1.488491 1.545796   0.01873101   0.01654115   0.29770495   0.16768277
    ## 9    9    1.347531 1.404836   0.01702294   0.01590984   0.31472789   0.18359262
    ## 10  10    1.291010 1.348314   0.01633805   0.01534867   0.33106594   0.19894129
    ## 11  11    1.165587 1.222891   0.01481825   0.01484362   0.34588418   0.21378492
    ## 
    ## $eigenvalues
    ## [1] 8.2 5.1
    ## 
    ## $site_vectors
    ## # A tibble: 200 × 15
    ##    field_key sample_key  Axis.1   Axis.2  Axis.3  Axis.4   Axis.5   Axis.6
    ##        <dbl> <chr>        <dbl>    <dbl>   <dbl>   <dbl>    <dbl>    <dbl>
    ##  1         1 1           0.186  -0.130   -0.0367 -0.0640  0.0754  -0.105  
    ##  2         1 2           0.223   0.0732  -0.0260 -0.0145  0.0412  -0.216  
    ##  3         1 4           0.196  -0.0473  -0.0133  0.0262 -0.00645 -0.130  
    ##  4         1 5           0.178  -0.0133  -0.121  -0.206  -0.0348   0.0215 
    ##  5         1 6           0.212   0.00541 -0.117  -0.173   0.0334  -0.169  
    ##  6         1 7           0.0497 -0.0479   0.0220 -0.109  -0.0426   0.0127 
    ##  7         1 9           0.119  -0.0863  -0.0587 -0.0990  0.0623  -0.0324 
    ##  8         1 10          0.187   0.0135  -0.0704 -0.0394 -0.0233  -0.0261 
    ##  9         2 1          -0.0542 -0.00306  0.326   0.0518 -0.0599  -0.00542
    ## 10         2 2          -0.0373  0.0504   0.101   0.0261 -0.106   -0.155  
    ## # ℹ 190 more rows
    ## # ℹ 7 more variables: Axis.7 <dbl>, Axis.8 <dbl>, Axis.9 <dbl>, Axis.10 <dbl>,
    ## #   field_name <chr>, region <chr>, yr_since <chr>
    ## 
    ## $broken_stick_plot

<img src="microbial_communities_files/figure-gfm/pcoa_its_samps-1.png" style="display: block; margin: auto;" />

    ## 
    ## $permanova
    ## Permutation test for adonis under reduced model
    ## Terms added sequentially (first to last)
    ## Blocks:  region 
    ## Plots: field_key, plot permutation: free
    ## Permutation: free
    ## Number of permutations: 1999
    ## 
    ## adonis2(formula = d ~ field_type, data = env_w, permutations = h)
    ##             Df SumOfSqs     R2      F Pr(>F)    
    ## field_type   2    5.775 0.0812 8.7045  5e-04 ***
    ## Residual   197   65.348 0.9188                  
    ## Total      199   71.123 1.0000                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Axis 1 explains 8.2% and axis 2 explains 5.1% of the variation in the
community data. Both axes are important based on the broken stick model,
in fact, the broken stick model shows that 10 axes are important in
explaining variation with this dataset. The relatively low percent
variation explained on axes 1 and 2 is partly due to the high number of
dimensions used when all samples from fields are included. The fidelity
of samples to fields was strong based on a permutation test when
restricting permutations to fields (=plots in `how()`) within regions
(=blocks in `how()`) $(R^2=0.08, p=5\times 10^{-4})$.

### PCoA with 18S gene, uninformed distance

“Uninformed” refers to Bray-Curtis or another distance informational
metric available in `vegdist()`. It is as opposed to UNIFRAC distance,
which is informed by phylogeny.

``` r
(pcoa_amf_bray <- pcoa_fun(distab$amf_bray, df_name = "18S gene, 97% OTU, Bray-Curtis distance"))
```

    ## $dataset
    ## [1] "18S gene, 97% OTU, Bray-Curtis distance"
    ## 
    ## $components_exceed_broken_stick
    ## [1] 5
    ## 
    ## $correction_note
    ## [1] "No correction was applied to the negative eigenvalues"
    ## 
    ## $values
    ##   Dim Eigenvalues Relative_eig Rel_corr_eig Broken_stick Cum_corr_eig
    ## 1   1   1.2073294   0.27431093   0.25119312   0.16236050    0.2511931
    ## 2   2   0.7872231   0.17886079   0.16522784   0.11888224    0.4164210
    ## 3   3   0.6177502   0.14035575   0.13054904   0.09714311    0.5469700
    ## 4   4   0.4064211   0.09234078   0.08730531   0.08265036    0.6342753
    ## 5   5   0.3182007   0.07229670   0.06925300   0.07178079    0.7035283
    ## 6   6   0.2210942   0.05023363   0.04938234   0.06308514    0.7529106
    ##   Cumul_br_stick
    ## 1      0.1623605
    ## 2      0.2812427
    ## 3      0.3783858
    ## 4      0.4610362
    ## 5      0.5328170
    ## 6      0.5959021
    ## 
    ## $eigenvalues
    ## [1] 27.4 17.9
    ## 
    ## $site_vectors
    ##    field_key       Axis.1      Axis.2      Axis.3       Axis.4      Axis.5
    ## 1          1  0.197120583  0.23189804 -0.19499867  0.098536786  0.04573574
    ## 2          2 -0.004325496 -0.27761806 -0.20923662  0.218157865 -0.12400202
    ## 3          3 -0.413622332  0.22857623  0.13880197 -0.012050192 -0.01795532
    ## 4          4  0.060292715 -0.03274167  0.27140624  0.168841843  0.03941375
    ## 5          5 -0.016689490 -0.07407059  0.16758204  0.353928271  0.02393383
    ## 6          6 -0.195793493  0.02617203 -0.12954454 -0.051276286 -0.38019527
    ## 7          7 -0.421677409  0.24766154 -0.18656484  0.234755985  0.01911424
    ## 8          8  0.118429182 -0.01726759  0.13599084  0.055245280 -0.07162789
    ## 9          9  0.212682581  0.24174496  0.20287360 -0.055120798 -0.07068451
    ## 10        10  0.253795367  0.16866001  0.08460577 -0.163236994 -0.04209945
    ## 11        11  0.142704032  0.11171729  0.15670335 -0.115242647 -0.10433560
    ## 12        12  0.111887037 -0.09328665 -0.15520266 -0.068859852 -0.04948825
    ## 13        13  0.158253564  0.03657408 -0.13760088 -0.093219604 -0.01308640
    ## 14        14  0.064724932 -0.04551044 -0.25753062 -0.059515262  0.01367401
    ## 15        15  0.320919153  0.25757011  0.09603121  0.046764494 -0.03100038
    ## 16        16 -0.410345267  0.08886275 -0.07057320 -0.116917120  0.05283496
    ## 17        17  0.028080671 -0.19507616 -0.07894365 -0.088160440  0.16343869
    ## 18        18 -0.097668245 -0.23806692  0.05328197 -0.156391984  0.03001277
    ## 19        19  0.035013626 -0.30006857 -0.06081734 -0.048024954 -0.10385080
    ## 20        20  0.258535328  0.16772225 -0.20421316  0.029946793  0.19754498
    ## 21        21  0.183905273 -0.06236805 -0.13580028  0.006839604  0.11706442
    ## 22        22 -0.032086419 -0.20082351  0.14149033 -0.035001380  0.10349579
    ## 23        23 -0.132982820 -0.13838964  0.14982304 -0.022048395  0.02644282
    ## 24        24 -0.434781243  0.11373456  0.04674017 -0.145325831  0.14288784
    ## 25        25  0.013628169 -0.24560600  0.17569592  0.017374818  0.03273204
    ##    region field_type yr_since
    ## 1      BM   restored       16
    ## 2      BM   restored        3
    ## 3      FG       corn        -
    ## 4      FG    remnant        +
    ## 5      FG   restored       15
    ## 6      FL       corn        -
    ## 7      FL       corn        -
    ## 8      FL    remnant        +
    ## 9      FL   restored       40
    ## 10     FL   restored       36
    ## 11     FL   restored       35
    ## 12     FL   restored       10
    ## 13     FL   restored       10
    ## 14     FL   restored       10
    ## 15     BM   restored       28
    ## 16     LP       corn        -
    ## 17     LP    remnant        +
    ## 18     LP   restored        4
    ## 19     LP   restored        4
    ## 20     BM    remnant        +
    ## 21     BM   restored       18
    ## 22     BM   restored        7
    ## 23     BM   restored        2
    ## 24     BM       corn        -
    ## 25     BM   restored       11
    ## 
    ## $broken_stick_plot

<img src="microbial_communities_files/figure-gfm/pcoa_amf_otu-1.png" style="display: block; margin: auto;" />

    ## 
    ## $permanova
    ## Permutation test for adonis under reduced model
    ## Terms added sequentially (first to last)
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1999
    ## 
    ## adonis2(formula = d ~ field_type, data = env, permutations = nperm, strata = region)
    ##            Df SumOfSqs      R2      F Pr(>F)   
    ## field_type  2   1.0913 0.24796 3.6268  0.003 **
    ## Residual   22   3.3100 0.75204                 
    ## Total      24   4.4013 1.00000                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Four axes are significant by a broken stick model, between them
explaining 68.6% of the variation in AMF among fields. It may be
worthwhile to examine structure on Axes 3 and 4 sometime. The most
substantial variation here is on the first axis (27.4%) with Axis 2
explaining 17.9% of the variation in AMF abundances. Testing the design
factor *field_type* (with *region* treated as a block using the `strata`
argument of `adonis2`) revealed a significant clustering
$(R^2=0.25, p=0.003)$.

Let’s view a plot with abundances of community subgroups inset.

``` r
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
```

``` r
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
```

<img src="microbial_communities_files/figure-gfm/amf_families_fig-1.png" style="display: block; margin: auto;" />

Community trajectories revealed in the ordination correlate with field
type. Corn fields stand well apart with AMF communities, with restored
and remnant fields clustering closer than we had seen with
ITS-identified fungi. Restoration age along Axis 1 follows a near-linear
progression in Blue Mounds fields; with Fermi, we see a weaker age
progression and instead a strong separation between “ring fields” and
switchgrass plots as before. Restored fields’ fidelity to remnants seems
stronger with AMF than we had seen with general fungi.

What’s becoming apparent here is that Axis 1 separates strongly on
*field_type* and years since restoration, and Axis 2 further separates
on years since restoration. A consistent signal of region isn’t obvious.

Let’s test the relationship between age and community axis scores with
restored fields only.

``` r
amf_resto_scores <-
    pcoa_amf_bray$site_vectors %>%
    filter(field_type == "restored") %>%
    mutate(yr_since = as.numeric(yr_since))
```

``` r
summary(lm(Axis.1 ~ yr_since,
           data = amf_resto_scores))
```

    ## 
    ## Call:
    ## lm(formula = Axis.1 ~ yr_since, data = amf_resto_scores)
    ## 
    ## Residuals:
    ##       Min        1Q    Median        3Q       Max 
    ## -0.113666 -0.075444  0.005346  0.069948  0.134159 
    ## 
    ## Coefficients:
    ##              Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept) -0.035169   0.035297  -0.996 0.335975    
    ## yr_since     0.007926   0.001789   4.429 0.000572 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.08675 on 14 degrees of freedom
    ## Multiple R-squared:  0.5836, Adjusted R-squared:  0.5538 
    ## F-statistic: 19.62 on 1 and 14 DF,  p-value: 0.0005717

``` r
summary(lm(Axis.2 ~ yr_since,
           data = amf_resto_scores))
```

    ## 
    ## Call:
    ## lm(formula = Axis.2 ~ yr_since, data = amf_resto_scores)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.14868 -0.06036 -0.03930  0.06613  0.26559 
    ## 
    ## Coefficients:
    ##              Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept) -0.236024   0.046669  -5.057 0.000175 ***
    ## yr_since     0.012646   0.002366   5.345 0.000103 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.1147 on 14 degrees of freedom
    ## Multiple R-squared:  0.6711, Adjusted R-squared:  0.6476 
    ## F-statistic: 28.57 on 1 and 14 DF,  p-value: 0.0001034

``` r
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
```

<img src="microbial_communities_files/figure-gfm/amf_yrs_scores_fig-1.png" style="display: block; margin: auto;" />

Both axes correlate significantly and strongly with years since
restoration. Axis 2 shows a stronger relationship
$(R^2_{Adj}=0.61, p<0.001)$, and Axis 1 shows a moderately strong
relationship $(R^2_{Adj}=0.46, p<0.005)$

### PCoA with 18S gene, UNIFRAC distance

``` r
(pcoa_amf_uni <- pcoa_fun(distab$amf_uni, df_name = "18S gene, 97% OTU, UNIFRAC distance", corr = "lingoes"))
```

    ## $dataset
    ## [1] "18S gene, 97% OTU, UNIFRAC distance"
    ## 
    ## $components_exceed_broken_stick
    ## [1] 3
    ## 
    ## $correction_note
    ## [1] "Lingoes correction applied to negative eigenvalues: D' = -0.5*D^2 - 0.00730538070686336 , except diagonal elements"
    ## 
    ## $values
    ##   Dim Eigenvalues   Corr_eig Rel_corr_eig Broken_stick Cum_corr_eig
    ## 1   1  0.08910476 0.09641014    0.2302806   0.16236050    0.2302806
    ## 2   2  0.05546845 0.06277383    0.1499385   0.11888224    0.3802192
    ## 3   3  0.03768300 0.04498838    0.1074571   0.09714311    0.4876762
    ## 4   4  0.02287671 0.03018210    0.0720915   0.08265036    0.5597677
    ##   Cum_br_stick
    ## 1    0.1623605
    ## 2    0.2812427
    ## 3    0.3783858
    ## 4    0.4610362
    ## 
    ## $eigenvalues
    ## [1] 23 15
    ## 
    ## $site_vectors
    ##    field_key       Axis.1       Axis.2       Axis.3 region field_type yr_since
    ## 1          1  0.005000811 -0.019022480 -0.075661370     BM   restored       16
    ## 2          2  0.184329497 -0.015254402 -0.007699769     BM   restored        3
    ## 3          3 -0.072487397  0.067424513  0.039009283     FG       corn        -
    ## 4          4 -0.041221010 -0.036068309  0.021185618     FG    remnant        +
    ## 5          5  0.024753793 -0.052063386  0.054757181     FG   restored       15
    ## 6          6  0.132027889  0.107077144  0.033575362     FL       corn        -
    ## 7          7 -0.042177309  0.100432373 -0.083348303     FL       corn        -
    ## 8          8  0.012829998 -0.029240145  0.019755980     FL    remnant        +
    ## 9          9 -0.048029289 -0.037910042  0.024674421     FL   restored       40
    ## 10        10 -0.047531002 -0.032394265  0.016179264     FL   restored       36
    ## 11        11 -0.032019394 -0.015349526  0.017878940     FL   restored       35
    ## 12        12  0.070369149 -0.022703786 -0.026097774     FL   restored       10
    ## 13        13  0.014601133 -0.009117052 -0.041202654     FL   restored       10
    ## 14        14  0.027018121 -0.003352651 -0.048179202     FL   restored       10
    ## 15        15 -0.040888003 -0.041071752 -0.022777348     BM   restored       28
    ## 16        16 -0.029033616  0.101034837 -0.007893224     LP       corn        -
    ## 17        17 -0.010295789 -0.010944908  0.003923329     LP    remnant        +
    ## 18        18  0.001881795  0.008404165  0.035231885     LP   restored        4
    ## 19        19  0.066396377 -0.002797405  0.021616717     LP   restored        4
    ## 20        20 -0.045321553 -0.044288430 -0.057552793     BM    remnant        +
    ## 21        21 -0.010605035 -0.038274923 -0.046884326     BM   restored       18
    ## 22        22 -0.024417886 -0.009405314  0.017526673     BM   restored        7
    ## 23        23 -0.029424159  0.005009133  0.025573301     BM   restored        2
    ## 24        24 -0.076762389  0.062360323  0.037484938     BM       corn        -
    ## 25        25  0.011005267 -0.032483714  0.048923871     BM   restored       11
    ## 
    ## $broken_stick_plot

<img src="microbial_communities_files/figure-gfm/pcoa_amf_uni-1.png" style="display: block; margin: auto;" />

    ## 
    ## $permanova
    ## Permutation test for adonis under reduced model
    ## Terms added sequentially (first to last)
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1999
    ## 
    ## adonis2(formula = d ~ field_type, data = env, permutations = nperm, strata = region)
    ##            Df SumOfSqs      R2      F Pr(>F)   
    ## field_type  2 0.054762 0.22505 3.1945  0.003 **
    ## Residual   22 0.188572 0.77495                 
    ## Total      24 0.243335 1.00000                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Three axes are significant by a broken stick model, between them
explaining 48.8% of the variation in AMF among fields. The most
substantial variation here is on the first axis (23%) with Axis 2
explaining 15% of the variation in AMF abundances. Testing the design
factor *field_type* (with *region* treated as a block using the `strata`
argument of `adonis2`) revealed a significant clustering
$(R^2=0.23, p=0.003)$.

Let’s view a plot with abundances of community subgroups inset.

``` r
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
```

``` r
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
```

<img src="microbial_communities_files/figure-gfm/amf_uni_families_fig-1.png" style="display: block; margin: auto;" />

Community trajectories revealed in the ordination separate cornfields
from everything else. Using UNIFRAC distance has really dissolved most
of what was apparent with the Bray-Curtis distance.  
Corn fields stand well apart with AMF communities, but no signal appears
for other field types or for years since restoration. I guess what this
shows is that for AMF, restored fields almost immediately resemble
remnants (but there must be some outlier taxa in Eric Rahnheim’s place).

What’s becoming apparent here is that Axis 1 separates strongly on
*field_type* and years since restoration, and Axis 2 further separates
on years since restoration. A consistent signal of region isn’t obvious.

Let’s test the relationship between age and community axis scores with
restored fields only. I don’t expect much.

``` r
amf_uni_resto_scores <-
    pcoa_amf_uni$site_vectors %>%
    filter(field_type == "restored") %>%
    mutate(yr_since = as.numeric(yr_since))
```

``` r
summary(lm(
    Axis.1 ~ yr_since,
    data = amf_uni_resto_scores
))
```

    ## 
    ## Call:
    ## lm(formula = Axis.1 ~ yr_since, data = amf_uni_resto_scores)
    ## 
    ## Residuals:
    ##       Min        1Q    Median        3Q       Max 
    ## -0.077622 -0.015330 -0.003244  0.011231  0.138891 
    ## 
    ## Coefficients:
    ##              Estimate Std. Error t value Pr(>|t|)  
    ## (Intercept)  0.053716   0.020091   2.674   0.0182 *
    ## yr_since    -0.002759   0.001019  -2.709   0.0170 *
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.04938 on 14 degrees of freedom
    ## Multiple R-squared:  0.3439, Adjusted R-squared:  0.297 
    ## F-statistic: 7.338 on 1 and 14 DF,  p-value: 0.01696

``` r
summary(lm(
    Axis.2 ~ yr_since,
    data = amf_uni_resto_scores
))
```

    ## 
    ## Call:
    ## lm(formula = Axis.2 ~ yr_since, data = amf_uni_resto_scores)
    ## 
    ## Residuals:
    ##       Min        1Q    Median        3Q       Max 
    ## -0.032677 -0.008332  0.002915  0.008418  0.020946 
    ## 
    ## Coefficients:
    ##               Estimate Std. Error t value Pr(>|t|)  
    ## (Intercept) -0.0067045  0.0060005  -1.117   0.2827  
    ## yr_since    -0.0008454  0.0003042  -2.779   0.0148 *
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.01475 on 14 degrees of freedom
    ## Multiple R-squared:  0.3555, Adjusted R-squared:  0.3095 
    ## F-statistic: 7.723 on 1 and 14 DF,  p-value: 0.01478

``` r
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
```

<img src="microbial_communities_files/figure-gfm/amf_uni_yrs_scores_fig-1.png" style="display: block; margin: auto;" />

Both axes correlate significantly but with less than moderate strength
with years since restoration. Axis 2 again shows a stronger relationship
$(R^2_{Adj}=0.31, p<0.05)$, and Axis 1 is close with
$(R^2_{Adj}=0.30, p<0.05)$
