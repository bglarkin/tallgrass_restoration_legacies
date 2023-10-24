Plant data: communities and traits
================
Beau Larkin

Last updated: 24 October, 2023

- [Description](#description)
- [Packages and libraries](#packages-and-libraries)
- [Functions](#functions)
  - [Cleanplot PCA](#cleanplot-pca)
  - [PCoA](#pcoa)
- [Data](#data)
  - [Data wrangling: traits data](#data-wrangling-traits-data)
    - [Abundance data (16 sites)](#abundance-data-16-sites)
    - [Presence data (20 sites)](#presence-data-20-sites)
- [Results](#results)
  - [Traits: sites with abundance
    data](#traits-sites-with-abundance-data)
  - [Traits: sites with presence data](#traits-sites-with-presence-data)
  - [Plant Communities](#plant-communities)

# Description

Plant data comprises two data sets. Mike Healy did quadrat surveys at
all sites except Fermi, recording plant abundance. At Fermi, we only
have relevé data with presence/absence. These data were provided by Mike
Miller. Wisconsin sites were surveyed in Aug-Sept 2016. Fermi sites were
surveyed in summer 2017

Plant metadata includes taxonomy and life history traits, and should
cover both plant data sets. With the abundance-data sites, trait data
are reported in percent cover. In the presence-data sites, traits are in
counts of species with that trait per field.

This script produces basic visualization and diagnostic views of the
plant data. Two traits matrices are output, one for sites with abundance
data (16 sites), and one for sites with presence data only (20 sites).

Plant data may be used in matrix form for ordination, etc. Abundance and
presence/absence matrices are available. The data may also be transposed
and summarized as abundance in taxonomic or life history classes. This
makes sense for the abundance data, but maybe makes less sense with the
presence/absence data because the summary would really be richness, not
abundance.

Plants data were obtained in the [TRY Plant Trait
Database](https://www.try-db.org/TryWeb/Home.php) ([Kattge et
al. 2010](https://besjournals.onlinelibrary.wiley.com/doi/10.1111/j.2041-210X.2010.00067.x),
[Kattge et
al. 2011](https://onlinelibrary.wiley.com/doi/10.1111/j.1365-2486.2011.02451.x))
in 2016.

# Packages and libraries

``` r
packages_needed = c("GGally", "tidyverse", "vegan", "colorspace", "ape")
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

# Functions

## Cleanplot PCA

Cleanplot PCA produces informative visualizations of PCA ordinations
[(Borcard et
al. 2018)](http://link.springer.com/10.1007/978-3-319-71404-2)

``` r
source(paste0(getwd(), "/supporting_files/cleanplot_pca.txt"))
```

## PCoA

``` r
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
```

# Data

Plant community data includes:

- Metadata, taxonomy and traits
- Abundance data, surveyed in 2016, sites limited to Wisconsin only
- Presence data, from Fermi in 2015, all other sites converted to
  presence data, does not include Fermi switchgrass or corn fields.

``` r
plant <- list(
    meta = read_csv(paste0(getwd(), "/clean_data/spe_plant_meta.csv"), show_col_types = FALSE) %>% 
        rename_with(tolower),
    ab   = read_csv(paste0(getwd(), "/clean_data/spe_plant_abund.csv"), show_col_types = FALSE) %>% 
        rename(field_name = SITE),
    pr   = read_csv(paste0(getwd(), "/clean_data/spe_plant_presence.csv"), show_col_types = FALSE) %>% 
        rename(field_name = SITE)
)
```

Metadata from sites, as in previous

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
```

``` r
sites_noc <- sites %>% 
    filter(field_type != "corn")
```

Distance matrices are needed for ordinations of the plant data.
Bray-Curtis distance is used for abundance data, the Jaccard similarity
is used for binary data.

``` r
distab = list(
    p_ab = vegdist(data.frame(plant$ab, row.names = 1), method = "bray"),
    p_pr = vegdist(data.frame(plant$pr, row.names = 1), method = "jac", binary = TRUE)
)
```

## Data wrangling: traits data

### Abundance data (16 sites)

The alignment of plant codes among files must be confirmed. Show any
mismatched codes.

``` r
which(!(colnames(plant$ab[, -1]) %in% plant$meta$code))
```

    ## integer(0)

``` r
which(!(colnames(plant$pr[, -1]) %in% plant$meta$code))
```

    ## integer(0)

The codes are aligned. Factor levels in metadata must be cleaned up and
simplified.

``` r
recode_group <- c(`annual/biennial` = "biennial", `annual/perennial` = "perennial", `C4-grass` = "C4_grass", `C3-grass` = "C3_grass")
filter_group <- c("bare", "litter", "unknown", "fern", "moss", "sedge", "Unknown", "shrub", "tree", "C3/C4", "rush")
```

Plant traits are summarized among the set of sites with plant abundance
data

``` r
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
```

### Presence data (20 sites)

``` r
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
```

# Results

## Traits: sites with abundance data

Run a PCA on chord-transformed traits data from sites with abundance
data, perform typical basic diagnostics. This should be done without
corn fields because they exert too strong a difference on everything
else.

``` r
p_ab_trait_ch <- decostand(data.frame(p_ab_trait %>% filter(field_name %in% sites_noc$field_name), row.names = 1), "normalize")
p_ab_trait_pca <- rda(p_ab_trait_ch)
p_ab_trait_pca %>% summary(., display = NULL)
```

    ## 
    ## Call:
    ## rda(X = p_ab_trait_ch) 
    ## 
    ## Partitioning of variance:
    ##               Inertia Proportion
    ## Total         0.07381          1
    ## Unconstrained 0.07381          1
    ## 
    ## Eigenvalues, and their contribution to the variance 
    ## 
    ## Importance of components:
    ##                           PC1     PC2      PC3      PC4      PC5       PC6
    ## Eigenvalue            0.04172 0.01537 0.008637 0.005015 0.001408 0.0008968
    ## Proportion Explained  0.56521 0.20826 0.117010 0.067942 0.019081 0.0121501
    ## Cumulative Proportion 0.56521 0.77347 0.890479 0.958420 0.977501 0.9896511
    ##                             PC7       PC8       PC9      PC10
    ## Eigenvalue            0.0004438 0.0001811 0.0001265 1.243e-05
    ## Proportion Explained  0.0060131 0.0024530 0.0017144 1.684e-04
    ## Cumulative Proportion 0.9956642 0.9981172 0.9998316 1.000e+00
    ## 
    ## Scaling 2 for species and site scores
    ## * Species are scaled proportional to eigenvalues
    ## * Sites are unscaled: weighted dispersion equal on all dimensions
    ## * General scaling constant of scores:

``` r
screeplot(p_ab_trait_pca, bstick = TRUE)
```

<img src="plant_files/figure-gfm/p_ab_trait_pca_scree-1.png" style="display: block; margin: auto;" />

``` r
cleanplot.pca(p_ab_trait_pca)
```

<img src="plant_files/figure-gfm/p_ab_trait_pca_cleanplot-1.png" style="display: block; margin: auto;" /><img src="plant_files/figure-gfm/p_ab_trait_pca_cleanplot-2.png" style="display: block; margin: auto;" />

Axis 1 & 2 explain 77% of the variation, and both eigenvalues exceed the
broken stick model. Traits forb, perennial, annual, nonnative, and C4
grass exceed the unit circle, suggesting a strong correlation with site
differences. Traits appear collinear, explore which ones produce high
VIF.

``` r
sort(diag(solve(cor(data.frame(p_ab_trait, row.names = 1)))), decreasing = TRUE) 
```

    ##   perennial      annual      native   nonnative    C4_grass        forb 
    ## 1091.599068 1070.069928  752.469656  655.604119  113.993616   93.683065 
    ##    biennial   shrubTree      legume    C3_grass 
    ##   49.116860    7.197170    4.623733    4.462052

Many are very high. Traits levels need not be mutually exclusive, but it
appears here that they are. Native opposes nonnative abundance,
perennial opposes annual, and forb opposes C4 grass. The last case may
be driven by Karla Ott’s field, but in other cases. What happens to VIF
when one of these opposing factor levels are removed?

``` r
sort(diag(solve(cor(data.frame(p_ab_trait, row.names = 1) %>% select(-annual, -nonnative, -C4_grass)))), decreasing = TRUE)
```

    ##    native perennial      forb  biennial shrubTree    legume  C3_grass 
    ##  9.315159  6.395496  5.055116  2.039153  1.801386  1.757270  1.412914

It’s clear that some are correlated and there is good reason to remove
them. It’s useful here because even after these are removed, they
suggest their opposite factor level. This will help inform forward
selection later. Export the traits matrix for sites with abundance data:

``` r
write_csv(p_ab_trait, paste0(getwd(), "/clean_data/plant_trait_abund.csv"))
```

## Traits: sites with presence data

Run a PCA on chord-transformed traits data from sites with abundance
data, perform typical basic diagnostics.

``` r
p_pr_trait_ch <- decostand(data.frame(p_pr_trait %>% filter(field_name %in% sites_noc$field_name), row.names = 1), "normalize")
p_pr_trait_pca <- rda(p_pr_trait_ch)
p_pr_trait_pca %>% summary(., display = NULL)
```

    ## 
    ## Call:
    ## rda(X = p_pr_trait_ch) 
    ## 
    ## Partitioning of variance:
    ##               Inertia Proportion
    ## Total         0.02141          1
    ## Unconstrained 0.02141          1
    ## 
    ## Eigenvalues, and their contribution to the variance 
    ## 
    ## Importance of components:
    ##                           PC1      PC2      PC3      PC4       PC5       PC6
    ## Eigenvalue            0.01349 0.003402 0.002393 0.001243 0.0003813 0.0002555
    ## Proportion Explained  0.63016 0.158910 0.111772 0.058062 0.0178111 0.0119346
    ## Cumulative Proportion 0.63016 0.789070 0.900842 0.958904 0.9767150 0.9886496
    ##                             PC7       PC8       PC9      PC10
    ## Eigenvalue            0.0001259 6.532e-05 0.0000496 2.145e-06
    ## Proportion Explained  0.0058828 3.051e-03 0.0023166 1.002e-04
    ## Cumulative Proportion 0.9945323 9.976e-01 0.9998998 1.000e+00
    ## 
    ## Scaling 2 for species and site scores
    ## * Species are scaled proportional to eigenvalues
    ## * Sites are unscaled: weighted dispersion equal on all dimensions
    ## * General scaling constant of scores:

``` r
screeplot(p_pr_trait_pca, bstick = TRUE)
```

<img src="plant_files/figure-gfm/p_pr_trait_pca_scree-1.png" style="display: block; margin: auto;" />

``` r
cleanplot.pca(p_pr_trait_pca)
```

<img src="plant_files/figure-gfm/p_pr_trait_pca_cleanplot-1.png" style="display: block; margin: auto;" /><img src="plant_files/figure-gfm/p_pr_trait_pca_cleanplot-2.png" style="display: block; margin: auto;" />

Axis 1 & 2 explain 79% of the variation. The axis 1 eigenvalue is the
only one which exceeds a broken stick model. Traits forb, biennial,
perennial, native, and nonnative exceed the unit circle, suggesting a
strong correlation with site differences. Traits appear collinear,
explore which ones produce high VIF.

``` r
sort(diag(solve(cor(data.frame(p_pr_trait, row.names = 1)))), decreasing = TRUE) 
```

    ##    perennial       native         forb    nonnative     biennial       annual 
    ## 12225.829278  9552.497611   760.415589   111.667352    82.714461    43.158554 
    ##     C3_grass       legume     C4_grass    shrubTree 
    ##    10.752644     9.665848     8.332982     7.784560

Many are very high, and lie in opposition by factor levels as in the
abundance data, but possibly a little less strong in terms of pure
linear correlation. Biennial probably describes forbs, and can easily be
discarded. It also looks like most perennials are also native?

``` r
sort(diag(solve(cor(data.frame(p_pr_trait, row.names = 1) %>% select(-nonnative, -biennial, -native)))), decreasing = TRUE)
```

    ##  perennial       forb   C4_grass     annual     legume  shrubTree   C3_grass 
    ## 185.188816 138.977387   6.865699   4.597058   3.778917   2.424485   2.070412

It’s clear that some are correlated and there is good reason to remove
them. It’s useful here because even after these are removed, they
suggest their opposite factor level. This will help inform forward
selection later. Export the traits matrix for sites with abundance data:

``` r
write_csv(p_pr_trait, paste0(getwd(), "/clean_data/plant_trait_presence.csv"))
```

## Plant Communities

In restored fields, plant communities don’t reflect natural community
assembly. Still, it’s useful to examine an ordination of sites to
develop an understanding of how they differ. Plant traits data are
probably more useful to looking at development of plant communities over
time after restoration. Traits may be filtered more than species, and
species’ occurrence may not be uniform across sites, though a species’
realized niche may include sites where it is not found. \### Sites with
abundance data An ordiation is run on plant abundance data using
`pcoa_fun()`.

``` r
(pcoa_ab <- pcoa_fun(plant$ab, distab$p_ab, corr="lingoes", df_name = "plant abundance data, 16 sites"))
```

    ## $dataset
    ## [1] "plant abundance data, 16 sites"
    ## 
    ## $components_exceed_broken_stick
    ## [1] 1
    ## 
    ## $correction_note
    ## [1] "Lingoes correction applied to negative eigenvalues: D' = -0.5*D^2 - 0.000809415888073743 , except diagonal elements"
    ## 
    ## $values
    ##   Dim Eigenvalues  Corr_eig Rel_corr_eig Broken_stick Cum_corr_eig Cum_br_stick
    ## 1   1   1.6884108 1.6892202    0.3128013    0.2322545    0.3128013    0.2322545
    ## 2   2   0.6812888 0.6820982    0.1263075    0.1608259    0.4391088    0.3930803
    ## 3   3   0.5707939 0.5716033    0.1058466    0.1251116    0.5449554    0.5181919
    ## 
    ## $eigenvalues
    ## [1] 31.3 12.6
    ## 
    ## $site_vectors
    ##    field_name      Axis.1      Axis.2 field_key region field_type yr_since
    ## 1       BBRP1  0.24056775 -0.16522640         1     BM   restored       16
    ## 2       ERRP1  0.08413967  0.41643560         2     BM   restored        3
    ## 3        FGC1 -0.66433300 -0.03290989         3     FG       corn        -
    ## 4      FGREM1  0.04431444 -0.01450857         4     FG    remnant        +
    ## 5       FGRP1  0.12276250 -0.15986548         5     FG   restored       15
    ## 6       KORP1  0.15461767 -0.26031738        15     BM   restored       28
    ## 7        LPC1 -0.66168712 -0.03298898        16     LP       corn        -
    ## 8      LPREM1  0.15732628  0.05593888        17     LP    remnant        +
    ## 9       LPRP1  0.20256709  0.11248606        18     LP   restored        4
    ## 10      LPRP2  0.22751545 -0.02626885        19     LP   restored        4
    ## 11     MBREM1  0.17430970 -0.33069721        20     BM    remnant        +
    ## 12      MBRP1  0.18213508 -0.20141260        21     BM   restored       18
    ## 13      MHRP1  0.18090418  0.03336249        22     BM   restored        7
    ## 14      MHRP2  0.04005768  0.39612929        23     BM   restored        2
    ## 15       PHC1 -0.67283397 -0.03367491        24     BM       corn        -
    ## 16      PHRP1  0.18763660  0.24351797        25     BM   restored       11
    ## 
    ## $broken_stick_plot

![](plant_files/figure-gfm/pcoa_ab-1.png)<!-- -->

    ## 
    ## $permanova
    ## Permutation test for adonis under reduced model
    ## Terms added sequentially (first to last)
    ## Blocks:  region 
    ## Plots: field_name, plot permutation: free
    ## Permutation: none
    ## Number of permutations: 1999
    ## 
    ## adonis2(formula = d ~ field_type, data = env_w, permutations = h)
    ##            Df SumOfSqs      R2      F Pr(>F)   
    ## field_type  2   1.9713 0.36586 3.7501  0.005 **
    ## Residual   13   3.4169 0.63414                 
    ## Total      15   5.3882 1.00000                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Axis 1 explains 31.3% of the variation and is the only eigenvalue that
exceeds a broken stick model. The most substantial variation here will
be on the first axis. Axis 2 explains 12.6% of the variation and was not
very close to the broken stick value. Testing the design factor
*field_type* (with *region* treated as a block using arguments to
`how()` revealed a significant clustering $(R^2=0.37,~p=0.005)$. Let’s
view a plot of these results.

``` r
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
```

<img src="plant_files/figure-gfm/pcoa_ab_fig-1.png" style="display: block; margin: auto;" />

``` r
#" ### Sites with presence data
```

An ordiation is run on plant presence data using `pcoa_fun()`. The
dataset includes 20 sites. This analysis isn’t appropriate because the
blocks are unbalanced (no cornfield data from Fermi), but it still shows
differences with plant data.

``` r
(pcoa_pr <- pcoa_fun(plant$pr, distab$p_pr, corr="none", df_name = "plant presence data, 20 sites"))
```

    ## $dataset
    ## [1] "plant presence data, 20 sites"
    ## 
    ## $components_exceed_broken_stick
    ## [1] 2
    ## 
    ## $correction_note
    ## [1] "There were no negative eigenvalues. No correction was applied"
    ## 
    ## $values
    ##   Dim Eigenvalues Relative_eig Broken_stick Cumul_eig Cumul_br_stick
    ## 1   1   1.4264339   0.19876306    0.1867231 0.1987631      0.1867231
    ## 2   2   1.2066484   0.16813758    0.1340916 0.3669006      0.3208147
    ## 3   3   0.6135197   0.08548946    0.1077758 0.4523901      0.4285905
    ## 
    ## $eigenvalues
    ## [1] 19.9 16.8
    ## 
    ## $site_vectors
    ##    field_name       Axis.1      Axis.2 field_key region field_type yr_since
    ## 1       BBRP1  0.007406404  0.24131308         1     BM   restored       16
    ## 2       ERRP1  0.067469134  0.21951909         2     BM   restored        3
    ## 3        FGC1  0.512236442 -0.22516150         3     FG       corn        -
    ## 4      FGREM1 -0.057759438 -0.02094116         4     FG    remnant        +
    ## 5       FGRP1 -0.035178611  0.12877987         5     FG   restored       15
    ## 6      FLREM1 -0.346670841 -0.38303529         8     FL    remnant        +
    ## 7       FLRP1 -0.367062459 -0.37699824         9     FL   restored       40
    ## 8       FLRP4 -0.366696539 -0.38481716        10     FL   restored       36
    ## 9       FLRP5 -0.321881936 -0.25314288        11     FL   restored       35
    ## 10      KORP1  0.054932725  0.16097174        15     BM   restored       28
    ## 11       LPC1  0.582345898 -0.29432253        16     LP       corn        -
    ## 12     LPREM1 -0.022447100  0.25575985        17     LP    remnant        +
    ## 13      LPRP1 -0.024723803  0.26041246        18     LP   restored        4
    ## 14      LPRP2 -0.099697392  0.22860489        19     LP   restored        4
    ## 15     MBREM1 -0.028171733  0.07734331        20     BM    remnant        +
    ## 16      MBRP1 -0.072344127  0.13855455        21     BM   restored       18
    ## 17      MHRP1  0.021730449  0.21161280        22     BM   restored        7
    ## 18      MHRP2  0.002954358  0.12383740        23     BM   restored        2
    ## 19       PHC1  0.547087359 -0.32089787        24     BM       corn        -
    ## 20      PHRP1 -0.053528790  0.21260759        25     BM   restored       11
    ## 
    ## $broken_stick_plot

![](plant_files/figure-gfm/pcoa_pr-1.png)<!-- -->

    ## 
    ## $permanova
    ## Permutation test for adonis under reduced model
    ## Terms added sequentially (first to last)
    ## Blocks:  region 
    ## Plots: field_name, plot permutation: free
    ## Permutation: none
    ## Number of permutations: 1999
    ## 
    ## adonis2(formula = d ~ field_type, data = env_w, permutations = h)
    ##            Df SumOfSqs    R2      F Pr(>F)    
    ## field_type  2   1.6722 0.233 2.5822  5e-04 ***
    ## Residual   17   5.5044 0.767                  
    ## Total      19   7.1766 1.000                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Axis 1 explains 19.9% of the variation and axis 2 explains 12.6% of the
variation. These two eigenvalues exceed the broken stick value. stick
value. Testing the design factor *field_type* (with *region* treated as
a block using arguments to `how()` revealed a significant clustering
$(R^2=0.23,~p=5\times 10^{-4})$. Let’s view a plot of these results.

``` r
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
```

<img src="plant_files/figure-gfm/pcoa_pr_fig-1.png" style="display: block; margin: auto;" />

The regional signal is most obvious here.
