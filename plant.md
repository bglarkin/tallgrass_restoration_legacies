Plant data: communities and traits
================
Beau Larkin

Last updated: 23 October, 2023

- [Description](#description)
- [Packages and libraries](#packages-and-libraries)
- [Functions](#functions)
- [Data](#data)
  - [Data wrangling: traits data](#data-wrangling-traits-data)
    - [Abundance data (16 sites)](#abundance-data-16-sites)
    - [Presence data (20 sites)](#presence-data-20-sites)
- [Results](#results)
  - [Visualization and output of abundance
    data](#visualization-and-output-of-abundance-data)
  - [Visualization and output of presence
    data](#visualization-and-output-of-presence-data)

# Description

Plant data comprises two separate data sets. Mike Healy did quadrat
surveys at all sites except Fermi, recording plant abundance. At Fermi,
we only have relevé data with presence/absence. These data were provided
by Mike Miller. Wisconsin sites were surveyed in Aug-Sept 2016. Fermi
sites were surveyed in summer 2017

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
packages_needed = c("GGally", "tidyverse", "vegan", "colorspace")
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

Cleanplot PCA produces informative visualizations of PCA ordinations
[(Borcard et
al. 2018)](http://link.springer.com/10.1007/978-3-319-71404-2)

``` r
source(paste0(getwd(), "/supporting_files/cleanplot_pca.txt"))
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

## Visualization and output of abundance data

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

## Visualization and output of presence data

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
write_csv(p_ab_trait, paste0(getwd(), "/clean_data/plant_trait_abund.csv"))
```
