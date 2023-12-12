Constrained and summary analysis
================
Beau Larkin

Last updated: 12 December, 2023

- [Description](#description)
- [Packages and libraries](#packages-and-libraries)
- [Functions](#functions)
  - [PCoA axes](#pcoa-axes)
  - [Constrained analyses](#constrained-analyses)
  - [Plot results](#plot-results)
- [Data](#data)
  - [Site metadata and experimental
    design](#site-metadata-and-experimental-design)
  - [Plant data](#plant-data)
    - [Traits](#traits)
    - [Plant communities](#plant-communities)
  - [Environmental data](#environmental-data)
  - [Microbial data](#microbial-data)
    - [Fungal communities](#fungal-communities)
    - [Species metadata](#species-metadata)
  - [Response data](#response-data)
    - [Fungal biomass](#fungal-biomass)
    - [Water stable aggregates](#water-stable-aggregates)
- [Analysis and Results](#analysis-and-results)
  - [Plant community axes](#plant-community-axes)
  - [Microbial communities with explanatory and
    covariables](#microbial-communities-with-explanatory-and-covariables)
    - [General fungal community (ITS sequence
      abundance)](#general-fungal-community-its-sequence-abundance)
    - [AMF community (18S sequence
      abundance)](#amf-community-18s-sequence-abundance)
  - [Effects of soil and plants](#effects-of-soil-and-plants)
  - [Response correlations](#response-correlations)

# Description

This script will combine and summarize the many single datasets
presented so far. Several multivariate analyses will be considered.
Symmetrical analyses will attempt to discern the relative contributions
of plant and soil data on site ordination. Asymmetric constrained
analyses will attempt to determine the most significant contributors of
fungal community difference.

Other data, like microbial biomass and soil water stable aggregates may
also be considered.

All explanatory data are provided as field-based averages. Using them to
constrain microbial community data will require alignment with
field-based summary data and cannot use subsamples of microbial data.

Note that the plant data are not available for all sites, leading to
some complication in the analysis. See further notes below.

# Packages and libraries

``` r
packages_needed = c("tidyverse", "vegan", "conflicted", "knitr", "GGally", "ape", "ade4", "GGally")
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

``` r
conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")
```

# Functions

## PCoA axes

Produce community axes for use in constrained analyses later.

``` r
pcoa_fun <- function(s, ft=c("restored"), rg, method="bray", binary=FALSE, corr="none") {
    d <- vegdist(
        data.frame(
            s %>% 
                filter(field_type %in% ft, region %in% rg) %>% 
                select(-field_type, -region) %>% 
                select(field_name, where(~ is.numeric(.) && sum(.) > 0)),
            row.names = 1),
        method = method,
        binary = binary)
    p <- pcoa(d, correction = corr)
    p_vals <- data.frame(p$values) %>% 
        rownames_to_column(var = "Dim") %>% 
        mutate(Dim = as.integer(Dim))
    p_vec <- data.frame(p$vectors)
    if(corr == "none" | ncol(p_vals) == 6) {
        p_ncomp <- with(p_vals, which(Relative_eig < Broken_stick)[1]-1)
    } else {
        p_ncomp <- with(p_vals, which(Rel_corr_eig < Broken_stick)[1]-1)
    }
    ncomp <- if(p_ncomp <= 2) {2} else {p_ncomp}
    # Ordination plot
    scores <- p_vec[, 1:2]
    # Output data
    output <- list(correction_note = p$note,
                   values          = p_vals[1:(ncomp+1), ],
                   site_vectors    = scores)
    return(output)
}
```

## Constrained analyses

This function performs a db-RDA on microbial data with soil/site
covariables and forward selects on plant communities, agricultural
nutrients, and soil carbon.

``` r
dbrda_fun <- function(s, pspe_pcoa="none", ft, rg) {
    fspe_bray <- vegdist(
        data.frame(
            s %>% 
                filter(field_type %in% ft, region %in% rg) %>% 
                select(-field_type, -region) %>% 
                select(field_name, where(~ is.numeric(.) && sum(.) > 0)),
            row.names = 1),
        method = "bray")
    if(is.data.frame(pspe_pcoa) == TRUE) {
        pspe_ax <- pspe_pcoa %>% 
            rownames_to_column("field_name") %>% 
            rename(plant1 = Axis.1, plant2 = Axis.2)
    }
    ptr_norm <- decostand(
        data.frame(
            ptr %>% 
                filter(field_type %in% ft, region %in% rg) %>% 
                select(-field_type, -region),
            row.names = 1),
        "normalize") %>% 
        rownames_to_column("field_name")
    soil_cov_z <- decostand(
        data.frame(
            soil %>% 
                filter(field_type %in% ft, region %in% rg) %>% 
                select(-field_type, -region, -OM, -NO3, -P, -K),
            row.names = 1),
        "standardize")
    soil_cov_sc <- data.frame(scores(rda(soil_cov_z), display = "sites")) %>% 
        rename(soil1 = PC1, soil2 = PC2) %>% 
        rownames_to_column("field_name")
    soil_expl_z <- decostand(
        data.frame(
            soil %>% 
                filter(field_type %in% ft, region %in% rg) %>% 
                select(field_name, OM, NO3, P, K, -field_type, -region),
            row.names = 1),
        "standardize") %>% 
        rownames_to_column("field_name")
    yr_z <- decostand(
        data.frame(
            sites %>% 
                filter(field_type %in% ft, region %in% rg) %>% 
                select(field_name, yr_since),
            row.names = 1),
        "standardize") %>% 
        rownames_to_column("field_name")
    # Create explanatory data frame and covariables matrix
    env <- if(is.data.frame(pspe_pcoa) == FALSE) {
        soil_cov_sc %>% 
            left_join(ptr_norm, by = join_by(field_name)) %>% 
            left_join(soil_expl_z, by = join_by(field_name)) %>% 
            left_join(yr_z, by = join_by(field_name)) %>% 
            column_to_rownames("field_name") %>% 
            drop_na()
    } else {
        soil_cov_sc %>% 
            left_join(pspe_ax, by = join_by(field_name)) %>% 
            left_join(soil_expl_z, by = join_by(field_name)) %>% 
            left_join(yr_z, by = join_by(field_name)) %>% 
            column_to_rownames("field_name") %>% 
            drop_na()
    }
    covars <- as.matrix(env[, 1:2])
    expl <- env[, 3:ncol(env)]
    # Forward select on explanatory data with covariables
    mod_null <- dbrda(fspe_bray ~ 1 + Condition(covars), data = expl, sqrt.dist = TRUE)
    mod_full <- dbrda(fspe_bray ~ . + Condition(covars), data = expl, sqrt.dist = TRUE)
    mod_step <- ordistep(mod_null, 
                         scope = formula(mod_full), 
                         direction = "forward", 
                         permutations = how(nperm = 1999), 
                         trace = FALSE)
    mod_r2   <- RsquareAdj(mod_step, permutations = 1999)
    mod_glax <- anova(mod_step, permutations = how(nperm = 1999))
    mod_inax <- anova(mod_step, by = "axis", permutations = how(nperm = 1999))
    # Produce plot data including borderline vars if possible
    mod_scor <- if(is.data.frame(pspe_pcoa) == FALSE) {
        scores(
            dbrda(fspe_bray ~ yr_since + forb + C4_grass + Condition(covars), data = expl, sqrt.dist = TRUE),
            choices = c(1,2),
            display = c("bp", "sites"), tidy = FALSE
        )
    } else {
        scores(
            dbrda(fspe_bray ~ yr_since + K + plant1 + OM + Condition(covars), data = expl, sqrt.dist = TRUE),
            choices = c(1,2),
            display = c("bp", "sites"), tidy = FALSE
        )
    }
    return(list(
        plot_data = mod_scor,
        select_mod = mod_step,
        R2 = mod_r2 %>% map(\(x) round(x, 2)),
        global_axis_test = mod_glax,
        individual_axis_test = mod_inax
    ))
} 
```

## Plot results

``` r
plot_dbrda <- function(site_sc, site_bp) {
    site_df <- site_sc %>% 
        data.frame() %>% 
        rownames_to_column(var = "field_name") %>% 
        left_join(sites, by = join_by(field_name))
    bp_df <- site_bp %>% 
        data.frame() %>% 
        rownames_to_column(var = "envvar") %>% 
        mutate(
            origin = 0,
            m = dbRDA2 / dbRDA1, 
            d = sqrt(dbRDA1^2 + dbRDA2^2), 
            dadd = sqrt((max(dbRDA1)-min(dbRDA2))^2 + (max(dbRDA2)-min(dbRDA2))^2)*0.1,
            labx = ((d+dadd)*cos(atan(m)))*(dbRDA1/abs(dbRDA1)), 
            laby = ((d+dadd)*sin(atan(m)))*(dbRDA1/abs(dbRDA1)))
    ggplot(site_df, aes(x = dbRDA1, y = dbRDA2)) +
        geom_text(aes(label = field_name)) +
        geom_segment(data = bp_df, 
                     aes(x = origin, xend = dbRDA1, y = origin, yend = dbRDA2), 
                     arrow = arrow(length = unit(3, "mm")),
                     color = "blue") +
        geom_text(data = bp_df, 
                  aes(x = labx, y = laby, label = envvar),
                  color = "blue") +
        theme_bw()
}
```

# Data

## Site metadata and experimental design

``` r
sites <-
    read_csv(paste0(getwd(), "/clean_data/sites.csv"), show_col_types = FALSE) %>%
    mutate(
        field_type = factor(
            field_type,
            ordered = TRUE,
            levels = c("corn", "restored", "remnant"))) %>%
    select(-lat, -long, -yr_restore, -yr_rank) %>% 
    arrange(field_key)
```

## Plant data

### Traits

Plant abundance data was only available at 16 sites (none at Fermi).
These were translated into percent cover of traits in `plant.R`. Site
metadata are joined to allow custom filtering of sites later.

``` r
ptr <- read_csv(paste0(getwd(), "/clean_data/plant_trait_abund.csv"), show_col_types = FALSE) %>% 
    left_join(sites %>% select(starts_with("field"), region), by = join_by("field_name")) %>% 
    left_join(read_csv(paste0(getwd(), "/clean_data/spe_plant_abund.csv"), show_col_types = FALSE) %>% 
                  rename(field_name = SITE) %>% select(field_name, BARESOIL, LITTER), by = join_by(field_name)) %>% 
    select(field_name, field_type, region, everything(), -field_key)
```

Plant releve data was available from four Fermi sites, but survey data
aren’t correctable between Fermi and other sites. Also, translating
counts of species to counts of traits isn’t appropriate. Trait data
isn’t included for the plant presence dataset.

### Plant communities

Both abundance and presence data are used. Plant site-species matrices
use field names for rows, site metadata will be joined.

``` r
pspe <- list(
    ab = read_csv(paste0(getwd(), "/clean_data/spe_plant_abund.csv"), show_col_types = FALSE),
    pr = read_csv(paste0(getwd(), "/clean_data/spe_plant_presence.csv"), show_col_types = FALSE)
) %>% map(function(x) x %>% 
              rename(field_name = SITE) %>% 
              left_join(sites %>% select(starts_with("field"), region), by = join_by("field_name")) %>% 
              select(field_name, field_type, region, everything(), -field_key, -BARESOIL, -LITTER))
```

## Environmental data

Use precipitation as proxy for soil moisture

``` r
rain = read_csv(paste0(getwd(), "/clean_data/site_precip_normal.csv"), show_col_types = FALSE)
```

Create subsets of soil environmental data to align with plant abundance
or presence sites

``` r
soil <- read_csv(paste0(getwd(), "/clean_data/soil.csv"), show_col_types = FALSE) %>% 
    left_join(rain, by = join_by("field_key")) %>% 
    filter(field_key %in% sites$field_key) %>% 
    left_join(sites %>% select(starts_with("field"), region), by = join_by("field_name", "field_key")) %>% 
    select(field_name, field_type, region, everything(), -field_key)
```

## Microbial data

### Fungal communities

Fungal species matrices must have field names instead of field keys to
align all datasets. Create subsets of fungal species matrices to align
with plant abundance or presence sites Sites-species tables contain
rarefied sequence abundances. This list includes composition summarized
in fields.

- Fermi switchgrass prairies are removed because they have no plant data
  associated.
- CSV files were produced in [process_data.R](process_data.md)

``` r
fspe <- list(
    its = read_csv(paste0(getwd(), "/clean_data/spe_ITS_rfy.csv"), show_col_types = FALSE),
    amf = read_csv(paste0(getwd(), "/clean_data/spe_18S_rfy.csv"), show_col_types = FALSE)
) %>% map(function(x) x %>% 
              left_join(sites %>% select(starts_with("field"), region), by = join_by("field_key")) %>% 
              filter(!(field_name %in% c("FLRSP1", "FLRSP2", "FLRSP3"))) %>% 
              select(field_name, field_type, region, everything(), -field_key))
```

### Species metadata

The OTUs and sequence abundances in these files matches the rarefied
data in `spe$` above. CSV files were produced in the [microbial
diversity script](microbial_diversity.md).

``` r
fspe_meta <- list(
    its = read_csv(paste0(getwd(), "/clean_data/speTaxa_ITS_rfy.csv"), show_col_types = FALSE),
    amf =  read_csv(paste0(getwd(), "/clean_data/speTaxa_18S_rfy.csv"), show_col_types = FALSE)
)
```

## Response data

### Fungal biomass

``` r
fb <- read_csv(paste0(getwd(), "/clean_data/plfa.csv"), show_col_types = FALSE) %>% 
    left_join(sites %>% select(starts_with("field"), region), by = join_by("field_key", "field_name")) %>% 
    select(field_name, field_type, region, everything(), -field_key, -starts_with("fa_"))
```

### Water stable aggregates

``` r
# Remove rows from old field sites (26 and 27)
wsa <- read_csv(paste0(getwd(), "/clean_data/wsa.csv"), show_col_types = FALSE)[-c(26:27), ] %>% 
    left_join(sites, by = "field_key")
```

# Analysis and Results

A great number of symmetric and asymmetric comparative and constrained
analyses have been attempted with these data. The best and simplest is a
partial db-RDA, executed with the custom function `dbrda_fun()`. The
strategy is to use most of the soil abiotic data and precipitation as
covariables to remove site-dependent affects on the microbial community.
These covariables are numerous, so they are simplified by being
normalized and then transformed into two PCA axes. Then, variables of
experimental interest are forward selected. These are the variables that
we expect to change as a result of the restoration, or affect the
microbial community directly. Agricultural nutrients, soil organic
matter, and plant community/traits data are in this group.

As before, the plant data aren’t consistent across all regions, so the
analysis is adjusted accordingly as described below. Restricting this
analysis to Blue Mounds only is probably the most defensible approach,
but with the covariables, other regions can be included too. It’s also
possible that a precise approach of the permutation scheme (with `how()`
from package `permute`) would properly handle the design of blocks and
replicates used here.

## Plant community axes

In this analysis, both plant community and traits data can be used. To
analyze all regions, community data based on presence/absence is the
only possibility. Let’s use `pcoa_fun()` to produce community axes for
the abundance and presence/absence plant data.

``` r
(pspe_pcoa_ab <- pcoa_fun(pspe$ab, rg = c("BM", "FG", "LP")))
```

    ## $correction_note
    ## [1] "There were no negative eigenvalues. No correction was applied"
    ## 
    ## $values
    ##   Dim Eigenvalues Relative_eig Broken_stick Cumul_eig Cumul_br_stick
    ## 1   1   0.6592570    0.2365141    0.3143298 0.2365141      0.3143298
    ## 2   2   0.5440994    0.1952003    0.2032187 0.4317144      0.5175485
    ## 3   3   0.4680283    0.1679092    0.1476631 0.5996236      0.6652116
    ## 
    ## $site_vectors
    ##            Axis.1      Axis.2
    ## BBRP1  0.20290126 -0.09287398
    ## ERRP1 -0.43044076 -0.07760580
    ## FGRP1  0.28666433 -0.07970389
    ## KORP1  0.26355718 -0.14999596
    ## LPRP1 -0.02010419  0.47968897
    ## LPRP2  0.12484815  0.38958210
    ## MBRP1  0.17920064 -0.21338125
    ## MHRP1  0.01611089 -0.20409549
    ## MHRP2 -0.45084157 -0.14870043
    ## PHRP1 -0.17189593  0.09708574

``` r
(pspe_pcoa_pr <- pcoa_fun(pspe$pr, rg = c("BM", "FG", "LP", "FL"), method = "jaccard", binary = TRUE))
```

    ## $correction_note
    ## [1] "There were no negative eigenvalues. No correction was applied"
    ## 
    ## $values
    ##   Dim Eigenvalues Relative_eig Broken_stick Cumul_eig Cumul_br_stick
    ## 1   1   0.9562197    0.2299358    0.2586009 0.2299358      0.2586009
    ## 2   2   0.5920732    0.1423719    0.1752676 0.3723077      0.4338684
    ## 3   3   0.4634133    0.1114339    0.1336009 0.4837416      0.5674693
    ## 
    ## $site_vectors
    ##            Axis.1       Axis.2
    ## BBRP1  0.20641293  0.094262023
    ## ERRP1  0.21439889 -0.274249184
    ## FGRP1  0.08005329  0.046027339
    ## FLRP1 -0.52071381 -0.016041532
    ## FLRP4 -0.51926721  0.004443971
    ## FLRP5 -0.41527736  0.051578203
    ## KORP1  0.16890218  0.327780631
    ## LPRP1  0.19555015  0.233535348
    ## LPRP2  0.11823091  0.253683185
    ## MBRP1  0.05335394 -0.139092020
    ## MHRP1  0.18966781 -0.305187933
    ## MHRP2  0.08054080 -0.389175199
    ## PHRP1  0.14814748  0.112435167

## Microbial communities with explanatory and covariables

Partial db-RDA

### General fungal community (ITS sequence abundance)

#### Blue Mounds with plant traits data

``` r
dbrda_fun(
    s = fspe$its,
    pspe_pcoa = "none",
    ft = c("restored"),
    rg = c("BM")
)[c(3, 4, 2)]
```

    ## $R2
    ## $R2$r.squared
    ## numeric(0)
    ## 
    ## $R2$adj.r.squared
    ## numeric(0)
    ## 
    ## 
    ## $global_axis_test
    ## No constrained component
    ## 
    ## Model: dbrda(formula = fspe_bray ~ 1 + Condition(covars), data = expl, sqrt.dist = TRUE)
    ##          Df SumOfSqs  F Pr(>F)
    ## Model     0   0.0000  0       
    ## Residual  6   1.4099          
    ## 
    ## $select_mod
    ## Call: dbrda(formula = fspe_bray ~ 1 + Condition(covars), data = expl,
    ## sqrt.dist = TRUE)
    ## 
    ##               Inertia Proportion Rank
    ## Total          2.1520     1.0000     
    ## Conditional    0.7420     0.3448    2
    ## Unconstrained  1.4099     0.6552    4
    ## Inertia is Bray distance 
    ## 
    ## Eigenvalues for unconstrained axes:
    ##   MDS1   MDS2   MDS3   MDS4 
    ## 0.4836 0.3676 0.3127 0.2460

No explanatory variables were selected, and the set of permutations was
less than the 1999 selected, suggesting that this is a pretty small
dataset.

#### Wisconsin sites with plant traits

``` r
(dbrda_wi_tr_its <-
        dbrda_fun(
            s = fspe$its,
            pspe_pcoa = "none",
            ft = c("restored"),
            rg = c("BM", "LP", "FG")
        ))[c(3, 4, 2)]
```

    ## $R2
    ## $R2$r.squared
    ## [1] 0.17
    ## 
    ## $R2$adj.r.squared
    ## [1] 0.09
    ## 
    ## 
    ## $global_axis_test
    ## Permutation test for dbrda under reduced model
    ## Permutation: free
    ## Number of permutations: 1999
    ## 
    ## Model: dbrda(formula = fspe_bray ~ Condition(covars) + yr_since, data = expl, sqrt.dist = TRUE)
    ##          Df SumOfSqs      F Pr(>F)  
    ## Model     1  0.54694 1.7384  0.017 *
    ## Residual  6  1.88772                
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## $select_mod
    ## Call: dbrda(formula = fspe_bray ~ Condition(covars) + yr_since, data =
    ## expl, sqrt.dist = TRUE)
    ## 
    ##               Inertia Proportion Rank
    ## Total          3.1951     1.0000     
    ## Conditional    0.7604     0.2380    2
    ## Constrained    0.5469     0.1712    1
    ## Unconstrained  1.8877     0.5908    6
    ## Inertia is Bray distance 
    ## 
    ## Eigenvalues for constrained axes:
    ## dbRDA1 
    ## 0.5469 
    ## 
    ## Eigenvalues for unconstrained axes:
    ##   MDS1   MDS2   MDS3   MDS4   MDS5   MDS6 
    ## 0.3876 0.3722 0.3477 0.3010 0.2563 0.2230

Global and individual axis test for axis 1 were significant. Years since
restoration was selected, and it explains 17% of the variation. Forb and
C4 grass are runners-up but appear highly correlated with years (not
shown). Let’s view a plot and include forb and C4 grass for
visualization purposes:

``` r
plot_dbrda(site_sc = dbrda_wi_tr_its$plot_data$sites,
           site_bp = dbrda_wi_tr_its$plot_data$biplot)
```

<img src="tgr_constrained_files/figure-gfm/plot_wi_tr_its-1.png" style="display: block; margin: auto;" />

``` r
write_csv(as_tibble(dbrda_wi_tr_its$plot_data$sites, rownames = "field_name"), "tgr_constrained_files/wi_tr_its_sitelocs.csv")
write_csv(as_tibble(dbrda_wi_tr_its$plot_data$biplot, rownames = "envvar"), "tgr_constrained_files/wi_tr_its_bp.csv")
```

#### Wisconsin sites with plant community axes

``` r
(
    dbrda_wi_ab_its <-
        dbrda_fun(
            s = fspe$its,
            pspe_pcoa = pspe_pcoa_ab$site_vectors,
            ft = c("restored"),
            rg = c("BM", "LP", "FG")
        )
)[c(3, 4, 2)]
```

    ## $R2
    ## $R2$r.squared
    ## [1] 0.17
    ## 
    ## $R2$adj.r.squared
    ## [1] 0.09
    ## 
    ## 
    ## $global_axis_test
    ## Permutation test for dbrda under reduced model
    ## Permutation: free
    ## Number of permutations: 1999
    ## 
    ## Model: dbrda(formula = fspe_bray ~ Condition(covars) + yr_since, data = expl, sqrt.dist = TRUE)
    ##          Df SumOfSqs      F Pr(>F)  
    ## Model     1  0.54694 1.7384 0.0115 *
    ## Residual  6  1.88772                
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## $select_mod
    ## Call: dbrda(formula = fspe_bray ~ Condition(covars) + yr_since, data =
    ## expl, sqrt.dist = TRUE)
    ## 
    ##               Inertia Proportion Rank
    ## Total          3.1951     1.0000     
    ## Conditional    0.7604     0.2380    2
    ## Constrained    0.5469     0.1712    1
    ## Unconstrained  1.8877     0.5908    6
    ## Inertia is Bray distance 
    ## 
    ## Eigenvalues for constrained axes:
    ## dbRDA1 
    ## 0.5469 
    ## 
    ## Eigenvalues for unconstrained axes:
    ##   MDS1   MDS2   MDS3   MDS4   MDS5   MDS6 
    ## 0.3876 0.3722 0.3477 0.3010 0.2563 0.2230

Global and individual axis test for axis 1 were significant. Years since
restoration was selected, and it explains 17% of the variation. Since
the plant community axes failed to contribute explanatory power, the
result of the test is identical to the previous one.

#### All regions with plant community axes

``` r
(
    dbrda_all_pr_its <-
        dbrda_fun(
            s = fspe$its,
            pspe_pcoa = pspe_pcoa_pr$site_vectors,
            ft = c("restored"),
            rg = c("BM", "LP", "FG", "FL")
        )
)[c(3, 4, 2)]
```

    ## $R2
    ## $R2$r.squared
    ## [1] 0.16
    ## 
    ## $R2$adj.r.squared
    ## [1] 0.11
    ## 
    ## 
    ## $global_axis_test
    ## Permutation test for dbrda under reduced model
    ## Permutation: free
    ## Number of permutations: 1999
    ## 
    ## Model: dbrda(formula = fspe_bray ~ Condition(covars) + yr_since, data = expl, sqrt.dist = TRUE)
    ##          Df SumOfSqs      F Pr(>F)    
    ## Model     1   0.7023 2.2343  5e-04 ***
    ## Residual  9   2.8289                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## $select_mod
    ## Call: dbrda(formula = fspe_bray ~ Condition(covars) + yr_since, data =
    ## expl, sqrt.dist = TRUE)
    ## 
    ##               Inertia Proportion Rank
    ## Total          4.3407     1.0000     
    ## Conditional    0.8095     0.1865    2
    ## Constrained    0.7023     0.1618    1
    ## Unconstrained  2.8289     0.6517    9
    ## Inertia is Bray distance 
    ## 
    ## Eigenvalues for constrained axes:
    ## dbRDA1 
    ## 0.7023 
    ## 
    ## Eigenvalues for unconstrained axes:
    ##   MDS1   MDS2   MDS3   MDS4   MDS5   MDS6   MDS7   MDS8   MDS9 
    ## 0.4449 0.3870 0.3617 0.3448 0.3054 0.2633 0.2533 0.2474 0.2210

Global and individual axis test for axis 1 were significant and strong.
Years since restoration was selected and it explains 16% of the
variation here. Potassium and plant axis 1 were runners up. It’s nice to
see that the addition of a very different plant community didn’t make
much of a difference, and that the conditional variation was 18%. I
expected more with how different Fermi soils are. Let’s view a plot and
add SOM, plant axis 1, and potassium for visualization purposes.

``` r
plot_dbrda(site_sc = dbrda_all_pr_its$plot_data$sites,
           site_bp = dbrda_all_pr_its$plot_data$biplot)
```

<img src="tgr_constrained_files/figure-gfm/plot_all_pr_its-1.png" style="display: block; margin: auto;" />

``` r
write_csv(as_tibble(dbrda_all_pr_its$plot_data$sites, rownames = "field_name"), "tgr_constrained_files/all_its_sitelocs.csv")
write_csv(as_tibble(dbrda_all_pr_its$plot_data$biplot, rownames = "envvar"), "tgr_constrained_files/all_its_bp.csv")
```

### AMF community (18S sequence abundance)

#### Blue Mounds with plant traits data

``` r
dbrda_fun(
    s = fspe$amf,
    pspe_pcoa = "none",
    ft = c("restored"),
    rg = c("BM")
)[c(3, 4, 2)]
```

    ## $R2
    ## $R2$r.squared
    ## numeric(0)
    ## 
    ## $R2$adj.r.squared
    ## numeric(0)
    ## 
    ## 
    ## $global_axis_test
    ## No constrained component
    ## 
    ## Model: dbrda(formula = fspe_bray ~ 1 + Condition(covars), data = expl, sqrt.dist = TRUE)
    ##          Df SumOfSqs  F Pr(>F)
    ## Model     0   0.0000  0       
    ## Residual  6   1.1161          
    ## 
    ## $select_mod
    ## Call: dbrda(formula = fspe_bray ~ 1 + Condition(covars), data = expl,
    ## sqrt.dist = TRUE)
    ## 
    ##               Inertia Proportion Rank
    ## Total          1.6960     1.0000     
    ## Conditional    0.5799     0.3419    2
    ## Unconstrained  1.1161     0.6581    4
    ## Inertia is Bray distance 
    ## 
    ## Eigenvalues for unconstrained axes:
    ##   MDS1   MDS2   MDS3   MDS4 
    ## 0.5240 0.2568 0.1960 0.1393

No explanatory variables were selected. Blue Mounds has too few sites
for the number of conditional and explanatory variables used, perhaps.

#### Wisconsin sites with plant traits data

``` r
(dbrda_wi_tr_amf <-
        dbrda_fun(
            s = fspe$amf,
            pspe_pcoa = "none",
            ft = c("restored"),
            rg = c("BM", "LP", "FG")
        ))[c(3, 4, 2)]
```

    ## $R2
    ## $R2$r.squared
    ## [1] 0.23
    ## 
    ## $R2$adj.r.squared
    ## [1] 0.19
    ## 
    ## 
    ## $global_axis_test
    ## Permutation test for dbrda under reduced model
    ## Permutation: free
    ## Number of permutations: 1999
    ## 
    ## Model: dbrda(formula = fspe_bray ~ Condition(covars) + yr_since, data = expl, sqrt.dist = TRUE)
    ##          Df SumOfSqs      F Pr(>F)   
    ## Model     1  0.57437 2.7177 0.0055 **
    ## Residual  6  1.26808                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## $select_mod
    ## Call: dbrda(formula = fspe_bray ~ Condition(covars) + yr_since, data =
    ## expl, sqrt.dist = TRUE)
    ## 
    ##               Inertia Proportion Rank
    ## Total          2.4978     1.0000     
    ## Conditional    0.6554     0.2624    2
    ## Constrained    0.5744     0.2299    1
    ## Unconstrained  1.2681     0.5077    6
    ## Inertia is Bray distance 
    ## 
    ## Eigenvalues for constrained axes:
    ## dbRDA1 
    ## 0.5744 
    ## 
    ## Eigenvalues for unconstrained axes:
    ##    MDS1    MDS2    MDS3    MDS4    MDS5    MDS6 
    ## 0.30525 0.27886 0.21898 0.18494 0.15579 0.12426

Global and individual axis test for axis 1 were significant in site rank
at p\<0.01. Years since restoration was the selected explanatory
variable, explaining 23% of the variation in communities. Forb and C4
grass were runners up and appear highly correlated with years since
restoration. Let’s view a plot and include forb and C4 grass for
visualization purposes:

``` r
plot_dbrda(site_sc = dbrda_wi_tr_amf$plot_data$sites,
           site_bp = dbrda_wi_tr_amf$plot_data$biplot)
```

<img src="tgr_constrained_files/figure-gfm/plot_wi_tr_amf-1.png" style="display: block; margin: auto;" />

#### Wisconsin sites with plant community axes

``` r
(
    dbrda_wi_ab_amf <-
        dbrda_fun(
            s = fspe$amf,
            pspe_pcoa = pspe_pcoa_ab$site_vectors,
            ft = c("restored"),
            rg = c("BM", "LP", "FG")
        )
)[c(3, 4, 2)]
```

    ## $R2
    ## $R2$r.squared
    ## [1] 0.23
    ## 
    ## $R2$adj.r.squared
    ## [1] 0.19
    ## 
    ## 
    ## $global_axis_test
    ## Permutation test for dbrda under reduced model
    ## Permutation: free
    ## Number of permutations: 1999
    ## 
    ## Model: dbrda(formula = fspe_bray ~ Condition(covars) + yr_since, data = expl, sqrt.dist = TRUE)
    ##          Df SumOfSqs      F Pr(>F)   
    ## Model     1  0.57437 2.7177  0.004 **
    ## Residual  6  1.26808                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## $select_mod
    ## Call: dbrda(formula = fspe_bray ~ Condition(covars) + yr_since, data =
    ## expl, sqrt.dist = TRUE)
    ## 
    ##               Inertia Proportion Rank
    ## Total          2.4978     1.0000     
    ## Conditional    0.6554     0.2624    2
    ## Constrained    0.5744     0.2299    1
    ## Unconstrained  1.2681     0.5077    6
    ## Inertia is Bray distance 
    ## 
    ## Eigenvalues for constrained axes:
    ## dbRDA1 
    ## 0.5744 
    ## 
    ## Eigenvalues for unconstrained axes:
    ##    MDS1    MDS2    MDS3    MDS4    MDS5    MDS6 
    ## 0.30525 0.27886 0.21898 0.18494 0.15579 0.12426

Global and individual axis test for axis 1 were significant. Years since
restoration was selected with the same strength as the previous test.

#### All regions with plant community axes

``` r
(
    dbrda_all_pr_amf <-
        dbrda_fun(
            s = fspe$amf,
            pspe_pcoa = pspe_pcoa_pr$site_vectors,
            ft = c("restored"),
            rg = c("BM", "LP", "FG", "FL")
        )
)[c(3, 4, 2)]
```

    ## $R2
    ## $R2$r.squared
    ## [1] 0.21
    ## 
    ## $R2$adj.r.squared
    ## [1] 0.18
    ## 
    ## 
    ## $global_axis_test
    ## Permutation test for dbrda under reduced model
    ## Permutation: free
    ## Number of permutations: 1999
    ## 
    ## Model: dbrda(formula = fspe_bray ~ Condition(covars) + yr_since, data = expl, sqrt.dist = TRUE)
    ##          Df SumOfSqs      F Pr(>F)    
    ## Model     1    0.695 3.2851  0.001 ***
    ## Residual  9    1.904                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## $select_mod
    ## Call: dbrda(formula = fspe_bray ~ Condition(covars) + yr_since, data =
    ## expl, sqrt.dist = TRUE)
    ## 
    ##               Inertia Proportion Rank
    ## Total          3.2951     1.0000     
    ## Conditional    0.6961     0.2112    2
    ## Constrained    0.6950     0.2109    1
    ## Unconstrained  1.9040     0.5778    9
    ## Inertia is Bray distance 
    ## 
    ## Eigenvalues for constrained axes:
    ## dbRDA1 
    ##  0.695 
    ## 
    ## Eigenvalues for unconstrained axes:
    ##   MDS1   MDS2   MDS3   MDS4   MDS5   MDS6   MDS7   MDS8   MDS9 
    ## 0.3596 0.2898 0.2732 0.2252 0.1779 0.1656 0.1561 0.1355 0.1211

Global and individual axis test for axis 1 were significant in site
rank. Years since restoration explained 21% of the variation and
potassium was a runner up. Let’s view a plot and add SOM, plant axis 1,
and potassium for visualization purposes.

``` r
plot_dbrda(site_sc = dbrda_all_pr_amf$plot_data$sites,
           site_bp = dbrda_all_pr_amf$plot_data$biplot)
```

<img src="tgr_constrained_files/figure-gfm/plot_all_pr_amf-1.png" style="display: block; margin: auto;" />

Microbial communities align with years since restoration across regions
and types (general fungi and amf).

## Effects of soil and plants

We know that restored plant communities differ among fields, and that
those differences change in a systematic way over time. Within our
chronosequence sites at Blue Mounds, what is the relative effect of soil
and plant data on soil microbes?

Variation partitioning will be used. Since forward selection failed to
find a parsimonious number of important explanatory variables, we’ll
just use the entire plant and soil datasets here to look for the
relative contribution of these weak effects.

With so many explanatory axes, the analysis failed with raw data, so
explanatory data will first be transformed into PCoA axes.

``` r
soil_pcoa <- pcoa_fun(soil, rg = c("BM"), corr = "lingoes")
vpdat_its <- list(
    Y = fspe$its %>% filter(region == "BM", field_type == "restored") %>% 
        select(-field_type, -region) %>% 
        data.frame(., row.names = 1),
    X1 = pspe_pcoa_ab$site_vectors[-c(3,5,6), ],
    X2 = soil_pcoa$site_vectors
)
vpdat_its_zcols <- vpdat_its %>% map(\(df) which(apply(df, 2, sum) == 0))
```

``` r
(vp_its <- varpart(vpdat_its$Y %>% select(-vpdat_its_zcols$Y), vpdat_its$X1, vpdat_its$X2))
```

    ## 
    ## Partition of variance in RDA 
    ## 
    ## Call: varpart(Y = vpdat_its$Y %>% select(-vpdat_its_zcols$Y), X =
    ## vpdat_its$X1, vpdat_its$X2)
    ## 
    ## Explanatory tables:
    ## X1:  vpdat_its$X1
    ## X2:  vpdat_its$X2 
    ## 
    ## No. of explanatory tables: 2 
    ## Total variation (SS): 264796663 
    ##             Variance: 44132777 
    ## No. of observations: 7 
    ## 
    ## Partition table:
    ##                      Df R.squared Adj.R.squared Testable
    ## [a+c] = X1            2   0.38627       0.07941     TRUE
    ## [b+c] = X2            2   0.34566       0.01849     TRUE
    ## [a+b+c] = X1+X2       4   0.70270       0.10810     TRUE
    ## Individual fractions                                    
    ## [a] = X1|X2           2                 0.08962     TRUE
    ## [b] = X2|X1           2                 0.02869     TRUE
    ## [c]                   0                -0.01021    FALSE
    ## [d] = Residuals                         0.89190    FALSE
    ## ---
    ## Use function 'rda' to test significance of fractions of interest

``` r
plot(vp_its, digits = 2, bg = c("tan", "palegreen"))
```

<img src="tgr_constrained_files/figure-gfm/varpart_its_plot-1.png" style="display: block; margin: auto;" />

With soil axes accounted for, plant axes explain 9% of the ITS fungal
community variation in Blue Mounds. The unique explanation made by soil
variables is about one-third as much. Neither is a powerful explanation.

## Response correlations

Fungal communities varied with years since restoration, C4 grass, and
forbs. How do these predictors affect microbial biomass and function?

``` r
func_vars <-
    sites %>%
    left_join(ptr %>% select(field_name, C4_grass, forb), by = join_by(field_name)) %>%
    left_join(fb %>% select(field_name, fungi, amf), by = join_by(field_name)) %>%
    left_join(wsa %>% select(field_name, wsa), by = join_by(field_name)) %>%
    rename(
        C4_grass_pct = C4_grass,
        forb_pct = forb,
        mass_fungi = fungi,
        mass_amf = amf
    ) %>%
    select(-field_key)
```

Plant traits data are only available in Wisconsin; try these first.

``` r
ggpairs(
    data = func_vars %>% filter(field_type == "restored", region != "FL"),
    columns = 4:9
) +
    theme_bw()
```

<img src="tgr_constrained_files/figure-gfm/pairs_wi-1.png" style="display: block; margin: auto;" />

In all restored sites, response correlations are still possible without
plant traits.

``` r
ggpairs(
    data = func_vars %>% filter(field_type == "restored") %>% 
        select(-C4_grass_pct, -forb_pct),
    columns = 4:7
) +
    theme_bw()
```

<img src="tgr_constrained_files/figure-gfm/pairs_all-1.png" style="display: block; margin: auto;" />

Intercorrelations among responses are obvious and not very useful. Years
since restoration correlates positively with fungal biomass, and with
sites from Fermi included, this relationship is fairly strong. Fermi
sites, with abundant SOM, are probably influencing this result more than
years are, though.

\*\* What if we correlate the constrained axes with responses?\*\* This
might help a little by reordering sites based on time since restoration
*and* community differences, reducing the leverage of sites from Fermi.
Let’s arrange the first axes from each dbRDA and join them with response
data.

``` r
axis_corr <- 
    list(
        wi_tr_its = dbrda_wi_tr_its$plot_data$sites,
        wi_ab_its = dbrda_wi_ab_its$plot_data$sites,
        all_pr_its = dbrda_all_pr_its$plot_data$sites,
        wi_tr_amf = dbrda_wi_tr_amf$plot_data$sites,
        wi_ab_amf = dbrda_wi_ab_amf$plot_data$sites,
        all_pr_amf = dbrda_all_pr_amf$plot_data$sites
    ) %>% 
    map( ~ .x %>%
             data.frame() %>% 
             rownames_to_column(var = "field_name") %>% 
             left_join(sites %>% select(field_name, yr_since), by = join_by(field_name)) %>% 
             left_join(fb %>% select(field_name, fungi, amf), by = join_by(field_name)) %>% 
             left_join(wsa %>% select(field_name, wsa), by = join_by(field_name)) %>% 
             rename(mass_fungi = fungi, mass_amf = amf) %>% 
             select(-dbRDA2)
    )
```

``` r
axis_corr_plot <-
    axis_corr %>%
    map( ~ .x %>%
            ggpairs(columns = 2:6) +
            theme_bw())
```

Pairs panels were investigated (not shown), and only fungal biomass
stood out, as before. Results of simple linear models might help
interpret the importance of fungal biomass. Unfortunately, the only
tests of this that showed significance were those that included Fermi,
which returns us to the confounded design and inability statistically to
handle this problem due to lack of reps.

This isn’t a satisfying analysis. The most defensible chronosequence,
the sites from Blue Mounds, failed to significantly relate to any
constraints, and no functional responses showed any meaningful
relationship here.
