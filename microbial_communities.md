Microbial data: community differences
================
Beau Larkin

Last updated: 12 March, 2023

- <a href="#description" id="toc-description">Description</a>
- <a href="#packages-and-libraries"
  id="toc-packages-and-libraries">Packages and libraries</a>
- <a href="#data" id="toc-data">Data</a>
  - <a href="#sites-species-tables"
    id="toc-sites-species-tables">Sites-species tables</a>
  - <a href="#species-metadata" id="toc-species-metadata">Species
    metadata</a>
  - <a href="#site-metadata" id="toc-site-metadata">Site metadata</a>
  - <a href="#distance-tables" id="toc-distance-tables">Distance tables</a>
  - <a href="#functions" id="toc-functions">Functions</a>
- <a href="#results" id="toc-results">Results</a>
  - <a href="#ordinations" id="toc-ordinations">Ordinations</a>
    - <a href="#pcoa-with-its-gene-otu-clusters"
      id="toc-pcoa-with-its-gene-otu-clusters">PCoA with ITS gene, OTU
      clusters</a>
    - <a href="#pcoa-with-18s-gene-uninformed-distance"
      id="toc-pcoa-with-18s-gene-uninformed-distance">PCoA with 18S gene,
      uninformed distance</a>
    - <a href="#pcoa-with-18s-gene-unifrac-distance"
      id="toc-pcoa-with-18s-gene-unifrac-distance">PCoA with 18S gene, UNIFRAC
      distance</a>

# Description

Microbial data include site-species tables derived from high-throughput
sequencing and clustering in QIIME by Lorinda Bullington and PLFA/NLFA
data which Ylva Lekberg did.

This presents basic visualizations of community differences among
sites/regions based on ITS data.

One goal here is to see whether choosing OTU or SV clusters presents
qualitatively different outcomes in ordinations. We will choose one
(OTUs) if they start to look similar, as they have so far.

Species distance matrices are resampled to the minimum number which
successfully amplified per field. This was done to equalize sampling
effort. This procedure can easily be undone in the [process_data
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

Samples-species tables with rarefied sequence abundances. CSV files were
produced in [process_data.R](process_data.md)

``` r
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

Needed for figure interpretation and permanova designs.

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

## Distance tables

Creating distance objects from the samples-species tables is done with
the typical process of `vegdist()` in vegan. Bray-Curtis or Ruzicka
distance are both appropriate methods for these data, but Bray-Curtis
has produced axes with better explanatory power (Ruzicka is used with
method=“jaccard”) With the 18S data, we can take advantage of
phylogenetic relationships in a UNIFRAC distance matrix. The UNIFRAC
distance was produced in QIIME II and needs some wrangling to conform to
the standards of a distance object in R. The following list contains
vegdist-produced distance objects for ITS and 18S, and it includes
UNIFRAC distance for 18S.

``` r
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
```

## Functions

A function handles the Principal Components Analysis (PCoA) diagnostics,
with outputs and figures saved to a list for later use.

``` r
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
```

# Results

## Ordinations

Bray-Curtis or Ruzicka distance are both appropriate, but Bray-Curtis
has produced axes with better explanatory power (Ruzicka is used with
method=“jaccard”)

In trial runs, no negative eigenvalues were observed (not shown). No
correction is needed for these ordinations.

### PCoA with ITS gene, OTU clusters

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
    ## 1   1   1.2529397   0.17860673   0.15733159 0.1786067      0.1573316
    ## 2   2   0.7693629   0.10967279   0.11566492 0.2882795      0.2729965
    ## 3   3   0.5885190   0.08389347   0.09483159 0.3721730      0.3678281
    ## 
    ## $eigenvalues
    ## [1] 17.9 11.0
    ## 
    ## $site_vectors
    ##    field_key      Axis.1        Axis.2 region field_type yr_since
    ## 1          1  0.22889326  0.0630662983     BM   restored       16
    ## 2          2 -0.11117308 -0.0281387494     BM   restored        3
    ## 3          3 -0.31826809 -0.0006234872     FG       corn        -
    ## 4          4  0.10862190  0.2780476298     FG    remnant        +
    ## 5          5 -0.04816594  0.2759278897     FG   restored       15
    ## 6          6 -0.28641178 -0.1210065337     FL       corn        -
    ## 7          7 -0.32608067  0.0063205971     FL       corn        -
    ## 8          8  0.07996342  0.1789179971     FL    remnant        +
    ## 9          9  0.22819377  0.2400639142     FL   restored       40
    ## 10        10  0.31740215  0.1223057714     FL   restored       36
    ## 11        11  0.21345934  0.2149288690     FL   restored       35
    ## 12        12  0.21750925 -0.2098882992     FL   restored       10
    ## 13        13  0.18055205 -0.1205774656     FL   restored       10
    ## 14        14  0.21775717 -0.2687214844     FL   restored       10
    ## 15        15  0.25833525  0.0387953241     BM   restored       28
    ## 16        16 -0.36081618 -0.1194989086     LP       corn        -
    ## 17        17  0.06563181 -0.1531959944     LP    remnant        +
    ## 18        18 -0.24227804 -0.0604847820     LP   restored        4
    ## 19        19 -0.13782033 -0.0896144556     LP   restored        4
    ## 20        20  0.23956962 -0.3320133587     BM    remnant        +
    ## 21        21  0.18553214 -0.2968120636     BM   restored       18
    ## 22        22 -0.11722721  0.1873339170     BM   restored        7
    ## 23        23 -0.22147390  0.0611701134     BM   restored        2
    ## 24        24 -0.34348248 -0.0275405697     BM       corn        -
    ## 25        25 -0.02822344  0.1612378312     BM   restored       11
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
    ## field_type  2   1.2007 0.17116 2.2715  0.001 ***
    ## Residual   22   5.8144 0.82884                  
    ## Total      24   7.0151 1.00000                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Axis 1 explains 17.9% of the variation and is the only eigenvalue that
exceeds a broken stick model. The most substantial variation here will
be on the first axis, although axis 2 explains 11% of the variation and
was very close to the broken stick value. Testing the design factor
*field_type* (with *region* treated as a block using the `strata`
argument of `adonis2`) revealed a significant clustering
$(R^2=0.17, p=0.001)$.

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
```

``` r
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
    ## -0.18318 -0.09727 -0.04444  0.10319  0.20917 
    ## 
    ## Coefficients:
    ##              Estimate Std. Error t value Pr(>|t|)   
    ## (Intercept) -0.104215   0.054894  -1.898  0.07844 . 
    ## yr_since     0.011280   0.002783   4.053  0.00119 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.1349 on 14 degrees of freedom
    ## Multiple R-squared:  0.5399, Adjusted R-squared:  0.507 
    ## F-statistic: 16.43 on 1 and 14 DF,  p-value: 0.001186

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

It’s probably better to do with with a new ordination of just
restoration sites and a constrained ordination with years and other
environmental variables.

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
    ## [1] 4
    ## 
    ## $correction_note
    ## [1] "No correction was applied to the negative eigenvalues"
    ## 
    ## $values
    ##   Dim Eigenvalues Relative_eig Rel_corr_eig Broken_stick Cum_corr_eig
    ## 1   1   1.2121589   0.27150740   0.24544804   0.16236050    0.2454480
    ## 2   2   0.8184553   0.18332305   0.16726204   0.11888224    0.4127101
    ## 3   3   0.5844705   0.13091359   0.12079476   0.09714311    0.5335048
    ## 4   4   0.4140203   0.09273501   0.08694488   0.08265036    0.6204497
    ## 5   5   0.3103469   0.06951355   0.06635627   0.07178079    0.6868060
    ##   Cumul_br_stick
    ## 1      0.1623605
    ## 2      0.2812427
    ## 3      0.3783858
    ## 4      0.4610362
    ## 5      0.5328170
    ## 
    ## $eigenvalues
    ## [1] 27.2 18.3
    ## 
    ## $site_vectors
    ##    field_key      Axis.1       Axis.2      Axis.3       Axis.4 region
    ## 1          1  0.19701412  0.252225472 -0.20165034 -0.086966501     BM
    ## 2          2  0.01751762 -0.270096774 -0.20822461 -0.241162872     BM
    ## 3          3 -0.42929415  0.219043681  0.12875132 -0.003013836     FG
    ## 4          4  0.05147254  0.039304784  0.25555480 -0.156876296     FG
    ## 5          5 -0.02489329 -0.041343487  0.18694096 -0.361772727     FG
    ## 6          6 -0.19102744 -0.009212571 -0.11841789  0.061284937     FL
    ## 7          7 -0.43919512  0.283347149 -0.19778043 -0.199267082     FL
    ## 8          8  0.10038590 -0.026270323  0.13794482 -0.029892972     FL
    ## 9          9  0.19155558  0.232678704  0.20434480  0.071483878     FL
    ## 10        10  0.24752162  0.127639437  0.07663480  0.186392886     FL
    ## 11        11  0.12616725  0.117527402  0.14284140  0.133722486     FL
    ## 12        12  0.14763130 -0.096342656 -0.14557491  0.022364688     FL
    ## 13        13  0.17082397  0.036008697 -0.15066866  0.091973080     FL
    ## 14        14  0.08070522 -0.060309769 -0.22613829  0.060949633     FL
    ## 15        15  0.31208752  0.301787002  0.08886715 -0.019619367     BM
    ## 16        16 -0.40535594  0.054165861 -0.09013882  0.138177087     LP
    ## 17        17  0.01726908 -0.217243396 -0.05418968  0.072707040     LP
    ## 18        18 -0.08510416 -0.254286674  0.05276762  0.150446905     LP
    ## 19        19  0.04761389 -0.304640046 -0.05819587  0.032269206     LP
    ## 20        20  0.24699196  0.173914281 -0.19353391 -0.030958072     BM
    ## 21        21  0.18486204 -0.058548123 -0.11440842 -0.024023002     BM
    ## 22        22 -0.01953250 -0.214479197  0.11427676  0.006353751     BM
    ## 23        23 -0.13023127 -0.132524065  0.15652450 -0.018805790     BM
    ## 24        24 -0.43607391  0.071504188  0.02810668  0.165022799     BM
    ## 25        25  0.02108816 -0.223849578  0.18536621 -0.020789858     BM
    ##    field_type yr_since
    ## 1    restored       16
    ## 2    restored        3
    ## 3        corn        -
    ## 4     remnant        +
    ## 5    restored       15
    ## 6        corn        -
    ## 7        corn        -
    ## 8     remnant        +
    ## 9    restored       40
    ## 10   restored       36
    ## 11   restored       35
    ## 12   restored       10
    ## 13   restored       10
    ## 14   restored       10
    ## 15   restored       28
    ## 16       corn        -
    ## 17    remnant        +
    ## 18   restored        4
    ## 19   restored        4
    ## 20    remnant        +
    ## 21   restored       18
    ## 22   restored        7
    ## 23   restored        2
    ## 24       corn        -
    ## 25   restored       11
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
    ## field_type  2   1.0873 0.24354 3.5414 0.0015 **
    ## Residual   22   3.3773 0.75646                 
    ## Total      24   4.4646 1.00000                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Four axes are significant by a broken stick model, between them
explaining 67.8% of the variation in AMF among fields. It may be
worthwhile to examine structure on Axes 3 and 4 sometime. The most
substantial variation here is on the first axis (27.2%) with Axis 2
explaining 18.3% of the variation in AMF abundances. Testing the design
factor *field_type* (with *region* treated as a block using the `strata`
argument of `adonis2`) revealed a significant clustering
$(R^2=0.24, p=0.002)$.

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
    ##      Min       1Q   Median       3Q      Max 
    ## -0.12806 -0.07851  0.01215  0.07969  0.13220 
    ## 
    ## Coefficients:
    ##              Estimate Std. Error t value Pr(>|t|)   
    ## (Intercept) -0.016172   0.037191  -0.435  0.67031   
    ## yr_since     0.007002   0.001886   3.714  0.00231 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.0914 on 14 degrees of freedom
    ## Multiple R-squared:  0.4963, Adjusted R-squared:  0.4603 
    ## F-statistic: 13.79 on 1 and 14 DF,  p-value: 0.002313

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
    ## -0.12978 -0.07912 -0.04487  0.05337  0.28352 
    ## 
    ## Coefficients:
    ##              Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept) -0.232185   0.049866  -4.656 0.000371 ***
    ## yr_since     0.012556   0.002528   4.967 0.000207 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.1226 on 14 degrees of freedom
    ## Multiple R-squared:  0.6379, Adjusted R-squared:  0.6121 
    ## F-statistic: 24.67 on 1 and 14 DF,  p-value: 0.000207

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

Try fitting years since restoration with envfit…
