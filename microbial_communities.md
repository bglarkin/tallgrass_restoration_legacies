Microbial data: community differences
================
Beau Larkin

Last updated: 11 September, 2024

- [Description](#description)
- [Packages and libraries](#packages-and-libraries)
  - [Functions](#functions)
- [Data](#data)
  - [Site metadata](#site-metadata)
  - [Sites-species tables](#sites-species-tables)
  - [Species metadata](#species-metadata)
  - [Distance tables](#distance-tables)
- [Results](#results)
  - [ITS gene, OTU clustering](#its-gene-otu-clustering)
    - [PCoA with Blue Mounds restored fields, all
      subsamples](#pcoa-with-blue-mounds-restored-fields-all-subsamples)
    - [PCoA with all fields and regions, all
      subsamples](#pcoa-with-all-fields-and-regions-all-subsamples)
    - [PCoA in Faville Grove, all
      subsamples](#pcoa-in-faville-grove-all-subsamples)
    - [PCoA in Fermilab, all
      subsamples](#pcoa-in-fermilab-all-subsamples)
    - [PCoA in Lake Petite Prairie, all
      subsamples](#pcoa-in-lake-petite-prairie-all-subsamples)
    - [PCoA ordination, all regions, all
      subsamples](#pcoa-ordination-all-regions-all-subsamples)
  - [18S gene, OTU clustering](#18s-gene-otu-clustering)
    - [PCoA with abundances summed in fields, Bray-Curtis
      distance](#pcoa-with-abundances-summed-in-fields-bray-curtis-distance)
    - [PCoA with abundances summed in fields, UNIFRAC
      distance](#pcoa-with-abundances-summed-in-fields-unifrac-distance)
    - [PCoA with Blue Mounds restored fields, all
      subsamples](#pcoa-with-blue-mounds-restored-fields-all-subsamples-1)
    - [PCoA with all fields and regions, all
      subsamples](#pcoa-with-all-fields-and-regions-all-subsamples-1)
    - [PCoA in Blue Mounds, all
      subsamples](#pcoa-in-blue-mounds-all-subsamples-1)
    - [PCoA in Faville Grove, all
      subsamples](#pcoa-in-faville-grove-all-subsamples-1)
    - [PCoA in Fermilab, all
      subsamples](#pcoa-in-fermilab-all-subsamples-1)
    - [PCoA in Lake Petite Prairie, all
      subsamples](#pcoa-in-lake-petite-prairie-all-subsamples-1)
    - [PCoA ordination, all regions, all
      subsamples](#pcoa-ordination-all-regions-all-subsamples-1)

# Description

Microbial data include site-species tables derived from high-throughput
sequencing and clustering in QIIME by Lorinda Bullington and PLFA/NLFA
data which Ylva Lekberg did.

This presents basic visualizations of community differences among
sites/regions based on ITS and 18S data.

During data processing, not all subsamples were retained. Some had
failed to amplify and others had very few sequences, leading to the
potential for a loss of information during rarefication. With the loss
of some subsamples, all fields were resampled to the same lower number
of samples. This was done to equalize sampling effort (from a
statistical perspective). This procedure can easily be undone in the
[process_data script](process_data.md). Whether 9, 8, or 7 subsamples
are retained, the interpretation of analyses presented here would be the
same (not shown).

Pairwise contrasts in multivariate analysis were accomplished with a
custom function adapted from [O’Leary et
al. 2021](https://link.springer.com/article/10.1007/s12237-021-00917-2).

# Packages and libraries

``` r
packages_needed = c("tidyverse", "vegan", "colorspace", "ape", "knitr", "gridExtra")
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

**Functions are stored** in a separate
[script](supporting_files/microbial_communities_functions.md) to reduce
clutter here and allow for easier editing.

``` r
source("supporting_files/microbial_communities_functions.R")
```

# Data

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
        yr_since = replace(yr_since, which(field_type == "remnant"), NA),
        yr_since = replace(yr_since, which(field_type == "corn"), NA)) %>%
    select(-lat, -long, -yr_restore, -yr_rank) %>% 
    arrange(field_key)
sites_resto_bm <- 
    sites %>% 
    filter(field_type == "restored",
           region == "BM") %>% 
    select(-field_name, -region) %>% 
    mutate(yr_since = as.numeric(yr_since))   
```

## Sites-species tables

Sites-species tables with rarefied sequence abundances. This list
includes composition summarized by fields or unsummarized (all samples).
It also includes subsets by region. All subsets have zero sum columns
removed.  
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
    its_samps_bm = read_csv(
        paste0(getwd(), "/clean_data/spe_ITS_rfy_samples.csv"),
        show_col_types = FALSE
    ) %>% 
        left_join(sites %>% select(field_key, region), by = join_by(field_key)) %>% 
        filter(region == "BM") %>% 
        select(-region),
    its_samps_fg = read_csv(
        paste0(getwd(), "/clean_data/spe_ITS_rfy_samples.csv"),
        show_col_types = FALSE
    ) %>% 
        left_join(sites %>% select(field_key, region), by = join_by(field_key)) %>% 
        filter(region == "FG") %>% 
        select(-region),
    its_samps_fl = read_csv(
        paste0(getwd(), "/clean_data/spe_ITS_rfy_samples.csv"),
        show_col_types = FALSE
    ) %>% 
        left_join(sites %>% select(field_key, region), by = join_by(field_key)) %>% 
        filter(region == "FL") %>% 
        select(-region),
    its_samps_lp = read_csv(
        paste0(getwd(), "/clean_data/spe_ITS_rfy_samples.csv"),
        show_col_types = FALSE
    ) %>% 
        left_join(sites %>% select(field_key, region), by = join_by(field_key)) %>% 
        filter(region == "LP") %>% 
        select(-region),
    amf = read_csv(
        paste0(getwd(), "/clean_data/spe_18S_rfy.csv"),
        show_col_types = FALSE
    ),
    amf_samps = read_csv(
        paste0(getwd(), "/clean_data/spe_18S_rfy_samples.csv"),
        show_col_types = FALSE
    ),
    amf_samps_bm = read_csv(
        paste0(getwd(), "/clean_data/spe_18S_rfy_samples.csv"),
        show_col_types = FALSE
    ) %>% 
        left_join(sites %>% select(field_key, region), by = join_by(field_key)) %>% 
        filter(region == "BM") %>% 
        select(-region),
    amf_samps_fg = read_csv(
        paste0(getwd(), "/clean_data/spe_18S_rfy_samples.csv"),
        show_col_types = FALSE
    ) %>% 
        left_join(sites %>% select(field_key, region), by = join_by(field_key)) %>% 
        filter(region == "FG") %>% 
        select(-region),
    amf_samps_fl = read_csv(
        paste0(getwd(), "/clean_data/spe_18S_rfy_samples.csv"),
        show_col_types = FALSE
    ) %>% 
        left_join(sites %>% select(field_key, region), by = join_by(field_key)) %>% 
        filter(region == "FL") %>% 
        select(-region),
    amf_samps_lp = read_csv(
        paste0(getwd(), "/clean_data/spe_18S_rfy_samples.csv"),
        show_col_types = FALSE
    ) %>% 
        left_join(sites %>% select(field_key, region), by = join_by(field_key)) %>% 
        filter(region == "LP") %>% 
        select(-region)
) %>% 
    map(. %>% select(where( ~ sum(.) != 0)))
```

## Species metadata

Needed to make inset figures showing most important categories of
species. The OTUs and sequence abundances in these files matches the
rarefied data in `spe$` above. CSV files were produced in the guild
taxonomy [script](microbial_guild_taxonomy.md).

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
- *gene_samps_region*: objects are distances matrices taken from
  rarefied data, subsetted to region, with zero sum columns removed.
  Samples in each field depend on the gene-based dataset, see above.

``` r
index <- "bray"
distab <- list(
    its       = vegdist(data.frame(spe$its, row.names = 1), method = index),
    its_samps = vegdist(
        data.frame(
            spe$its_samps %>% 
                mutate(field_sample = paste(field_key, sample, sep = "_")) %>% 
                column_to_rownames(var = "field_sample") %>% 
                select(-field_key, -sample)
        ) %>% select(where(~ sum(.) > 0)), method = index),
    its_resto_bm = vegdist(
        data.frame(
            spe$its %>% 
                filter(field_key %in% sites_resto_bm$field_key), 
            row.names = 1
        ) %>% select(where(~ sum(.) > 0)), method = index),
    its_resto_samps_bm = vegdist(
        data.frame(
            spe$its_samps %>% 
                filter(field_key %in% sites_resto_bm$field_key) %>% 
                mutate(field_sample = paste(field_key, sample, sep = "_")) %>% 
                column_to_rownames(var = "field_sample") %>% 
                select(-field_key, -sample)
        ) %>% select(where(~ sum(.) > 0)), method = index),
    its_samps_bm = vegdist(
        data.frame(
            spe$its_samps_bm %>% 
                mutate(field_sample = paste(field_key, sample, sep = "_")) %>% 
                column_to_rownames(var = "field_sample") %>% 
                select(-field_key, -sample)
        ), method = index), # zero sum columns were already removed in the spe list
    its_samps_fg = vegdist(
        data.frame(
            spe$its_samps_fg %>% 
                mutate(field_sample = paste(field_key, sample, sep = "_")) %>% 
                column_to_rownames(var = "field_sample") %>% 
                select(-field_key, -sample)
        ), method = index), # zero sum columns were already removed in the spe list
    its_samps_fl = vegdist(
        data.frame(
            spe$its_samps_fl %>% 
                mutate(field_sample = paste(field_key, sample, sep = "_")) %>% 
                column_to_rownames(var = "field_sample") %>% 
                select(-field_key, -sample)
        ), method = index), # zero sum columns were already removed in the spe list
    its_samps_lp = vegdist(
        data.frame(
            spe$its_samps_lp %>% 
                mutate(field_sample = paste(field_key, sample, sep = "_")) %>% 
                column_to_rownames(var = "field_sample") %>% 
                select(-field_key, -sample)
        ), method = index), # zero sum columns were already removed in the spe list
    amf_bray  = vegdist(data.frame(spe$amf, row.names = 1), method = index),
    amf_samps = vegdist(
        data.frame(
            spe$amf_samps %>% 
                mutate(field_sample = paste(field_key, sample, sep = "_")) %>% 
                column_to_rownames(var = "field_sample") %>% 
                select(-field_key, -sample)
        ) %>% select(where(~ sum(.) > 0)), method = index),
    amf_resto_bm = vegdist(
        data.frame(
            spe$amf %>% 
                filter(field_key %in% sites_resto_bm$field_key), 
            row.names = 1
        ) %>% select(where(~ sum(.) > 0)), method = index),
    amf_resto_samps_bm = vegdist(
        data.frame(
            spe$amf_samps %>% 
                filter(field_key %in% sites_resto_bm$field_key) %>% 
                mutate(field_sample = paste(field_key, sample, sep = "_")) %>% 
                column_to_rownames(var = "field_sample") %>% 
                select(-field_key, -sample)
        ) %>% select(where(~ sum(.) > 0)), method = index),
    amf_samps_bm = vegdist(
        data.frame(
            spe$amf_samps_bm %>% 
                mutate(field_sample = paste(field_key, sample, sep = "_")) %>% 
                column_to_rownames(var = "field_sample") %>% 
                select(-field_key, -sample)
        ), method = index), # zero sum columns were already removed in the spe list
    amf_samps_fg = vegdist(
        data.frame(
            spe$amf_samps_fg %>% 
                mutate(field_sample = paste(field_key, sample, sep = "_")) %>% 
                column_to_rownames(var = "field_sample") %>% 
                select(-field_key, -sample)
        ), method = index), # zero sum columns were already removed in the spe list
    amf_samps_fl = vegdist(
        data.frame(
            spe$amf_samps_fl %>% 
                mutate(field_sample = paste(field_key, sample, sep = "_")) %>% 
                column_to_rownames(var = "field_sample") %>% 
                select(-field_key, -sample)
        ), method = index), # zero sum columns were already removed in the spe list
    amf_samps_lp = vegdist(
        data.frame(
            spe$amf_samps_lp %>% 
                mutate(field_sample = paste(field_key, sample, sep = "_")) %>% 
                column_to_rownames(var = "field_sample") %>% 
                select(-field_key, -sample)
        ), method = index), # zero sum columns were already removed in the spe list
    amf_uni   = sites %>%
        select(field_name, field_key) %>%
        left_join(
            read_delim(paste0(getwd(), "/otu_tables/18S/18S_weighted_Unifrac.tsv"), show_col_types = FALSE),
            by = join_by(field_name)) %>%
        select(field_key, everything(),-field_name) %>%
        data.frame(row.names = 1) %>%
        as.dist()
) 
```

# Results

#### Ordinations

Bray-Curtis, Morisita-Horn, or Ruzicka distance are appropriate, but
Bray-Curtis has produced axes with better explanatory power.

## ITS gene, OTU clustering

In trial runs, no negative eigenvalues were observed (not shown). No
\### PCoA with abundances summed in fields correction is needed for
these ordinations.

``` r
(pcoa_its <- pcoa_fun(spe$its, distab$its, adonis_index = "bray", df_name = "ITS gene, 97% OTU"))
```

    ## 'nperm' >= set of all permutations: complete enumeration.

    ## Set of permutations < 'minperm'. Generating entire set.

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
    ## 1   1   1.2802276   0.18661535   0.15733159 0.1866154      0.1573316
    ## 2   2   0.7864201   0.11463436   0.11566492 0.3012497      0.2729965
    ## 3   3   0.5929217   0.08642861   0.09483159 0.3876783      0.3678281
    ## 
    ## $eigenvalues
    ## [1] 18.7 11.5
    ## 
    ## $site_vectors
    ##    field_key      Axis.1       Axis.2 region field_type yr_since
    ## 1          1  0.22567685 -0.053150223     BM   restored       16
    ## 2          2 -0.10332809  0.007845206     BM   restored        3
    ## 3          3 -0.33074971  0.025086045     FG       corn       NA
    ## 4          4  0.10090403 -0.301897155     FG    remnant       NA
    ## 5          5 -0.05949033 -0.285299292     FG   restored       15
    ## 6          6 -0.29465345  0.130140919     FL       corn       NA
    ## 7          7 -0.32606680  0.027555097     FL       corn       NA
    ## 8          8  0.09372057 -0.222121172     FL    remnant       NA
    ## 9          9  0.24857486 -0.192591329     FL   restored       40
    ## 10        10  0.32222241 -0.081876654     FL   restored       36
    ## 11        11  0.22000364 -0.168616366     FL   restored       35
    ## 12        12  0.21559326  0.238376330     FL   restored       10
    ## 13        13  0.18477609  0.115943948     FL   restored       10
    ## 14        14  0.19591933  0.254293029     FL   restored       10
    ## 15        15  0.24591083 -0.054154050     BM   restored       28
    ## 16        16 -0.37640561  0.120798839     LP       corn       NA
    ## 17        17  0.07471032  0.182194215     LP    remnant       NA
    ## 18        18 -0.24345537  0.065408981     LP   restored        4
    ## 19        19 -0.14663486  0.091447342     LP   restored        4
    ## 20        20  0.24492434  0.308576799     BM    remnant       NA
    ## 21        21  0.19239945  0.278486280     BM   restored       18
    ## 22        22 -0.11048797 -0.215563334     BM   restored        7
    ## 23        23 -0.22011657 -0.089520293     BM   restored        2
    ## 24        24 -0.33782242  0.029967438     BM       corn       NA
    ## 25        25 -0.01612480 -0.211330601     BM   restored       11
    ## 
    ## $broken_stick_plot

<img src="microbial_communities_files/figure-gfm/pcoa_its_otu-1.png" style="display: block; margin: auto;" />

    ## 
    ## $permanova
    ## Permutation test for adonis under reduced model
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1999
    ## 
    ## adonis2(formula = d ~ field_type, data = env, permutations = nperm, method = adonis_index, add = if (corr == "none") FALSE else "lingoes", strata = region)
    ##          Df SumOfSqs      R2      F Pr(>F)    
    ## Model     2   1.2184 0.17761 2.3756  5e-04 ***
    ## Residual 22   5.6418 0.82239                  
    ## Total    24   6.8602 1.00000                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## $pairwise_contrasts
    ## 
    ## 
    ## group1     group2        R2   F_value   df1   df2     p_value   p_value_adj
    ## ---------  --------  ------  --------  ----  ----  ----------  ------------
    ## restored   corn       0.156     3.523     1    19   0.0015000        0.0045
    ## restored   remnant    0.055     1.056     1    18   0.1455000        0.1455
    ## corn       remnant    0.289     2.850     1     7   0.0416667        0.0625

Axis 1 explains 18.7% of the variation and is the only eigenvalue that
exceeds a broken stick model. The most substantial variation here will
be on the first axis, although axis 2 explains 11.5% of the variation
and was very close to the broken stick value. Testing the design factor
*field_type* (with *region* treated as a block using the `strata`
argument of `adonis2`) revealed a significant clustering
$(R^2=0.18,~p=5\times 10^{-4})$.

Let’s view a plot with abundances of community subgroups inset.

``` r
pcoa_its$ord <-
    ggplot(pcoa_its$site_vectors, aes(x = Axis.1, y = Axis.2)) +
    geom_point(aes(fill = field_type, shape = region), size = 8) +
    geom_text(aes(label = yr_since), size = 4) +
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
        caption = "Text in icons for restored fields indicates years since restoration."
    ) +
    lims(y = c(-0.35,0.44)) +
    theme_bw() +
    guides(fill = guide_legend(override.aes = list(shape = 21)))
pcoa_its$inset <-
    spe_meta$its %>%
    filter(primary_lifestyle %in% c("plant_pathogen", "soil_saprotroph")) %>%
    mutate(field_type = factor(field_type, ordered = TRUE, levels = c("corn", "restored", "remnant"))) %>%
    group_by(primary_lifestyle, field_type, field_name) %>%
    summarize(sum_seq_abund = sum(seq_abund), .groups = "drop_last") %>% 
    summarize(avg_seq_abund = mean(sum_seq_abund), .groups = "drop") %>%
    ggplot(aes(x = primary_lifestyle, y = avg_seq_abund, fill = field_type)) +
    geom_col(position = "dodge") +
    labs(y = "Seq. abund. (avg)") +
    scale_fill_discrete_qualitative(palette = "Harmonic") +
    scale_x_discrete(label = c("plnt path", "soil sapr")) +
    # coord_flip() +
    theme_classic() +
    theme(legend.position = "none", axis.title.x = element_blank())
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

    ## Warning: Removed 9 rows containing missing values or values outside the scale range
    ## (`geom_text()`).

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

The most appropriate way to look at communities vs. field age is with
the Blue Mounds restored fields. The function `pcoa_its_samps_bm()` will
take care of this. Field age will be fitted to the ordination and tested
using `envfit()`.

### PCoA with Blue Mounds restored fields, all subsamples

In trial runs, no negative eigenvalues were observed (not shown). No
correction is needed for these ordinations.

``` r
(pcoa_its_resto_samps_bm <- pcoa_samps_bm_fun(spe$its_samps, 
                                        distab$its_resto_samps_bm, 
                                        sites_resto_bm, 
                                        adonis_index = "bray", 
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
    ## 1   1    2.143437   0.11526409   0.08352022 0.1152641     0.08352022
    ## 2   2    1.424838   0.07662119   0.06533840 0.1918853     0.14885863
    ## 3   3    1.036238   0.05572406   0.05624749 0.2476093     0.20510612
    ## 
    ## $eigenvalues
    ## [1] 11.5  7.7
    ## 
    ## $site_vectors
    ## # A tibble: 56 × 5
    ##    field_key sample_key  Axis.1  Axis.2 yr_since
    ##        <dbl> <chr>        <dbl>   <dbl>    <dbl>
    ##  1         1 1          -0.188  -0.200        16
    ##  2         1 2          -0.212  -0.0549       16
    ##  3         1 4          -0.182  -0.0656       16
    ##  4         1 5          -0.281  -0.146        16
    ##  5         1 6          -0.328  -0.162        16
    ##  6         1 7          -0.108  -0.0827       16
    ##  7         1 9          -0.209  -0.150        16
    ##  8         1 10         -0.241  -0.0307       16
    ##  9         2 1           0.248   0.216         3
    ## 10         2 2           0.0876  0.0830        3
    ## # ℹ 46 more rows
    ## 
    ## $broken_stick_plot

<img src="microbial_communities_files/figure-gfm/pcoa_its_resto_samps_bm-1.png" style="display: block; margin: auto;" />

    ## 
    ## $permanova
    ## Permutation test for adonis under reduced model
    ## Permutation: free
    ## Number of permutations: 1999
    ## 
    ## adonis2(formula = d ~ field_key, data = env_w, permutations = nperm, method = adonis_index)
    ##          Df SumOfSqs      R2      F Pr(>F)    
    ## Model     1   0.7819 0.04205 2.3702  5e-04 ***
    ## Residual 54  17.8140 0.95795                  
    ## Total    55  18.5959 1.00000                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## $vector_fit
    ## 
    ## ***VECTORS
    ## 
    ##             Axis.1    Axis.2     r2 Pr(>r)  
    ## yr_since -0.996470  0.083919 0.7272  0.022 *
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Plots: field_key, plot permutation: free
    ## Permutation: none
    ## Number of permutations: 5039
    ## 
    ## 
    ## 
    ## $vector_fit_scores
    ##              Axis.1     Axis.2
    ## yr_since -0.8497566 0.07156277

Axis 1 explains 11.5% and axis 2 explains 7.7% of the variation in the
community data. Both axes are important based on the broken stick model.
Indeed, the first four axes are borderline important. The relatively low
percent variation explained is partly due to the high number of
dimensions used when all samples from fields are included. The fidelity
of samples to fields was significant based on a permutation test
$(R^2=0.04,~p=5\times 10^{-4})$. In this case, the partial $R^2$ shows
the proportion of sum of squares from the total. It is a low number here
because so much unexplained variation exists, resulting in a high sum of
squares that is outside the assignment of subsamples to fields.

Years since restoration has a moderately strong correlation with
communities and was significant with a permutation test where samples
were constrained within fields to account for lack of independence \#’
$(R^2=0.73,~p=0.02)$.

Let’s view an ordination plot with hulls around subsamples and a fitted
vector for field age overlaid.

``` r
centroid_its_bm <- aggregate(cbind(Axis.1, Axis.2) ~ field_key, data = pcoa_its_resto_samps_bm$site_vectors, mean) %>% 
    left_join(sites %>% select(field_key, yr_since), by = join_by(field_key))
hull_its_bm <- pcoa_its_resto_samps_bm$site_vectors %>% 
    group_by(field_key) %>% 
    slice(chull(Axis.1, Axis.2))
```

``` r
ggplot(pcoa_its_resto_samps_bm$site_vectors, aes(x = Axis.1, y = Axis.2)) +
    geom_point(fill = "#5CBD92", shape = 21) +
    geom_polygon(data = hull_its_bm, aes(group = as.character(field_key)), fill = "#5CBD92", alpha = 0.3) +
    geom_point(data = centroid_its_bm, fill = "#5CBD92", size = 8, shape = 21) +
    geom_text(data = centroid_its_bm, aes(label = yr_since)) +
    geom_segment(aes(x = 0, 
                     y = 0, 
                     xend = pcoa_its_resto_samps_bm$vector_fit_scores[1] * 0.4, 
                     yend = pcoa_its_resto_samps_bm$vector_fit_scores[2] * 0.4),
                 color = "blue", 
                 arrow = arrow(length = unit(3, "mm"))) +
    labs(
        x = paste0("Axis 1 (", pcoa_its_resto_samps_bm$eigenvalues[1], "%)"),
        y = paste0("Axis 2 (", pcoa_its_resto_samps_bm$eigenvalues[2], "%)"),
        title = paste0(
            "PCoA Ordination (",
            pcoa_its_resto_samps_bm$dataset,
            ")"
        ),
        caption = "Text indicates years since restoration.\nYears since restoration significant at p<0.05."
    ) +
    theme_bw() +
    theme(legend.position = "none")
```

    ## Warning in geom_segment(aes(x = 0, y = 0, xend = pcoa_its_resto_samps_bm$vector_fit_scores[1] * : All aesthetics have length 1, but the data has 56 rows.
    ## ℹ Please consider using `annotate()` or provide this layer with data containing
    ##   a single row.

<img src="microbial_communities_files/figure-gfm/its_samps_bm_fig-1.png" style="display: block; margin: auto;" />

### PCoA with all fields and regions, all subsamples

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
                                  adonis_index = "bray", 
                                  df_name = "ITS gene, 97% OTU"))
```

    ## 'nperm' >= set of all permutations: complete enumeration.

    ## Set of permutations < 'minperm'. Generating entire set.

    ## $dataset
    ## [1] "ITS gene, 97% OTU"
    ## 
    ## $components_exceed_broken_stick
    ## [1] 11
    ## 
    ## $correction_note
    ## [1] "Lingoes correction applied to negative eigenvalues: D' = -0.5*D^2 - 0.0582972612582258 , except diagonal elements"
    ## 
    ## $values
    ##    Dim Eigenvalues Corr_eig Rel_corr_eig Broken_stick Cum_corr_eig Cum_br_stick
    ## 1    1    6.695143 6.753440   0.08163976   0.02963639   0.08163976   0.02963639
    ## 2    2    4.127910 4.186207   0.05060547   0.02458589   0.13224523   0.05422228
    ## 3    3    3.319685 3.377982   0.04083514   0.02206064   0.17308037   0.07628292
    ## 4    4    2.624358 2.682656   0.03242960   0.02037713   0.20550998   0.09666005
    ## 5    5    2.134346 2.192643   0.02650603   0.01911451   0.23201601   0.11577456
    ## 6    6    2.040445 2.098743   0.02537090   0.01810441   0.25738691   0.13387896
    ## 7    7    1.695953 1.754250   0.02120646   0.01726266   0.27859337   0.15114162
    ## 8    8    1.482016 1.540313   0.01862026   0.01654115   0.29721364   0.16768277
    ## 9    9    1.342005 1.400302   0.01692772   0.01590984   0.31414136   0.18359262
    ## 10  10    1.298941 1.357238   0.01640714   0.01534867   0.33054850   0.19894129
    ## 11  11    1.170003 1.228301   0.01484846   0.01484362   0.34539696   0.21378492
    ## 12  12    1.106225 1.164522   0.01407747   0.01438449   0.35947442   0.22816940
    ## 
    ## $eigenvalues
    ## [1] 8.2 5.1
    ## 
    ## $site_vectors
    ## # A tibble: 200 × 17
    ##    field_key sample_key  Axis.1   Axis.2   Axis.3  Axis.4  Axis.5   Axis.6
    ##        <dbl> <chr>        <dbl>    <dbl>    <dbl>   <dbl>   <dbl>    <dbl>
    ##  1         1 1          -0.181  -0.129   -0.0310  -0.0711  0.0695 -0.102  
    ##  2         1 2          -0.222   0.0710  -0.0228  -0.0157  0.0133 -0.217  
    ##  3         1 4          -0.192  -0.0437  -0.00762  0.0210 -0.0157 -0.124  
    ##  4         1 5          -0.182  -0.0190  -0.128   -0.207  -0.0270  0.0333 
    ##  5         1 6          -0.207   0.00685 -0.120   -0.174   0.0201 -0.167  
    ##  6         1 7          -0.0507 -0.0490   0.0194  -0.107  -0.0310  0.0261 
    ##  7         1 9          -0.121  -0.0848  -0.0596  -0.0955  0.0630 -0.0284 
    ##  8         1 10         -0.187   0.0149  -0.0701  -0.0417 -0.0225 -0.0223 
    ##  9         2 1           0.0564 -0.00475  0.331    0.0520 -0.0671  0.00535
    ## 10         2 2           0.0429  0.0510   0.0898   0.0130 -0.123  -0.130  
    ## # ℹ 190 more rows
    ## # ℹ 9 more variables: Axis.7 <dbl>, Axis.8 <dbl>, Axis.9 <dbl>, Axis.10 <dbl>,
    ## #   Axis.11 <dbl>, field_name <chr>, region <chr>, field_type <ord>,
    ## #   yr_since <dbl>
    ## 
    ## $broken_stick_plot

<img src="microbial_communities_files/figure-gfm/pcoa_its_samps-1.png" style="display: block; margin: auto;" />

    ## 
    ## $permanova
    ## Permutation test for adonis under reduced model
    ## Blocks:  region 
    ## Plots: field_key, plot permutation: free
    ## Permutation: none
    ## Number of permutations: 1999
    ## 
    ## adonis2(formula = d ~ field_type, data = env_w, permutations = gl_perm_design, method = adonis_index)
    ##           Df SumOfSqs      R2      F Pr(>F)    
    ## Model      2    5.750 0.08084 8.6636  0.001 ***
    ## Residual 197   65.372 0.91916                  
    ## Total    199   71.121 1.00000                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## $pairwise_contrasts
    ##     group1  group2    R2 F_value df1 df2   p_value p_value_adj
    ## 1 restored    corn 0.066  11.742   1 166 0.0015000   0.0045000
    ## 2 restored remnant 0.017   2.795   1 158 0.2620000   0.2620000
    ## 3     corn remnant 0.134  10.850   1  70 0.1428571   0.2142857
    ## 
    ## $format
    ## [1] "pandoc"

``` r
write_delim(pcoa_its_samps$permanova, "microbial_communities_files/pcoa_its_samps_permanova.txt")
write_delim(pcoa_its_samps$pairwise_contrasts %>% mutate(across(starts_with("p_value"), ~ round(.x, 3))), "microbial_communities_files/pcoa_its_samps_pairwise.txt")
```

Axis 1 explains 8.2% and axis 2 explains 5.1% of the variation in the
community data. Both axes are important based on the broken stick model,
in fact, the broken stick model shows that 11 axes are important in
explaining variation with this dataset. The relatively low percent
variation explained on axes 1 and 2 is partly due to the high number of
dimensions used when all samples from fields are included. The fidelity
of samples to fields was strong based on a permutation test when
restricting permutations to fields (=plots in `how()`) within regions
(=blocks in `how()`) $(R^2=0.08,~p=0.001)$.

Let’s view an ordination plot with hulls around subsamples.

``` r
centroid_its <- aggregate(cbind(Axis.1, Axis.2) ~ field_key, data = pcoa_its_samps$site_vectors, mean) %>% 
    left_join(sites %>% select(field_key, yr_since, field_type, region), by = join_by(field_key))
hull_its <- pcoa_its_samps$site_vectors %>% 
    group_by(field_key) %>% 
    slice(chull(Axis.1, Axis.2))
```

``` r
its_samps_fig <- 
    ggplot(pcoa_its_samps$site_vectors, aes(x = -1 * (Axis.1), y = Axis.2)) +
    geom_vline(xintercept = 0, linewidth = 0.1) +
    geom_hline(yintercept = 0, linewidth = 0.1) +
    geom_point(aes(fill = field_type), shape = 21, alpha = 0.8, color = "gray10") +
    geom_polygon(data = hull_its, aes(group = field_key, fill = field_type), alpha = 0.3) +
    geom_point(data = centroid_its, aes(fill = field_type, shape = region), size = 6) +
    geom_text(data = centroid_its, aes(label = yr_since), size = 3) +
    labs(
        x = paste0("Axis 1 (", pcoa_its_samps$eigenvalues[1], "%)"),
        y = paste0("Axis 2 (", pcoa_its_samps$eigenvalues[2], "%)"),
        title = paste0(
            "PCoA Ordination (",
            pcoa_its_samps$dataset,
            ")"
        ),
        caption = "Text indicates years since restoration."
    ) +
    lims(y = c(-0.35, 0.48)) +
    scale_fill_discrete_qualitative(name = "Field Type", palette = "Harmonic") +
    scale_shape_manual(name = "Region", values = c(21, 22, 23, 24)) +
    theme_bw() +
    guides(fill = guide_legend(override.aes = list(shape = 21)))
```

``` r
(its_samps_guilds_fig <- 
    its_samps_fig +
    annotation_custom(
        ggplotGrob(
            pcoa_its$inset + 
                theme(
            plot.background = element_rect(colour = "black", fill = "gray90"), 
            axis.title.y = element_text(size = 8)
        )),
        xmin = -0.40,
        xmax = -0.05,
        ymin = 0.20,
        ymax = 0.48
    ))
```

    ## Warning: Removed 9 rows containing missing values or values outside the scale range
    ## (`geom_text()`).

<img src="microbial_communities_files/figure-gfm/its_samps_guilds_fig-1.png" style="display: block; margin: auto;" />

#### PCoA in Blue Mounds, all subsamples

This is as above with the diagnostics and permutation tests. Pairwise
contrasts among field types should be ignored here because there is no
replication.

``` r
(pcoa_its_samps_bm <- pcoa_samps_fun(
    s = spe$its_samps_bm,
    d = distab$its_samps_bm,
    env = sites %>% filter(region == "BM"),
    corr = "none",
    df_name = "Blue Mounds, ITS gene, 97% OTU"
))
```

    ## $dataset
    ## [1] "Blue Mounds, ITS gene, 97% OTU"
    ## 
    ## $components_exceed_broken_stick
    ## [1] 4
    ## 
    ## $correction_note
    ## [1] "There were no negative eigenvalues. No correction was applied"
    ## 
    ## $values
    ##   Dim Eigenvalues Relative_eig Broken_stick Cumul_eig Cumul_br_stick
    ## 1   1   2.9575981   0.11566017   0.06826650 0.1156602      0.0682665
    ## 2   2   1.8106716   0.07080833   0.05418199 0.1864685      0.1224485
    ## 3   3   1.6165316   0.06321627   0.04713974 0.2496848      0.1695882
    ## 4   4   1.1148355   0.04359688   0.04244490 0.2932817      0.2120331
    ## 5   5   0.9949093   0.03890704   0.03892377 0.3321887      0.2509569
    ## 
    ## $eigenvalues
    ## [1] 11.6  7.1
    ## 
    ## $site_vectors
    ## # A tibble: 72 × 10
    ##    field_key sample_key  Axis.1   Axis.2  Axis.3  Axis.4 field_name region
    ##        <dbl> <chr>        <dbl>    <dbl>   <dbl>   <dbl> <chr>      <chr> 
    ##  1         1 1           0.0796  0.127   -0.294  -0.0408 BBRP1      BM    
    ##  2         1 2           0.157   0.0948  -0.161   0.0207 BBRP1      BM    
    ##  3         1 4           0.0879  0.131   -0.201   0.0850 BBRP1      BM    
    ##  4         1 5           0.212  -0.00571 -0.220  -0.238  BBRP1      BM    
    ##  5         1 6           0.194   0.0325  -0.310  -0.0374 BBRP1      BM    
    ##  6         1 7           0.0343  0.0566  -0.146  -0.151  BBRP1      BM    
    ##  7         1 9           0.101   0.0611  -0.252  -0.0863 BBRP1      BM    
    ##  8         1 10          0.165   0.0729  -0.162  -0.0572 BBRP1      BM    
    ##  9         2 1          -0.191   0.243    0.215   0.0802 ERRP1      BM    
    ## 10         2 2          -0.106   0.0241   0.0643  0.0464 ERRP1      BM    
    ## # ℹ 62 more rows
    ## # ℹ 2 more variables: field_type <ord>, yr_since <dbl>
    ## 
    ## $broken_stick_plot

![](microbial_communities_files/figure-gfm/pcoa_its_samps_bm-1.png)<!-- -->

    ## 
    ## $permanova
    ## Permutation test for adonis under reduced model
    ## Blocks:  region 
    ## Plots: field_key, plot permutation: free
    ## Permutation: none
    ## Number of permutations: 1999
    ## 
    ## adonis2(formula = d ~ field_type, data = env_w, permutations = gl_perm_design, method = adonis_index)
    ##          Df SumOfSqs     R2      F Pr(>F)  
    ## Model     2   3.0456 0.1191 4.6645  0.033 *
    ## Residual 69  22.5259 0.8809                
    ## Total    71  25.5714 1.0000                
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## $pairwise_contrasts
    ##     group1  group2    R2 F_value df1 df2 p_value p_value_adj
    ## 1 restored remnant 0.068   4.509   1  62   0.135      0.2025
    ## 2 restored    corn 0.068   4.550   1  62   0.124      0.2025
    ## 3  remnant    corn 0.299   5.978   1  14   1.000      1.0000
    ## 
    ## $format
    ## [1] "pandoc"

Field type remains significant.

### PCoA in Faville Grove, all subsamples

This is as above with the diagnostics and permutation tests. Pairwise
contrasts among field types should be ignored here because there is no
replication.

``` r
(pcoa_its_samps_fg <- pcoa_samps_fun(
    s = spe$its_samps_fg,
    d = distab$its_samps_fg,
    env = sites %>% filter(region == "FG"),
    corr = "none",
    df_name = "Faville Grove, ITS gene, 97% OTU"
))
```

    ## $dataset
    ## [1] "Faville Grove, ITS gene, 97% OTU"
    ## 
    ## $components_exceed_broken_stick
    ## [1] 1
    ## 
    ## $correction_note
    ## [1] "There were no negative eigenvalues. No correction was applied"
    ## 
    ## $values
    ##   Dim Eigenvalues Relative_eig Broken_stick Cumul_eig Cumul_br_stick
    ## 1   1   1.9224683   0.26517808   0.16236050 0.2651781      0.1623605
    ## 2   2   0.8536510   0.11774943   0.11888224 0.3829275      0.2812427
    ## 3   3   0.4475412   0.06173215   0.09714311 0.4446597      0.3783858
    ## 
    ## $eigenvalues
    ## [1] 26.5 11.8
    ## 
    ## $site_vectors
    ## # A tibble: 24 × 8
    ##    field_key sample_key Axis.1   Axis.2 field_name region field_type yr_since
    ##        <dbl> <chr>       <dbl>    <dbl> <chr>      <chr>  <ord>         <dbl>
    ##  1         3 1          -0.367  0.00593 FGC1       FG     corn             NA
    ##  2         3 2          -0.413  0.0783  FGC1       FG     corn             NA
    ##  3         3 3          -0.430  0.0551  FGC1       FG     corn             NA
    ##  4         3 5          -0.330  0.0109  FGC1       FG     corn             NA
    ##  5         3 6          -0.410  0.0383  FGC1       FG     corn             NA
    ##  6         3 7          -0.414 -0.00224 FGC1       FG     corn             NA
    ##  7         3 9          -0.414  0.0689  FGC1       FG     corn             NA
    ##  8         3 10         -0.381 -0.00980 FGC1       FG     corn             NA
    ##  9         4 1           0.260  0.260   FGREM1     FG     remnant          NA
    ## 10         4 2           0.290  0.247   FGREM1     FG     remnant          NA
    ## # ℹ 14 more rows
    ## 
    ## $broken_stick_plot

![](microbial_communities_files/figure-gfm/pcoa_its_samps_fg-1.png)<!-- -->

    ## 
    ## $permanova
    ## Permutation test for adonis under reduced model
    ## Blocks:  region 
    ## Plots: field_key, plot permutation: free
    ## Permutation: none
    ## Number of permutations: 5
    ## 
    ## adonis2(formula = d ~ field_type, data = env_w, permutations = gl_perm_design, method = adonis_index)
    ##          Df SumOfSqs      R2      F Pr(>F)
    ## Model     2   2.7090 0.37367 6.2644      1
    ## Residual 21   4.5407 0.62633              
    ## Total    23   7.2497 1.00000              
    ## 
    ## $pairwise_contrasts
    ##    group1   group2    R2 F_value df1 df2 p_value p_value_adj
    ## 1    corn  remnant 0.341   7.247   1  14       1           1
    ## 2    corn restored 0.341   7.232   1  14       1           1
    ## 3 remnant restored 0.226   4.093   1  14       1           1
    ## 
    ## $format
    ## [1] "pandoc"

Field type is not significant here.

### PCoA in Fermilab, all subsamples

This is as above with the diagnostics and permutation tests. Pairwise
contrasts among field types should be ignored here because there is no
replication.

``` r
(pcoa_its_samps_fl <- pcoa_samps_fun(
    s = spe$its_samps_fl,
    d = distab$its_samps_fl,
    env = sites %>% filter(region == "FL"),
    corr = "lingoes",
    df_name = "Fermilab, ITS gene, 97% OTU"
))
```

    ## $dataset
    ## [1] "Fermilab, ITS gene, 97% OTU"
    ## 
    ## $components_exceed_broken_stick
    ## [1] 2
    ## 
    ## $correction_note
    ## [1] "Lingoes correction applied to negative eigenvalues: D' = -0.5*D^2 - 0.0150971138378135 , except diagonal elements"
    ## 
    ## $values
    ##   Dim Eigenvalues Corr_eig Rel_corr_eig Broken_stick Cum_corr_eig Cum_br_stick
    ## 1   1    3.583733 3.598830   0.14566937   0.06904053    0.1456694   0.06904053
    ## 2   2    2.543095 2.558192   0.10354759   0.05475481    0.2492170   0.12379534
    ## 3   3    1.069182 1.084279   0.04388821   0.04761195    0.2931052   0.17140729
    ## 
    ## $eigenvalues
    ## [1] 14.6 10.4
    ## 
    ## $site_vectors
    ## # A tibble: 72 × 8
    ##    field_key sample_key Axis.1  Axis.2 field_name region field_type yr_since
    ##        <dbl> <chr>       <dbl>   <dbl> <chr>      <chr>  <ord>         <dbl>
    ##  1         6 1          -0.424  0.0237 FLC1       FL     corn             NA
    ##  2         6 2          -0.398 -0.0106 FLC1       FL     corn             NA
    ##  3         6 4          -0.427  0.0382 FLC1       FL     corn             NA
    ##  4         6 5          -0.421  0.0356 FLC1       FL     corn             NA
    ##  5         6 6          -0.304 -0.0238 FLC1       FL     corn             NA
    ##  6         6 7          -0.383  0.0397 FLC1       FL     corn             NA
    ##  7         6 9          -0.338  0.0550 FLC1       FL     corn             NA
    ##  8         6 10         -0.344 -0.0191 FLC1       FL     corn             NA
    ##  9         7 1          -0.425  0.105  FLC2       FL     corn             NA
    ## 10         7 3          -0.401  0.108  FLC2       FL     corn             NA
    ## # ℹ 62 more rows
    ## 
    ## $broken_stick_plot

![](microbial_communities_files/figure-gfm/pcoa_its_samps_fl-1.png)<!-- -->

    ## 
    ## $permanova
    ## Permutation test for adonis under reduced model
    ## Blocks:  region 
    ## Plots: field_key, plot permutation: free
    ## Permutation: none
    ## Number of permutations: 1999
    ## 
    ## adonis2(formula = d ~ field_type, data = env_w, permutations = gl_perm_design, method = adonis_index)
    ##          Df SumOfSqs      R2      F Pr(>F)  
    ## Model     2   4.1615 0.17609 7.3733 0.0155 *
    ## Residual 69  19.4720 0.82391                
    ## Total    71  23.6336 1.00000                
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## $pairwise_contrasts
    ##    group1   group2    R2 F_value df1 df2   p_value p_value_adj
    ## 1    corn  remnant 0.244   7.102   1  22 0.3333333       0.500
    ## 2    corn restored 0.159  11.688   1  62 0.0380000       0.114
    ## 3 remnant restored 0.046   2.601   1  54 0.5820000       0.582
    ## 
    ## $format
    ## [1] "pandoc"

Field type is again significant by permutation test.

### PCoA in Lake Petite Prairie, all subsamples

This is as above with the diagnostics and permutation tests. Pairwise
contrasts among field types should be ignored here because there is no
replication.

``` r
(pcoa_its_samps_lp <- pcoa_samps_fun(
    s = spe$its_samps_lp,
    d = distab$its_samps_lp,
    env = sites %>% filter(region == "LP"),
    corr = "none",
    df_name = "Lake Petite Prairie, ITS gene, 97% OTU"
))
```

    ## $dataset
    ## [1] "Lake Petite Prairie, ITS gene, 97% OTU"
    ## 
    ## $components_exceed_broken_stick
    ## [1] 3
    ## 
    ## $correction_note
    ## [1] "There were no negative eigenvalues. No correction was applied"
    ## 
    ## $values
    ##   Dim Eigenvalues Relative_eig Broken_stick Cumul_eig Cumul_br_stick
    ## 1   1   1.3420307   0.15635442   0.12991114 0.1563544      0.1299111
    ## 2   2   0.9397633   0.10948791   0.09765307 0.2658423      0.2275642
    ## 3   3   0.7198034   0.08386131   0.08152404 0.3497036      0.3090882
    ## 4   4   0.4956706   0.05774852   0.07077135 0.4074522      0.3798596
    ## 
    ## $eigenvalues
    ## [1] 15.6 10.9
    ## 
    ## $site_vectors
    ## # A tibble: 32 × 9
    ##    field_key sample_key Axis.1 Axis.2  Axis.3 field_name region field_type
    ##        <dbl> <chr>       <dbl>  <dbl>   <dbl> <chr>      <chr>  <ord>     
    ##  1        16 1          -0.240 0.0836 -0.170  LPC1       LP     corn      
    ##  2        16 3          -0.230 0.117  -0.170  LPC1       LP     corn      
    ##  3        16 5          -0.262 0.0883 -0.174  LPC1       LP     corn      
    ##  4        16 6          -0.269 0.0936 -0.137  LPC1       LP     corn      
    ##  5        16 7          -0.250 0.0691 -0.140  LPC1       LP     corn      
    ##  6        16 8          -0.204 0.0322 -0.228  LPC1       LP     corn      
    ##  7        16 9          -0.255 0.0797 -0.117  LPC1       LP     corn      
    ##  8        16 10         -0.191 0.0510 -0.136  LPC1       LP     corn      
    ##  9        17 1           0.238 0.145   0.101  LPREM1     LP     remnant   
    ## 10        17 2           0.286 0.205   0.0277 LPREM1     LP     remnant   
    ## # ℹ 22 more rows
    ## # ℹ 1 more variable: yr_since <dbl>
    ## 
    ## $broken_stick_plot

![](microbial_communities_files/figure-gfm/pcoa_its_samps_lp-1.png)<!-- -->

    ## 
    ## $permanova
    ## Permutation test for adonis under reduced model
    ## Blocks:  region 
    ## Plots: field_key, plot permutation: free
    ## Permutation: none
    ## Number of permutations: 23
    ## 
    ## adonis2(formula = d ~ field_type, data = env_w, permutations = gl_perm_design, method = adonis_index)
    ##          Df SumOfSqs      R2      F Pr(>F)
    ## Model     2   1.9352 0.22546 4.2207 0.1667
    ## Residual 29   6.6481 0.77454              
    ## Total    31   8.5833 1.00000              
    ## 
    ## $pairwise_contrasts
    ##    group1   group2    R2 F_value df1 df2   p_value p_value_adj
    ## 1    corn  remnant 0.265   5.035   1  14 1.0000000         1.0
    ## 2    corn restored 0.162   4.251   1  22 0.3333333         0.5
    ## 3 remnant restored 0.144   3.714   1  22 0.3333333         0.5
    ## 
    ## $format
    ## [1] "pandoc"

Let’s view an ordination plot with hulls around subsamples for each
indidual region.

### PCoA ordination, all regions, all subsamples

``` r
pcoa_its_site_vectors <- bind_rows(
    list(
        `Blue Mounds`   = pcoa_its_samps_bm$site_vectors,
        `Faville Grove` = pcoa_its_samps_fg$site_vectors,
        `Fermilab`      = pcoa_its_samps_fl$site_vectors,
        `Lake Petite`   = pcoa_its_samps_lp$site_vectors
    ),
    .id = "place"
)
pcoa_its_eigenvalues <- bind_rows(
    list(
        `Blue Mounds`   = pcoa_its_samps_bm$eigenvalues,
        `Faville Grove` = pcoa_its_samps_fg$eigenvalues,
        `Fermilab`      = pcoa_its_samps_fl$eigenvalues,
        `Lake Petite`   = pcoa_its_samps_lp$eigenvalues
    ),
    .id = "place"
) %>% 
    mutate(axis = c(1,2)) %>% 
    pivot_longer(cols = 1:4, names_to = "place", values_to = "eigenvalue") %>% 
    select(place, axis, eigenvalue) %>% 
    arrange(place, axis) %>% 
    pivot_wider(names_from = axis, names_prefix = "axis_", values_from = eigenvalue)
centroid_regions_its <- aggregate(cbind(Axis.1, Axis.2) ~ place + field_key, data = pcoa_its_site_vectors, mean) %>% 
    left_join(sites %>% select(field_key, yr_since, field_type, region), by = join_by(field_key))
hull_regions_its <- pcoa_its_site_vectors %>% 
    group_by(place, field_key) %>% 
    slice(chull(Axis.1, Axis.2))
```

``` r
(its_samps_regions_fig <- 
    ggplot(pcoa_its_site_vectors, aes(x = Axis.1, y = Axis.2)) +
    facet_wrap(vars(place), scales = "free") +
    geom_vline(xintercept = 0, linewidth = 0.1) +
    geom_hline(yintercept = 0, linewidth = 0.1) +
    geom_point(aes(fill = field_type), shape = 21, alpha = 0.8, color = "gray10") +
    geom_polygon(data = hull_regions_its, aes(group = field_key, fill = field_type), alpha = 0.3) +
    geom_point(data = centroid_regions_its, aes(fill = field_type, shape = region), size = 5) +
    geom_text(data = centroid_regions_its, aes(label = yr_since), size = 2.5) +
    labs(
        x = paste0("Axis 1"),
        y = paste0("Axis 2"),
        caption = "ITS gene. Text indicates years since restoration."
    ) +
    scale_fill_discrete_qualitative(name = "Field Type", palette = "Harmonic") +
    scale_shape_manual(name = "Region", values = c(21, 22, 23, 24)) +
    theme_bw() +
    guides(fill = guide_legend(override.aes = list(shape = 21))))
```

<img src="microbial_communities_files/figure-gfm/its_samps_regions_fig-1.png" style="display: block; margin: auto;" />

The eigenvalues are shown below:

``` r
kable(pcoa_its_eigenvalues, format = "pandoc")
```

| place         | axis_1 | axis_2 |
|:--------------|-------:|-------:|
| Blue Mounds   |   11.6 |    7.1 |
| Faville Grove |   26.5 |   11.8 |
| Fermilab      |   14.6 |   10.4 |
| Lake Petite   |   15.6 |   10.9 |

``` r
write_csv(pcoa_its_eigenvalues, file = "microbial_communities_files/pcoa_its_eig.csv")
```

Let’s view and save a plot that shows all the data together and broken
out by regions.

``` r
grid.arrange(
    its_samps_guilds_fig + labs(caption = "") + theme(plot.title = element_blank()), 
    its_samps_regions_fig + labs(caption = "") + theme(legend.position = "none"), 
    ncol = 1,
    heights = c(1.1,0.9)
    )
```

<img src="microbial_communities_files/figure-gfm/its_samps_unified_fig-1.png" style="display: block; margin: auto;" />

Then, we’ll follow up with panels showing trends with the most abundant
guilds.

``` r
spe_meta$its %>%
    filter(primary_lifestyle %in% c("plant_pathogen", "soil_saprotroph")) %>%
    mutate(field_type = factor(field_type, ordered = TRUE,
                               levels = c("corn", "restored", "remnant")),
           pl_labs = case_match(primary_lifestyle, "plant_pathogen" ~ "Plant Pathogens", "soil_saprotroph" ~ "Soil Saprotrophs")) %>%
    group_by(region, primary_lifestyle, pl_labs, field_type, field_name) %>%
    summarize(sum_seq_abund = sum(seq_abund), .groups = "drop_last") %>% 
    summarize(avg_seq_abund = mean(sum_seq_abund), .groups = "drop") %>%
    ggplot(aes(x = region, y = avg_seq_abund, fill = field_type)) +
    facet_wrap(vars(pl_labs), scales = "free_y") +
    geom_col(position = position_dodge(width = 0.9), color = "black", linewidth = 0.2) +
    labs(y = "Sequence abundance (avg)") +
    scale_fill_discrete_qualitative(name = "Field Type", palette = "Harmonic") +
    theme_bw() +
    theme(axis.title.x = element_blank())
```

<img src="microbial_communities_files/figure-gfm/its_guilds_regions_fig-1.png" style="display: block; margin: auto;" />

## 18S gene, OTU clustering

### PCoA with abundances summed in fields, Bray-Curtis distance

No negative eigenvalues produced, no correction applied.

``` r
(pcoa_amf_bray <- pcoa_fun(s = spe$amf, d = distab$amf_bray, adonis_index = "bray", df_name = "18S gene, 97% OTU, Bray-Curtis distance"))
```

    ## 'nperm' >= set of all permutations: complete enumeration.

    ## Set of permutations < 'minperm'. Generating entire set.

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
    ## 1   1   1.2083007   0.27452361   0.25081040   0.16236050    0.2508104
    ## 2   2   0.7848755   0.17832223   0.16440577   0.11888224    0.4152162
    ## 3   3   0.6169072   0.14016014   0.13012996   0.09714311    0.5453461
    ## 4   4   0.4105128   0.09326774   0.08801289   0.08265036    0.6333590
    ## 5   5   0.3160729   0.07181116   0.06874137   0.07178079    0.7021004
    ## 6   6   0.2172378   0.04935602   0.04857297   0.06308514    0.7506734
    ##   Cumul_br_stick
    ## 1      0.1623605
    ## 2      0.2812427
    ## 3      0.3783858
    ## 4      0.4610362
    ## 5      0.5328170
    ## 6      0.5959021
    ## 
    ## $eigenvalues
    ## [1] 27.5 17.8
    ## 
    ## $site_vectors
    ##    field_key       Axis.1      Axis.2      Axis.3       Axis.4       Axis.5
    ## 1          1  0.201278427  0.23072057 -0.19719994  0.097269137  0.037999587
    ## 2          2 -0.004545958 -0.27676224 -0.20820926  0.222988195 -0.126619060
    ## 3          3 -0.408632815  0.22940385  0.14063974 -0.014968305 -0.024042782
    ## 4          4  0.061297074 -0.02901905  0.27182252  0.168233793  0.044593735
    ## 5          5 -0.013670972 -0.07529477  0.16303671  0.355419689  0.024264322
    ## 6          6 -0.197333247  0.02461424 -0.12738083 -0.052601613 -0.379110692
    ## 7          7 -0.417892984  0.25283387 -0.18515725  0.237554856  0.019913563
    ## 8          8  0.124573031 -0.02002805  0.13544871  0.055883949 -0.066914540
    ## 9          9  0.215002820  0.23523651  0.19953453 -0.052826500 -0.070689478
    ## 10        10  0.250784693  0.16257541  0.08704855 -0.169570802 -0.046771518
    ## 11        11  0.144346094  0.10953107  0.15761634 -0.111349517 -0.106150078
    ## 12        12  0.107749457 -0.09685873 -0.15563029 -0.071331704 -0.041140249
    ## 13        13  0.154957066  0.04276735 -0.13583811 -0.094210977 -0.004730908
    ## 14        14  0.067303808 -0.04551409 -0.25994362 -0.060938126  0.010641690
    ## 15        15  0.323952864  0.25819324  0.10189862  0.046253376 -0.027746545
    ## 16        16 -0.411897333  0.09279246 -0.06816752 -0.117768671  0.053929245
    ## 17        17  0.023716745 -0.19198574 -0.07857517 -0.087265530  0.163438882
    ## 18        18 -0.105709183 -0.23176808  0.05583277 -0.157707698  0.027491130
    ## 19        19  0.036509898 -0.30154080 -0.06430800 -0.048186493 -0.103969506
    ## 20        20  0.261020753  0.16750376 -0.20426953  0.027949293  0.193624545
    ## 21        21  0.183771133 -0.06518390 -0.13553273  0.003577215  0.115712859
    ## 22        22 -0.039411148 -0.20148357  0.13675957 -0.036617485  0.105101657
    ## 23        23 -0.129270093 -0.13982768  0.15487600 -0.016050167  0.020909240
    ## 24        24 -0.435592254  0.11638846  0.04234255 -0.138252863  0.144879218
    ## 25        25  0.007692123 -0.24729409  0.17335564  0.014516947  0.035385683
    ##    region field_type yr_since
    ## 1      BM   restored       16
    ## 2      BM   restored        3
    ## 3      FG       corn       NA
    ## 4      FG    remnant       NA
    ## 5      FG   restored       15
    ## 6      FL       corn       NA
    ## 7      FL       corn       NA
    ## 8      FL    remnant       NA
    ## 9      FL   restored       40
    ## 10     FL   restored       36
    ## 11     FL   restored       35
    ## 12     FL   restored       10
    ## 13     FL   restored       10
    ## 14     FL   restored       10
    ## 15     BM   restored       28
    ## 16     LP       corn       NA
    ## 17     LP    remnant       NA
    ## 18     LP   restored        4
    ## 19     LP   restored        4
    ## 20     BM    remnant       NA
    ## 21     BM   restored       18
    ## 22     BM   restored        7
    ## 23     BM   restored        2
    ## 24     BM       corn       NA
    ## 25     BM   restored       11
    ## 
    ## $broken_stick_plot

<img src="microbial_communities_files/figure-gfm/pcoa_amf_otu-1.png" style="display: block; margin: auto;" />

    ## 
    ## $permanova
    ## Permutation test for adonis under reduced model
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1999
    ## 
    ## adonis2(formula = d ~ field_type, data = env, permutations = nperm, method = adonis_index, add = if (corr == "none") FALSE else "lingoes", strata = region)
    ##          Df SumOfSqs    R2      F Pr(>F)   
    ## Model     2   1.0916 0.248 3.6276 0.0015 **
    ## Residual 22   3.3099 0.752                 
    ## Total    24   4.4014 1.000                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## $pairwise_contrasts
    ## 
    ## 
    ## group1     group2        R2   F_value   df1   df2     p_value   p_value_adj
    ## ---------  --------  ------  --------  ----  ----  ----------  ------------
    ## restored   corn       0.253     6.432     1    19   0.0015000        0.0045
    ## restored   remnant    0.023     0.423     1    18   0.9650000        0.9650
    ## corn       remnant    0.383     4.339     1     7   0.0416667        0.0625

Four axes are significant by a broken stick model, between them
explaining 68.6% of the variation in AMF among fields. It may be
worthwhile to examine structure on Axes 3 and 4 sometime. The most
substantial variation here is on the first axis (27.5%) with Axis 2
explaining 17.8% of the variation in AMF abundances. Testing the design
factor *field_type* (with *region* treated as a block using the `strata`
argument of `adonis2`) revealed a significant clustering
$(R^2=0.25,~p=0.002)$.

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
    group_by(family, field_type, field_name) %>% 
    summarize(sum_seq_abund = sum(seq_abund), .groups = "drop_last") %>% 
    summarize(avg_seq_abund = mean(sum_seq_abund), .groups = "drop") %>%
    ggplot(aes(x = family, y = avg_seq_abund)) +
    geom_col(position = "dodge", aes(fill = field_type)) +
    labs(y = "Seq. abund. (avg)") +
    scale_fill_discrete_qualitative(palette = "Harmonic") +
    scale_x_discrete(label = function(x) abbreviate(x, minlength = 6)) +
    # coord_flip() +
    theme_classic() +
    theme(legend.position = "none", axis.title.x = element_blank())
```

``` r
(amf_families_fig <- 
    pcoa_amf_bray$ord +
    annotation_custom(
        ggplotGrob(
            pcoa_amf_bray$inset + 
                theme(
            plot.background = element_rect(colour = "black", fill = "gray90")
        )),
        xmin = -0.63,
        xmax = -0.2,
        ymin = -0.32,
        ymax = -0.10
    ))
```

    ## Warning: Removed 9 rows containing missing values or values outside the scale range
    ## (`geom_text()`).

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

### PCoA with abundances summed in fields, UNIFRAC distance

``` r
(pcoa_amf_uni <- pcoa_fun(s = spe$amf, d = distab$amf_uni, df_name = "18S gene, 97% OTU, UNIFRAC distance", corr = "lingoes"))
```

    ## 'nperm' >= set of all permutations: complete enumeration.

    ## Set of permutations < 'minperm'. Generating entire set.

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
    ## 3          3 -0.072487397  0.067424513  0.039009283     FG       corn       NA
    ## 4          4 -0.041221010 -0.036068309  0.021185618     FG    remnant       NA
    ## 5          5  0.024753793 -0.052063386  0.054757181     FG   restored       15
    ## 6          6  0.132027889  0.107077144  0.033575362     FL       corn       NA
    ## 7          7 -0.042177309  0.100432373 -0.083348303     FL       corn       NA
    ## 8          8  0.012829998 -0.029240145  0.019755980     FL    remnant       NA
    ## 9          9 -0.048029289 -0.037910042  0.024674421     FL   restored       40
    ## 10        10 -0.047531002 -0.032394265  0.016179264     FL   restored       36
    ## 11        11 -0.032019394 -0.015349526  0.017878940     FL   restored       35
    ## 12        12  0.070369149 -0.022703786 -0.026097774     FL   restored       10
    ## 13        13  0.014601133 -0.009117052 -0.041202654     FL   restored       10
    ## 14        14  0.027018121 -0.003352651 -0.048179202     FL   restored       10
    ## 15        15 -0.040888003 -0.041071752 -0.022777348     BM   restored       28
    ## 16        16 -0.029033616  0.101034837 -0.007893224     LP       corn       NA
    ## 17        17 -0.010295789 -0.010944908  0.003923329     LP    remnant       NA
    ## 18        18  0.001881795  0.008404165  0.035231885     LP   restored        4
    ## 19        19  0.066396377 -0.002797405  0.021616717     LP   restored        4
    ## 20        20 -0.045321553 -0.044288430 -0.057552793     BM    remnant       NA
    ## 21        21 -0.010605035 -0.038274923 -0.046884326     BM   restored       18
    ## 22        22 -0.024417886 -0.009405314  0.017526673     BM   restored        7
    ## 23        23 -0.029424159  0.005009133  0.025573301     BM   restored        2
    ## 24        24 -0.076762389  0.062360323  0.037484938     BM       corn       NA
    ## 25        25  0.011005267 -0.032483714  0.048923871     BM   restored       11
    ## 
    ## $broken_stick_plot

<img src="microbial_communities_files/figure-gfm/pcoa_amf_uni-1.png" style="display: block; margin: auto;" />

    ## 
    ## $permanova
    ## Permutation test for adonis under reduced model
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 1999
    ## 
    ## adonis2(formula = d ~ field_type, data = env, permutations = nperm, method = adonis_index, add = if (corr == "none") FALSE else "lingoes", strata = region)
    ##          Df SumOfSqs     R2      F Pr(>F)   
    ## Model     2  0.06937 0.1657 2.1847  0.006 **
    ## Residual 22  0.34929 0.8343                 
    ## Total    24  0.41866 1.0000                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## $pairwise_contrasts
    ## 
    ## 
    ## group1     group2        R2   F_value   df1   df2     p_value   p_value_adj
    ## ---------  --------  ------  --------  ----  ----  ----------  ------------
    ## restored   corn       0.237     5.892     1    19   0.0010000        0.0030
    ## restored   remnant    0.026     0.476     1    18   0.9680000        0.9680
    ## corn       remnant    0.383     4.339     1     7   0.0416667        0.0625

Three axes are significant by a broken stick model, between them
explaining 48.8% of the variation in AMF among fields. The most
substantial variation here is on the first axis (23%) with Axis 2
explaining 15% of the variation in AMF abundances. Testing the design
factor *field_type* (with *region* treated as a block using the `strata`
argument of `adonis2`) revealed a significant clustering
$(R^2=0.17,~p=0.006)$.

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

    ## Warning: Removed 9 rows containing missing values or values outside the scale range
    ## (`geom_text()`).

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
$(R^2_{Adj}=0.31,~p<0.05)$, and Axis 1 is close with
$(R^2_{Adj}=0.30,~p<0.05)$

Correlating age with axis scores isn’t appropriate because the axis
scores were produced with corn and remnant fields included. A better way
is to look at the Blue Mounds restored fields only. For now, we’ll
return to Bray-Curtis distance.

### PCoA with Blue Mounds restored fields, all subsamples

**Bray-Curtis distance used**. A Lingoes correction was applied to the
negative eigenvalues.

``` r
(pcoa_amf_resto_samps_bm <- pcoa_samps_bm_fun(spe$amf_samps, 
                                        distab$amf_resto_samps_bm, 
                                        sites_resto_bm, 
                                        corr="lingoes",
                                        df_name="BM restored, 18S gene, 97% OTU, BC dist."))
```

    ## Set of permutations < 'minperm'. Generating entire set.
    ## Set of permutations < 'minperm'. Generating entire set.

    ## $dataset
    ## [1] "BM restored, 18S gene, 97% OTU, BC dist."
    ## 
    ## $components_exceed_broken_stick
    ## [1] 5
    ## 
    ## $correction_note
    ## [1] "Lingoes correction applied to negative eigenvalues: D' = -0.5*D^2 - 0.119261876193705 , except diagonal elements"
    ## 
    ## $values
    ##   Dim Eigenvalues  Corr_eig Rel_corr_eig Broken_stick Cum_corr_eig Cum_br_stick
    ## 1   1   2.7073026 2.8265645   0.15970588   0.09442476    0.1597059   0.09442476
    ## 2   2   2.0569441 2.1762060   0.12295948   0.07314817    0.2826654   0.16757293
    ## 3   3   1.0088495 1.1281114   0.06374028   0.06250987    0.3464056   0.23008280
    ## 4   4   0.9040475 1.0233094   0.05781879   0.05541767    0.4042244   0.28550047
    ## 5   5   0.8490533 0.9683152   0.05471152   0.05009852    0.4589359   0.33559899
    ## 6   6   0.6152107 0.7344726   0.04149900   0.04584320    0.5004349   0.38144219
    ## 
    ## $eigenvalues
    ## [1] 16.0 12.3
    ## 
    ## $site_vectors
    ## # A tibble: 49 × 8
    ##    field_key sample_key  Axis.1  Axis.2  Axis.3  Axis.4  Axis.5 yr_since
    ##        <dbl> <chr>        <dbl>   <dbl>   <dbl>   <dbl>   <dbl>    <dbl>
    ##  1         1 1          -0.132  -0.0584 -0.134  -0.0808 -0.128        16
    ##  2         1 2          -0.160  -0.295  -0.102  -0.157  -0.140        16
    ##  3         1 4          -0.453   0.191  -0.121  -0.154   0.0811       16
    ##  4         1 5          -0.0740 -0.338  -0.300  -0.126   0.136        16
    ##  5         1 7          -0.272  -0.332  -0.0561 -0.147  -0.0969       16
    ##  6         1 8          -0.356   0.0832 -0.0336  0.163  -0.140        16
    ##  7         1 10         -0.376  -0.193  -0.0231 -0.124  -0.0660       16
    ##  8         2 1           0.0679 -0.293   0.317   0.122   0.109         3
    ##  9         2 3           0.147  -0.158   0.0775  0.162  -0.0474        3
    ## 10         2 5           0.0251 -0.416   0.205  -0.0773  0.106         3
    ## # ℹ 39 more rows
    ## 
    ## $broken_stick_plot

<img src="microbial_communities_files/figure-gfm/pcoa_amf_resto_samps_bm-1.png" style="display: block; margin: auto;" />

    ## 
    ## $permanova
    ## Permutation test for adonis under reduced model
    ## Permutation: free
    ## Number of permutations: 1999
    ## 
    ## adonis2(formula = d ~ field_key, data = env_w, permutations = nperm, method = adonis_index)
    ##          Df SumOfSqs      R2      F Pr(>F)    
    ## Model     1    1.413 0.11801 6.2883  5e-04 ***
    ## Residual 47   10.561 0.88199                  
    ## Total    48   11.974 1.00000                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## $vector_fit
    ## 
    ## ***VECTORS
    ## 
    ##            Axis.1   Axis.2   Axis.3   Axis.4   Axis.5     r2 Pr(>r)  
    ## yr_since -0.85500  0.39634 -0.18249  0.18450  0.21108 0.7779 0.0205 *
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Plots: field_key, plot permutation: free
    ## Permutation: none
    ## Number of permutations: 5039
    ## 
    ## 
    ## 
    ## $vector_fit_scores
    ##             Axis.1    Axis.2     Axis.3    Axis.4    Axis.5
    ## yr_since -0.754091 0.3495607 -0.1609538 0.1627285 0.1861653

Axis 1 explains 16% and axis 2 explains 12.3% of the variation in the
community data. Both axes are important based on the broken stick model
(5 relative corrected eigenvalues exceed the broken stick model). The
relatively low percent variation explained is partly due to the high
number of dimensions used when all samples from fields are included. The
fidelity of samples to fields was significant based on a permutation
test $(R^2=0.12,~p=5\times 10^{-4})$. In this case, the partial $R^2$
shows the proportion of sum of squares from the total. It is a low
number here because so much unexplained variation exists, resulting in a
high sum of squares that is outside the assignment of subsamples to
fields.

Years since restoration has a moderately strong correlation with
communities and was significant with a permutation test where samples
were constrained within fields to account for lack of independence
$(R^2=0.78,~p=0.02)$.

Let’s view an ordination plot with hulls around subsamples and a fitted
vector for field age overlaid.

``` r
centroid_amf_bm <- aggregate(cbind(Axis.1, Axis.2) ~ field_key, data = pcoa_amf_resto_samps_bm$site_vectors, mean) %>% 
    left_join(sites %>% select(field_key, yr_since), by = join_by(field_key))
hull_amf_bm <- pcoa_amf_resto_samps_bm$site_vectors %>% 
    group_by(field_key) %>% 
    slice(chull(Axis.1, Axis.2))
```

``` r
ggplot(pcoa_amf_resto_samps_bm$site_vectors, aes(x = Axis.1, y = Axis.2)) +
    geom_point(fill = "#5CBD92", shape = 21) +
    geom_polygon(data = hull_amf_bm, aes(group = as.character(field_key)), fill = "#5CBD92", alpha = 0.3) +
    geom_point(data = centroid_amf_bm, fill = "#5CBD92", size = 8, shape = 21) +
    geom_text(data = centroid_amf_bm, aes(label = yr_since)) +
    geom_segment(aes(x = 0, 
                     y = 0, 
                     xend = pcoa_amf_resto_samps_bm$vector_fit_scores[1] * 0.6, 
                     yend = pcoa_amf_resto_samps_bm$vector_fit_scores[2] * 0.6),
                 color = "blue", 
                 arrow = arrow(length = unit(3, "mm"))) +
    labs(
        x = paste0("Axis 1 (", pcoa_amf_resto_samps_bm$eigenvalues[1], "%)"),
        y = paste0("Axis 2 (", pcoa_amf_resto_samps_bm$eigenvalues[2], "%)"),
        title = paste0(
            "PCoA Ordination (",
            pcoa_amf_resto_samps_bm$dataset,
            ")"
        ),
        caption = "Text indicates years since restoration.\nYears since restoration significant at p<0.05."
    ) +
    theme_bw() +
    theme(legend.position = "none")
```

    ## Warning in geom_segment(aes(x = 0, y = 0, xend = pcoa_amf_resto_samps_bm$vector_fit_scores[1] * : All aesthetics have length 1, but the data has 49 rows.
    ## ℹ Please consider using `annotate()` or provide this layer with data containing
    ##   a single row.

<img src="microbial_communities_files/figure-gfm/amf_samps_bm_fig-1.png" style="display: block; margin: auto;" />

### PCoA with all fields and regions, all subsamples

**Bray-Curtis distance used.** This leverages the information from all
subsamples. Modifications to `how()` from package
[permute](https://cran.r-project.org/package=permute) allow for the more
complex design.

Negative eigenvalues were produced in trial runs (not shown). A Lingoes
correction was applied.

``` r
(pcoa_amf_samps <- pcoa_samps_fun(spe$amf_samps, 
                                  distab$amf_samps, 
                                  corr="lingoes", 
                                  df_name = "18S gene, 97% OTU"))
```

    ## 'nperm' >= set of all permutations: complete enumeration.

    ## Set of permutations < 'minperm'. Generating entire set.

    ## $dataset
    ## [1] "18S gene, 97% OTU"
    ## 
    ## $components_exceed_broken_stick
    ## [1] 10
    ## 
    ## $correction_note
    ## [1] "Lingoes correction applied to negative eigenvalues: D' = -0.5*D^2 - 0.364033325702051 , except diagonal elements"
    ## 
    ## $values
    ##    Dim Eigenvalues Corr_eig Rel_corr_eig Broken_stick Cum_corr_eig Cum_br_stick
    ## 1    1    7.676267 8.040301   0.07133171   0.03314101   0.07133171   0.03314101
    ## 2    2    5.335630 5.699663   0.05056611   0.02736066   0.12189782   0.06050167
    ## 3    3    5.273113 5.637146   0.05001148   0.02447049   0.17190930   0.08497216
    ## 4    4    4.262246 4.626279   0.04104329   0.02254371   0.21295259   0.10751587
    ## 5    5    3.074543 3.438577   0.03050627   0.02109862   0.24345886   0.12861449
    ## 6    6    2.473389 2.837422   0.02517296   0.01994255   0.26863182   0.14855704
    ## 7    7    2.286809 2.650843   0.02351767   0.01897916   0.29214949   0.16753620
    ## 8    8    2.072857 2.436891   0.02161954   0.01815340   0.31376903   0.18568960
    ## 9    9    1.761787 2.125820   0.01885979   0.01743085   0.33262882   0.20312045
    ## 10  10    1.607071 1.971104   0.01748719   0.01678859   0.35011601   0.21990904
    ## 11  11    1.406162 1.770195   0.01570477   0.01621056   0.36582078   0.23611960
    ## 
    ## $eigenvalues
    ## [1] 7.1 5.1
    ## 
    ## $site_vectors
    ## # A tibble: 175 × 16
    ##    field_key sample_key  Axis.1  Axis.2   Axis.3   Axis.4   Axis.5   Axis.6
    ##        <dbl> <chr>        <dbl>   <dbl>    <dbl>    <dbl>    <dbl>    <dbl>
    ##  1         1 1          -0.0904 -0.0316  0.0467  -0.0900  -0.178   -0.0274 
    ##  2         1 2          -0.0936 -0.297   0.115    0.0144  -0.119   -0.150  
    ##  3         1 4          -0.336   0.0638  0.274   -0.168   -0.128    0.0594 
    ##  4         1 5           0.0945 -0.341   0.00930 -0.245   -0.0715  -0.0994 
    ##  5         1 7          -0.177  -0.322   0.146    0.00950 -0.116   -0.00601
    ##  6         1 8          -0.353   0.0356  0.117   -0.0988  -0.0266   0.0307 
    ##  7         1 10         -0.318  -0.222   0.151   -0.00647 -0.176    0.00379
    ##  8         2 1          -0.0394 -0.281  -0.221    0.154    0.219    0.178  
    ##  9         2 3          -0.0273 -0.0965 -0.248    0.172   -0.0396   0.0450 
    ## 10         2 5          -0.0450 -0.440  -0.144    0.180    0.00349  0.0903 
    ## # ℹ 165 more rows
    ## # ℹ 8 more variables: Axis.7 <dbl>, Axis.8 <dbl>, Axis.9 <dbl>, Axis.10 <dbl>,
    ## #   field_name <chr>, region <chr>, field_type <ord>, yr_since <dbl>
    ## 
    ## $broken_stick_plot

<img src="microbial_communities_files/figure-gfm/pcoa_amf_samps-1.png" style="display: block; margin: auto;" />

    ## 
    ## $permanova
    ## Permutation test for adonis under reduced model
    ## Blocks:  region 
    ## Plots: field_key, plot permutation: free
    ## Permutation: none
    ## Number of permutations: 1999
    ## 
    ## adonis2(formula = d ~ field_type, data = env_w, permutations = gl_perm_design, method = adonis_index)
    ##           Df SumOfSqs      R2      F Pr(>F)    
    ## Model      2    5.437 0.11012 10.642  0.001 ***
    ## Residual 172   43.938 0.88988                  
    ## Total    174   49.375 1.00000                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## $pairwise_contrasts
    ##     group1  group2    R2 F_value df1 df2   p_value p_value_adj
    ## 1 restored    corn 0.057   8.781   1 145 0.0010000   0.0030000
    ## 2 restored remnant 0.008   1.087   1 138 0.9960000   0.9960000
    ## 3     corn remnant 0.113   7.771   1  61 0.1020408   0.1530612
    ## 
    ## $format
    ## [1] "pandoc"

``` r
write_delim(pcoa_amf_samps$permanova %>% round(., 3), "microbial_communities_files/pcoa_amf_samps_permanova.txt")
write_delim(pcoa_amf_samps$pairwise_contrasts %>% mutate(across(starts_with("p_value"), ~ round(.x, 3))), "microbial_communities_files/pcoa_amf_samps_pairwise.txt")
```

Axis 1 explains 7.1% and axis 2 explains 5.1% of the variation in the
community data. Both axes are important based on the broken stick model,
in fact, the broken stick model shows that 10 axes are important in
explaining variation with this dataset. The relatively low percent
variation explained on axes 1 and 2 is partly due to the high number of
dimensions used when all samples from fields are included. The fidelity
of samples to fields was strong based on a permutation test when
restricting permutations to fields (=plots in `how()`) within regions
(=blocks in `how()`) $(R^2=0.11,~p=0.001)$.

Let’s view an ordination plot with hulls around subsamples.

``` r
centroid_amf <- aggregate(cbind(Axis.1, Axis.2) ~ field_key, data = pcoa_amf_samps$site_vectors, mean) %>% 
    left_join(sites %>% select(field_key, yr_since, field_type, region), by = join_by(field_key))
hull_amf <- pcoa_amf_samps$site_vectors %>% 
    group_by(field_key) %>% 
    slice(chull(Axis.1, Axis.2))
```

``` r
amf_samps_fig <- 
    ggplot(pcoa_amf_samps$site_vectors, aes(x = Axis.1, y = Axis.2)) +
    geom_vline(xintercept = 0, linewidth = 0.1) +
    geom_hline(yintercept = 0, linewidth = 0.1) +
    geom_point(aes(fill = field_type), shape = 21, alpha = 0.8, color = "gray10") +
    geom_polygon(data = hull_amf, aes(group = field_key, fill = field_type), alpha = 0.3) +
    geom_point(data = centroid_amf, aes(fill = field_type, shape = region), size = 6) +
    geom_text(data = centroid_amf, aes(label = yr_since), size = 3) +
    labs(
        x = paste0("Axis 1 (", pcoa_amf_samps$eigenvalues[1], "%)"),
        y = paste0("Axis 2 (", pcoa_amf_samps$eigenvalues[2], "%)"),
        title = paste0(
            "PCoA Ordination (",
            pcoa_amf_samps$dataset,
            ")"
        ),
        caption = "Text indicates years since restoration."
    ) +
    lims(y = c(-0.60,0.34)) +
    scale_fill_discrete_qualitative(name = "Field Type", palette = "Harmonic") +
    scale_shape_manual(name = "Region", values = c(21, 22, 23, 24)) +
    theme_bw() +
    guides(fill = guide_legend(override.aes = list(shape = 21)))
```

``` r
(amf_samps_families_fig <- 
        amf_samps_fig +
        annotation_custom(
            ggplotGrob(
                pcoa_amf_bray$inset + 
                    theme(
                        plot.background = element_rect(colour = "black", fill = "gray90"), 
                        axis.title.y = element_text(size = 8)
                    )),
            xmin = -0.52,
            xmax = -0.10,
            ymin = -0.62,
            ymax = -0.34
        ))
```

    ## Warning: Removed 9 rows containing missing values or values outside the scale range
    ## (`geom_text()`).

<img src="microbial_communities_files/figure-gfm/amf_samps_families_fig-1.png" style="display: block; margin: auto;" />

### PCoA in Blue Mounds, all subsamples

This is as above with the diagnostics and permutation tests. Pairwise
contrasts among field types should be ignored here because there is no
replication.

``` r
(pcoa_amf_samps_bm <- pcoa_samps_fun(
    s = spe$amf_samps_bm,
    d = distab$amf_samps_bm,
    env = sites %>% filter(region == "BM"),
    corr = "lingoes",
    df_name = "Blue Mounds, 18S gene, 97% OTU"
))
```

    ## $dataset
    ## [1] "Blue Mounds, 18S gene, 97% OTU"
    ## 
    ## $components_exceed_broken_stick
    ## [1] 5
    ## 
    ## $correction_note
    ## [1] "Lingoes correction applied to negative eigenvalues: D' = -0.5*D^2 - 0.150085589665598 , except diagonal elements"
    ## 
    ## $values
    ##   Dim Eigenvalues  Corr_eig Rel_corr_eig Broken_stick Cum_corr_eig Cum_br_stick
    ## 1   1   3.8108399 3.9609255   0.15319463   0.07698793    0.1531946   0.07698793
    ## 2   2   2.3572107 2.5072963   0.09697338   0.06059449    0.2501680   0.13758242
    ## 3   3   1.8924637 2.0425493   0.07899861   0.05239777    0.3291666   0.18998019
    ## 4   4   1.1684699 1.3185555   0.05099708   0.04693329    0.3801637   0.23691348
    ## 5   5   0.9657084 1.1157940   0.04315498   0.04283493    0.4233187   0.27974840
    ## 6   6   0.7690392 0.9191248   0.03554851   0.03955624    0.4588672   0.31930464
    ## 
    ## $eigenvalues
    ## [1] 15.3  9.7
    ## 
    ## $site_vectors
    ## # A tibble: 63 × 11
    ##    field_key sample_key  Axis.1   Axis.2   Axis.3  Axis.4    Axis.5 field_name
    ##        <dbl> <chr>        <dbl>    <dbl>    <dbl>   <dbl>     <dbl> <chr>     
    ##  1         1 1          -0.115   0.00421 -0.0978  -0.0433  0.100    BBRP1     
    ##  2         1 2          -0.196  -0.173   -0.247   -0.119   0.107    BBRP1     
    ##  3         1 4          -0.363   0.287    0.00732 -0.130  -0.0930   BBRP1     
    ##  4         1 5          -0.0857 -0.107   -0.429    0.0821 -0.131    BBRP1     
    ##  5         1 7          -0.289  -0.189   -0.233   -0.148   0.0302   BBRP1     
    ##  6         1 8          -0.340   0.112    0.110    0.0592  0.112    BBRP1     
    ##  7         1 10         -0.380  -0.0993  -0.0970  -0.157   0.000142 BBRP1     
    ##  8         2 1           0.0234 -0.396    0.0981  -0.0453 -0.176    ERRP1     
    ##  9         2 3           0.0771 -0.273    0.0835   0.105   0.0428   ERRP1     
    ## 10         2 5          -0.0476 -0.460   -0.0759  -0.142  -0.158    ERRP1     
    ## # ℹ 53 more rows
    ## # ℹ 3 more variables: region <chr>, field_type <ord>, yr_since <dbl>
    ## 
    ## $broken_stick_plot

![](microbial_communities_files/figure-gfm/pcoa_amf_samps_bm-1.png)<!-- -->

    ## 
    ## $permanova
    ## Permutation test for adonis under reduced model
    ## Blocks:  region 
    ## Plots: field_key, plot permutation: free
    ## Permutation: none
    ## Number of permutations: 1999
    ## 
    ## adonis2(formula = d ~ field_type, data = env_w, permutations = gl_perm_design, method = adonis_index)
    ##          Df SumOfSqs     R2      F Pr(>F)  
    ## Model     2   2.8251 0.1707 6.1751 0.0925 .
    ## Residual 60  13.7251 0.8293                
    ## Total    62  16.5502 1.0000                
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## $pairwise_contrasts
    ##     group1  group2    R2 F_value df1 df2 p_value p_value_adj
    ## 1 restored remnant 0.053   3.021   1  54   0.389      0.5835
    ## 2 restored    corn 0.085   5.037   1  54   0.123      0.3690
    ## 3  remnant    corn 0.469  10.617   1  12   1.000      1.0000
    ## 
    ## $format
    ## [1] "pandoc"

Field type trends significant. Four axes significant.

### PCoA in Faville Grove, all subsamples

This is as above with the diagnostics and permutation tests. Pairwise
contrasts among field types should be ignored here because there is no
replication.

``` r
(pcoa_amf_samps_fg <- pcoa_samps_fun(
    s = spe$amf_samps_fg,
    d = distab$amf_samps_fg,
    env = sites %>% filter(region == "FG"),
    corr = "lingoes",
    df_name = "Faville Grove, 18S gene, 97% OTU"
))
```

    ## $dataset
    ## [1] "Faville Grove, 18S gene, 97% OTU"
    ## 
    ## $components_exceed_broken_stick
    ## [1] 3
    ## 
    ## $correction_note
    ## [1] "Lingoes correction applied to negative eigenvalues: D' = -0.5*D^2 - 0.0667727147236192 , except diagonal elements"
    ## 
    ## $values
    ##   Dim Eigenvalues  Corr_eig Rel_corr_eig Broken_stick Cum_corr_eig Cum_br_stick
    ## 1   1   1.8546304 1.9214032   0.28640834   0.18672314    0.2864083    0.1867231
    ## 2   2   0.9966864 1.0634591   0.15852142   0.13409156    0.4449298    0.3208147
    ## 3   3   0.6888312 0.7556039   0.11263189   0.10777577    0.5575616    0.4285905
    ## 4   4   0.4456716 0.5124443   0.07638601   0.09023191    0.6339477    0.5188224
    ## 
    ## $eigenvalues
    ## [1] 28.6 15.9
    ## 
    ## $site_vectors
    ## # A tibble: 21 × 9
    ##    field_key sample_key  Axis.1  Axis.2   Axis.3 field_name region field_type
    ##        <dbl> <chr>        <dbl>   <dbl>    <dbl> <chr>      <chr>  <ord>     
    ##  1         3 1          -0.319  -0.0878  0.246   FGC1       FG     corn      
    ##  2         3 2          -0.420   0.0464  0.201   FGC1       FG     corn      
    ##  3         3 4          -0.395  -0.138   0.0322  FGC1       FG     corn      
    ##  4         3 6          -0.375  -0.129  -0.145   FGC1       FG     corn      
    ##  5         3 7          -0.406  -0.0791 -0.00955 FGC1       FG     corn      
    ##  6         3 8          -0.317  -0.119  -0.0741  FGC1       FG     corn      
    ##  7         3 9          -0.402  -0.0936 -0.0739  FGC1       FG     corn      
    ##  8         4 1          -0.153   0.177  -0.275   FGREM1     FG     remnant   
    ##  9         4 2           0.0963  0.422   0.0928  FGREM1     FG     remnant   
    ## 10         4 3           0.239   0.338   0.251   FGREM1     FG     remnant   
    ## # ℹ 11 more rows
    ## # ℹ 1 more variable: yr_since <dbl>
    ## 
    ## $broken_stick_plot

![](microbial_communities_files/figure-gfm/pcoa_amf_samps_fg-1.png)<!-- -->

    ## 
    ## $permanova
    ## Permutation test for adonis under reduced model
    ## Blocks:  region 
    ## Plots: field_key, plot permutation: free
    ## Permutation: none
    ## Number of permutations: 5
    ## 
    ## adonis2(formula = d ~ field_type, data = env_w, permutations = gl_perm_design, method = adonis_index)
    ##          Df SumOfSqs      R2      F Pr(>F)
    ## Model     2   2.5767 0.47955 8.2929      1
    ## Residual 18   2.7964 0.52045              
    ## Total    20   5.3732 1.00000              
    ## 
    ## $pairwise_contrasts
    ##    group1   group2    R2 F_value df1 df2 p_value p_value_adj
    ## 1    corn  remnant 0.347   6.374   1  12       1           1
    ## 2    corn restored 0.451   9.857   1  12       1           1
    ## 3 remnant restored 0.328   5.863   1  12       1           1
    ## 
    ## $format
    ## [1] "pandoc"

Field type is not significant here. Three significant axes.

### PCoA in Fermilab, all subsamples

This is as above with the diagnostics and permutation tests. Pairwise
contrasts among field types should be ignored here because there is no
replication.

``` r
(pcoa_amf_samps_fl <- pcoa_samps_fun(
    s = spe$amf_samps_fl,
    d = distab$amf_samps_fl,
    env = sites %>% filter(region == "FL"),
    corr = "lingoes",
    df_name = "Fermilab, 18S gene, 97% OTU"
))
```

    ## $dataset
    ## [1] "Fermilab, 18S gene, 97% OTU"
    ## 
    ## $components_exceed_broken_stick
    ## [1] 6
    ## 
    ## $correction_note
    ## [1] "Lingoes correction applied to negative eigenvalues: D' = -0.5*D^2 - 0.17165314762112 , except diagonal elements"
    ## 
    ## $values
    ##   Dim Eigenvalues  Corr_eig Rel_corr_eig Broken_stick Cum_corr_eig Cum_br_stick
    ## 1   1   3.2625003 3.4341535   0.12298932   0.07698793    0.1229893   0.07698793
    ## 2   2   2.2829425 2.4545957   0.08790785   0.06059449    0.2108972   0.13758242
    ## 3   3   1.7309910 1.9026442   0.06814049   0.05239777    0.2790377   0.18998019
    ## 4   4   1.3529051 1.5245582   0.05459988   0.04693329    0.3336376   0.23691348
    ## 5   5   1.1886826 1.3603358   0.04871849   0.04283493    0.3823560   0.27974840
    ## 6   6   0.9394285 1.1110816   0.03979181   0.03955624    0.4221478   0.31930464
    ## 7   7   0.8120811 0.9837343   0.03523104   0.03682400    0.4573789   0.35612864
    ## 
    ## $eigenvalues
    ## [1] 12.3  8.8
    ## 
    ## $site_vectors
    ## # A tibble: 63 × 12
    ##    field_key sample_key Axis.1  Axis.2   Axis.3   Axis.4  Axis.5   Axis.6
    ##        <dbl> <chr>       <dbl>   <dbl>    <dbl>    <dbl>   <dbl>    <dbl>
    ##  1         6 2          -0.360 -0.0408  0.159    0.172    0.362   0.00206
    ##  2         6 4          -0.296  0.120  -0.176    0.280    0.181   0.0518 
    ##  3         6 5          -0.247 -0.0736 -0.383    0.298   -0.0463 -0.0424 
    ##  4         6 6          -0.208  0.108   0.00242  0.189    0.314   0.100  
    ##  5         6 7          -0.266 -0.119  -0.284    0.336    0.0996 -0.0770 
    ##  6         6 8          -0.267  0.0156  0.0106  -0.113   -0.0266 -0.0589 
    ##  7         6 9          -0.418 -0.117   0.0802  -0.00529  0.117  -0.00970
    ##  8         7 2          -0.271 -0.401  -0.256   -0.121   -0.108   0.0959 
    ##  9         7 3          -0.332 -0.325  -0.0808  -0.238   -0.122   0.229  
    ## 10         7 4          -0.379 -0.199   0.0903  -0.0298  -0.103  -0.175  
    ## # ℹ 53 more rows
    ## # ℹ 4 more variables: field_name <chr>, region <chr>, field_type <ord>,
    ## #   yr_since <dbl>
    ## 
    ## $broken_stick_plot

![](microbial_communities_files/figure-gfm/pcoa_amf_samps_fl-1.png)<!-- -->

    ## 
    ## $permanova
    ## Permutation test for adonis under reduced model
    ## Blocks:  region 
    ## Plots: field_key, plot permutation: free
    ## Permutation: none
    ## Number of permutations: 1999
    ## 
    ## adonis2(formula = d ~ field_type, data = env_w, permutations = gl_perm_design, method = adonis_index)
    ##          Df SumOfSqs      R2      F Pr(>F)  
    ## Model     2   2.8438 0.16458 5.9099  0.045 *
    ## Residual 60  14.4360 0.83542                
    ## Total    62  17.2799 1.00000                
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## $pairwise_contrasts
    ##    group1   group2    R2 F_value df1 df2   p_value p_value_adj
    ## 1    corn  remnant 0.194   4.574   1  19 0.3333333      0.5000
    ## 2    corn restored 0.103   6.222   1  54 0.0390000      0.1170
    ## 3 remnant restored 0.033   1.607   1  47 0.7095000      0.7095
    ## 
    ## $format
    ## [1] "pandoc"

Field type is again significant by permutation test. Six axes are
significant.

### PCoA in Lake Petite Prairie, all subsamples

This is as above with the diagnostics and permutation tests. Pairwise
contrasts among field types should be ignored here because there is no
replication.

``` r
(pcoa_amf_samps_lp <- pcoa_samps_fun(
    s = spe$amf_samps_lp,
    d = distab$amf_samps_lp,
    env = sites %>% filter(region == "LP"),
    corr = "lingoes",
    df_name = "Lake Petite Prairie, 18S gene, 97% OTU"
))
```

    ## $dataset
    ## [1] "Lake Petite Prairie, 18S gene, 97% OTU"
    ## 
    ## $components_exceed_broken_stick
    ## [1] 3
    ## 
    ## $correction_note
    ## [1] "Lingoes correction applied to negative eigenvalues: D' = -0.5*D^2 - 0.0749881482986054 , except diagonal elements"
    ## 
    ## $values
    ##   Dim Eigenvalues  Corr_eig Rel_corr_eig Broken_stick Cum_corr_eig Cum_br_stick
    ## 1   1   1.5851379 1.6601261   0.20920798   0.14824691    0.2092080    0.1482469
    ## 2   2   1.0792307 1.1542189   0.14545389   0.10978537    0.3546619    0.2580323
    ## 3   3   0.6889139 0.7639021   0.09626643   0.09055460    0.4509283    0.3485869
    ## 4   4   0.5418262 0.6168143   0.07773053   0.07773409    0.5286588    0.4263210
    ## 
    ## $eigenvalues
    ## [1] 20.9 14.5
    ## 
    ## $site_vectors
    ## # A tibble: 28 × 9
    ##    field_key sample_key  Axis.1  Axis.2  Axis.3 field_name region field_type
    ##        <dbl> <chr>        <dbl>   <dbl>   <dbl> <chr>      <chr>  <ord>     
    ##  1        16 1          -0.363   0.138   0.0318 LPC1       LP     corn      
    ##  2        16 2          -0.139   0.222   0.0813 LPC1       LP     corn      
    ##  3        16 4          -0.379   0.131   0.0340 LPC1       LP     corn      
    ##  4        16 5          -0.275   0.105  -0.0396 LPC1       LP     corn      
    ##  5        16 6          -0.347   0.174   0.0623 LPC1       LP     corn      
    ##  6        16 7          -0.298   0.188   0.0541 LPC1       LP     corn      
    ##  7        16 8          -0.323   0.168   0.0336 LPC1       LP     corn      
    ##  8        17 1          -0.0821 -0.0824  0.138  LPREM1     LP     remnant   
    ##  9        17 2           0.216   0.472  -0.164  LPREM1     LP     remnant   
    ## 10        17 4           0.0223 -0.179  -0.0203 LPREM1     LP     remnant   
    ## # ℹ 18 more rows
    ## # ℹ 1 more variable: yr_since <dbl>
    ## 
    ## $broken_stick_plot

![](microbial_communities_files/figure-gfm/pcoa_amf_samps_lp-1.png)<!-- -->

    ## 
    ## $permanova
    ## Permutation test for adonis under reduced model
    ## Blocks:  region 
    ## Plots: field_key, plot permutation: free
    ## Permutation: none
    ## Number of permutations: 23
    ## 
    ## adonis2(formula = d ~ field_type, data = env_w, permutations = gl_perm_design, method = adonis_index)
    ##          Df SumOfSqs     R2      F Pr(>F)
    ## Model     2   1.5634 0.2645 4.4953    0.5
    ## Residual 25   4.3472 0.7355              
    ## Total    27   5.9106 1.0000              
    ## 
    ## $pairwise_contrasts
    ##    group1   group2    R2 F_value df1 df2   p_value p_value_adj
    ## 1    corn  remnant 0.292   4.941   1  12 1.0000000           1
    ## 2    corn restored 0.255   6.510   1  19 0.3333333           1
    ## 3 remnant restored 0.082   1.708   1  19 1.0000000           1
    ## 
    ## $format
    ## [1] "pandoc"

Field type not significant with three important axes.

Let’s view an ordination plot with hulls around subsamples for each
indidual region.

### PCoA ordination, all regions, all subsamples

``` r
pcoa_amf_site_vectors <- bind_rows(
    list(
        `Blue Mounds`   = pcoa_amf_samps_bm$site_vectors,
        `Faville Grove` = pcoa_amf_samps_fg$site_vectors,
        `Fermilab`      = pcoa_amf_samps_fl$site_vectors,
        `Lake Petite`   = pcoa_amf_samps_lp$site_vectors
    ),
    .id = "place"
)
pcoa_amf_eigenvalues <- bind_rows(
    list(
        `Blue Mounds`   = pcoa_amf_samps_bm$eigenvalues,
        `Faville Grove` = pcoa_amf_samps_fg$eigenvalues,
        `Fermilab`      = pcoa_amf_samps_fl$eigenvalues,
        `Lake Petite`   = pcoa_amf_samps_lp$eigenvalues
    ),
    .id = "place"
) %>% 
    mutate(axis = c(1,2)) %>% 
    pivot_longer(cols = 1:4, names_to = "place", values_to = "eigenvalue") %>% 
    select(place, axis, eigenvalue) %>% 
    arrange(place, axis) %>% 
    pivot_wider(names_from = axis, names_prefix = "axis_", values_from = eigenvalue)
centroid_regions_amf <- aggregate(cbind(Axis.1, Axis.2) ~ place + field_key, data = pcoa_amf_site_vectors, mean) %>% 
    left_join(sites %>% select(field_key, yr_since, field_type, region), by = join_by(field_key))
hull_regions_amf <- pcoa_amf_site_vectors %>% 
    group_by(place, field_key) %>% 
    slice(chull(Axis.1, Axis.2))
```

``` r
(amf_samps_regions_fig <- 
    ggplot(pcoa_amf_site_vectors, aes(x = Axis.1, y = Axis.2)) +
    facet_wrap(vars(place), scales = "free") +
    geom_vline(xintercept = 0, linewidth = 0.1) +
    geom_hline(yintercept = 0, linewidth = 0.1) +
    geom_point(aes(fill = field_type), shape = 21, alpha = 0.8, color = "gray10") +
    geom_polygon(data = hull_regions_amf, aes(group = field_key, fill = field_type), alpha = 0.3) +
    geom_point(data = centroid_regions_amf, aes(fill = field_type, shape = region), size = 5) +
    geom_text(data = centroid_regions_amf, aes(label = yr_since), size = 2.5) +
    labs(
        x = paste0("Axis 1"),
        y = paste0("Axis 2"),
        caption = "18S gene (AMF). Text indicates years since restoration."
    ) +
    scale_fill_discrete_qualitative(name = "Field Type", palette = "Harmonic") +
    scale_shape_manual(name = "Region", values = c(21, 22, 23, 24)) +
    theme_bw() +
    guides(fill = guide_legend(override.aes = list(shape = 21))))
```

<img src="microbial_communities_files/figure-gfm/amf_samps_regions_fig-1.png" style="display: block; margin: auto;" />

The eigenvalues are shown below:

``` r
kable(pcoa_amf_eigenvalues, format = "pandoc")
```

| place         | axis_1 | axis_2 |
|:--------------|-------:|-------:|
| Blue Mounds   |   15.3 |    9.7 |
| Faville Grove |   28.6 |   15.9 |
| Fermilab      |   12.3 |    8.8 |
| Lake Petite   |   20.9 |   14.5 |

``` r
write_csv(pcoa_amf_eigenvalues, file = "microbial_communities_files/pcoa_amf_eig.csv")
```

Let’s view and save a plot that shows all the data together and broken
out by regions.

``` r
grid.arrange(
    amf_samps_families_fig + labs(caption = "") + theme(plot.title = element_blank()), 
    amf_samps_regions_fig + labs(caption = "") + theme(legend.position = "none"), 
    ncol = 1,
    heights = c(1.1,0.9)
)
```

<img src="microbial_communities_files/figure-gfm/amf_samps_unified_fig-1.png" style="display: block; margin: auto;" />

Then, we’ll follow up with panels showing trends with the most abundant
guilds.

``` r
spe_meta$amf %>%
    filter(family %in% c("Claroideoglomeraceae", "Paraglomeraceae", "Diversisporaceae", "Gigasporaceae")) %>% 
    mutate(field_type = factor(field_type, ordered = TRUE, 
                               levels = c("corn", "restored", "remnant"))) %>%
    group_by(region, family, field_type, field_name) %>%
    summarize(sum_seq_abund = sum(seq_abund), .groups = "drop_last") %>% 
    summarize(avg_seq_abund = mean(sum_seq_abund), .groups = "drop") %>%
    ggplot(aes(x = region, y = avg_seq_abund, fill = field_type)) +
    facet_wrap(vars(family), scales = "free_y") +
    geom_col(position = "dodge") +
    labs(y = "Sequence abundance (avg)") +
    scale_fill_discrete_qualitative(name = "Field Type", palette = "Harmonic") +
    theme_bw() +
    theme(axis.title.x = element_blank())
```

<img src="microbial_communities_files/figure-gfm/amf_guilds_regions_fig-1.png" style="display: block; margin: auto;" />
