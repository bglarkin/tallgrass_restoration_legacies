Microbial data: community differences
================
Beau Larkin

Last updated: 10 September, 2024

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
    ## 1   1    1.296359   0.18839901   0.15733159 0.1883990      0.1573316
    ## 2   2    0.791862   0.11508078   0.11566492 0.3034798      0.2729965
    ## 3   3    0.591334   0.08593818   0.09483159 0.3894180      0.3678281
    ## 
    ## $eigenvalues
    ## [1] 18.8 11.5
    ## 
    ## $site_vectors
    ##    field_key      Axis.1       Axis.2 region field_type yr_since
    ## 1          1  0.22586722 -0.052356698     BM   restored       16
    ## 2          2 -0.10660664  0.012741658     BM   restored        3
    ## 3          3 -0.33278089  0.004021646     FG       corn       NA
    ## 4          4  0.10685336 -0.302955832     FG    remnant       NA
    ## 5          5 -0.05350523 -0.292028919     FG   restored       15
    ## 6          6 -0.29754073  0.116424095     FL       corn       NA
    ## 7          7 -0.32889451  0.012010862     FL       corn       NA
    ## 8          8  0.15570034 -0.225578056     FL    remnant       NA
    ## 9          9  0.25001099 -0.188920864     FL   restored       40
    ## 10        10  0.32102356 -0.072299858     FL   restored       36
    ## 11        11  0.21968673 -0.165329590     FL   restored       35
    ## 12        12  0.20801552  0.242040565     FL   restored       10
    ## 13        13  0.18036935  0.124368237     FL   restored       10
    ## 14        14  0.18886639  0.256707526     FL   restored       10
    ## 15        15  0.24405130 -0.050732509     BM   restored       28
    ## 16        16 -0.38142478  0.111641388     LP       corn       NA
    ## 17        17  0.06697260  0.187011960     LP    remnant       NA
    ## 18        18 -0.24900969  0.072302073     LP   restored        4
    ## 19        19 -0.14762866  0.092191002     LP   restored        4
    ## 20        20  0.23557693  0.314408550     BM    remnant       NA
    ## 21        21  0.18422378  0.286002819     BM   restored       18
    ## 22        22 -0.11000262 -0.210284192     BM   restored        7
    ## 23        23 -0.22255414 -0.084904366     BM   restored        2
    ## 24        24 -0.34370145  0.020532440     BM       corn       NA
    ## 25        25 -0.01356873 -0.207013936     BM   restored       11
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
    ## Model     2   1.2362 0.17966 2.4091  5e-04 ***
    ## Residual 22   5.6447 0.82034                  
    ## Total    24   6.8809 1.00000                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## $pairwise_contrasts
    ## 
    ## 
    ## group1     group2        R2   F_value   df1   df2     p_value   p_value_adj
    ## ---------  --------  ------  --------  ----  ----  ----------  ------------
    ## restored   corn       0.157     3.532     1    19   0.0010000        0.0030
    ## restored   remnant    0.057     1.090     1    18   0.1000000        0.1000
    ## corn       remnant    0.296     2.940     1     7   0.0416667        0.0625

Axis 1 explains 18.8% of the variation and is the only eigenvalue that
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
    ## 1   1    2.141085   0.11526197   0.08352022 0.1152620     0.08352022
    ## 2   2    1.427693   0.07685761   0.06533840 0.1921196     0.14885863
    ## 3   3    1.029600   0.05542691   0.05624749 0.2475465     0.20510612
    ## 
    ## $eigenvalues
    ## [1] 11.5  7.7
    ## 
    ## $site_vectors
    ## # A tibble: 56 × 5
    ##    field_key sample_key  Axis.1  Axis.2 yr_since
    ##        <dbl> <chr>        <dbl>   <dbl>    <dbl>
    ##  1         1 1          -0.192  -0.195        16
    ##  2         1 2          -0.214  -0.0466       16
    ##  3         1 4          -0.184  -0.0670       16
    ##  4         1 5          -0.278  -0.141        16
    ##  5         1 6          -0.324  -0.163        16
    ##  6         1 7          -0.100  -0.0791       16
    ##  7         1 9          -0.211  -0.152        16
    ##  8         1 10         -0.245  -0.0307       16
    ##  9         2 1           0.244   0.221         3
    ## 10         2 2           0.0916  0.0800        3
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
    ## Model     1   0.7814 0.04207 2.3714  5e-04 ***
    ## Residual 54  17.7944 0.95793                  
    ## Total    55  18.5758 1.00000                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## $vector_fit
    ## 
    ## ***VECTORS
    ## 
    ##             Axis.1    Axis.2     r2 Pr(>r)  
    ## yr_since -0.997500  0.070667 0.7273  0.024 *
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
    ## yr_since -0.8506727 0.06026507

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
    ## [1] 10
    ## 
    ## $correction_note
    ## [1] "Lingoes correction applied to negative eigenvalues: D' = -0.5*D^2 - 0.0586364536544847 , except diagonal elements"
    ## 
    ## $values
    ##    Dim Eigenvalues Corr_eig Rel_corr_eig Broken_stick Cum_corr_eig Cum_br_stick
    ## 1    1    6.650014 6.708650   0.08109678   0.02963639   0.08109678   0.02963639
    ## 2    2    4.129312 4.187948   0.05062556   0.02458589   0.13172234   0.05422228
    ## 3    3    3.328468 3.387105   0.04094464   0.02206064   0.17266698   0.07628292
    ## 4    4    2.646543 2.705180   0.03270127   0.02037713   0.20536825   0.09666005
    ## 5    5    2.154798 2.213435   0.02675686   0.01911451   0.23212511   0.11577456
    ## 6    6    2.024480 2.083117   0.02518153   0.01810441   0.25730664   0.13387896
    ## 7    7    1.695887 1.754523   0.02120936   0.01726266   0.27851600   0.15114162
    ## 8    8    1.478246 1.536883   0.01857844   0.01654115   0.29709444   0.16768277
    ## 9    9    1.332238 1.390874   0.01681343   0.01590984   0.31390787   0.18359262
    ## 10  10    1.304042 1.362678   0.01647259   0.01534867   0.33038046   0.19894129
    ## 11  11    1.162318 1.220955   0.01475937   0.01484362   0.34513984   0.21378492
    ## 
    ## $eigenvalues
    ## [1] 8.1 5.1
    ## 
    ## $site_vectors
    ## # A tibble: 200 × 16
    ##    field_key sample_key  Axis.1   Axis.2  Axis.3  Axis.4  Axis.5  Axis.6  Axis.7
    ##        <dbl> <chr>        <dbl>    <dbl>   <dbl>   <dbl>   <dbl>   <dbl>   <dbl>
    ##  1         1 1           0.178  -0.125   -0.0387 -0.0709  0.0742 -0.109  -0.154 
    ##  2         1 2           0.219   0.0781  -0.0188 -0.0121  0.0119 -0.217  -0.0249
    ##  3         1 4           0.189  -0.0433  -0.0129  0.0214 -0.0201 -0.123  -0.0628
    ##  4         1 5           0.175  -0.0120  -0.126  -0.202  -0.0190  0.0235 -0.107 
    ##  5         1 6           0.205   0.0126  -0.118  -0.169   0.0227 -0.177  -0.268 
    ##  6         1 7           0.0467 -0.0510   0.0236 -0.109  -0.0366  0.0205 -0.125 
    ##  7         1 9           0.117  -0.0848  -0.0643 -0.0979  0.0631 -0.0364 -0.130 
    ##  8         1 10          0.186   0.0157  -0.0734 -0.0480 -0.0242 -0.0235 -0.0929
    ##  9         2 1          -0.0528 -0.00984  0.329   0.0503 -0.0784  0.0143  0.0496
    ## 10         2 2          -0.0460  0.0476   0.0864  0.0210 -0.134  -0.136  -0.0177
    ## # ℹ 190 more rows
    ## # ℹ 7 more variables: Axis.8 <dbl>, Axis.9 <dbl>, Axis.10 <dbl>,
    ## #   field_name <chr>, region <chr>, field_type <ord>, yr_since <dbl>
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
    ##           Df SumOfSqs      R2     F Pr(>F)    
    ## Model      2    5.863 0.08251 8.858  5e-04 ***
    ## Residual 197   65.193 0.91749                 
    ## Total    199   71.055 1.00000                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## $pairwise_contrasts
    ##     group1  group2    R2 F_value df1 df2   p_value p_value_adj
    ## 1 restored    corn 0.067  11.831   1 166 0.0005000      0.0015
    ## 2 restored remnant 0.018   2.942   1 158 0.1715000      0.1715
    ## 3     corn remnant 0.139  11.338   1  70 0.1428571      0.1715
    ## 
    ## $format
    ## [1] "pandoc"

``` r
write_delim(pcoa_its_samps$permanova, "microbial_communities_files/pcoa_its_samps_permanova.txt")
write_delim(pcoa_its_samps$pairwise_contrasts %>% mutate(across(starts_with("p_value"), ~ round(.x, 3))), "microbial_communities_files/pcoa_its_samps_pairwise.txt")
```

Axis 1 explains 8.1% and axis 2 explains 5.1% of the variation in the
community data. Both axes are important based on the broken stick model,
in fact, the broken stick model shows that 10 axes are important in
explaining variation with this dataset. The relatively low percent
variation explained on axes 1 and 2 is partly due to the high number of
dimensions used when all samples from fields are included. The fidelity
of samples to fields was strong based on a permutation test when
restricting permutations to fields (=plots in `how()`) within regions
(=blocks in `how()`) $(R^2=0.08,~p=5\times 10^{-4})$.

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
    ggplot(pcoa_its_samps$site_vectors, aes(x = Axis.1, y = Axis.2)) +
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
    ## 1   1   2.9643344   0.11598498   0.06826650 0.1159850      0.0682665
    ## 2   2   1.8319853   0.07167976   0.05418199 0.1876647      0.1224485
    ## 3   3   1.6118996   0.06306851   0.04713974 0.2507332      0.1695882
    ## 4   4   1.1131856   0.04355541   0.04244490 0.2942887      0.2120331
    ## 5   5   0.9934423   0.03887024   0.03892377 0.3331589      0.2509569
    ## 
    ## $eigenvalues
    ## [1] 11.6  7.2
    ## 
    ## $site_vectors
    ## # A tibble: 72 × 10
    ##    field_key sample_key  Axis.1   Axis.2  Axis.3  Axis.4 field_name region
    ##        <dbl> <chr>        <dbl>    <dbl>   <dbl>   <dbl> <chr>      <chr> 
    ##  1         1 1          -0.0849  0.131   -0.290   0.0297 BBRP1      BM    
    ##  2         1 2          -0.159   0.0973  -0.154  -0.0364 BBRP1      BM    
    ##  3         1 4          -0.0893  0.132   -0.200  -0.0937 BBRP1      BM    
    ##  4         1 5          -0.207  -0.00291 -0.219   0.220  BBRP1      BM    
    ##  5         1 6          -0.193   0.0393  -0.307   0.0346 BBRP1      BM    
    ##  6         1 7          -0.0286  0.0602  -0.139   0.142  BBRP1      BM    
    ##  7         1 9          -0.103   0.0615  -0.253   0.0722 BBRP1      BM    
    ##  8         1 10         -0.167   0.0759  -0.162   0.0434 BBRP1      BM    
    ##  9         2 1           0.189   0.232    0.223  -0.0882 ERRP1      BM    
    ## 10         2 2           0.112   0.0186   0.0634 -0.0607 ERRP1      BM    
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
    ##          Df SumOfSqs      R2      F Pr(>F)  
    ## Model     2   3.0705 0.12014 4.7107 0.0265 *
    ## Residual 69  22.4875 0.87986                
    ## Total    71  25.5579 1.00000                
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## $pairwise_contrasts
    ##     group1  group2    R2 F_value df1 df2 p_value p_value_adj
    ## 1 restored remnant 0.068   4.537   1  62  0.1245     0.18675
    ## 2 restored    corn 0.069   4.609   1  62  0.1125     0.18675
    ## 3  remnant    corn 0.302   6.049   1  14  1.0000     1.00000
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
    ## 1   1   1.9261547   0.26589636   0.16236050 0.2658964      0.1623605
    ## 2   2   0.8595206   0.11865267   0.11888224 0.3845490      0.2812427
    ## 3   3   0.4437595   0.06125886   0.09714311 0.4458079      0.3783858
    ## 
    ## $eigenvalues
    ## [1] 26.6 11.9
    ## 
    ## $site_vectors
    ## # A tibble: 24 × 8
    ##    field_key sample_key Axis.1   Axis.2 field_name region field_type yr_since
    ##        <dbl> <chr>       <dbl>    <dbl> <chr>      <chr>  <ord>         <dbl>
    ##  1         3 1          -0.367 -0.00792 FGC1       FG     corn             NA
    ##  2         3 2          -0.407 -0.0850  FGC1       FG     corn             NA
    ##  3         3 3          -0.429 -0.0579  FGC1       FG     corn             NA
    ##  4         3 5          -0.331 -0.00673 FGC1       FG     corn             NA
    ##  5         3 6          -0.409 -0.0442  FGC1       FG     corn             NA
    ##  6         3 7          -0.420  0.00875 FGC1       FG     corn             NA
    ##  7         3 9          -0.414 -0.0761  FGC1       FG     corn             NA
    ##  8         3 10         -0.382  0.0102  FGC1       FG     corn             NA
    ##  9         4 1           0.261 -0.268   FGREM1     FG     remnant          NA
    ## 10         4 2           0.294 -0.248   FGREM1     FG     remnant          NA
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
    ## Model     2   2.7135 0.37458 6.2887      1
    ## Residual 21   4.5305 0.62542              
    ## Total    23   7.2440 1.00000              
    ## 
    ## $pairwise_contrasts
    ##    group1   group2    R2 F_value df1 df2 p_value p_value_adj
    ## 1    corn  remnant 0.343   7.298   1  14       1           1
    ## 2    corn restored 0.341   7.246   1  14       1           1
    ## 3 remnant restored 0.227   4.106   1  14       1           1
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
    ## [1] "Lingoes correction applied to negative eigenvalues: D' = -0.5*D^2 - 0.0178526765158346 , except diagonal elements"
    ## 
    ## $values
    ##   Dim Eigenvalues Corr_eig Rel_corr_eig Broken_stick Cum_corr_eig Cum_br_stick
    ## 1   1    3.553783 3.571636   0.14403460   0.06904053    0.1440346   0.06904053
    ## 2   2    2.544772 2.562624   0.10334384   0.05475481    0.2473784   0.12379534
    ## 3   3    1.073956 1.091809   0.04402976   0.04761195    0.2914082   0.17140729
    ## 
    ## $eigenvalues
    ## [1] 14.4 10.3
    ## 
    ## $site_vectors
    ## # A tibble: 72 × 8
    ##    field_key sample_key Axis.1  Axis.2 field_name region field_type yr_since
    ##        <dbl> <chr>       <dbl>   <dbl> <chr>      <chr>  <ord>         <dbl>
    ##  1         6 1          -0.435  0.0237 FLC1       FL     corn             NA
    ##  2         6 2          -0.400 -0.0113 FLC1       FL     corn             NA
    ##  3         6 4          -0.430  0.0373 FLC1       FL     corn             NA
    ##  4         6 5          -0.422  0.0325 FLC1       FL     corn             NA
    ##  5         6 6          -0.298 -0.0221 FLC1       FL     corn             NA
    ##  6         6 7          -0.384  0.0401 FLC1       FL     corn             NA
    ##  7         6 9          -0.345  0.0574 FLC1       FL     corn             NA
    ##  8         6 10         -0.345 -0.0209 FLC1       FL     corn             NA
    ##  9         7 1          -0.432  0.106  FLC2       FL     corn             NA
    ## 10         7 3          -0.404  0.119  FLC2       FL     corn             NA
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
    ## Model     2   4.2493 0.18059 7.6036 0.0105 *
    ## Residual 69  19.2803 0.81941                
    ## Total    71  23.5295 1.00000                
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## $pairwise_contrasts
    ##    group1   group2    R2 F_value df1 df2   p_value p_value_adj
    ## 1    corn  remnant 0.268   8.063   1  22 0.3333333       0.426
    ## 2    corn restored 0.157  11.565   1  62 0.0420000       0.126
    ## 3 remnant restored 0.049   2.787   1  54 0.4260000       0.426
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
    ## 1   1   1.3451849   0.15690313   0.12991114 0.1569031      0.1299111
    ## 2   2   0.9382020   0.10943241   0.09765307 0.2663355      0.2275642
    ## 3   3   0.7176612   0.08370841   0.08152404 0.3500440      0.3090882
    ## 4   4   0.4858044   0.05666450   0.07077135 0.4067085      0.3798596
    ## 
    ## $eigenvalues
    ## [1] 15.7 10.9
    ## 
    ## $site_vectors
    ## # A tibble: 32 × 9
    ##    field_key sample_key Axis.1 Axis.2  Axis.3 field_name region field_type
    ##        <dbl> <chr>       <dbl>  <dbl>   <dbl> <chr>      <chr>  <ord>     
    ##  1        16 1          -0.236 0.0935 -0.184  LPC1       LP     corn      
    ##  2        16 3          -0.225 0.122  -0.171  LPC1       LP     corn      
    ##  3        16 5          -0.259 0.0900 -0.177  LPC1       LP     corn      
    ##  4        16 6          -0.266 0.102  -0.141  LPC1       LP     corn      
    ##  5        16 7          -0.248 0.0775 -0.143  LPC1       LP     corn      
    ##  6        16 8          -0.200 0.0340 -0.220  LPC1       LP     corn      
    ##  7        16 9          -0.249 0.0783 -0.119  LPC1       LP     corn      
    ##  8        16 10         -0.185 0.0531 -0.136  LPC1       LP     corn      
    ##  9        17 1           0.246 0.144   0.0949 LPREM1     LP     remnant   
    ## 10        17 2           0.291 0.204   0.0188 LPREM1     LP     remnant   
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
    ## Model     2   1.9325 0.22541 4.2195 0.1667
    ## Residual 29   6.6409 0.77459              
    ## Total    31   8.5733 1.00000              
    ## 
    ## $pairwise_contrasts
    ##    group1   group2    R2 F_value df1 df2   p_value p_value_adj
    ## 1    corn  remnant 0.265   5.044   1  14 1.0000000         1.0
    ## 2    corn restored 0.160   4.195   1  22 0.3333333         0.5
    ## 3 remnant restored 0.146   3.747   1  22 0.3333333         0.5
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
| Blue Mounds   |   11.6 |    7.2 |
| Faville Grove |   26.6 |   11.9 |
| Fermilab      |   14.4 |   10.3 |
| Lake Petite   |   15.7 |   10.9 |

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
    ## [1] 4
    ## 
    ## $correction_note
    ## [1] "No correction was applied to the negative eigenvalues"
    ## 
    ## $values
    ##   Dim Eigenvalues Relative_eig Rel_corr_eig Broken_stick Cum_corr_eig
    ## 1   1   1.2051279   0.27393921   0.25251062   0.16236050    0.2525106
    ## 2   2   0.7864744   0.17877454   0.16612548   0.11888224    0.4186361
    ## 3   3   0.6176689   0.14040314   0.13129409   0.09714311    0.5499302
    ## 4   4   0.4131360   0.09391049   0.08909067   0.08265036    0.6390208
    ## 5   5   0.3148359   0.07156577   0.06880739   0.07178079    0.7078282
    ##   Cumul_br_stick
    ## 1      0.1623605
    ## 2      0.2812427
    ## 3      0.3783858
    ## 4      0.4610362
    ## 5      0.5328170
    ## 
    ## $eigenvalues
    ## [1] 27.4 17.9
    ## 
    ## $site_vectors
    ##    field_key       Axis.1      Axis.2      Axis.3       Axis.4 region
    ## 1          1  0.199672424  0.23294856 -0.19435525  0.101605877     BM
    ## 2          2 -0.006972709 -0.27311674 -0.21225603  0.222666251     BM
    ## 3          3 -0.409557123  0.22782817  0.14309643 -0.017079209     FG
    ## 4          4  0.061067854 -0.03577972  0.26966828  0.166220414     FG
    ## 5          5 -0.015576502 -0.07462533  0.15848369  0.362307735     FG
    ## 6          6 -0.193624666  0.02690118 -0.12686722 -0.063421371     FL
    ## 7          7 -0.417659472  0.25561762 -0.18022720  0.232188041     FL
    ## 8          8  0.117713688 -0.02366959  0.13386643  0.051013897     FL
    ## 9          9  0.215731740  0.23400233  0.20182423 -0.054448210     FL
    ## 10        10  0.252574926  0.16138609  0.08836220 -0.166339112     FL
    ## 11        11  0.142663499  0.11252861  0.16308489 -0.113195823     FL
    ## 12        12  0.112122513 -0.09486125 -0.15492559 -0.074325476     FL
    ## 13        13  0.159474041  0.03725657 -0.13819717 -0.093736838     FL
    ## 14        14  0.065325625 -0.03703295 -0.26243051 -0.060578968     FL
    ## 15        15  0.323952309  0.25547118  0.10569566  0.043329383     BM
    ## 16        16 -0.410562025  0.09338195 -0.07190035 -0.113509999     LP
    ## 17        17  0.022585426 -0.19164738 -0.07475391 -0.092801325     LP
    ## 18        18 -0.105238210 -0.23766850  0.04831640 -0.155537215     LP
    ## 19        19  0.034199005 -0.29972627 -0.06664255 -0.050421780     LP
    ## 20        20  0.260512117  0.17085475 -0.20274993  0.034477037     BM
    ## 21        21  0.183288132 -0.06359058 -0.13695878  0.006401245     BM
    ## 22        22 -0.040023707 -0.20371269  0.14146691 -0.038384195     BM
    ## 23        23 -0.131375912 -0.14191328  0.15001266 -0.006479241     BM
    ## 24        24 -0.433585515  0.11699196  0.04599047 -0.137543457     BM
    ## 25        25  0.013292544 -0.24782468  0.17239626  0.017592336     BM
    ##    field_type yr_since
    ## 1    restored       16
    ## 2    restored        3
    ## 3        corn       NA
    ## 4     remnant       NA
    ## 5    restored       15
    ## 6        corn       NA
    ## 7        corn       NA
    ## 8     remnant       NA
    ## 9    restored       40
    ## 10   restored       36
    ## 11   restored       35
    ## 12   restored       10
    ## 13   restored       10
    ## 14   restored       10
    ## 15   restored       28
    ## 16       corn       NA
    ## 17    remnant       NA
    ## 18   restored        4
    ## 19   restored        4
    ## 20    remnant       NA
    ## 21   restored       18
    ## 22   restored        7
    ## 23   restored        2
    ## 24       corn       NA
    ## 25   restored       11
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
    ##          Df SumOfSqs     R2      F Pr(>F)   
    ## Model     2   1.0853 0.2467 3.6024  0.003 **
    ## Residual 22   3.3140 0.7533                 
    ## Total    24   4.3993 1.0000                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## $pairwise_contrasts
    ## 
    ## 
    ## group1     group2        R2   F_value   df1   df2     p_value   p_value_adj
    ## ---------  --------  ------  --------  ----  ----  ----------  ------------
    ## restored   corn       0.252     6.410     1    19   0.0010000        0.0030
    ## restored   remnant    0.022     0.409     1    18   0.9755000        0.9755
    ## corn       remnant    0.381     4.307     1     7   0.0416667        0.0625

Four axes are significant by a broken stick model, between them
explaining 68.7% of the variation in AMF among fields. It may be
worthwhile to examine structure on Axes 3 and 4 sometime. The most
substantial variation here is on the first axis (27.4%) with Axis 2
explaining 17.9% of the variation in AMF abundances. Testing the design
factor *field_type* (with *region* treated as a block using the `strata`
argument of `adonis2`) revealed a significant clustering
$(R^2=0.25,~p=0.003)$.

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
    ## Model     2  0.06937 0.1657 2.1847 0.0025 **
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
    ## restored   corn       0.239     5.961     1    19   0.0010000        0.0030
    ## restored   remnant    0.025     0.459     1    18   0.9725000        0.9725
    ## corn       remnant    0.381     4.307     1     7   0.0416667        0.0625

Three axes are significant by a broken stick model, between them
explaining 48.8% of the variation in AMF among fields. The most
substantial variation here is on the first axis (23%) with Axis 2
explaining 15% of the variation in AMF abundances. Testing the design
factor *field_type* (with *region* treated as a block using the `strata`
argument of `adonis2`) revealed a significant clustering
$(R^2=0.17,~p=0.002)$.

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
    ## [1] "Lingoes correction applied to negative eigenvalues: D' = -0.5*D^2 - 0.113202541256844 , except diagonal elements"
    ## 
    ## $values
    ##   Dim Eigenvalues Corr_eig Rel_corr_eig Broken_stick Cum_corr_eig Cum_br_stick
    ## 1   1   2.7063211 2.819524   0.16208229   0.09442476    0.1620823   0.09442476
    ## 2   2   2.0640256 2.177228   0.12515948   0.07314817    0.2872418   0.16757293
    ## 3   3   0.9798528 1.093055   0.06283505   0.06250987    0.3500768   0.23008280
    ## 4   4   0.8946780 1.007881   0.05793872   0.05541767    0.4080155   0.28550047
    ## 5   5   0.8284124 0.941615   0.05412939   0.05009852    0.4621449   0.33559899
    ## 6   6   0.6123964 0.725599   0.04171156   0.04584320    0.5038565   0.38144219
    ## 
    ## $eigenvalues
    ## [1] 16.2 12.5
    ## 
    ## $site_vectors
    ## # A tibble: 49 × 8
    ##    field_key sample_key  Axis.1  Axis.2  Axis.3  Axis.4   Axis.5 yr_since
    ##        <dbl> <chr>        <dbl>   <dbl>   <dbl>   <dbl>    <dbl>    <dbl>
    ##  1         1 1          -0.125  -0.0631 -0.126  -0.103  -0.0672        16
    ##  2         1 2          -0.160  -0.308  -0.0892 -0.192  -0.0538        16
    ##  3         1 4          -0.471   0.192  -0.0884 -0.133   0.134         16
    ##  4         1 5          -0.0809 -0.332  -0.292  -0.0513  0.210         16
    ##  5         1 7          -0.263  -0.342  -0.0794 -0.159  -0.0236        16
    ##  6         1 8          -0.357   0.0900 -0.0588  0.104  -0.190         16
    ##  7         1 10         -0.390  -0.169  -0.0453 -0.141  -0.00111       16
    ##  8         2 1           0.0595 -0.302   0.330   0.155   0.0356         3
    ##  9         2 3           0.165  -0.138   0.0779  0.127  -0.0976         3
    ## 10         2 5           0.0289 -0.420   0.196  -0.0360  0.0851         3
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
    ## Model     1   1.4272 0.11932 6.3676  5e-04 ***
    ## Residual 47  10.5347 0.88068                  
    ## Total    48  11.9619 1.00000                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## $vector_fit
    ## 
    ## ***VECTORS
    ## 
    ##            Axis.1   Axis.2   Axis.3   Axis.4   Axis.5     r2 Pr(>r)  
    ## yr_since -0.85611  0.39339 -0.17736  0.25939  0.11654 0.7686  0.018 *
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Plots: field_key, plot permutation: free
    ## Permutation: none
    ## Number of permutations: 5039
    ## 
    ## 
    ## 
    ## $vector_fit_scores
    ##              Axis.1    Axis.2     Axis.3    Axis.4    Axis.5
    ## yr_since -0.7505256 0.3448762 -0.1554887 0.2274027 0.1021696

Axis 1 explains 16.2% and axis 2 explains 12.5% of the variation in the
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
$(R^2=0.77,~p=0.02)$.

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
    ## [1] "Lingoes correction applied to negative eigenvalues: D' = -0.5*D^2 - 0.359448791519202 , except diagonal elements"
    ## 
    ## $values
    ##    Dim Eigenvalues Corr_eig Rel_corr_eig Broken_stick Cum_corr_eig Cum_br_stick
    ## 1    1    7.684857 8.044306   0.07189480   0.03314101    0.0718948   0.03314101
    ## 2    2    5.400271 5.759720   0.05147664   0.02736066    0.1233714   0.06050167
    ## 3    3    5.246161 5.605610   0.05009931   0.02447049    0.1734707   0.08497216
    ## 4    4    4.219016 4.578465   0.04091936   0.02254371    0.2143901   0.10751587
    ## 5    5    3.127413 3.486862   0.03116331   0.02109862    0.2455534   0.12861449
    ## 6    6    2.470697 2.830146   0.02529401   0.01994255    0.2708474   0.14855704
    ## 7    7    2.233959 2.593407   0.02317820   0.01897916    0.2940256   0.16753620
    ## 8    8    2.056173 2.415622   0.02158926   0.01815340    0.3156149   0.18568960
    ## 9    9    1.761203 2.120652   0.01895301   0.01743085    0.3345679   0.20312045
    ## 10  10    1.608313 1.967761   0.01758658   0.01678859    0.3521545   0.21990904
    ## 11  11    1.395898 1.755347   0.01568815   0.01621056    0.3678426   0.23611960
    ## 
    ## $eigenvalues
    ## [1] 7.2 5.1
    ## 
    ## $site_vectors
    ## # A tibble: 175 × 16
    ##    field_key sample_key  Axis.1   Axis.2  Axis.3   Axis.4  Axis.5   Axis.6
    ##        <dbl> <chr>        <dbl>    <dbl>   <dbl>    <dbl>   <dbl>    <dbl>
    ##  1         1 1          -0.0987 -0.0415   0.0210 -0.0715   0.195  -0.0407 
    ##  2         1 2          -0.101  -0.336    0.0285  0.0428   0.0839 -0.136  
    ##  3         1 4          -0.363  -0.0130   0.278  -0.149    0.139   0.0579 
    ##  4         1 5           0.0849 -0.338   -0.0842 -0.230    0.0925 -0.0926 
    ##  5         1 7          -0.180  -0.358    0.0462  0.0252   0.115  -0.00544
    ##  6         1 8          -0.350   0.00496  0.117  -0.111    0.0157  0.0251 
    ##  7         1 10         -0.335  -0.226    0.0989  0.00318  0.189   0.00736
    ##  8         2 1          -0.0339 -0.210   -0.278   0.141   -0.237   0.204  
    ##  9         2 3          -0.0152 -0.00195 -0.266   0.190    0.0256  0.0396 
    ## 10         2 5          -0.0474 -0.376   -0.251   0.182   -0.0182  0.106  
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
    ##           Df SumOfSqs     R2      F Pr(>F)    
    ## Model      2    5.393 0.1093 10.553  5e-04 ***
    ## Residual 172   43.953 0.8907                  
    ## Total    174   49.346 1.0000                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## $pairwise_contrasts
    ##     group1  group2    R2 F_value df1 df2   p_value p_value_adj
    ## 1 restored    corn 0.057   8.800   1 145 0.0030000   0.0090000
    ## 2 restored remnant 0.008   1.060   1 138 0.9975000   0.9975000
    ## 3     corn remnant 0.110   7.528   1  61 0.1020408   0.1530612
    ## 
    ## $format
    ## [1] "pandoc"

``` r
write_delim(pcoa_amf_samps$permanova %>% round(., 3), "microbial_communities_files/pcoa_amf_samps_permanova.txt")
write_delim(pcoa_amf_samps$pairwise_contrasts %>% mutate(across(starts_with("p_value"), ~ round(.x, 3))), "microbial_communities_files/pcoa_amf_samps_pairwise.txt")
```

Axis 1 explains 7.2% and axis 2 explains 5.1% of the variation in the
community data. Both axes are important based on the broken stick model,
in fact, the broken stick model shows that 10 axes are important in
explaining variation with this dataset. The relatively low percent
variation explained on axes 1 and 2 is partly due to the high number of
dimensions used when all samples from fields are included. The fidelity
of samples to fields was strong based on a permutation test when
restricting permutations to fields (=plots in `how()`) within regions
(=blocks in `how()`) $(R^2=0.11,~p=5\times 10^{-4})$.

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
    ## [1] "Lingoes correction applied to negative eigenvalues: D' = -0.5*D^2 - 0.156410552447331 , except diagonal elements"
    ## 
    ## $values
    ##   Dim Eigenvalues  Corr_eig Rel_corr_eig Broken_stick Cum_corr_eig Cum_br_stick
    ## 1   1   3.8237409 3.9801515   0.15170440   0.07698793    0.1517044   0.07698793
    ## 2   2   2.3505034 2.5069139   0.09555161   0.06059449    0.2472560   0.13758242
    ## 3   3   1.8578334 2.0142439   0.07677338   0.05239777    0.3240294   0.18998019
    ## 4   4   1.1653433 1.3217539   0.05037896   0.04693329    0.3744083   0.23691348
    ## 5   5   0.9699419 1.1263525   0.04293119   0.04283493    0.4173395   0.27974840
    ## 6   6   0.7828558 0.9392663   0.03580035   0.03955624    0.4531399   0.31930464
    ## 
    ## $eigenvalues
    ## [1] 15.2  9.6
    ## 
    ## $site_vectors
    ## # A tibble: 63 × 11
    ##    field_key sample_key  Axis.1  Axis.2   Axis.3  Axis.4    Axis.5 field_name
    ##        <dbl> <chr>        <dbl>   <dbl>    <dbl>   <dbl>     <dbl> <chr>     
    ##  1         1 1          -0.111  -0.0179 -0.0915  -0.0420  0.0837   BBRP1     
    ##  2         1 2          -0.192  -0.207  -0.235   -0.137   0.0887   BBRP1     
    ##  3         1 4          -0.383   0.279   0.00243 -0.146  -0.0774   BBRP1     
    ##  4         1 5          -0.0880 -0.128  -0.425    0.0990 -0.132    BBRP1     
    ##  5         1 7          -0.285  -0.217  -0.235   -0.131   0.0233   BBRP1     
    ##  6         1 8          -0.337   0.127   0.0975   0.0506  0.0975   BBRP1     
    ##  7         1 10         -0.389  -0.0837 -0.0889  -0.154  -0.000683 BBRP1     
    ##  8         2 1           0.0218 -0.385   0.115   -0.0473 -0.225    ERRP1     
    ##  9         2 3           0.0983 -0.257   0.121    0.101   0.0225   ERRP1     
    ## 10         2 5          -0.0450 -0.463  -0.0524  -0.123  -0.139    ERRP1     
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
    ##          Df SumOfSqs      R2      F Pr(>F)  
    ## Model     2   2.8326 0.17127 6.2001 0.0935 .
    ## Residual 60  13.7061 0.82873                
    ## Total    62  16.5388 1.00000                
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## $pairwise_contrasts
    ##     group1  group2    R2 F_value df1 df2 p_value p_value_adj
    ## 1 restored remnant 0.053   3.011   1  54  0.3675     0.55125
    ## 2 restored    corn 0.085   5.044   1  54  0.1085     0.32550
    ## 3  remnant    corn 0.479  11.011   1  12  1.0000     1.00000
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
    ## [1] "Lingoes correction applied to negative eigenvalues: D' = -0.5*D^2 - 0.0743328751373954 , except diagonal elements"
    ## 
    ## $values
    ##   Dim Eigenvalues  Corr_eig Rel_corr_eig Broken_stick Cum_corr_eig Cum_br_stick
    ## 1   1   1.8298566 1.9041895   0.28142752   0.18672314    0.2814275    0.1867231
    ## 2   2   0.9680232 1.0423561   0.15405383   0.13409156    0.4354814    0.3208147
    ## 3   3   0.6689469 0.7432798   0.10985219   0.10777577    0.5453336    0.4285905
    ## 4   4   0.4340836 0.5084165   0.07514084   0.09023191    0.6204744    0.5188224
    ## 
    ## $eigenvalues
    ## [1] 28.1 15.4
    ## 
    ## $site_vectors
    ## # A tibble: 21 × 9
    ##    field_key sample_key  Axis.1  Axis.2    Axis.3 field_name region field_type
    ##        <dbl> <chr>        <dbl>   <dbl>     <dbl> <chr>      <chr>  <ord>     
    ##  1         3 1          -0.317  -0.0804  0.235    FGC1       FG     corn      
    ##  2         3 2          -0.427   0.0424  0.196    FGC1       FG     corn      
    ##  3         3 4          -0.387  -0.134   0.0197   FGC1       FG     corn      
    ##  4         3 6          -0.362  -0.137  -0.131    FGC1       FG     corn      
    ##  5         3 7          -0.415  -0.0523 -0.000978 FGC1       FG     corn      
    ##  6         3 8          -0.321  -0.141  -0.0691   FGC1       FG     corn      
    ##  7         3 9          -0.388  -0.105  -0.0529   FGC1       FG     corn      
    ##  8         4 1          -0.146   0.153  -0.284    FGREM1     FG     remnant   
    ##  9         4 2           0.0797  0.414   0.0829   FGREM1     FG     remnant   
    ## 10         4 3           0.246   0.356   0.230    FGREM1     FG     remnant   
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
    ## Model     2   2.5246 0.47819 8.2476      1
    ## Residual 18   2.7549 0.52181              
    ## Total    20   5.2795 1.00000              
    ## 
    ## $pairwise_contrasts
    ##    group1   group2    R2 F_value df1 df2 p_value p_value_adj
    ## 1    corn  remnant 0.347   6.372   1  12       1           1
    ## 2    corn restored 0.446   9.645   1  12       1           1
    ## 3 remnant restored 0.325   5.781   1  12       1           1
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
    ## [1] "Lingoes correction applied to negative eigenvalues: D' = -0.5*D^2 - 0.172590559881542 , except diagonal elements"
    ## 
    ## $values
    ##   Dim Eigenvalues  Corr_eig Rel_corr_eig Broken_stick Cum_corr_eig Cum_br_stick
    ## 1   1   3.2451623 3.4177529   0.12198989   0.07698793    0.1219899   0.07698793
    ## 2   2   2.2888502 2.4614408   0.08785623   0.06059449    0.2098461   0.13758242
    ## 3   3   1.7087876 1.8813782   0.06715205   0.05239777    0.2769982   0.18998019
    ## 4   4   1.3684528 1.5410434   0.05500448   0.04693329    0.3320027   0.23691348
    ## 5   5   1.1391586 1.3117492   0.04682028   0.04283493    0.3788229   0.27974840
    ## 6   6   0.9516749 1.1242654   0.04012842   0.03955624    0.4189514   0.31930464
    ## 7   7   0.7767970 0.9493875   0.03388650   0.03682400    0.4528379   0.35612864
    ## 
    ## $eigenvalues
    ## [1] 12.2  8.8
    ## 
    ## $site_vectors
    ## # A tibble: 63 × 12
    ##    field_key sample_key Axis.1  Axis.2   Axis.3  Axis.4   Axis.5  Axis.6
    ##        <dbl> <chr>       <dbl>   <dbl>    <dbl>   <dbl>    <dbl>   <dbl>
    ##  1         6 2          -0.368  0.0266 -0.171    0.191   0.363    0.0205
    ##  2         6 4          -0.311 -0.104   0.153    0.295   0.165    0.0232
    ##  3         6 5          -0.262  0.115   0.342    0.320  -0.0550  -0.0843
    ##  4         6 6          -0.205 -0.126   0.00526  0.197   0.293    0.0804
    ##  5         6 7          -0.277  0.129   0.245    0.355   0.0968  -0.0958
    ##  6         6 8          -0.249 -0.0271  0.0196  -0.121  -0.00904 -0.0457
    ##  7         6 9          -0.420  0.108  -0.107   -0.0196  0.110   -0.0299
    ##  8         7 2          -0.267  0.416   0.249   -0.107  -0.110    0.0839
    ##  9         7 3          -0.315  0.349   0.105   -0.213  -0.111    0.246 
    ## 10         7 4          -0.378  0.194  -0.104   -0.0811 -0.108   -0.170 
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
    ## Model     2   2.8426 0.16416 5.8921  0.044 *
    ## Residual 60  14.4734 0.83584                
    ## Total    62  17.3161 1.00000                
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## $pairwise_contrasts
    ##    group1   group2    R2 F_value df1 df2   p_value p_value_adj
    ## 1    corn  remnant 0.189   4.418   1  19 0.3333333       0.500
    ## 2    corn restored 0.103   6.212   1  54 0.0370000       0.111
    ## 3 remnant restored 0.034   1.648   1  47 0.7180000       0.718
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
    ## [1] "Lingoes correction applied to negative eigenvalues: D' = -0.5*D^2 - 0.0745347378374308 , except diagonal elements"
    ## 
    ## $values
    ##   Dim Eigenvalues  Corr_eig Rel_corr_eig Broken_stick Cum_corr_eig Cum_br_stick
    ## 1   1   1.5763205 1.6508552   0.20834650   0.14824691    0.2083465    0.1482469
    ## 2   2   1.1341999 1.2087346   0.15254859   0.10978537    0.3608951    0.2580323
    ## 3   3   0.6917901 0.7663248   0.09671417   0.09055460    0.4576093    0.3485869
    ## 4   4   0.5127171 0.5872518   0.07411423   0.07773409    0.5317235    0.4263210
    ## 
    ## $eigenvalues
    ## [1] 20.8 15.3
    ## 
    ## $site_vectors
    ## # A tibble: 28 × 9
    ##    field_key sample_key  Axis.1 Axis.2  Axis.3 field_name region field_type
    ##        <dbl> <chr>        <dbl>  <dbl>   <dbl> <chr>      <chr>  <ord>     
    ##  1        16 1          -0.356   0.152  0.0142 LPC1       LP     corn      
    ##  2        16 2          -0.121   0.224  0.100  LPC1       LP     corn      
    ##  3        16 4          -0.372   0.145  0.0278 LPC1       LP     corn      
    ##  4        16 5          -0.263   0.116 -0.0308 LPC1       LP     corn      
    ##  5        16 6          -0.354   0.193  0.0585 LPC1       LP     corn      
    ##  6        16 7          -0.291   0.224  0.0458 LPC1       LP     corn      
    ##  7        16 8          -0.308   0.188  0.0312 LPC1       LP     corn      
    ##  8        17 1          -0.0831 -0.104  0.125  LPREM1     LP     remnant   
    ##  9        17 2           0.233   0.468 -0.178  LPREM1     LP     remnant   
    ## 10        17 4           0.0281 -0.210  0.0128 LPREM1     LP     remnant   
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
    ##          Df SumOfSqs      R2      F Pr(>F)
    ## Model     2   1.5381 0.26021 4.3966    0.5
    ## Residual 25   4.3730 0.73979              
    ## Total    27   5.9112 1.00000              
    ## 
    ## $pairwise_contrasts
    ##    group1   group2    R2 F_value df1 df2   p_value p_value_adj
    ## 1    corn  remnant 0.286   4.814   1  12 1.0000000           1
    ## 2    corn restored 0.256   6.541   1  19 0.3333333           1
    ## 3 remnant restored 0.079   1.634   1  19 1.0000000           1
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
| Blue Mounds   |   15.2 |    9.6 |
| Faville Grove |   28.1 |   15.4 |
| Fermilab      |   12.2 |    8.8 |
| Lake Petite   |   20.8 |   15.3 |

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
