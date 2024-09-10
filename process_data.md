Database assembly: species data
================
Beau Larkin

Last updated: 10 September, 2024

- [Description](#description)
  - [Workflow](#workflow)
  - [ITS data (all fungi)](#its-data-all-fungi)
  - [18S data (mycorrhizae)](#18s-data-mycorrhizae)
  - [Desired outcome](#desired-outcome)
- [Resources](#resources)
  - [Packages and libraries](#packages-and-libraries)
  - [Functions](#functions)
- [Load and process data](#load-and-process-data)
  - [Import files](#import-files)
  - [ETL using `etl()`](#etl-using-etl)
  - [Post-processing 18S data](#post-processing-18s-data)

# Description

Microbial sequence abundances were produced by Lorinda Bullington in
QIIME2. ETL must be performed on output text files to allow downstream
analysis. Iterative steps are needed to find out the optimal number of
samples to keep from each field to ensure equal sampling effort and
adequate representation of diversity.

## Workflow

1.  The script `process_data.R` is run first. A few samples failed to
    amplify, resulting in some fields characterized by 9 samples and
    others by 10. Also, one ITS sample was ambiguously assigned to site
    and was removed. To balance sampling effort across fields, the top 9
    samples by sequence abundance are chosen from each field. **Assign
    “pre” to the argument `process_step` in the etl function so that
    files are created in the /clean_data/pre/… directory.**
2.  Next, `microbial_diagnostics_pre.R` is run to investigate sequencing
    depth in samples and species accumulation in fields. A few samples
    are known to have low sequence abundance (an order of magnitude
    lower than the maximum), and the consequence of rarefying to this
    small depth must be known. A new cutoff for sequence depth, and
    definition of further samples which must be cut, is recommended.
3.  Then, `process_data.R` is run again, this time with the number of
    samples retained per field set to the levels recommended in
    `microbial_diagnostics_pre.R`. As of 2023-10-11, the recommended
    number of samples to keep from all fields is **8 from the ITS
    dataset** and **7 from the 18S dataset.**
    - If downstream analyses have been completed, then it’s likely that
      the `process_data.R` script has been left at this step.
    - **Assign “post” to the argument `process_step` in the etl function
      so that files are created in the /clean_data/… directory.**
4.  Finally, `microbial_diagnostics_post.R` is run. It is very similar
    to the “\_pre” script, but a different file is used so that the two
    may be compared.

## ITS data (all fungi)

Sequence abundances in 97% similar OTUs in individual samples form the
base data. The abundances are raw (not rarefied).

- ITS taxonomy are included in a separate file.

## 18S data (mycorrhizae)

Sequence abundance in 97% similar OTUs in individual samples. The
abundances are raw (not rarefied).

- 18S taxonomy are included in a separate file.
- A unifrac distance matrix will be created and included after sample
  selection and sequence depth rarefaction.

## Desired outcome

For each raw table, species OTU codes must be aligned with short, unique
keys, and then species tables must be transposed into sites-species
matrices. Some fields only retained nine samples. To correct for survey
effort, the nine samples from each field with the greatest total
sequence abundance will be chosen, and sequences summed for each OTU
within fields. Rarefaction of sequencing depth to the minimum total
number of sequences will be applied to summed OTUs.

For each taxonomy table, taxonomy strings must be parsed and unnecessary
characters removed. A function is used to streamline the pipeline and to
reduce errors. [Fungal
traits](https://link.springer.com/article/10.1007/s13225-020-00466-2)
data will be joined with the ITS taxonomy

For all tables, short and unique rownames must be created to allow for
easy joining of species and metadata tables.

For all sequences, zero-abundance and singleton OTUs must be removed
after OTUs have been rarefied and summed within fields.

For the 18S data, a second table is needed to produce a UNIFRAC distance
matrix. The table must have OTUs in rows with OTU ids.

# Resources

## Packages and libraries

``` r
packages_needed = c("tidyverse", "vegan", "knitr")
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

*NOTE:* Due to rounding, a very few OTUs are retained or lost (\<1%)
when function `Rarefy()` is rerun. These different outcomes change
nothing about how results would be interpreted, but they do change axis
limits and other trivial parameters that cause headaches later. Be
advised that downstream changes will be needed if `Rarefy()` is rerun
and new files are created by `write_csv()` in the steps at the end of
this function.

``` r
etl <- function(spe, taxa, samps, traits=NULL, varname, gene, cluster_type, colname_prefix, folder, process_step) {
    
    # Variable definitions
    # spe             = Dataframe or tibble with QIIME2 sequence abundances output, 
    #                   OTUs in rows and samples in columns.
    # taxa            = Dataframe or tibble with QIIME2 taxonomy outputs; OTUs in
    #                   rows and metadata in columns. 
    # samps           = Samples to keep from each field
    # traits          = Additional dataframe of traits or guilds.
    # varname         = An unique key will be created to replace the cumbersome cluster 
    #                   hash which is produced by QIIME2. Varname, a string, begins the
    #                   column name for this new, short key. Unquoted. Example: otu_num
    # gene            = Gene region targeted in sequencing, e.g.: "ITS", "18S", "16S". 
    #                   Quoted. Used to select() column names, so must match text in 
    #                   column names. Also used to create distinct file names.
    # cluster_type    = Clustering algorithm output, e.g.: "otu", "sv". Quoted.
    #                   Used to create simple cluster IDs. 
    # colname_prefix  = Existing, verbose prefix to text of OTU column names. The function removes
    #                   this prefix to make OTU names more concise. 
    # folder          = The function creates output files in the working directory by
    #                   default. To use a subfolder, use this variable. Quoted. 
    #                   Include the "/" before the folder name. If no folder 
    #                   name is desired, use "".
    # process_step    = logical toggle to indicate whether this is the pre or post step in 
    #                   data processing. 
    
    set.seed <- 397
    
    varname <- enquo(varname)
    
    data <- spe %>% left_join(taxa, by = join_by(`#OTU ID`))
    
    # Produce metadata for ITS or 18S data, write to file
    if(gene == "ITS") {
        meta <-
            data %>%
            mutate(!!varname := paste0(cluster_type, "_", row_number())) %>%
            select(!starts_with(gene)) %>%
            rename(otu_ID = `#OTU ID`) %>%
            select(!!varname, everything()) %>%
            separate_wider_delim(
                taxonomy,
                delim = ";",
                names = c(
                    "kingdom",
                    "phylum",
                    "class",
                    "order",
                    "family",
                    "genus",
                    "species"
                ),
                cols_remove = TRUE,
                too_few = "align_start"
            ) %>%
            mutate(kingdom = str_sub(kingdom, 4, nchar(kingdom)),
                   phylum  = str_sub(phylum,  4, nchar(phylum)),
                   class   = str_sub(class,   4, nchar(class)),
                   order   = str_sub(order,   4, nchar(order)),
                   family  = str_sub(family,  4, nchar(family)),
                   genus   = str_sub(genus,   4, nchar(genus)),
                   species = str_sub(species, 4, nchar(species))) %>%
            left_join(traits, by = join_by(phylum, class, order, family, genus)) %>%
            select(-kingdom, -Confidence)
    } else {
        meta <-
            data %>%
            mutate(!!varname := paste0(cluster_type, "_", row_number())) %>%
            select(!starts_with(gene)) %>%
            rename(otu_ID = `#OTU ID`) %>%
            select(!!varname, everything()) %>%
            separate(taxonomy,
                     c("class", "order", "family", "genus", "taxon", "accession"),
                     sep = ";", remove = TRUE, fill = "right") %>%
            select(-Confidence)
    }

    spe_t <-
        data.frame(
            data %>%
                mutate(!!varname := paste0(cluster_type, "_", row_number())) %>%
                select(!!varname, starts_with(gene)),
            row.names = 1
        ) %>%
        t() %>%
        as.data.frame() %>%
        rownames_to_column() %>%
        mutate(rowname = str_remove(rowname, colname_prefix)) %>%
        separate_wider_delim(cols = rowname, delim = "_", names = c("field_key", "sample"))
    
    # Display minimum number of samples in a field
    min_samples <- 
        spe_t %>%
        group_by(field_key) %>% 
        summarize(n = n(), .groups = "drop") %>% 
        pull(n) %>% 
        min()
    
    # Display number of samples in all fields
    samples_fields <-
        spe_t %>%
        group_by(field_key) %>%
        summarize(n = n(), .groups = "drop") %>%
        mutate(field_key = as.numeric(field_key)) %>% 
        left_join(sites, by = join_by(field_key)) %>%
        select(field_key, field_name, region, n) %>% 
        arrange(field_key) %>% 
        kable(format = "pandoc", caption = "Number of samples available in each field")
    
    # Raw (not rarefied) sequence abundances, top n samples, write to file
    spe_topn <- 
        spe_t %>%
        mutate(sum = rowSums(across(starts_with(cluster_type)))) %>%
        group_by(field_key) %>% 
        slice_max(sum, n=samps) %>%
        select(-sum) %>% 
        arrange(field_key, sample)
    # Remove zero abundance columns
    strip_cols1 <- which(apply(spe_topn[, -c(1,2)], 2, sum) == 0)
    spe_samps_raw <- 
        {
            if (length(strip_cols1) == 0) data.frame(spe_topn)
            else data.frame(spe_topn[, -strip_cols1])
        } %>% 
        mutate(field_key = as.numeric(field_key),
               sample = as.numeric(sample)) %>% 
        arrange(field_key, sample) %>% 
        as_tibble()
    
    # Rarefied sequence abundances, top n samples, write to file
    spe_samps_raw_df <- 
        spe_samps_raw %>% 
        mutate(field_sample = paste(field_key, sample, sep = "_")) %>% 
        select(field_sample, everything(), -field_key, -sample) %>% 
        data.frame(., row.names = 1)
    depth_spe_samps_rfy <- min(rowSums(spe_samps_raw_df))
    spe_samps_rrfd <- rrarefy(spe_samps_raw_df, depth_spe_samps_rfy)
    # Remove zero abundance columns
    strip_cols2 <- which(apply(spe_samps_rrfd, 2, sum) == 0)
    spe_samps_rfy <- 
        {
            if (length(strip_cols2) == 0) data.frame(spe_samps_rrfd)
            else data.frame(spe_samps_rrfd[, -strip_cols2])
        } %>% 
        rownames_to_column(var = "field_sample") %>%
        separate_wider_delim(cols = field_sample, delim = "_", names = c("field_key", "sample")) %>% 
        mutate(field_key = as.numeric(field_key),
               sample = as.numeric(sample)) %>% 
        arrange(field_key, sample) %>% 
        as_tibble()
    
    # Produce summaries of raw sequence data for each field, from top n samples, write to file
    spe_raw_sum <- 
        spe_samps_raw %>% 
        group_by(field_key) %>% 
        summarize(across(starts_with(cluster_type), ~ sum(.x)), .groups = "drop") %>%
        mutate(field_key = as.numeric(field_key)) %>% 
        arrange(field_key)
    # Remove zero abundance columns
    strip_cols3 <- which(apply(spe_raw_sum, 2, sum) <= 0)
    spe_raw <- 
        if(length(strip_cols3) == 0) {
            spe_raw_sum
        } else {
            spe_raw_sum[, -strip_cols3]
        }
    
    # Rarefy summed raw sequence data for each field, from top n samples, write to file
    spe_raw_df <- data.frame(spe_raw, row.names = 1)
    depth_spe_rfy <- min(rowSums(spe_raw_df))
    spe_rrfd <- rrarefy(spe_raw_df, depth_spe_rfy)
    # Remove zero abundance columns
    strip_cols4 <- which(apply(spe_rrfd, 2, sum) <= 0)
    spe_rfy <- 
        {
            if (length(strip_cols4) == 0) data.frame(spe_rrfd) 
            else data.frame(spe_rrfd[, -strip_cols4]) 
        } %>% 
        rownames_to_column(var = "field_key") %>%
        mutate(field_key = as.numeric(field_key)) %>% 
        arrange(field_key) %>% 
        as_tibble()
    
    if (process_step == "pre") {
        write_csv(meta, paste0(getwd(), folder, "/pre/spe_", gene, "_metadata.csv"))
        write_csv(spe_samps_raw, paste0(getwd(), folder, "/pre/spe_", gene, "_raw_samples.csv"))
        write_csv(spe_samps_rfy, paste0(getwd(), folder, "/pre/spe_", gene, "_rfy_samples.csv"))
        write_csv(spe_raw, paste0(getwd(), folder, "/pre/spe_", gene, "_raw.csv"))
        write_csv(spe_rfy, paste0(getwd(), folder, "/pre/spe_", gene, "_rfy.csv"))
    } else {
        write_csv(meta, paste0(getwd(), folder, "/spe_", gene, "_metadata.csv"))
        write_csv(spe_samps_raw, paste0(getwd(), folder, "/spe_", gene, "_raw_samples.csv"))
        write_csv(spe_samps_rfy, paste0(getwd(), folder, "/spe_", gene, "_rfy_samples.csv"))
        write_csv(spe_raw, paste0(getwd(), folder, "/spe_", gene, "_raw.csv"))
        write_csv(spe_rfy, paste0(getwd(), folder, "/spe_", gene, "_rfy.csv"))
    }
    
    out <- list(
        min_samples         = min_samples,
        samples_retained    = samps,
        samples_fields      = samples_fields,
        spe_meta            = meta,
        spe_samps_raw       = spe_samps_raw,
        depth_spe_samps_rfy = depth_spe_samps_rfy,
        spe_samps_rfy       = spe_samps_rfy,
        spe_raw             = spe_raw,
        depth_spe_rfy       = depth_spe_rfy,
        spe_rfy             = spe_rfy
    )
    
    return(out)
    
}
```

# Load and process data

## Import files

``` r
its_otu  <- read_delim(paste0(getwd(), "/otu_tables/ITS/ITS_otu_raw.txt"), 
                       show_col_types = FALSE)
its_taxa <- read_delim(paste0(getwd(), "/otu_tables/ITS/ITS_otu_taxonomy.txt"), 
                       show_col_types = FALSE)
# The 18S OTU file contains an unknown site label in the last column; remove it
amf_otu  <- read_delim(paste0(getwd(), "/otu_tables/18S/18S_otu_raw.txt"), 
                       show_col_types = FALSE) %>% select(-last_col())
amf_taxa <- read_delim(paste0(getwd(), "/otu_tables/18S/18S_otu_taxonomy.txt"), 
                       show_col_types = FALSE)
traits   <- read_csv(paste0(getwd(),   "/otu_tables/2023-02-23_fungal_traits.csv"), 
                     show_col_types = FALSE) %>% select(phylum:primary_lifestyle)
# Site metadata
```

``` r
sites    <- read_csv(paste0(getwd(), "/clean_data/sites.csv"), show_col_types = FALSE)
```

## ETL using `etl()`

Schema:
`process_qiime(spe, taxa, samps, traits=NULL, varname, gene, cluster_type, colname_prefix, folder)`
**Note:** If doing the first step of this workflow, retain 9 samples per
field from each dataset and proceed to further diagnostics in
`microbial_diagnostics_pre.R`. **Note:** If doing the third step of this
workflow, Retain 8 samples from ITS and 7 samples from 18S per field
based on results from `microbial_diagnostics_pre.R`. Proceed to
`microbial_diagnostics_post.R` for final exploration of the datasets.

``` r
its <-
    etl(
        spe = its_otu,
        taxa = its_taxa,
        samps = 8,
        traits = traits,
        varname = otu_num,
        gene = "ITS",
        cluster_type = "otu",
        colname_prefix = "ITS_TGP_",
        folder = "/clean_data",
        process_step = "post"
    )
its
```

    ## $min_samples
    ## [1] 9
    ## 
    ## $samples_retained
    ## [1] 8
    ## 
    ## $samples_fields
    ## 
    ## 
    ## Table: Number of samples available in each field
    ## 
    ##  field_key  field_name   region     n
    ## ----------  -----------  -------  ---
    ##          1  BBRP1        BM        10
    ##          2  ERRP1        BM        10
    ##          3  FGC1         FG        10
    ##          4  FGREM1       FG        10
    ##          5  FGRP1        FG        10
    ##          6  FLC1         FL        10
    ##          7  FLC2         FL        10
    ##          8  FLREM1       FL         9
    ##          9  FLRP1        FL         9
    ##         10  FLRP4        FL        10
    ##         11  FLRP5        FL        10
    ##         12  FLRSP1       FL        10
    ##         13  FLRSP2       FL        10
    ##         14  FLRSP3       FL         9
    ##         15  KORP1        BM        10
    ##         16  LPC1         LP        10
    ##         17  LPREM1       LP        10
    ##         18  LPRP1        LP        10
    ##         19  LPRP2        LP        10
    ##         20  MBREM1       BM        10
    ##         21  MBRP1        BM        10
    ##         22  MHRP1        BM        10
    ##         23  MHRP2        BM        10
    ##         24  PHC1         BM        10
    ##         25  PHRP1        BM        10
    ## 
    ## $spe_meta
    ## # A tibble: 3,175 × 9
    ##    otu_num otu_ID      phylum class order family genus species primary_lifestyle
    ##    <chr>   <chr>       <chr>  <chr> <chr> <chr>  <chr> <chr>   <chr>            
    ##  1 otu_1   352d386293… Ascom… Sord… Hypo… Nectr… Fusa… Fusari… plant_pathogen   
    ##  2 otu_2   a78342f18e… Morti… Mort… Mort… Morti… Mort… Mortie… soil_saprotroph  
    ##  3 otu_3   dabfbac17a… Ascom… Sord… Glom… Plect… Gibe… <NA>    plant_pathogen   
    ##  4 otu_4   bdee5c3cbc… Ascom… Euro… Chae… Herpo… unid… uniden… <NA>             
    ##  5 otu_5   73514f6e23… Ascom… Sord… Hypo… Nectr… <NA>  <NA>    <NA>             
    ##  6 otu_6   c92d481c3f… Ascom… Doth… Pleo… Didym… <NA>  <NA>    <NA>             
    ##  7 otu_7   204f1bd97d… Ascom… Doth… Pleo… Peric… Peri… <NA>    plant_pathogen   
    ##  8 otu_8   0ab6be0adc… Ascom… Euro… Chae… Herpo… unid… uniden… <NA>             
    ##  9 otu_9   3c7865fa95… Basid… Trem… Cyst… Mraki… Taus… Tauson… soil_saprotroph  
    ## 10 otu_10  caa87147e4… Ascom… <NA>  <NA>  <NA>   <NA>  <NA>    <NA>             
    ## # ℹ 3,165 more rows
    ## 
    ## $spe_samps_raw
    ## # A tibble: 200 × 3,029
    ##    field_key sample otu_1 otu_2 otu_3 otu_4 otu_5 otu_6 otu_7 otu_8 otu_9 otu_10
    ##        <dbl>  <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>  <dbl>
    ##  1         1      1   170   241     8     0    45   469  1502     0     0    108
    ##  2         1      2    72  1320     0   231    28    42   502    10     0     80
    ##  3         1      4   102   842     3     0    19    80   239     0     0     31
    ##  4         1      5   153     0     4     0    27   203   498   274     0     98
    ##  5         1      6   117     0    13  2286    24   152   862     0     0     33
    ##  6         1      7   416    11    25     6    48   225   112     0     0    154
    ##  7         1      9   115    34    32     0    33   565   217     0     0     82
    ##  8         1     10    85   137    15    54    33   118   294     0     0    118
    ##  9         2      1   721   352   199     0    93     0    84     0     0    680
    ## 10         2      2    90   212    42     0    51     0    43     0     0    233
    ## # ℹ 190 more rows
    ## # ℹ 3,017 more variables: otu_11 <dbl>, otu_12 <dbl>, otu_13 <dbl>,
    ## #   otu_14 <dbl>, otu_15 <dbl>, otu_16 <dbl>, otu_17 <dbl>, otu_18 <dbl>,
    ## #   otu_19 <dbl>, otu_20 <dbl>, otu_21 <dbl>, otu_22 <dbl>, otu_23 <dbl>,
    ## #   otu_24 <dbl>, otu_25 <dbl>, otu_26 <dbl>, otu_27 <dbl>, otu_28 <dbl>,
    ## #   otu_29 <dbl>, otu_30 <dbl>, otu_31 <dbl>, otu_32 <dbl>, otu_33 <dbl>,
    ## #   otu_34 <dbl>, otu_35 <dbl>, otu_36 <dbl>, otu_37 <dbl>, otu_38 <dbl>, …
    ## 
    ## $depth_spe_samps_rfy
    ## [1] 5321
    ## 
    ## $spe_samps_rfy
    ## # A tibble: 200 × 2,855
    ##    field_key sample otu_1 otu_2 otu_3 otu_4 otu_5 otu_6 otu_7 otu_8 otu_9 otu_10
    ##        <dbl>  <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>  <dbl>
    ##  1         1      1   104   126     4     0    29   268   851     0     0     58
    ##  2         1      2    59  1043     0   186    26    27   400     9     0     68
    ##  3         1      4    87   627     2     0    12    57   181     0     0     23
    ##  4         1      5   115     0     4     0    17   142   321   182     0     64
    ##  5         1      6    88     0    11  1689    17   111   625     0     0     28
    ##  6         1      7   330     7    19     6    38   171    86     0     0    122
    ##  7         1      9   103    28    26     0    32   485   194     0     0     72
    ##  8         1     10    63   104    10    46    23    86   230     0     0     94
    ##  9         2      1   500   237   132     0    65     0    56     0     0    447
    ## 10         2      2    51   144    31     0    30     0    23     0     0    126
    ## # ℹ 190 more rows
    ## # ℹ 2,843 more variables: otu_11 <dbl>, otu_12 <dbl>, otu_13 <dbl>,
    ## #   otu_14 <dbl>, otu_15 <dbl>, otu_16 <dbl>, otu_17 <dbl>, otu_18 <dbl>,
    ## #   otu_19 <dbl>, otu_20 <dbl>, otu_21 <dbl>, otu_22 <dbl>, otu_23 <dbl>,
    ## #   otu_24 <dbl>, otu_25 <dbl>, otu_26 <dbl>, otu_27 <dbl>, otu_28 <dbl>,
    ## #   otu_29 <dbl>, otu_30 <dbl>, otu_31 <dbl>, otu_32 <dbl>, otu_33 <dbl>,
    ## #   otu_34 <dbl>, otu_35 <dbl>, otu_36 <dbl>, otu_37 <dbl>, otu_38 <dbl>, …
    ## 
    ## $spe_raw
    ## # A tibble: 25 × 2,890
    ##    field_key otu_1 otu_2 otu_3 otu_4 otu_5 otu_6 otu_7 otu_8 otu_9 otu_10 otu_11
    ##        <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>  <dbl>  <dbl>
    ##  1         1  1230  2585   100  2577   257  1854  4226   284     0    704      0
    ##  2         2  4082  2526  1231    92   855   131   632     0    11   2344    140
    ##  3         3   492     7   489     0   808   387   379     0    30      0    692
    ##  4         4  2351     0     5     0  1779  2874  2036     0    85    917      0
    ##  5         5   663   191   648     0  1798  3766  6229     0   101   1306    196
    ##  6         6  1717   490   282     0   774  1033  1548     0  4099      6   3816
    ##  7         7  1831  1922  1503     9   902  3530   126     0 11026      0    410
    ##  8         8  1721   659   462  1395  1173  1773  2188     0     0   1513     60
    ##  9         9  2835  1581  1414  1143   306  1704  1397    27     0    754     79
    ## 10        10  1001  1033  1007   540   266   938  1573  2037     0    318     34
    ## # ℹ 15 more rows
    ## # ℹ 2,878 more variables: otu_12 <dbl>, otu_13 <dbl>, otu_14 <dbl>,
    ## #   otu_15 <dbl>, otu_16 <dbl>, otu_17 <dbl>, otu_18 <dbl>, otu_19 <dbl>,
    ## #   otu_20 <dbl>, otu_21 <dbl>, otu_22 <dbl>, otu_23 <dbl>, otu_24 <dbl>,
    ## #   otu_25 <dbl>, otu_26 <dbl>, otu_27 <dbl>, otu_28 <dbl>, otu_29 <dbl>,
    ## #   otu_30 <dbl>, otu_31 <dbl>, otu_32 <dbl>, otu_33 <dbl>, otu_34 <dbl>,
    ## #   otu_35 <dbl>, otu_36 <dbl>, otu_37 <dbl>, otu_38 <dbl>, otu_39 <dbl>, …
    ## 
    ## $depth_spe_rfy
    ## [1] 58170
    ## 
    ## $spe_rfy
    ## # A tibble: 25 × 2,886
    ##    field_key otu_1 otu_2 otu_3 otu_4 otu_5 otu_6 otu_7 otu_8 otu_9 otu_10 otu_11
    ##        <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>  <dbl>  <dbl>
    ##  1         1  1230  2585   100  2577   257  1854  4226   284     0    704      0
    ##  2         2  3614  2198  1072    80   766   116   565     0     9   2071    125
    ##  3         3   444     6   433     0   706   338   330     0    27      0    614
    ##  4         4  2241     0     5     0  1686  2736  1918     0    79    866      0
    ##  5         5   535   149   526     0  1479  3089  5070     0    86   1074    162
    ##  6         6  1369   383   214     0   614   803  1205     0  3276      6   2958
    ##  7         7  1521  1622  1270     7   754  2984   104     0  9232      0    346
    ##  8         8  1566   595   431  1260  1057  1615  1998     0     0   1397     57
    ##  9         9  1968  1088  1016   814   210  1213   988    17     0    538     55
    ## 10        10   844   870   844   435   228   778  1335  1700     0    271     24
    ## # ℹ 15 more rows
    ## # ℹ 2,874 more variables: otu_12 <dbl>, otu_13 <dbl>, otu_14 <dbl>,
    ## #   otu_15 <dbl>, otu_16 <dbl>, otu_17 <dbl>, otu_18 <dbl>, otu_19 <dbl>,
    ## #   otu_20 <dbl>, otu_21 <dbl>, otu_22 <dbl>, otu_23 <dbl>, otu_24 <dbl>,
    ## #   otu_25 <dbl>, otu_26 <dbl>, otu_27 <dbl>, otu_28 <dbl>, otu_29 <dbl>,
    ## #   otu_30 <dbl>, otu_31 <dbl>, otu_32 <dbl>, otu_33 <dbl>, otu_34 <dbl>,
    ## #   otu_35 <dbl>, otu_36 <dbl>, otu_37 <dbl>, otu_38 <dbl>, otu_39 <dbl>, …

``` r
amf <-
    etl(
        spe = amf_otu,
        taxa = amf_taxa,
        samps = 7,
        varname = otu_num,
        gene = "18S",
        cluster_type = "otu",
        colname_prefix = "X18S_TGP_",
        folder = "/clean_data",
        process_step = "post"
    )
amf
```

    ## $min_samples
    ## [1] 9
    ## 
    ## $samples_retained
    ## [1] 7
    ## 
    ## $samples_fields
    ## 
    ## 
    ## Table: Number of samples available in each field
    ## 
    ##  field_key  field_name   region     n
    ## ----------  -----------  -------  ---
    ##          1  BBRP1        BM         9
    ##          2  ERRP1        BM        10
    ##          3  FGC1         FG        10
    ##          4  FGREM1       FG        10
    ##          5  FGRP1        FG        10
    ##          6  FLC1         FL        10
    ##          7  FLC2         FL        10
    ##          8  FLREM1       FL        10
    ##          9  FLRP1        FL         9
    ##         10  FLRP4        FL        10
    ##         11  FLRP5        FL        10
    ##         12  FLRSP1       FL        10
    ##         13  FLRSP2       FL        10
    ##         14  FLRSP3       FL         9
    ##         15  KORP1        BM        10
    ##         16  LPC1         LP        10
    ##         17  LPREM1       LP        10
    ##         18  LPRP1        LP        10
    ##         19  LPRP2        LP        10
    ##         20  MBREM1       BM        10
    ##         21  MBRP1        BM         9
    ##         22  MHRP1        BM        10
    ##         23  MHRP2        BM        10
    ##         24  PHC1         BM        10
    ##         25  PHRP1        BM        10
    ## 
    ## $spe_meta
    ## # A tibble: 152 × 8
    ##    otu_num otu_ID                       class order family genus taxon accession
    ##    <chr>   <chr>                        <chr> <chr> <chr>  <chr> <chr> <chr>    
    ##  1 otu_1   320f3edc7b48ba5691766ccc71b… Glom… Glom… Glome… Glom… Glom… VTX00212 
    ##  2 otu_2   97cefc055a2fa1b8caf5b81635b… Glom… Glom… Glome… Glom… Glom… VTX00135 
    ##  3 otu_3   2c47b9acb976f06611bddfe10f9… Para… Para… Parag… Para… <NA>  <NA>     
    ##  4 otu_4   7487be824e0acea7ee8c22b94f1… Glom… Glom… Glome… Glom… Glom… VTX00195 
    ##  5 otu_5   e8cbf85e2f78be8ff21c0c3e0f5… Glom… Glom… Glome… Glom… <NA>  <NA>     
    ##  6 otu_6   87ff74cc0f55f12dc4056a4a8a2… Glom… Glom… Claro… Clar… <NA>  <NA>     
    ##  7 otu_7   f10a50b1220da2d3e275f9772b3… Glom… Glom… Glome… Glom… Glom… VTX00222 
    ##  8 otu_8   0d508e08bf3048aa561e7c9d96e… Glom… Glom… Glome… Glom… Glom… VTX00315 
    ##  9 otu_9   4a4251fdd8c94240584c5d4ff6a… Glom… Glom… Glome… Glom… Glom… VTX00214 
    ## 10 otu_10  210c71717f87c8e4912431ed98b… Glom… Glom… Claro… Clar… Clar… VTX00056 
    ## # ℹ 142 more rows
    ## 
    ## $spe_samps_raw
    ## # A tibble: 175 × 151
    ##    field_key sample otu_1 otu_2 otu_3 otu_4 otu_5 otu_6 otu_7 otu_8 otu_9 otu_10
    ##        <dbl>  <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>  <dbl>
    ##  1         1      1   120   195    51     0   307   120    88   143     0     18
    ##  2         1      2   102   195   378   267  1006   275   120   171    61      0
    ##  3         1      4    17   279    30     0   395    16   184   894     0      0
    ##  4         1      5     0     0   486    30   566   223   550    33     0     16
    ##  5         1      7    55   111   275     0   516   110    21   148     0      0
    ##  6         1      8   130   505   129     0   113    26    55   758    26    110
    ##  7         1     10   144   359   193     0   648    30    16   151     0     28
    ##  8         2      1  1251     0  2407     0   275   278    30     0   293     37
    ##  9         2      3   454    25    84     0   265   541    95     0   288    352
    ## 10         2      5   393     0  2247     0  1381    26     0     0   294    109
    ## # ℹ 165 more rows
    ## # ℹ 139 more variables: otu_11 <dbl>, otu_12 <dbl>, otu_13 <dbl>, otu_14 <dbl>,
    ## #   otu_15 <dbl>, otu_16 <dbl>, otu_17 <dbl>, otu_18 <dbl>, otu_19 <dbl>,
    ## #   otu_20 <dbl>, otu_21 <dbl>, otu_22 <dbl>, otu_23 <dbl>, otu_24 <dbl>,
    ## #   otu_25 <dbl>, otu_26 <dbl>, otu_27 <dbl>, otu_28 <dbl>, otu_29 <dbl>,
    ## #   otu_30 <dbl>, otu_31 <dbl>, otu_32 <dbl>, otu_33 <dbl>, otu_34 <dbl>,
    ## #   otu_35 <dbl>, otu_36 <dbl>, otu_37 <dbl>, otu_38 <dbl>, otu_39 <dbl>, …
    ## 
    ## $depth_spe_samps_rfy
    ## [1] 1364
    ## 
    ## $spe_samps_rfy
    ## # A tibble: 175 × 148
    ##    field_key sample otu_1 otu_2 otu_3 otu_4 otu_5 otu_6 otu_7 otu_8 otu_9 otu_10
    ##        <dbl>  <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>  <dbl>
    ##  1         1      1    94   145    34     0   241    93    55   109     0     15
    ##  2         1      2    43    50   138    90   313    82    30    53    19      0
    ##  3         1      4     5    92    10     0   142     2    61   305     0      0
    ##  4         1      5     0     0   189     8   221    83   233    13     0      5
    ##  5         1      7    32    63   163     0   322    58    11    86     0      0
    ##  6         1      8    65   225    62     0    49     9    22   329    15     47
    ##  7         1     10    68   189    86     0   297    16     7    76     0     11
    ##  8         2      1   310     0   584     0    64    77     5     0    58     11
    ##  9         2      3   135     6    27     0    82   166    33     0   103    113
    ## 10         2      5    99     0   524     0   288     7     0     0    68     21
    ## # ℹ 165 more rows
    ## # ℹ 136 more variables: otu_11 <dbl>, otu_12 <dbl>, otu_13 <dbl>, otu_14 <dbl>,
    ## #   otu_15 <dbl>, otu_16 <dbl>, otu_17 <dbl>, otu_18 <dbl>, otu_19 <dbl>,
    ## #   otu_20 <dbl>, otu_21 <dbl>, otu_22 <dbl>, otu_23 <dbl>, otu_24 <dbl>,
    ## #   otu_25 <dbl>, otu_26 <dbl>, otu_27 <dbl>, otu_28 <dbl>, otu_29 <dbl>,
    ## #   otu_30 <dbl>, otu_31 <dbl>, otu_32 <dbl>, otu_33 <dbl>, otu_34 <dbl>,
    ## #   otu_35 <dbl>, otu_36 <dbl>, otu_37 <dbl>, otu_38 <dbl>, otu_39 <dbl>, …
    ## 
    ## $spe_raw
    ## # A tibble: 25 × 147
    ##    field_key otu_1 otu_2 otu_3 otu_4 otu_5 otu_6 otu_7 otu_8 otu_9 otu_10 otu_11
    ##        <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>  <dbl>  <dbl>
    ##  1         1   568  1644  1542   297  3551   800  1034  2298    87    172   1485
    ##  2         2  4373    25  9072     0  4175  1415   462     0  1400   1265      0
    ##  3         3    44   350    31  3454  1386   595  1974     0   570     25      0
    ##  4         4  3402  3525     0  1223  1124  1298   466  1810  5069   2188    580
    ##  5         5   755  2901    27    58  4195  1675   779    18  2256   3249    123
    ##  6         6   725   125  6648  4679  1610  1895  2840     5   725    149      6
    ##  7         7    16     0  1677  1294  5984  1262    67   112   166      0      0
    ##  8         8  2221  2617   454   624  1688  1706  2341  1161  1126   2343   3033
    ##  9         9   725  2699   211  3035   487   489  1616  4418  1811   1397   4372
    ## 10        10  1850  1234   165  2431   192   708  1627  4483   228    737   1832
    ## # ℹ 15 more rows
    ## # ℹ 135 more variables: otu_12 <dbl>, otu_13 <dbl>, otu_14 <dbl>, otu_15 <dbl>,
    ## #   otu_16 <dbl>, otu_17 <dbl>, otu_18 <dbl>, otu_19 <dbl>, otu_20 <dbl>,
    ## #   otu_21 <dbl>, otu_22 <dbl>, otu_23 <dbl>, otu_24 <dbl>, otu_25 <dbl>,
    ## #   otu_26 <dbl>, otu_27 <dbl>, otu_28 <dbl>, otu_29 <dbl>, otu_30 <dbl>,
    ## #   otu_31 <dbl>, otu_32 <dbl>, otu_33 <dbl>, otu_34 <dbl>, otu_35 <dbl>,
    ## #   otu_36 <dbl>, otu_37 <dbl>, otu_38 <dbl>, otu_39 <dbl>, otu_40 <dbl>, …
    ## 
    ## $depth_spe_rfy
    ## [1] 17975
    ## 
    ## $spe_rfy
    ## # A tibble: 25 × 146
    ##    field_key otu_1 otu_2 otu_3 otu_4 otu_5 otu_6 otu_7 otu_8 otu_9 otu_10 otu_11
    ##        <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>  <dbl>  <dbl>
    ##  1         1   472  1359  1274   245  2958   673   856  1922    73    145   1263
    ##  2         2  2044    12  4261     0  2010   710   229     0   662    617      0
    ##  3         3    22   210    24  2108   889   370  1259     0   366     16      0
    ##  4         4  1869  1913     0   660   600   696   260   972  2770   1185    309
    ##  5         5   371  1346     8    23  2009   782   390     7  1079   1527     64
    ##  6         6   373    69  3432  2431   811   983  1491     4   378     83      4
    ##  7         7    13     0  1097   867  4066   841    45    68   120      0      0
    ##  8         8  1091  1277   200   309   874   871  1208   577   551   1124   1515
    ##  9         9   420  1544   122  1822   293   278   939  2580  1056    808   2507
    ## 10        10  1338   886   121  1724   141   509  1171  3238   160    534   1327
    ## # ℹ 15 more rows
    ## # ℹ 134 more variables: otu_12 <dbl>, otu_13 <dbl>, otu_14 <dbl>, otu_15 <dbl>,
    ## #   otu_16 <dbl>, otu_17 <dbl>, otu_18 <dbl>, otu_19 <dbl>, otu_20 <dbl>,
    ## #   otu_21 <dbl>, otu_22 <dbl>, otu_23 <dbl>, otu_24 <dbl>, otu_25 <dbl>,
    ## #   otu_26 <dbl>, otu_27 <dbl>, otu_28 <dbl>, otu_29 <dbl>, otu_30 <dbl>,
    ## #   otu_31 <dbl>, otu_32 <dbl>, otu_33 <dbl>, otu_34 <dbl>, otu_35 <dbl>,
    ## #   otu_36 <dbl>, otu_37 <dbl>, otu_38 <dbl>, otu_39 <dbl>, otu_40 <dbl>, …

## Post-processing 18S data

To produce a UNIFRAC distance table, the trimmed table `amf$spe_rfy`
must be imported back into QIIME, where sequence data can be used with
abundances to create a phylogeny and distance matrix. The data frame
must be transposed and use legal column names (i.e., non-numeric).

Site metadata is used to create better column names.

The resultant distance matrix will be imported again when needed for
multivariate analysis and ordination.

``` r
amf_export <-
    data.frame(
        sites %>%
            select(field_key, field_name) %>%
            left_join(amf$spe_rfy, by = join_by(field_key)) %>%
            select(-field_key),
        row.names = 1
    ) %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column(var = "otu_num") %>%
    left_join(amf$spe_meta %>% select(otu_num, otu_ID), by = join_by(otu_num)) %>%
    select(otu_ID, everything(), -otu_num)
write_tsv(amf_export, paste0(getwd(), "/otu_tables/18S/spe_18S_rfy_export.tsv"))
```
