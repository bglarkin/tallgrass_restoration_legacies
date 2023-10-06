Database assembly: species data
================
Beau Larkin

Last updated: 06 October, 2023

- [Description](#description)
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
analysis.

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
packages_needed = c("tidyverse", "vegan")
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
etl <- function(spe, taxa, samps, traits=NULL, varname, gene, cluster_type, colname_prefix, folder) {
    
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
        if(length(strip_cols1) == 0) {
            data.frame(spe_topn)
        } else {
            data.frame(spe_topn[, -strip_cols1])
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
        if(length(strip_cols2) == 0) {
            data.frame(spe_samps_rrfd)
        } else {
            spe_samps_rfy <- data.frame(spe_samps_rrfd[, -strip_cols2])
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
        if(length(strip_cols4) == 0) {
            data.frame(spe_rrfd)
        } else {
            data.frame(spe_rrfd[, -strip_cols4])
        } %>% 
            rownames_to_column(var = "field_key") %>%
            mutate(field_key = as.numeric(field_key)) %>% 
            arrange(field_key) %>% 
            as_tibble()
    
    write_csv(meta, paste0(getwd(), folder, "/spe_", gene, "_metadata.csv"))
    write_csv(spe_samps_raw, paste0(getwd(), folder, "/spe_", gene, "_raw_samples.csv"))
    write_csv(spe_samps_rfy, paste0(getwd(), folder, "/spe_", gene, "_rfy_samples.csv"))
    write_csv(spe_raw, paste0(getwd(), folder, "/spe_", gene, "_raw.csv"))
    write_csv(spe_rfy, paste0(getwd(), folder, "/spe_", gene, "_rfy.csv"))
    
    out <- list(
        min_samples         = min_samples,
        spe_meta            = as_tibble(meta),
        spe_samps_raw       = as_tibble(spe_samps_raw),
        depth_spe_samps_rfy = depth_spe_samps_rfy,
        spe_samps_rfy       = as_tibble(spe_samps_rfy),
        spe_raw             = as_tibble(spe_raw),
        depth_spe_rfy       = depth_spe_rfy,
        spe_rfy             = as_tibble(spe_rfy)
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
```

## ETL using `etl()`

Schema:
`process_qiime(spe, taxa, samps, traits=NULL, varname, gene, cluster_type, colname_prefix, folder)`

``` r
its <-
    etl(
        spe = its_otu,
        taxa = its_taxa,
        samps = 9,
        traits = traits,
        varname = otu_num,
        gene = "ITS",
        cluster_type = "otu",
        colname_prefix = "ITS_TGP_",
        folder = "/clean_data"
    )
its
```

    ## $min_samples
    ## [1] 9
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
    ## # A tibble: 225 × 3,127
    ##    field_key sample otu_1 otu_2 otu_3 otu_4 otu_5 otu_6 otu_7 otu_8 otu_9 otu_10
    ##        <dbl>  <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>  <dbl>
    ##  1         1      1   170   241     8     0    45   469  1502     0     0    108
    ##  2         1      2    72  1320     0   231    28    42   502    10     0     80
    ##  3         1      3   144     0     0   426    16   159   166   798     0     34
    ##  4         1      4   102   842     3     0    19    80   239     0     0     31
    ##  5         1      5   153     0     4     0    27   203   498   274     0     98
    ##  6         1      6   117     0    13  2286    24   152   862     0     0     33
    ##  7         1      7   416    11    25     6    48   225   112     0     0    154
    ##  8         1      9   115    34    32     0    33   565   217     0     0     82
    ##  9         1     10    85   137    15    54    33   118   294     0     0    118
    ## 10         2      1   721   352   199     0    93     0    84     0     0    680
    ## # ℹ 215 more rows
    ## # ℹ 3,115 more variables: otu_11 <dbl>, otu_12 <dbl>, otu_13 <dbl>,
    ## #   otu_14 <dbl>, otu_15 <dbl>, otu_16 <dbl>, otu_17 <dbl>, otu_18 <dbl>,
    ## #   otu_19 <dbl>, otu_20 <dbl>, otu_21 <dbl>, otu_22 <dbl>, otu_23 <dbl>,
    ## #   otu_24 <dbl>, otu_25 <dbl>, otu_26 <dbl>, otu_27 <dbl>, otu_28 <dbl>,
    ## #   otu_29 <dbl>, otu_30 <dbl>, otu_31 <dbl>, otu_32 <dbl>, otu_33 <dbl>,
    ## #   otu_34 <dbl>, otu_35 <dbl>, otu_36 <dbl>, otu_37 <dbl>, otu_38 <dbl>, …
    ## 
    ## $depth_spe_samps_rfy
    ## [1] 1629
    ## 
    ## $spe_samps_rfy
    ## # A tibble: 225 × 2,765
    ##    field_key sample otu_1 otu_2 otu_3 otu_4 otu_5 otu_6 otu_7 otu_8 otu_9 otu_10
    ##        <dbl>  <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>  <dbl>
    ##  1         1      1    30    37     3     0     9    83   251     0     0     23
    ##  2         1      2    19   307     0    61     6     9   124     3     0     23
    ##  3         1      3    42     0     0   125     7    50    44   210     0      9
    ##  4         1      4    26   191     0     0     3    18    59     0     0      5
    ##  5         1      5    37     0     3     0     8    40   104    50     0     18
    ##  6         1      6    29     0     3   490     7    37   207     0     0      8
    ##  7         1      7   106     5     5     1    14    50    32     0     0     34
    ##  8         1      9    30     6     8     0     9   158    70     0     0     25
    ##  9         1     10    24    29     2    10     8    28    79     0     0     31
    ## 10         2      1   160    61    40     0     9     0    21     0     0    140
    ## # ℹ 215 more rows
    ## # ℹ 2,753 more variables: otu_11 <dbl>, otu_12 <dbl>, otu_13 <dbl>,
    ## #   otu_14 <dbl>, otu_15 <dbl>, otu_16 <dbl>, otu_17 <dbl>, otu_18 <dbl>,
    ## #   otu_19 <dbl>, otu_20 <dbl>, otu_21 <dbl>, otu_22 <dbl>, otu_23 <dbl>,
    ## #   otu_24 <dbl>, otu_25 <dbl>, otu_26 <dbl>, otu_27 <dbl>, otu_28 <dbl>,
    ## #   otu_29 <dbl>, otu_30 <dbl>, otu_31 <dbl>, otu_32 <dbl>, otu_33 <dbl>,
    ## #   otu_34 <dbl>, otu_35 <dbl>, otu_36 <dbl>, otu_37 <dbl>, otu_38 <dbl>, …
    ## 
    ## $spe_raw
    ## # A tibble: 25 × 3,078
    ##    field_key otu_1 otu_2 otu_3 otu_4 otu_5 otu_6 otu_7 otu_8 otu_9 otu_10 otu_11
    ##        <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>  <dbl>  <dbl>
    ##  1         1  1374  2585   100  3003   273  2013  4392  1082     0    738      0
    ##  2         2  4709  2575  1247    92   986   131   762     0    11   2520    235
    ##  3         3   593     7   567     0   854   435   413     0    30      0    716
    ##  4         4  2456    31     5     0  2233  3000  2161     0    85    917      0
    ##  5         5   745   191   665     0  2012  3948  6380     0   101   1381    231
    ##  6         6  1746   490   282     0   795  1044  1553     0  4117      6   3906
    ##  7         7  1914  1999  1639     9  1008  3690   134     0 12627      0    434
    ##  8         8  2198   785  1067  1395  1244  1815  2211     0    38   2189    137
    ##  9         9  2954  1581  1446  1143   306  1765  1397    27     0    777     79
    ## 10        10  1187  1170  1108   586   318  1087  1762  2037     0    403     34
    ## # ℹ 15 more rows
    ## # ℹ 3,066 more variables: otu_12 <dbl>, otu_13 <dbl>, otu_14 <dbl>,
    ## #   otu_15 <dbl>, otu_16 <dbl>, otu_17 <dbl>, otu_18 <dbl>, otu_19 <dbl>,
    ## #   otu_20 <dbl>, otu_21 <dbl>, otu_22 <dbl>, otu_23 <dbl>, otu_24 <dbl>,
    ## #   otu_25 <dbl>, otu_26 <dbl>, otu_27 <dbl>, otu_28 <dbl>, otu_29 <dbl>,
    ## #   otu_30 <dbl>, otu_31 <dbl>, otu_32 <dbl>, otu_33 <dbl>, otu_34 <dbl>,
    ## #   otu_35 <dbl>, otu_36 <dbl>, otu_37 <dbl>, otu_38 <dbl>, otu_39 <dbl>, …
    ## 
    ## $depth_spe_rfy
    ## [1] 64141
    ## 
    ## $spe_rfy
    ## # A tibble: 25 × 3,074
    ##    field_key otu_1 otu_2 otu_3 otu_4 otu_5 otu_6 otu_7 otu_8 otu_9 otu_10 otu_11
    ##        <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>  <dbl>  <dbl>
    ##  1         1  1374  2585   100  3003   273  2013  4392  1082     0    738      0
    ##  2         2  4178  2293  1091    81   888   117   681     0    10   2228    198
    ##  3         3   525     7   498     0   747   393   369     0    29      0    639
    ##  4         4  2308    28     5     0  2104  2818  2045     0    83    871      0
    ##  5         5   614   155   535     0  1676  3303  5305     0    89   1157    194
    ##  6         6  1418   407   215     0   635   843  1256     0  3369      6   3188
    ##  7         7  1615  1688  1391     8   858  3112   113     0 10647      0    372
    ##  8         8  1885   688   913  1201  1103  1578  1932     0    33   1917    120
    ##  9         9  2240  1199  1099   853   235  1329  1060    22     0    606     66
    ## 10        10  1019  1006   932   507   273   920  1495  1746     0    356     26
    ## # ℹ 15 more rows
    ## # ℹ 3,062 more variables: otu_12 <dbl>, otu_13 <dbl>, otu_14 <dbl>,
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
        samps = 9,
        varname = otu_num,
        gene = "18S",
        cluster_type = "otu",
        colname_prefix = "X18S_TGP_",
        folder = "/clean_data"
    )
```

    ## Warning in rrarefy(spe_samps_raw_df, depth_spe_samps_rfy): function should be
    ## used for observed counts, but smallest count is 2

    ## Warning in rrarefy(spe_raw_df, depth_spe_rfy): function should be used for
    ## observed counts, but smallest count is 2

``` r
amf
```

    ## $min_samples
    ## [1] 9
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
    ## # A tibble: 225 × 154
    ##    field_key sample otu_1 otu_2 otu_3 otu_4 otu_5 otu_6 otu_7 otu_8 otu_9 otu_10
    ##    <chr>     <chr>  <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>  <dbl>
    ##  1 1         1        120   195    51     0   307   120    88   143     0     18
    ##  2 1         10       144   359   193     0   648    30    16   151     0     28
    ##  3 1         2        102   195   378   267  1006   275   120   171    61      0
    ##  4 1         4         17   279    30     0   395    16   184   894     0      0
    ##  5 1         5          0     0   486    30   566   223   550    33     0     16
    ##  6 1         6          0     0    44     0    39     0     0     0     0      0
    ##  7 1         7         55   111   275     0   516   110    21   148     0      0
    ##  8 1         8        130   505   129     0   113    26    55   758    26    110
    ##  9 1         9         23    35   118     0    29     0    60    84     0      0
    ## 10 10        1        399   256    38     0    10   190    59   599     0    321
    ## # ℹ 215 more rows
    ## # ℹ 142 more variables: otu_11 <dbl>, otu_12 <dbl>, otu_13 <dbl>, otu_14 <dbl>,
    ## #   otu_15 <dbl>, otu_16 <dbl>, otu_17 <dbl>, otu_18 <dbl>, otu_19 <dbl>,
    ## #   otu_20 <dbl>, otu_21 <dbl>, otu_22 <dbl>, otu_23 <dbl>, otu_24 <dbl>,
    ## #   otu_25 <dbl>, otu_26 <dbl>, otu_27 <dbl>, otu_28 <dbl>, otu_29 <dbl>,
    ## #   otu_30 <dbl>, otu_31 <dbl>, otu_32 <dbl>, otu_33 <dbl>, otu_34 <dbl>,
    ## #   otu_35 <dbl>, otu_36 <dbl>, otu_37 <dbl>, otu_38 <dbl>, otu_39 <dbl>, …
    ## 
    ## $depth_spe_samps_rfy
    ## [1] 163
    ## 
    ## $spe_samps_rfy
    ## # A tibble: 225 × 139
    ##    field_key sample otu_1 otu_2 otu_3 otu_4 otu_5 otu_6 otu_7 otu_8 otu_9 otu_10
    ##        <dbl>  <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>  <dbl>
    ##  1         1      1    13    19     5     0    23     8    10    13     0      0
    ##  2         1      2     4    10    17     2    37     8     6     8     6      0
    ##  3         1      4     0    12     3     0    17     2     8    36     0      0
    ##  4         1      5     0     0    19     1    23    10    29     4     0      0
    ##  5         1      6     0     0    44     0    39     0     0     0     0      0
    ##  6         1      7     5    10    27     0    21     9     1    10     0      0
    ##  7         1      8     4    30     7     0    13     0     2    43     0      9
    ##  8         1      9     4    10    32     0    10     0    17    24     0      0
    ##  9         1     10     6    20    17     0    34     1     3     9     0      0
    ## 10         2      1    35     0    81     0     5    11     0     0     6      1
    ## # ℹ 215 more rows
    ## # ℹ 127 more variables: otu_11 <dbl>, otu_12 <dbl>, otu_13 <dbl>, otu_14 <dbl>,
    ## #   otu_15 <dbl>, otu_16 <dbl>, otu_17 <dbl>, otu_18 <dbl>, otu_19 <dbl>,
    ## #   otu_20 <dbl>, otu_21 <dbl>, otu_22 <dbl>, otu_23 <dbl>, otu_24 <dbl>,
    ## #   otu_25 <dbl>, otu_26 <dbl>, otu_27 <dbl>, otu_28 <dbl>, otu_29 <dbl>,
    ## #   otu_30 <dbl>, otu_31 <dbl>, otu_32 <dbl>, otu_33 <dbl>, otu_34 <dbl>,
    ## #   otu_35 <dbl>, otu_36 <dbl>, otu_37 <dbl>, otu_38 <dbl>, otu_39 <dbl>, …
    ## 
    ## $spe_raw
    ## # A tibble: 25 × 153
    ##    field_key otu_1 otu_2 otu_3 otu_4 otu_5 otu_6 otu_7 otu_8 otu_9 otu_10 otu_11
    ##        <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>  <dbl>  <dbl>
    ##  1         1   591  1679  1704   297  3619   800  1094  2382    87    172   1485
    ##  2         2  5102    72 11494     0  5482  1552   690     0  1575   1469      0
    ##  3         3    44   370    31  4403  1677   759  2001     0   651     25      0
    ##  4         4  3960  4657     0  1223  1407  1376   707  2436  5684   2194    614
    ##  5         5  1058  3181    37    58  4505  1760  1078    18  2655   3749    127
    ##  6         6   742   134  6900  6012  1723  2200  3213     5   725    149      6
    ##  7         7    16     0  1860  1817  6286  1326   214   112   166      0      0
    ##  8         8  2504  3105   498  1181  2049  1896  2504  1475  1176   2419   3883
    ##  9         9  1246  2938   225  3035   521   635  1701  4710  1811   1587   5032
    ## 10        10  2255  1547   169  2511   282   737  1840  5165   265    783   2017
    ## # ℹ 15 more rows
    ## # ℹ 141 more variables: otu_12 <dbl>, otu_13 <dbl>, otu_14 <dbl>, otu_15 <dbl>,
    ## #   otu_16 <dbl>, otu_17 <dbl>, otu_18 <dbl>, otu_19 <dbl>, otu_20 <dbl>,
    ## #   otu_21 <dbl>, otu_22 <dbl>, otu_23 <dbl>, otu_24 <dbl>, otu_25 <dbl>,
    ## #   otu_26 <dbl>, otu_27 <dbl>, otu_28 <dbl>, otu_29 <dbl>, otu_30 <dbl>,
    ## #   otu_31 <dbl>, otu_32 <dbl>, otu_33 <dbl>, otu_34 <dbl>, otu_35 <dbl>,
    ## #   otu_36 <dbl>, otu_37 <dbl>, otu_38 <dbl>, otu_39 <dbl>, otu_40 <dbl>, …
    ## 
    ## $depth_spe_rfy
    ## [1] 19545
    ## 
    ## $spe_rfy
    ## # A tibble: 25 × 150
    ##    field_key otu_1 otu_2 otu_3 otu_4 otu_5 otu_6 otu_7 otu_8 otu_9 otu_10 otu_11
    ##        <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>  <dbl>  <dbl>
    ##  1         1   505  1472  1480   260  3177   691   956  2077    78    150   1280
    ##  2         2  2173    33  4964     0  2344   675   312     0   692    645      0
    ##  3         3    29   202    17  2577   947   426  1196     0   405     17      0
    ##  4         4  1968  2330     0   610   701   681   355  1223  2862   1100    306
    ##  5         5   434  1382    21    25  1998   792   473     4  1193   1655     55
    ##  6         6   374    71  3470  2939   838  1083  1605     2   332     76      2
    ##  7         7    12     0  1229  1153  4040   839   150    73   104      0      0
    ##  8         8  1154  1428   225   569   953   857  1166   701   557   1169   1785
    ##  9         9   701  1632   130  1737   306   343   966  2682  1081    886   2893
    ## 10        10  1508  1058   123  1756   209   520  1244  3570   178    536   1398
    ## # ℹ 15 more rows
    ## # ℹ 138 more variables: otu_12 <dbl>, otu_13 <dbl>, otu_14 <dbl>, otu_15 <dbl>,
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
sites <- read_csv(paste0(getwd(), "/clean_data/sites.csv"), show_col_types = FALSE)
```

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
