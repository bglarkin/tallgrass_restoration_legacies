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

*NOTE:* all `write_csv()` steps have been commented out as of 2023-03-14
to prevent overwriting existing files. This is because function
`Rarefy()` produces inconsistent results. Due to rounding, a very few
OTUs are retained or lost (\<1%) when this function is rerun. These
different outcomes change nothing about how results would be
interpreted, but they do change axis limits and other trivial parameters
that cause headaches later.

``` r
etl <- function(spe, taxa, samps, traits=NULL, varname, gene, cluster_type, colname_prefix, folder) {
    
    # Variable definitions
    # spe             = Dataframe or tibble with QIIME2 sequence abundances output, 
    #                   OTUs in rows and samples in columns.
    # taxa            = Dataframe or tibble with QIIME2 taxonomy outputs; OTUs in
    #                   rows and metadata in columns.   
    # samps           = Samples to keep from each field, default=6
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
    ## # A tibble: 225 × 2,764
    ##    field_key sample otu_1 otu_2 otu_3 otu_4 otu_5 otu_6 otu_7 otu_8 otu_9 otu_10
    ##        <dbl>  <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>  <dbl>
    ##  1         1      1    33    48     1     0     8    74   249     0     0     13
    ##  2         1      2    21   309     0    55     8    16   119     2     0     21
    ##  3         1      3    35     0     0   137     6    39    46   188     0     12
    ##  4         1      4    20   202     2     0     3    23    62     0     0      5
    ##  5         1      5    38     0     0     0     6    35   108    51     0     21
    ##  6         1      6    32     0     1   527     6    29   183     0     0      5
    ##  7         1      7   103     3     8     0    15    56    31     0     0     36
    ##  8         1      9    28     9     7     0     8   140    63     0     0     20
    ##  9         1     10    19    22     1    15     6    24    68     0     0     21
    ## 10         2      1   139    62    43     0    20     0    22     0     0    143
    ## # ℹ 215 more rows
    ## # ℹ 2,752 more variables: otu_11 <dbl>, otu_12 <dbl>, otu_13 <dbl>,
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
    ##  2         2  4178  2272  1124    74   880   111   678     0    11   2223    202
    ##  3         3   522     7   502     0   739   389   372     0    28      0    641
    ##  4         4  2305    31     5     0  2105  2819  2041     0    79    868      0
    ##  5         5   626   162   561     0  1659  3331  5332     0    86   1144    195
    ##  6         6  1380   385   240     0   652   845  1273     0  3381      6   3190
    ##  7         7  1608  1660  1408     9   828  3088   114     0 10698      0    364
    ##  8         8  1903   664   926  1223  1070  1569  1913     0    33   1896    115
    ##  9         9  2295  1185  1098   872   234  1385  1079    18     0    600     60
    ## 10        10  1014  1012   938   493   285   920  1489  1744     0    339     30
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
    ##     field_key sample otu_1 otu_2 otu_3 otu_4 otu_5 otu_6 otu_7 otu_8 otu_9
    ## 1           1      1   120   195    51     0   307   120    88   143     0
    ## 2           1     10   144   359   193     0   648    30    16   151     0
    ## 3           1      2   102   195   378   267  1006   275   120   171    61
    ## 4           1      4    17   279    30     0   395    16   184   894     0
    ## 5           1      5     0     0   486    30   566   223   550    33     0
    ## 6           1      6     0     0    44     0    39     0     0     0     0
    ## 7           1      7    55   111   275     0   516   110    21   148     0
    ## 8           1      8   130   505   129     0   113    26    55   758    26
    ## 9           1      9    23    35   118     0    29     0    60    84     0
    ## 10         10      1   399   256    38     0    10   190    59   599     0
    ## 11         10     10    21    42     0   835    96    59     0   139     0
    ## 12         10      2    80   203     4     9    71    29     0   598    37
    ## 13         10      3   235   195     0   101    14   125   118   948     0
    ## 14         10      4     0   231     0   190    17    23   278  1425   193
    ## 15         10      5   745    21    94  1238    39   182   130     0    35
    ## 16         10      6   325   110     0    71    19     0   213    84     0
    ## 17         10      7    92   165    33    67    16    55  1037   625     0
    ## 18         10      9   358   324     0     0     0    74     5   747     0
    ## 19         11      1   278    96     0   169   395   131   849     0   279
    ## 20         11     10    84   308    34     0     0    78   302   329    66
    ## 21         11      3    28   805    11     0    37   103   304    61    14
    ## 22         11      4    41   278    16  1112    18   169    80   573   100
    ## 23         11      5   356   373     0    24    20   128     0   896   195
    ## 24         11      6    11    91   368     0     0   129    17   394   660
    ## 25         11      7    51   682   114   731    76   235   666   281   256
    ## 26         11      8   442   290    58   857    17    21   361   243    84
    ## 27         11      9    20   286    66     0    49    81   344   448    60
    ## 28         12      1   125   453   625     0    14   229   383    78   221
    ## 29         12     10   185    63    73     0    66     3   239   226     0
    ## 30         12      3    26   159   112     0   240   164   723     0     8
    ## 31         12      4   145    12    72     0     9    33   130    12     0
    ## 32         12      5   116   351   125    25   611   423   544    28    71
    ## 33         12      6  1228   385   362     0   287     0     0   737    86
    ## 34         12      7   320     0  1192     0   260   105    36   111   120
    ## 35         12      8   425   192   128     0   274    54   856    97     5
    ## 36         12      9   214  1105   127     0    11   199   294     0     0
    ## 37         13      1   597   420   107     0   804   776   287   348    18
    ## 38         13     10   231   218   253     0   114   114    16   211     0
    ## 39         13      2     0   233    17     0   121   212    85     0   107
    ## 40         13      3     4   205   200   222   192   243   288   213    73
    ## 41         13      4    33   364   226   747   508   228   307   244     0
    ## 42         13      5   103  1293    18    31   106   484   424   318     0
    ## 43         13      7    44     0  1116     0    80    42   230   630     0
    ## 44         13      8   840    38   396     0    56    12   555   339     0
    ## 45         13      9   453   223    54     0   164   137   369   227     0
    ## 46         14      1   111   305   293    46   864   158   583     0     0
    ## 47         14      2   334   415   115     0   186   165   460    35     0
    ## 48         14      3    85   514   134     0   363   149   270     0   452
    ## 49         14      4    80   127   394     0   230   134   223     0     0
    ## 50         14      5   135   321   586     0   585   514   509     0   107
    ## 51         14      6   150     0    27    21     0   148   328     9    78
    ## 52         14      7   203     0   512     0   167   103   444  1006     0
    ## 53         14      8   780     0   617    71   505   209   338     0    22
    ## 54         14      9   153     0   332     0   643     7   115     0     0
    ## 55         15      1   103   285     0     0     0     0   118   585    25
    ## 56         15     10     0   251     0     0    33    35   174   450    82
    ## 57         15      2   320   273     0    75    63     0    61   120    73
    ## 58         15      3   246   433   101   628   233    66     0   747     5
    ## 59         15      4     0   189     0     0     0    46    20  1035     0
    ## 60         15      5   559   403   772     0    25     0    35   486     0
    ## 61         15      6   501  1423    98     0    42    87     9  3297     0
    ## 62         15      7   248   271     0     0   111    15   184  2058   171
    ## 63         15      9     9   480     0    12   320     0   165    31     0
    ## 64         16      1     0     0    31    16   244   224   741     0     0
    ## 65         16     10    12     0    14   398   382    27   230     0     0
    ## 66         16      2     0     0   370   129   183     0   339     0     0
    ## 67         16      3     0     0   180     0   264     0   110     0     0
    ## 68         16      4     0     0    28    29   227   126   563     0    66
    ## 69         16      5     8     0   423  1097   119   273  1073     0   222
    ## 70         16      6     0     0   318     0   739    78  2138    33   133
    ## 71         16      7     0     0   292   133   763    16   993     0    20
    ## 72         16      8     0     0    92    66   238   194   459     0    16
    ## 73         17      1   434   201    71   199   135     7  1199     0    18
    ## 74         17     10     5    78    68     0  1469    55    63     0    20
    ## 75         17      2   106     0   135     0   502   379    62     0     0
    ## 76         17      3    94   312   113     0   145    76   341     0     0
    ## 77         17      4   467   817     0   166   221   187   520     0     0
    ## 78         17      5  1517     0     0     0   172   189   448     0     0
    ## 79         17      6   340   246     0   232   173    63  1222     0   599
    ## 80         17      7   359   239    18    27   104    20   148   136     0
    ## 81         17      9    21     0   273   253   111    60    39     0     0
    ## 82         18      1   586     0    27   199   203   248   508     0   270
    ## 83         18     10    56    31     0     0    46   428   449     0   245
    ## 84         18      2     0     9     0  1298   144   336   135   372   686
    ## 85         18      3   321     0    82     0    47   165   110   207   409
    ## 86         18      4   592     0     8    22   173    61   597     0   797
    ## 87         18      5   318     0    17     0   341   120  1031     0   497
    ## 88         18      6     0     0     0    75   157   923   983    19   791
    ## 89         18      7  1396     0     0     0   146   378   930     0   196
    ## 90         18      8     0   396    78     0   237   232   832     0   271
    ## 91         19      1   252     9   221    28   148    25   724     0    57
    ## 92         19     10   674     0   356     0    33     0   482     0   142
    ## 93         19      2   618   799   113    30     0   328   260    13   370
    ## 94         19      3  1577    60   580   119   127   129   191     0     0
    ## 95         19      4   714     0   789     0    46    56    18     0     0
    ## 96         19      5     0     0  1088     0   123   177   499     0   835
    ## 97         19      6     0     0   695     0   131   440   102     0    43
    ## 98         19      7     0     0   909     0     0   213    16     0     0
    ## 99         19      8   322     0   604     0   133   327   317     0   202
    ## 100         2      1  1251     0  2407     0   275   278    30     0   293
    ## 101         2      2   505    31  1597     0  1078    52   151     0    20
    ## 102         2      3   454    25    84     0   265   541    95     0   288
    ## 103         2      4   224    16   825     0   229    85    77     0   155
    ## 104         2      5   393     0  2247     0  1381    26     0     0   294
    ## 105         2      6     0     0  1454     0   695   272    20     0   223
    ## 106         2      7  1621     0  1102     0   708    50     0     0    76
    ## 107         2      8   352     0  1252     0   125   109   317     0   226
    ## 108         2      9   302     0   526     0   726   139     0     0     0
    ## 109        20      1     0   153   187     0   636    15   152   747     0
    ## 110        20      2   119   514     0     0   240   112    33   456     0
    ## 111        20      3     0    55     0     0   110    27    71    91     0
    ## 112        20      4   202   189     0     0   166   337   297   397     0
    ## 113        20      5   259    54    38     0   434    30   205   394     0
    ## 114        20      6     0    21   100     0    23    32   157   138     0
    ## 115        20      7     0   151     0     0    13     6   202   379     0
    ## 116        20      8     0   247    54     0   467     0   111   331     0
    ## 117        20      9     0     7    15     0     5     4    63    20     0
    ## 118        21      1    42   274     0     0    69    64   225   205    23
    ## 119        21     10   315   351     0     0    94     0   176   408   104
    ## 120        21      2   383   135   180     0    85    57   140   325     0
    ## 121        21      3    14     0    13     0   134     8   174     0     0
    ## 122        21      4   230    62   183     0   357    46    61   230   699
    ## 123        21      5   168     0   188     0   537    71   160     0     0
    ## 124        21      7   196   203   122     0   142    86   271   521     0
    ## 125        21      8    55     0   106     0    21     0     0     0     0
    ## 126        21      9     0    32    87     0   121   280     0   208     0
    ## 127        22      1   186    22     0   103     0     0   528     0   282
    ## 128        22     10   467   591   187   332   477    87   634   812   907
    ## 129        22      2     0   123     0    34   691     0   951     0   479
    ## 130        22      3   112   195     0     0   230     0   464     0   137
    ## 131        22      5  1377   374     0     0  1667   178    93     0   490
    ## 132        22      6   805   153     0    51    16    83   521     0   828
    ## 133        22      7   425   174     0    70    91    50   403     0   125
    ## 134        22      8   362    50     0   108   792     5   365     0   243
    ## 135        22      9   157   121     0   345   111    21   148     0   243
    ## 136        23      1    46   220     0    65   319    25   409     0   677
    ## 137        23     10   125    46     0     0   415    69   271     0   313
    ## 138        23      2   303   227     0   301   244   145    59     0    74
    ## 139        23      3   646   104     0     0   107     0   341     0    20
    ## 140        23      4   696    69     0   721   351   181     0     0   951
    ## 141        23      5    62    11     0  1337   665   522   191     0   633
    ## 142        23      6   891     0    38   385   394   618   850     0   377
    ## 143        23      7   230    99    66   143   970   327   206     0   175
    ## 144        23      9   147    92     0   181    59   188   546     0   439
    ## 145        24      1     5     0     0     0  1008    59   749     0     0
    ## 146        24     10     0    86     0    77   251     0  1031     0     0
    ## 147        24      2     0     0     0    97  1637    84   227     0   249
    ## 148        24      3     0     0     0    46    81    15   476     0    37
    ## 149        24      4     0     0     0   337     0    97   696     0    38
    ## 150        24      5     0     0     0  1232    27    96   323     0     6
    ## 151        24      6     0     0     0    78   161   339  1533    13    44
    ## 152        24      8     0     0     0     0    51    16  1025     0   139
    ## 153        24      9     0     0     0   249   323    67   578     0     4
    ## 154        25      1   201   302     9    18    77   153   545     0   263
    ## 155        25     10   164     0    35     0   682   112   235     0   326
    ## 156        25      2   918   188    92    25   196    45   233     0   250
    ## 157        25      3   676   530   105   148    16   422    44     0   287
    ## 158        25      4   333     0    73     0   158    45   221     0  1356
    ## 159        25      5   517   644     0   100    47     0   683     0   432
    ## 160        25      6   760   356   179   531    39   194   307     0   571
    ## 161        25      8   524   626     0     0   183    76     0     0    14
    ## 162        25      9   277   301    17     0   132    10    73     0   581
    ## 163         3      1    25     0     0  1206  1059    71    84     0     0
    ## 164         3     10     0    20     0    42   130   158     0     0    81
    ## 165         3      2    15     0     0  2064    52   159    15     0     0
    ## 166         3      3     0     0     0   907   161     6    27     0     0
    ## 167         3      4     4     0     0   171   174     4   203     0    31
    ## 168         3      6     0   306    10     0     0   272   599     0     0
    ## 169         3      7     0    31     0     0     9    28   540     0   290
    ## 170         3      8     0     8    21     7    47    21    48     0   249
    ## 171         3      9     0     5     0     6    45    40   485     0     0
    ## 172         4      1   344   443     0    55    19    49     0   293   257
    ## 173         4      2   433    77     0     0    68   268     0   238  1352
    ## 174         4      3   307    35     0    84   122   141     0    97  1244
    ## 175         4      4   909    45     0     0    15    76    50   194   612
    ## 176         4      5   191     0     0  1065   210   392   262     0  1106
    ## 177         4      6   668  1216     0    19   536    88   154   672   379
    ## 178         4      7   550  1709     0     0   154   284     0   316   119
    ## 179         4      8   122   365     0     0   175     0   149     8   531
    ## 180         4      9   436   767     0     0   108    78    92   618    84
    ## 181         5     10   158     0     0     0   104    22   107     0   321
    ## 182         5      2    99   187     0     0   642   162   100     0   491
    ## 183         5      3     0   488     0     0    33    91   124     0   464
    ## 184         5      4    29    40     0     0   341   228   115     0    60
    ## 185         5      5   110   923    17     0   670   283    65    18   331
    ## 186         5      6    19   498    10     0   684    95     0     0   419
    ## 187         5      7    64   386     0     0  1594   569   108     0   129
    ## 188         5      8   434   379     0    58   231   247   267     0   362
    ## 189         5      9   145   280    10     0   206    63   192     0    78
    ## 190         6     10    17     9    99  1333    12    18   373     0     0
    ## 191         6      2     0     0   657     0    80    63     0     0   418
    ## 192         6      3     0     0   153     0   101   287     0     0     0
    ## 193         6      4     0     0  1812   243    18    19   708     5     0
    ## 194         6      5     0     0   392  1835    30     0   798     0     8
    ## 195         6      6   399     0  1857    56    19   286   498     0    73
    ## 196         6      7     0     0  1014  2316    84   260   516     0     0
    ## 197         6      8   326   125   539   229  1236   337   262     0   226
    ## 198         6      9     0     0   377     0   143   930    58     0     0
    ## 199         7     10     0     0    18   282   244     0   147     0     0
    ## 200         7      2     0     0     0   574    71   105    15     0    35
    ## 201         7      3    16     0    55   183   125    77     0    12     0
    ## 202         7      4     0     0   263   177  1143    39     0   100     0
    ## 203         7      5     0     0   556     0     0     0     0     0     0
    ## 204         7      6     0     0     0   360   813    48     0     0     0
    ## 205         7      7     0     0   673     0  3331    14    52     0    32
    ## 206         7      8     0     0   130     0   501   979     0     0    99
    ## 207         7      9     0     0   165   241    58    64     0     0     0
    ## 208         8     10    80   138    22   458    87    63     0   131     0
    ## 209         8      2     0   455     4    47    17    24   478   235   106
    ## 210         8      3   139   616    29     0    31   203   113   242   185
    ## 211         8      4   170   554    66     0    11    30    53     0   133
    ## 212         8      5  1045    16     0     0   127   658   914     0   231
    ## 213         8      6   491   128   271     0   425   574     0   684   102
    ## 214         8      7   234   509    84   467   872   173   500     0   310
    ## 215         8      8   203   350    22    99   274   127   163   183    50
    ## 216         8      9   142   339     0   110   205    44   283     0    59
    ## 217         9      1    64   440    76     7    64   156     0   719     0
    ## 218         9     10   234   131     0  1706    21   111   260   489    22
    ## 219         9      2   106   338     0     0     0    97   362    88   262
    ## 220         9      3   225   415     0    23    15     0   313   629    18
    ## 221         9      5     0   368    73   102   326    57   261  1818     0
    ## 222         9      6     0   893    62   652    61    31   397     0  1509
    ## 223         9      7    96   114     0   545     0    37    23   675     0
    ## 224         9      8   282    79     5     0    34   146    43   266     0
    ## 225         9      9   239   160     9     0     0     0    42    26     0
    ##     otu_10 otu_11 otu_12 otu_13 otu_14 otu_15 otu_16 otu_17 otu_18 otu_19
    ## 1       18      0      0     69      0      0     13    111      0      0
    ## 2       28    361      0      4      0     74      0     12      0      0
    ## 3        0      0      0     20      0    250      0     47      0      0
    ## 4        0   1021      0      0      0      0      0     47      0    221
    ## 5       16      0      0      0      0     14      6    104     12      0
    ## 6        0      0      0      0      0      0      0      0      0      0
    ## 7        0    103      0     55      0    164      0     12      0      0
    ## 8      110      0      8      0      0     39      0     11      0      0
    ## 9        0      0      0      0      0     28     17      0      0      0
    ## 10     321      0    470     86     16      0      0     63      0      0
    ## 11       7    281    103      0      0      0      0      0      0      0
    ## 12      14      0     38      0      8      0      0      0      0      0
    ## 13      48    887     28      0      0      8      0     34      0      0
    ## 14      48      0    214      0      0      0      0    113      8      0
    ## 15     190    153     64      0      0      0      0     84      0      0
    ## 16      32    185      0      0      0     59      0     63      0      0
    ## 17      57      0    108     64      0     27      0    258     16      0
    ## 18      66    511     16      8      0     76      0      0      0     13
    ## 19       0      0    204      0     49     12      7    348      0      0
    ## 20     253      4     47     63      6      0      0    113      0      0
    ## 21     201     54    527     37     64      0      0    127     37      0
    ## 22      88      0      0     99      0     35      0     18      5      0
    ## 23     100      0      0     63      0      0      0      0      0    106
    ## 24      73   1304     29     41     13     55     27      6      0      0
    ## 25      13   1309     32    213      0      0      5    235      0      0
    ## 26     295      0      6      0    149      0      0    274     37      0
    ## 27      57      0      0    480     38     14      0    186      0      0
    ## 28     213      0     43      0      0    133      0    121      0      0
    ## 29       7      0      0      0      0    275     17     67      0      0
    ## 30      18      0      0      0      0      0     29    203      0      0
    ## 31       3      0     38      0      0      4      0     49     28      0
    ## 32       0      0     19      0      0    116     17    236      0      0
    ## 33      23    674      0      0      0    201      0      0      0      0
    ## 34      97      0      8      0      0     58    385     52     85      0
    ## 35     192      0      0      0      0      0      0    268     67      0
    ## 36     363      0      0      0      0     11      8    227      0      0
    ## 37       0      0    100      0      0      0      0    196    592      0
    ## 38     177    103    139      0      0      0      0     24      0      0
    ## 39     114   1343    388     30      0      0      0     51      0      0
    ## 40     584    144      0      0      0      0      0     83      0      0
    ## 41       0     56      0      0      0      8      0    127      0      0
    ## 42       0    175     24     11      0      0      0    337      0      0
    ## 43      40      0     44      0      0      0      8     15      0      0
    ## 44       7      0      0      0      0      0      0    120      0      0
    ## 45      13    242    164      0      0      0      0    109      0      0
    ## 46       0      0      4      0      0      0     59    198      0      0
    ## 47       4      0      0      7      0      0     61    143      0      0
    ## 48       0      0     12      0      0      0     42    130      0      0
    ## 49     107      0     75      0      0     80    104     90      0      0
    ## 50       4      0     65      0      0     19     24    226      0      0
    ## 51      37      0      0      0      0      0     15    115      0      0
    ## 52       0      0      0      0     39      0     59    160      0      0
    ## 53      30      0     12    379      0    945    134     87      0      0
    ## 54       0      0      0      0      0      0     39     18      0      0
    ## 55       0   2205     64     21      0      0      0     35      0    115
    ## 56      62    901    170     30      0      0      0     79      5    161
    ## 57       0   1075    107     14      0      0      0      0     77      0
    ## 58     123      0     79     17      0     26      0      0      0    124
    ## 59      28    754    139      0      0     51      0      0      0      0
    ## 60     102   1596      0      0      0     26    231     88      0      0
    ## 61       0    410     12     27      0      0      0      0    129      0
    ## 62      18   2570    111      0      0      0      0     35    180      0
    ## 63      57    764    202      0     29      0      0     54      0      0
    ## 64       0      0     16    562      0      0      0    261      0      0
    ## 65       0      0     26    666      0      0     36    102     11      0
    ## 66       0      0      0    203      0      0     21    170      0      0
    ## 67       0      0      0   1235      0      0      0     36      0      0
    ## 68      11      0      9    620      0      0      0    253      0      0
    ## 69       0      0      6    779      0      0      0    386      0      0
    ## 70       0      0     77   1455      0      0     37    467      0      0
    ## 71       0      0     72    444      0      0    109    665      0      0
    ## 72       8      0     14    410      0     14      0    245      0      0
    ## 73      50      0    130      0      0     56      4    420      0      0
    ## 74     404      0     12      0      0    199      0      0      0      0
    ## 75      38      0      7      0      0      0      0     10      0      0
    ## 76      90      0    198      0      0      0      0    113      0      0
    ## 77       0      0    757    109      0      0      0    240    238      0
    ## 78     130      0    693     16      0      0      0    183    482      0
    ## 79     323      0    223     59      0      0    487    584    259      0
    ## 80      99      0    338      0      0     21     15     97      0      0
    ## 81     228      0    290      0      0    357     33     19      0      0
    ## 82       0      0    652     35      0    855    238    312      0      0
    ## 83     784      0    377      0      0    132     86    163    263      0
    ## 84       9      0    817      0      0      0     82    107      0      0
    ## 85      13      0    442    329      0    751    115     81      0      0
    ## 86      70      0    374     76      0     92    180    346      0      0
    ## 87      38      0    586    208      0     34    442    477      0      0
    ## 88       0      0    401     65      0     22    260    438      0      0
    ## 89       0      0     37      0      0      0     11    301    570      0
    ## 90       0      0    108      0      0    447    219    419      0      0
    ## 91       0      0    395     38      0    538      0    340     99      0
    ## 92       0      0     39      0      0    298      0    220    252      0
    ## 93     464      0    314    115      0    114    368    169    105      0
    ## 94     207      0    133     59      0    209     14     79      0      0
    ## 95     570      0    112      0      0    306      0     12      0      0
    ## 96       0      0   1316      0      0    107      0    180      0      0
    ## 97       7      0      0      0      0    439      0     97      0      0
    ## 98       0      0      0      0      0     52     36     26      0      0
    ## 99     476      0    373     26      0    586     37    260      0      0
    ## 100     37      0      0     24      0    125    191     17    314      0
    ## 101    119      0      0     16      0      0    123     55     16      0
    ## 102    352      0     69     49      0    438    479     73     19      0
    ## 103     85      0      0      0      0    837    159     68    119      0
    ## 104    109      0      0      0      0    722     39     16    293      0
    ## 105      8      0    224      7      0    298     66     16    352      0
    ## 106    365      0     44      0      0   1651    335      0     47      0
    ## 107    316      0    128     12      0    154    150    277    158      0
    ## 108     78      0     39     27      0      0     34      0      0      0
    ## 109     43      0      0      0      0    563      0     26      0      0
    ## 110      0    179      0      0      0     87      0     26      0      0
    ## 111      0    200      0      0      0    103      0     23      0      0
    ## 112      0     16      0      0      0     25      0     45      0      0
    ## 113      0      0      0      0      0     67      0     37      0      0
    ## 114      0      0      0      0      0      0      0     18      0      0
    ## 115     51     80      0      0      0     85      0     50      0      0
    ## 116      0      0      0      0      0      0      0     23      0      0
    ## 117     19      0      0      0      0      9      0     32      0      0
    ## 118     19      0      0      0      0     13      0     62      0      0
    ## 119    151      7     53     17      0      0     11     73      0      0
    ## 120     56      0     52      0      0     53      0     53      0      0
    ## 121    471      0     37      0      0     73      0     49      0      0
    ## 122     77      0      0      0      0      0      0      0      0      0
    ## 123     15    246     20      0      0    313     33     56      0      0
    ## 124      0      0     39     33      0      0    274    151      0      0
    ## 125      0      0      0      0      0      0     18      0      0      0
    ## 126     50      0     16     21      0    109      0      0      0      0
    ## 127     14      0    582     62     71      0      0    288     82      0
    ## 128    291     21    663      0    175      0      0    517     35     74
    ## 129     66      0    826      0      0      0    313    390    205    180
    ## 130     43      0    307     30     91      0     40    179      0      0
    ## 131    129      0    209    309    422      0      0     41     70      6
    ## 132    131      0   1259    101      0      0      0    221    250      0
    ## 133    395      0    372      0      0      0     23    201     16      0
    ## 134    301      0     91     38      0      0      0    186     15      0
    ## 135    678      0     57    279      0      0      0     77     22      0
    ## 136     33      0    159    256      0      0      0    216     51      0
    ## 137    715      0      8     83      0      0      0    134     57     38
    ## 138    117      0    354     77      0      0      0     17     86     11
    ## 139     61      0    664    172      0      0      0    128      0    309
    ## 140     53      0    267    182      0      0      0      0      0      0
    ## 141      0      0    474    103      0      0      0    139    370      0
    ## 142    767      0    249     81      0      0     88    444      8     73
    ## 143    127      0     25    214      0      0      0    117     55     65
    ## 144     95      0    152    106      0      0      0    267      0     58
    ## 145      0      0     34    137      0      0      0    249      0      0
    ## 146     17      0     21    257      0      0      0    632     25      0
    ## 147     40      0      0    345      0      0      0     86      0      0
    ## 148      8      0      0    867      0      0      0    198      0    216
    ## 149      6      0      0    985      0      0      0    199      0     59
    ## 150      0      0      0   1400      0      0      0    119      0      0
    ## 151     25      0     53    581      0      0      0    658      9      0
    ## 152     16      0     72    502      0     22      0    496      0      0
    ## 153      0      0      0     85      0      0      0    206      0      0
    ## 154    111      0    128      0      0    101     73    363    383      0
    ## 155   1609      0    270     78      0     21    227    139     68      0
    ## 156     97      0    443     10      0      0     51     75    115      0
    ## 157    305      0    193     11      0    286     75     34      0      0
    ## 158    428      0    351     82      0    612    183    103     53      0
    ## 159    856      0    431    166      0     22     87    375    453    143
    ## 160    791      0    963     57      0      0     57    161    230      0
    ## 161     18      0     23     35      0      0     70      0   1207      0
    ## 162    211      0    297     34      0    249    105     42     87    138
    ## 163      0      0      0    191      0      0      0     11      0      0
    ## 164      0      0      0    447      0      0      0      0      0    160
    ## 165      0      0      0    495      0      0      0      0      0      0
    ## 166      0      0      0    261      0      0      0     43      0      0
    ## 167      0      0     14    787      0      0      0     37      0     25
    ## 168      5      0      0    664      0      0      0    186      0    498
    ## 169      0      0      0    918      0      0      0    144      0     16
    ## 170     20      0      0    723      0      0      0     37      0    551
    ## 171      0      0      0    348      0      0      0    200      0    250
    ## 172    226    149     20    359      0      0      0     15      0     27
    ## 173     52    416     12     94      0      0      0      0      0      0
    ## 174     87      0    302     19      0      0      0      0     92     49
    ## 175    875     15    447    188      0      0      0     30     13      0
    ## 176     27      0     80    277     41      0      0    102      0      3
    ## 177    921      0      0     15      0      0      0     59    427    280
    ## 178      0      0      0    531      0      0      0      0     33      5
    ## 179      6      0     93     22      0      0      0     43     44      0
    ## 180      0     34    119    286      0      0      0     41      0     41
    ## 181      0      4    176     30    648      0      0     45     18     24
    ## 182    107    123    256      0    810      0      0     27      0      0
    ## 183    205      0     97      0    800      0     44     45      0     50
    ## 184    628      0    129      0    346      0      0     82     33     29
    ## 185    191      0     49      0    722      0     44     24     67    828
    ## 186    946      0    195     88   1361      0      9     21    129    347
    ## 187    488      0    132    177    788      0      0     45    393    204
    ## 188    684      0    429     38    781      0     37    219     15     42
    ## 189    500      0    181      0    803      0    110     74     17      0
    ## 190      0      0     38     29      0      0      0     98      0      0
    ## 191      0      0      0     53      0     13     36      0      0      0
    ## 192      0      0      0    883      0      0     50      0      0      0
    ## 193      0      0      0      0      0      0      0    283      0      0
    ## 194      0      0      0     70      0      0      0    159      0      0
    ## 195      0      0    297     49      0      0     86    216      0      0
    ## 196     30      6      0     66      0      0      0    195      0      0
    ## 197    119      0     33    361      9      0      6    139     64     72
    ## 198      0      0      0      0      0      0     15     11      0     40
    ## 199      0      0      0     49     54      0      0     65      0    553
    ## 200      0      0      0    278      0      0      7      7      0      0
    ## 201      0      0      0   1027      0      0      0      0      0      0
    ## 202      0      0      0     14      0      0      0      7      0      0
    ## 203      0      0      0   1625      0      0      0     62      0      0
    ## 204      0      0      0   1093      0      0      0      0      0      0
    ## 205      0      0      0     41      0      0      0     28      0      0
    ## 206      0      0      0    252      0      0      0      0      0      0
    ## 207      0      0      0    356     54      0      0      0      0      0
    ## 208     49    503    155    104      9    101      0      0      0      0
    ## 209     75    122    169    124      0      3      0    178     64    162
    ## 210      0      0     85      0     11     33      0     38     35     39
    ## 211     89     93    113     35    184     85    542     28     12      0
    ## 212      0      0    194    119      0      0      0    318    219      0
    ## 213     65   1535    315      0      0    146     72      0      0      0
    ## 214   1753      3    162     24    142      0    277    261      0      0
    ## 215     27    347     62     23     10    210     28     88     11      0
    ## 216    361   1280      0     32     66      0      0    119     43      0
    ## 217    323   1737     34     14     23    345      0      0      0      0
    ## 218     47    115    101     75     63      0      0     79      0      0
    ## 219    362     64     54     55     58      0      0    111      0      0
    ## 220     49      0     38     84     69      0      0    144      0      0
    ## 221    199   2084      0      0    104      0      0    171      0      0
    ## 222    400    372     12    196      0      0      0    226      0      0
    ## 223     17      0      0     77      0      0      0      0      0      0
    ## 224    126    629      0    114      0      0      0      9      0      0
    ## 225     64     31    457    135      0      0      0     19     19      0
    ##     otu_20 otu_21 otu_22 otu_23 otu_24 otu_25 otu_26 otu_27 otu_28 otu_29
    ## 1       47      0      0      0     76     33      0      0      0     64
    ## 2        0      0      0      0      0      0      0      0     64      0
    ## 3       17      0      0      0      0     19      0      0     79     14
    ## 4        0      0      0      0      0      0      0      0      0     49
    ## 5        0      0      0     94      0      0      0      0    202     71
    ## 6        0      0      0     25      0      0      0      0      0      0
    ## 7       31      0      0     15     21      0      0      0     49      0
    ## 8        0      0      0      0      0      8      0      0      0     23
    ## 9       16      0      0     34      0      0      0      0      7      0
    ## 10      26      0      0      0     91      0      0      0      0     49
    ## 11       0      0     47      0      0      0      0      0      0     11
    ## 12       0      0      0      0      0     16      0      0      0     15
    ## 13       0      0      0      0      0     53      0      0      0     42
    ## 14       0      0      0      0      0      0      0      0      0     76
    ## 15       0      0      0    371    381    307      0      0      0     44
    ## 16       0      0      0      0      0      4      0      0      0     60
    ## 17      84      0      0      0      0      0      0      0      0    189
    ## 18      12      0      0      0     49      0      0      0      0      0
    ## 19       0      0      0      0      0     55      0      0      0    191
    ## 20      12      0      0      0      0      0      0      0      0     36
    ## 21      29     12      0      0      0     18      0      0      0    108
    ## 22      46      0      9     11      0     41      0      0      0     19
    ## 23      35      0    147      0      0     16      0      0      0      0
    ## 24       9      0      0      0      0     16    158      0      0      0
    ## 25     121      0    289      0     63     10      0      0      0    154
    ## 26       0      0      0      0      0      0      0      0      0     14
    ## 27     193      0      0      0      0     93      0      0      0     96
    ## 28       0      0      0      0     37      0      0      0      0     71
    ## 29       0      0      0    123    146      0      0      0      0     41
    ## 30       8      0      0    212    116      0      0      0      0    151
    ## 31       0      0      0     44     14      0      0      0      0      0
    ## 32       0      0      0      5      0      0      0      0      0    112
    ## 33       0      0      0     78    396     47      0      0      0      0
    ## 34      19      0      0    123    183     11      0      0      0      0
    ## 35       0      0      0     62      0      0      0      0      0    146
    ## 36       0      0      0     80    142      0      0      0      0     90
    ## 37       0      0      0      0    286      0      0      0      0     81
    ## 38       0      0      0    243    276     70      0      0      0     27
    ## 39      20      0      0      0     34      0      0      0      0     24
    ## 40       0      0      0      0     18      0      0      0      0     61
    ## 41       0      0     46    140      0      0      0      0      0    101
    ## 42       0      0      0      0    972      0      0      0      0    141
    ## 43       0      0      0    715      0      0      0      0      0     20
    ## 44       0      0      0     11    647      0      0      0      0    115
    ## 45       0      0      0      0    895      0      0      0      0     81
    ## 46       0      0      0    475      0     63      0      0      0    126
    ## 47       0      0      0    222      0     28      0      0      0     99
    ## 48       0      0      0    320      0      0      0      0      0     34
    ## 49       0      0      0    219    435      0      0      0      0     41
    ## 50       0      0      0      0      0      0      0      0      0     96
    ## 51       0      0      0    297      0      4      0      0      0     96
    ## 52       0      0      0    869   1213      0      0      0      0     89
    ## 53     313      0      0    444      0      0      0      0      0     46
    ## 54       0      0      0    402    440      0      0      0      0     21
    ## 55      11      0    223      0      0    117      0      0      0     11
    ## 56      11    307      0      0      0    178      0      0      0     43
    ## 57      17     45      0      0    231     77      0      0      0      9
    ## 58      23     52      0      0      0      0      0      0      0      0
    ## 59       0      0      0      0      0     69      0      0      0      0
    ## 60       0      0      0     19      0     10      0      0      0     46
    ## 61       0      0      0      0     12     83      0      0     28      0
    ## 62       0      0      0      0      0     30      0      0      0     48
    ## 63       0      0      0      0      0     31      0      0      0     51
    ## 64     407      0      0      0      0      0      0      0      0    122
    ## 65     154     31    144      0      0      0    125      0      0     35
    ## 66      56      7    764      0      0      0      0      0      0     27
    ## 67     259      0      0      0      0      0      0      0      0      0
    ## 68     136     70      0      0      0      0      0      0      0    142
    ## 69     387    112     70      0      0      0      0      0      0    137
    ## 70     364    183      0      0      0      0     10      0      0    313
    ## 71     118    233      0      0      0      0      6      0      0    187
    ## 72     127     31      0      0      0      0      0      0      0     87
    ## 73       0      0      0      0    130      0      0      0      0    268
    ## 74       0      0      0    110    124      0      0      0    293     23
    ## 75       0      0      0    484      0     17      0      0    731      0
    ## 76       0      0      0      0      0     43      0      0      0     76
    ## 77      85      0      0      0     88     99      0      0      0    117
    ## 78      12    260      0      0    261    157      0      0      0    118
    ## 79      66      0      0      0    262     48      0      0      0    365
    ## 80       0     23      0      0     11     23     25      0      0      0
    ## 81       0      0      0    268      0      0      0      0    693     29
    ## 82       7      0      0     49     29      0      0      0      0    154
    ## 83       0    106      0      0      0    142      0      0      0     59
    ## 84       0      0      0      0    106    187      0      0      0     11
    ## 85     138      0      0      0      0      0      0      0      0     29
    ## 86      44      0      0      0    113      0      5      0      0    143
    ## 87      99     87      0      0    250    218      0      0      0    286
    ## 88      17     34      0      0      0     41      0      0      0    240
    ## 89       0      0      0      0      0      0      0      0      0    146
    ## 90       0     77      0     13     27    352    137      0      0    230
    ## 91      15      0      0      4    207     51    266      0      0    163
    ## 92       0      0      0      0    137      0      0      0      0    102
    ## 93      39      0      0      0    121     59      0      0      0     74
    ## 94      18      0      0      0     11      0      0      0      0     29
    ## 95       0      0      0    273      0      0      0      0      0      0
    ## 96       0      0      0     10      0      0      0      0      0     87
    ## 97       0      0      0      0      0      0      0      0      0     43
    ## 98       0      0      0     76      0      0      0      0      0     33
    ## 99      16      0      0    169     54      0      0      0      0    104
    ## 100     14      0      0      0     33     59      0      0      0     37
    ## 101     15      0      0      8      0     37     42      0      0     34
    ## 102     27      0      0      0    175     32      0      0      0     59
    ## 103      0      0      0      0    192     25      0      0      0     25
    ## 104      0      0      0      0    105    230    113      0      0      0
    ## 105     15      0      0      0    130    288      0      0      0      0
    ## 106      0      0      0      0    320    122    111      0      0      0
    ## 107     11      0      0     59     98      0    154      0      0     82
    ## 108     27      0      0      0     44    226    120      0      0      0
    ## 109      0      0      0    246     55      0      0      0    191      0
    ## 110      0      0      0      0    161      0      0      0    333     39
    ## 111      0      0      0      0     75      0      0      0    142     31
    ## 112      0      0      0    139     40      0      0      0    254     91
    ## 113      0      0      0      0     22    109      0      0     45     43
    ## 114      0      0      0    236     85      0      0      0     63     27
    ## 115      0      0      0    197     72      0      0      0     31      0
    ## 116      0      0      0     68     49     89      0      0    115      0
    ## 117      0      0      0    202     40      0      0      0     18     19
    ## 118      0      0      0      0      0      0      0      0     99     60
    ## 119      0      0      0     14     26     15      0      0     50     42
    ## 120      0      0      0    307    167      0      0      0    454     40
    ## 121      0      0      0    160      0      0      0      0     27     23
    ## 122      0      0      0     32    379     26      0      0    932      0
    ## 123      0      0      0    386     98     24      0      0    239     74
    ## 124     43      0    275    169     79      0      0      0    346     32
    ## 125      0      0      0     94     40      0      0      0      0     14
    ## 126     10      0      0     53     22      0      0      0    226      0
    ## 127     40    149    416      0    172      0    312      0      0    172
    ## 128      0      0     13      0     69    133    345      0      0     28
    ## 129      0      0      4      0     66     17    274      0      0    255
    ## 130      9      0      0      0     30     40      0      0      0     93
    ## 131    183      0      0      0    127    135    120      0      0     55
    ## 132     26    550    139      0    172      6    180      0      0    165
    ## 133      0      0      0      0    290     53     47      0      0    126
    ## 134     36      0      0      0      0    108     45      0      0     79
    ## 135    173    181      0      0    236     17     58      0      0      0
    ## 136    124      0      0      0     23     30    428      0      0     59
    ## 137     62    606     14      0      0      0    163      0      0     51
    ## 138     38    343      0      0     16    215    656      0      0     33
    ## 139     98    134      0      0     77    168    358      0      0    151
    ## 140    202      3      0      0     28     10    240      0      0     34
    ## 141     68      0      0      0    187    189    583      0      0     93
    ## 142     55      0     11      0    193     33    824      0      0    192
    ## 143    165      0     53      0     32     33     43      0      0    103
    ## 144     48      0    215      0      0      0      0      0      0    120
    ## 145     75     23      0      0     32      0      8      0      0    139
    ## 146     26      0      0      0      0      0    154      0      0    228
    ## 147     89    448      0      0      0      0    119      0      0     78
    ## 148    315     31      9      0      0      0      7      0      0     95
    ## 149    376     27      0      0     95     11      0      0      0     82
    ## 150    343    380     82      0      0     10   1436      0      0     63
    ## 151    148   1677      0      0      0      9    217      0      0    242
    ## 152    242    165      0      0      0      0      0      0      8    218
    ## 153     48     22      0      0      0      0     21      0      0     31
    ## 154      0    666      0      0      0    227      0      0      0    120
    ## 155     17     61      0      0      9    215      0      0      0     76
    ## 156      6    150    204      0      0     63      0      0      0     43
    ## 157     14      0    561      0      0     24      0      0      0      0
    ## 158     34      0      0      0      0    161      0      0      0    103
    ## 159     96     12      0      0      0     43      0      0      0    218
    ## 160     42      0     92      0      5    384      0      0      0    119
    ## 161      8      0    255      0     24      0      0      0      0      0
    ## 162     21      0     68      0      0     39    185      0      0     23
    ## 163     94   1199      0      0      0      7      0      0      0      0
    ## 164    193    551      0      0      0      0      0      0      0      0
    ## 165    239    302      0      0      0      0      0      0      0     12
    ## 166     33   1116      0      0      0      0      0      0      0     42
    ## 167    268    717      0      0      0      0      0      0      0     83
    ## 168    348   1014     25     16      0      0      0      0      0    244
    ## 169    416   2627     35      0      0      0      0      0      0     86
    ## 170    171    434     22      0      0      0     12      0      0     32
    ## 171    248    305      0      0      0      0      0      0      0    158
    ## 172    192    890      0      0      0      0      0      0      0     26
    ## 173     42    195      0      0      0     23      0      0      0      0
    ## 174     10     55      0      0     11     17      0      0      0      0
    ## 175     97      0      0      0      0      0      0      0      0      0
    ## 176    172      0     18      0      0     21      0      0      0     49
    ## 177      0      0      0      0     10      0      0      0      0     13
    ## 178    268      0      0      0      0      0      0      0      0      0
    ## 179     19      0    597      0      0      0      0      0      0     34
    ## 180    135      0     27      0     57     67      0      0      0     22
    ## 181      0    179     84      0      0      0      0      0      0     18
    ## 182      0      0    956      0      0     19      0      0      0     19
    ## 183      0     96    432      0      0      0      0      0      0     24
    ## 184      0     55    185      0     40     18      0      0      0     90
    ## 185      0     76    779      0      0      0      0      0      0     22
    ## 186     23      6    755      0      0     68      0      0      0      0
    ## 187    133    713      0      0      0      0      0      0      0     14
    ## 188     14      5     94      0    257     45     54      0      0    136
    ## 189      0    104     33      0      0     18      0      0      0     70
    ## 190     36     10      0      0      0      0      0      0      0     51
    ## 191      0      0      0    452     15      0      0    770      0      0
    ## 192    218      0      0      0      0      0      0    201      0      0
    ## 193      0      0    107      0      0      0     19    242      0     59
    ## 194      0     57      0      0      0      0      0    349      0     79
    ## 195      8      0   1207     67      0     19      0    409      0     27
    ## 196     15      0    241      0      0      0      0    983      0     91
    ## 197    256      0     90      0     46     30     63    141      0    133
    ## 198     27      0     34      0      0      0      0   1000      0      0
    ## 199     44      0      0      0      0      0      0      0      0     34
    ## 200    135     72     46      0      0      0      0      0      0      0
    ## 201    649     23      0      0      0      0      0      0      0      0
    ## 202     24     66      0      0      0      0      0     23      0      0
    ## 203    952      0      0      0      0      0      0      0      0      0
    ## 204    786    833      0      0      0      0      0      0      0      0
    ## 205     58     27     78      0      0      0      0      0      0     19
    ## 206     44      0      0      0      0      0      0     32      0      0
    ## 207    150    103      0      0      0      0      0      0      0      0
    ## 208     66      0     21      0      0     20      0      0      0     11
    ## 209     67     99    471      0     32     77     46      0      0    108
    ## 210      0      0    277     34      0      0     22      0      0     24
    ## 211     18      0    500      0      0    105      0      0      0     11
    ## 212     27     18      0      0      0     11      0      0      0    168
    ## 213      5      0      0      0      0     14      0      0      0      0
    ## 214      0    901    799      0     88     33      0      0      0    291
    ## 215     12    539      0      0     40      0      0      0      0     32
    ## 216     15      0    109      0      0    183      0      0      0    129
    ## 217      0      0      0      0      0     65      0      0      0      0
    ## 218     46      0      0      0     33      0      0      0      0     68
    ## 219     17      0     76      0      0      0      0      0      0    139
    ## 220     34      0     21      0      0      0      0      0      0     94
    ## 221      0      0      0      0      0      0      0      0      0     77
    ## 222     94      0      0      0     26      0      0      0      0      0
    ## 223     25      0     32      0     24      0      0      0      0      0
    ## 224     62      0      0      0     12     40      0      0      0     22
    ## 225    121      0      0      0      0     44      0      0      0     43
    ##     otu_30 otu_31 otu_32 otu_33 otu_34 otu_35 otu_36 otu_37 otu_38 otu_39
    ## 1        0      0     16      0     21      0      0      0      0      0
    ## 2        0      0      0      0     31      0      0      0      0      0
    ## 3        0      0      0      0      0     22      0     38      0      0
    ## 4        0      0      0      0    117      0      0    304      0      0
    ## 5        0      0      0      0      0      0      0      0      0      0
    ## 6        0      0      0      0      0      0      0      0      0      0
    ## 7       32      0     33      0      0      0      0      0      0      0
    ## 8        0      0      0      0    119      5      0      0    177      0
    ## 9       12      0     49      0     30      0      0      0      0      0
    ## 10       0      0      0      0     93      0      0     77      0      0
    ## 11       0      0      0      0      0      0      0      0      0      0
    ## 12       0      0      0      0     46      0      0      0      0      0
    ## 13       0      0      0      0    125      0      0     25      0      0
    ## 14       0      0      0      0    178      0      0     47      0      0
    ## 15       0      0      0      0      0     10      0    233      0      0
    ## 16       0      0      0      0     12      0      0    229      0      0
    ## 17       0      0      0      0     14      0      0    321      0      0
    ## 18       0      0      0      0     96      0      0    267      0      0
    ## 19       0      0      0      0      0     54      0      0      0      0
    ## 20      19      0      0      0      0     25      0      0      0      0
    ## 21       0      0      0      0      0      0      0      0      0      0
    ## 22       0      0      0      0     87     23     19      0     61     40
    ## 23      58      0      0      0     45     21     17    348    123      0
    ## 24       0      0      0      0      0    105     27      0      0     31
    ## 25       0      0      0      0     34     20      0      0     26      0
    ## 26     143      0      0      0     13      9     18      0      0     10
    ## 27       0      0      0      0      0      0      0     80    154      0
    ## 28      12      0      0      0      8     43      0      0      0      0
    ## 29      89      0      0      0     42      0      0    337    172      0
    ## 30     179      0      0      0      0      0      0      0      0      0
    ## 31      49      0      0      0      0      0      0      0     18      0
    ## 32     340      0      0      0      7      0      0      0    389      0
    ## 33     446      0      0      0     61     24     20      0     21      0
    ## 34     241      0      0      0     14     27      0      0      0      0
    ## 35     435      0      0      0     18      0      0      0    218      0
    ## 36     209      0      0      0      0      0      0      0      0      0
    ## 37     120      0      0      0     87      0      0    181     58      0
    ## 38      37      0      0      0     41      0      0    360     36      0
    ## 39     107      0      0      0      0     25      0     69      0      0
    ## 40      17      0      0      0     44     12      0    212      0      0
    ## 41     276      0      0      0     59      0      0    219     75      0
    ## 42       0      0      0      0      0      0      0     97     82      0
    ## 43      28      0      0      0    173      0      0      0     52      0
    ## 44       0      0      0      0     25      0      0    355     79      0
    ## 45       0      0      0      0     51      0      0      0    106      0
    ## 46      22      0      0      0      0      0      0      0     94      0
    ## 47       0      0      0      0      0      0      0     91    200      0
    ## 48       0      0      0      0      0     63      0    442     41      0
    ## 49       0      0      0      0      0      0      0     39      0      0
    ## 50     188      0      0      0      0    157      0    391     10      0
    ## 51      17      0      0      0      0     14      0      0      7      0
    ## 52       0      0      0      0    216      0      0     21     75      0
    ## 53       5      0      0      0      0     11      0    268     17      0
    ## 54       0      0      0      0      0      0      0      0     90      0
    ## 55       0      0      0      0    115      0      0    114      0      0
    ## 56      12      0      0      0     82      0      0     48      0      0
    ## 57       0      0      0      0      0      0     37      0      2     14
    ## 58       0      0      0      0    151      0      0    247      0      0
    ## 59       0      0      0      0     86      0      0     45      0      0
    ## 60       0      0      0      0     77      0      0      0      0      0
    ## 61       0      0      0      0      0      0      0      8    462      0
    ## 62       0      0      0      0    334      0     13    749      0     16
    ## 63       0      0      0      0      8      0      0      0      0      0
    ## 64       0      0      0     57      0      0      0      0      3      0
    ## 65       0      0      0     66      0      0      0      0      0      0
    ## 66       0      0      0    144      0      0      0      0    464      0
    ## 67       0      0      0     14      0      0      0      0     37      0
    ## 68       0      0      0    113      0      0      0      0    148      0
    ## 69       0     29      0    309      0     75      0      0    205      0
    ## 70       0   1469      0    232      0     39      0      0      0      0
    ## 71       0    450      0     30      0      0      0      0    172      0
    ## 72       0      0      0     47      0      0      0      0    801      0
    ## 73      87      0     30      0      0      0      0      0      0      0
    ## 74      12      0      0      0      0      0      0     46     37      0
    ## 75       0      0      0      0      0      0      0    177     54      0
    ## 76     194      0      0      0      0      0      0    427      0      0
    ## 77      29      0      0      0      0      0      0      0     21      0
    ## 78       7      0      0      0      0      0      0      0      0      0
    ## 79       0      0      8      0      0    278      0      0      0      0
    ## 80      30      0      0      0     41      0      0      0      0      0
    ## 81       4      0      0      0      0      0      0    257      0      0
    ## 82     200      0    195     51      0     57      0      0     24      0
    ## 83      28      0      0     40      0     18      0      0    195      0
    ## 84       0      0      0     59      0    156      0      0      0      0
    ## 85     198      0     56     47      0     51      0      0      0      0
    ## 86      31      0    260    133      0     85      0      0     47      0
    ## 87       0      0    100     48      0     87      0      0     27      0
    ## 88     115      0     94     13      0    197      0      0     52      0
    ## 89      38      0     81    194      0     56      0      0    363      0
    ## 90      13      0    206      0      0     66      0      0    153      0
    ## 91       0      0      0     21      0     14      0      0     35      0
    ## 92       0      0      0     50      0      0      0      0     29      0
    ## 93       0      0      0     30      0     86      0      0      0      0
    ## 94       0      0      0      0      0      0      0      0      0      0
    ## 95       0      0      0     56      0      0      0      0      0      0
    ## 96       0      0      0     24      0    348      0      0     38      0
    ## 97       0      0      0     13      0      0      0      0     22      0
    ## 98       0      0      0      0      0      0      0      0      0      0
    ## 99       0      0      0     19      0      0      0      0     24      0
    ## 100      0      0      0      0      0     29      0      0      0      0
    ## 101      0      0      0      0      0      0      0      0      0      0
    ## 102      0      0    441      0      0     26     74      0      0     66
    ## 103      0      0    239      0      0     29     25      0      0      0
    ## 104      0      0      0      0      0     73     35      0      0     20
    ## 105      0      0      0      0      0     39    132      0     40     40
    ## 106      0      0     87      0      6      0     22      0      0     18
    ## 107      0      0      0      0      0     23      0      0      0      0
    ## 108      0      0    696      0      0      0    506      0      0    238
    ## 109      0      0      0      0    151      0      0      8      0      0
    ## 110     66      0      0      0    100      0      0     22     52      0
    ## 111      7      0      0      0      0      0      0      0      0      0
    ## 112      0      0      0      0     82      0      0     78      0      0
    ## 113      0      0      0      0     73      0      0     55      0      0
    ## 114      0      0      0      0     67      0      0     12      0      0
    ## 115      0      0      0      0    101      0      0      0     96      0
    ## 116      0      0      0      0     51      0      0    110      0      0
    ## 117      0      0      0      0      0      0      0     12      0      0
    ## 118      0      0      0      0     20      0      0      0      0      0
    ## 119     28      0      0      0     33     38      0      0      0      0
    ## 120     19      0      0      0      0      0      0      0      0      0
    ## 121      2      0      0      6      0      0      0      0      0      0
    ## 122     50      0      0      0     27    104      0      0      0      0
    ## 123    117      0      0      0      0      0      0      0     19      0
    ## 124     70      0      0      0     94      0      0      0     98      0
    ## 125      0      0      0      0      0      0      0      0      0      0
    ## 126    129      0      0      0      0      0      0      0      0      0
    ## 127      0      8      0      0      0      0      0      0     74      0
    ## 128      0      0      0      0    106      0     42      0     42     49
    ## 129    293      0     29      0      0      0     18      0      0     19
    ## 130      0    174      0      0      0     18      0      0     20      0
    ## 131      0     64      0      0      0      0      0      0      0      0
    ## 132      0      0     41     32      0      0      0      0     18      0
    ## 133      0     41      0      0      0      0      0      0    101      0
    ## 134      0      0      0      0      0      0     36      0     14     27
    ## 135      0    148      0      0      0      0     19      0     26     14
    ## 136     13    149      0      0      0      0     37      0      0     22
    ## 137      0     79      0     34      0      0     15      0      2     11
    ## 138      0    200      0      8      0      0      0      0      0      0
    ## 139     19     39      0      0      0      0      0      0      0      0
    ## 140     29     36      6     53      0      0    122      0      5    201
    ## 141     24      0      0      0      0     52     42      0     26     15
    ## 142      0      0      0     27      0     75      0      0     47      0
    ## 143     12     26      0      0      0      8     21      0    284     15
    ## 144     50    318      0      0      0      0      7      0      0      0
    ## 145      0      0      0      0      0      0      0      0      0      0
    ## 146      0      0     47     87      0      0      0      0      0      0
    ## 147      0      0      0      0      0     40      0      0     22      6
    ## 148      7      0      0      0      0      0      0      0      0      0
    ## 149      0      0      0      0      0      0      0      0    241      0
    ## 150      0      0      0      0      0      0      0      0     15      0
    ## 151      0      0      0      0      0      0      0      0      0      0
    ## 152      0      0      0     86      0      0      0      0     19      0
    ## 153      0      0      0      0      0      0      0      0      0      0
    ## 154      6      0     25     17      0     56      0      0     48      0
    ## 155      0      0      0      0      0     34      0      0      9      0
    ## 156    161      0      0     10      0     42      0      0      0      0
    ## 157    120      0     39    112      0     45      0      0     34      0
    ## 158      0      0      0     31      0     50      0      0      0      0
    ## 159     40      0      0      0      0     59      0      0    119      0
    ## 160    185      0     15     72      0    106      0      0      0      0
    ## 161      0      0      0     12      0      0      0      0     27      0
    ## 162      0      0      0     24      0      0      0      0     29      0
    ## 163      0      0      0      0      0      0      0      0      0      0
    ## 164      0      0      0     20      0      0      0      0      0      0
    ## 165      0      0      0      0      0      0      0      0      0      0
    ## 166      0      0      0      0      0      0      0      0      0      0
    ## 167      0      0      0     51      0      0    440      0      0    336
    ## 168      0      3      0    258      0      0      0      0    171      0
    ## 169      0    205      0     97      0     15     13      0      0      0
    ## 170      0      0      0     18      0      0    281      0      0      0
    ## 171      0      0      0      0      0      0      0      0      0      0
    ## 172      0      0      0      0     78      0      0      0      0      0
    ## 173      0     26      0      0      0     87      0      0      0      0
    ## 174      0    227      0      0      0      0      0      0      0      0
    ## 175     46      0      0      0      0      0      0      0     21      0
    ## 176    124      0      0      0      0      0      0      0      0      0
    ## 177    139      0      0      0     47      0      0      0      0      0
    ## 178     29      0      0      0     47      0      0      0      0      0
    ## 179      0     15      0      0      0     28      0      0      0      0
    ## 180      0      0      0      0     34      0      0      0     16      0
    ## 181      0     92      0      0      0     10    392      0      0    432
    ## 182      0     57      0      8      0      0      0      0      2      5
    ## 183      0     73      0      0      0      0    186      0      0     77
    ## 184      0     20      0      0      0      0    179      0     16    113
    ## 185      0     42      0      0      0      0     69      0     27     87
    ## 186      0     56      0      9      0     21    345      0     34    317
    ## 187      0    103      0     10      0     17    352      0      0    329
    ## 188      0     91      0      0      0      0     88      0    141     69
    ## 189      0     11      0      0      0      0    340      0     10    145
    ## 190      0      8      0     62      0      0      0      0      0      0
    ## 191      0      0      0     58      0     77      0      0      0      0
    ## 192      0      0      0      0      0      0      0      0      0      0
    ## 193      0      0      0    199      0      0      0      0    146      0
    ## 194      0     31      0    278      0      0     39      0    762      0
    ## 195      0      0      0    295      0      0      0      0     29      0
    ## 196      0      0      0    166      0      0      0      0    246      0
    ## 197     15     26      0     29      0      0      0      0    336     24
    ## 198      0      0      0     91      0      0      0      0    118      0
    ## 199      0      0      0     25      0      0     12      0      0      0
    ## 200      0      0      0      4      0      0     53      0    288      0
    ## 201      0      0      0      0      0      0      0      0    549      0
    ## 202      0      0      0      0     22      0      0      0     21      0
    ## 203      0      0      0      0      0      0      0      0      0      0
    ## 204      0      0      0      0      0      0      0      0      0      0
    ## 205      0      0      0    387      0     12      0      0    229      0
    ## 206      0      0      0     10      0     34      0      0     16      0
    ## 207      0      0      0      0      0      0      0      0    294      0
    ## 208      0      0      0      0      0      0    113      0      7     24
    ## 209    257      0      0      0      0      0      0      0     12      0
    ## 210    117      0      0      0     70     43      0      0      0      0
    ## 211    141      0      0      0      0     23      0      0     25      0
    ## 212    180      0   1310     10      0     70      0      0    271      0
    ## 213      0      0      0      0     12     47      0      0      0      0
    ## 214     31      0      0     47      0     56    464      0     12      0
    ## 215     17     16      0     18      0      0      0      0      0      0
    ## 216      0      0      0      0      0      0     26      0      0      0
    ## 217      0      0      0      0     91      0     18     13     43      0
    ## 218      0      0      0      0     31      0      0      0      0      0
    ## 219      0      0      0      0     11      0    207      0     30     59
    ## 220     55      0      0      0      0      0      0      0      0      0
    ## 221      0      0      0      0    271      0      0    571      0      0
    ## 222      0      0      0      0      0      0      0      0      0      0
    ## 223    134      0      0      0     82      0      0      0      0      0
    ## 224      0      0      0      0     67      0      0      0     26      0
    ## 225      0      0      0      0      0      0      0      0      0      0
    ##     otu_40 otu_41 otu_42 otu_43 otu_44 otu_45 otu_46 otu_47 otu_48 otu_49
    ## 1        0      0      0      0      0     24      0      0      0    108
    ## 2        0      0      0      0      0      0      0      0    178    304
    ## 3        0      0      0      0      0     16      0      0    910    107
    ## 4        0      0      0      0     14      0      0      0      0     15
    ## 5        0      0      0      0      0      0      0      0    405     93
    ## 6        0      0      0      0      0      0      0      0     39      0
    ## 7        0      0      0      0      0     17      0      0    128    131
    ## 8        0      0      0      0     30      0      0      0    162     18
    ## 9        0      0      0      0      0      0      0      0     45      0
    ## 10     136      0      0      0      4     27      0      0      0      0
    ## 11       0      0      0      0      0      0      0      0      0    788
    ## 12      47      0      0      0      0      0      0      0      0     89
    ## 13      92      0      0      0      0      0      0      0      0      0
    ## 14       0      0      0      0      0      0      0     31      0      0
    ## 15     168      0      0      0      0      0      0      0      0    219
    ## 16       0      0      0      0      0      0      0      0      0    242
    ## 17      18      0      0      0      0     11      0      0      0    130
    ## 18      19      0      0      0      0      0      0      0      0      0
    ## 19       0      0      0      0     32      0      0     50      0      0
    ## 20      51      0      0      0      0     16      0     76      0      0
    ## 21      94      0      0      0      0      6      0      0      0    113
    ## 22      25      0      0      0      0     40      0      0      0      0
    ## 23      71      0      0      0      0     11      0      0      0      0
    ## 24      30      0      0      0      0      6      0      0      0      0
    ## 25     221      0      0      0      0     39      0      0      0     10
    ## 26     154      0      0      0      0      0      0     75      0     42
    ## 27      37      0      0      0      0     98      0      0      0     37
    ## 28      11      0      0      0      0      0      0      0      0      0
    ## 29       0      0      0      0      0      0      0      0      0      0
    ## 30     102      0      0      0      0      6      0      0      0      0
    ## 31      16      0      0      0      0      0      0      0      0      0
    ## 32       0      0      0      0      0      0      0      0      0      0
    ## 33     100      0      0      0     19      0      0      0      0      0
    ## 34       0      0      0      0      0      0      0      0      0      0
    ## 35     236      0      0      0      0      0     19      0      0      0
    ## 36       0      0      0      0      0      0      0      0      0      0
    ## 37      77      0      0      0      0      0      0      0      0      0
    ## 38       9      0      0      0      0      0      0      0      0      0
    ## 39      11      0      0      0      0      0      0      0      0      0
    ## 40       7      0      0      0      0      0      0      0      0      0
    ## 41       0      0      0      0      0      0      0      0      0      0
    ## 42      71      0      0      0     32      0      0      0      0      0
    ## 43       0      0      0      0      0      0      0      0    469      0
    ## 44       0      0      0      0      0      0      0      0      0      0
    ## 45      20      0      0      0      0      0      0      0      0      0
    ## 46      65      0      0      0      0      0      0      0      0      0
    ## 47       0      0      0      0      0      0      0      0      0      0
    ## 48      19      0      0      0      0      0      0      0      0      0
    ## 49       0      0      0      0      0      0      0      0      0      0
    ## 50       0      0      0      0      0      0      0      0      0      0
    ## 51      31      0      0      0      0      0      0      0      0      0
    ## 52      69      0      0      0      0      0      0      0      0      0
    ## 53       0      0      0      0      0     66      0      0      0      0
    ## 54      46      0      0      0      0      0      0      0      0      0
    ## 55       0      0      0      0      0     11      0      0      0      0
    ## 56       0      0      0      0      0     20      0      0      0      0
    ## 57       0      0     46      0     29      0      0      0      0     26
    ## 58       0      0      0      0      0      0      0      0      0      0
    ## 59       0      0     10      0      0      0      0      0      0      0
    ## 60       0      0      0      0      0      0      0      0      0    148
    ## 61       0      0      0      0     44      0      0      0      0    279
    ## 62       0      0     50      0      0      0      0      0      0      0
    ## 63       0      0     12      0     14      0      0      0      0      0
    ## 64       0    131      0     63      0      0      0     40      0      0
    ## 65       0      0      0      0      0    176      0      4      0      0
    ## 66       0      0     80      0      0     43      0    740      0      0
    ## 67       0      0      0      0      0    301      0      0      0      0
    ## 68       0      0     28     11      0    158      0      0      0      0
    ## 69       0      7      0      0      0    156      0      5      0      0
    ## 70       0      0      0      0      0    380      0      0      0     17
    ## 71       0      0      0      0      0    120      0     95      0      0
    ## 72       0      0      0      0      0     52      0      0      0      0
    ## 73       0      0      0      0      0      0      0      0      0      0
    ## 74       0      0      0      0      0      0      0      0      0      0
    ## 75       0      0     28      0      0      0      0      0      0      0
    ## 76       0      0      0      0      0      0      0      0      0      0
    ## 77       0      0    160      0     10     37      0      0      0      0
    ## 78       0      0      0      0    489      0      0      0      0      0
    ## 79       0      0      0      0     31      0      0      0      0      0
    ## 80       0      0      0      0    137      0      0      0      0      0
    ## 81       0      0      0      0      0      0      0      0      0      0
    ## 82       0      0      0      0      0     12      0      0      0      0
    ## 83       0      0      0      0      0      3      0      0      0      0
    ## 84       0      0      0      0      0      0      0      0      0      0
    ## 85       0      0      0      0      0     71      0      0      0      0
    ## 86       0      0     43      0      0     36      0      0      0      0
    ## 87       0      0      0      0      0     60      0      0      0      0
    ## 88       0      0      0      0    856     16      0      0      0      0
    ## 89       0      0      0      0     38      0      0      0      0      0
    ## 90       0      0      0      0      0      0      0      0      0      0
    ## 91       0      0     81      0     67      8      0      0      0      0
    ## 92       0      0      0      0      0      0      0      0      0      0
    ## 93       0     15      0      0      0     23      0      0      0      0
    ## 94       0      0      0      0      0     10      0      0      0      0
    ## 95       0      0      0      0      0      0      0      0      0      0
    ## 96       0      0      0      0      0      0      0      0      0      0
    ## 97       0      0      0      0      0      0      0      0      0      0
    ## 98       0      0      0      0      0      0      0      0      0      0
    ## 99       0      0      0      0      0     16      0      0      0      0
    ## 100      0      0      0      0      0      6      0      0      0      0
    ## 101      0      0      0      0      0      6      0      0      0      0
    ## 102      0      0     44      0      0     10      0     50      0      0
    ## 103      0      0      0      0      0      0      0      0      0      0
    ## 104      0      0      0      0      5      0      0      0      0      0
    ## 105      0      0     70      0      0      0    916      0      0      0
    ## 106      0      0     58      0     30      0      0      0      0      0
    ## 107      0      0      0      0      0      0   1030      0      0      0
    ## 108      0      0     12      0      0      4      0    264      0      0
    ## 109      0      0      0      0      0      0      0      0      0      8
    ## 110      0      0      0      0      0      0      0      0      0      0
    ## 111      0      0      0      0      0      0      0      0    272     26
    ## 112      0      0      0      0      0      0      0      0    187     58
    ## 113      0      0      0      0      0      0      0      0      0    120
    ## 114      0      0      0      0      0      0      0      0    320      0
    ## 115      0      0      0      0      0      0      0      0    207      0
    ## 116      0      0      0      0      0      0      0      0    482     10
    ## 117      0      0      0      0      0      0      0      0      0      0
    ## 118      0      0      0      0      0      0      0      0      0      0
    ## 119      0      0      0      0      0     13      0      0    129      0
    ## 120      0      0      0      0      0      0      0      0      0      0
    ## 121      0      0      0      0      0      0      0      0      0      0
    ## 122      0      0      0      0      0      0      0      0      0      0
    ## 123      0      0    314      0      0      0      0      0      0      0
    ## 124      0      0      0      0      0     16      0      0      0      0
    ## 125      0      0      0      0      0      0      0      0      0      0
    ## 126      0      0      0      0      0      0      0      0      0      0
    ## 127      0      0     47      0      0     21      0     97      0      0
    ## 128     17      0      0      0     68      0      0     38      0      0
    ## 129      0      0     96      0    167      0     75      0      0      0
    ## 130      0      0     15      0     48      0      0    225      0      0
    ## 131      0      0     74      0      0     83      0     80      0      0
    ## 132      0      0     83      0    150     18      0    126      0      0
    ## 133      0      0    305      0      0      0     78      0      0      0
    ## 134      0      0      0      0     77     13      0      0      0      0
    ## 135      0      0     31      0      6     93      0     10      0      0
    ## 136      0      0    307      0     34     71      0     17     55      0
    ## 137      0      0     49      0      3     31      0     31      0      0
    ## 138      0      0     44      0      0     11      0    188      0      0
    ## 139      0      0      0      0      0     41      0     31      0      0
    ## 140      0      0     66      0      0     49      0    140      0      0
    ## 141      0      0     56      0      0     35      0     15      0      0
    ## 142      0      0      0      0      0     24      0     14      0      0
    ## 143      0      0      0      0      0     51      0     47      0      0
    ## 144      0      0    170      0      6     29      0     53      0      0
    ## 145      0      3      0      0      0     30      0     39      0      0
    ## 146      0      0     23      0      0     39      0      0      0      0
    ## 147      0      0      0      0      0    112      0     22      0      0
    ## 148      0     35      0      0      0    190      0     93      0      0
    ## 149      0     73      0      0      0    171      0      7      0      0
    ## 150      0      0      0      0      0    330      0    331      0      0
    ## 151      0      0      0      0      0    160      0     60      0      0
    ## 152      0     59      0      0      0     45      0     66      0      0
    ## 153      0      0      0     33      0     74      0      0      0      0
    ## 154      0      0      0      0     14      0      0    381      0      0
    ## 155      0      0      0      0      0     31      0    205      0      0
    ## 156      0      0      0      0      0      0      0      0      0      0
    ## 157      0      0      0      0      0      8      0      0      0      0
    ## 158      0      0      0      0      0     14      0    157      0      0
    ## 159      0      0      0      0      0     25      0    307      0      0
    ## 160      0      0      0      0      0     13      0    277      0      0
    ## 161      0      0      0      0      0     11      0    112      0      0
    ## 162      0      0      0      0      0     11      0     11      0      0
    ## 163      0     43      0      0      0     19      0      0      0      0
    ## 164      0    250      0      0      0      0      0      0      0      0
    ## 165      0    111      0      0      0     59      0      0      0      0
    ## 166      0     13      0      0      0     46      0      0      0      0
    ## 167      0     23      0      0      0    197      0      0      0      0
    ## 168      0    127      0      0      0     70      0      0      0      0
    ## 169      0    171      0      0      0    135      0      0      0      0
    ## 170      0     12      0      0      0    154      0      0      0      0
    ## 171      0     63      0      0      0     28      0      0      0      0
    ## 172      0      0      0      0      0     93      0      0      0      0
    ## 173      0      0      0      0      7     28      0      0      0      0
    ## 174      0      0      0      0     34      6      0      0      0      0
    ## 175      0      0      0      0     15     29      0      0      0      0
    ## 176      0      0      0      0      0     50      0      0      0    226
    ## 177      0      0      0      0      0      0      0      0      0      0
    ## 178      0      0      0      0      0    116      0      0      0      0
    ## 179      0      0      0      0      4      9      0      0      0      0
    ## 180      0      0      0      0      0     49      0      0      0      0
    ## 181      0      0      0      0      0      0      0      0      0      0
    ## 182      0      0      0      0     17      0      0     31      0      0
    ## 183      0      0      0      0      2      0      0     17      0      0
    ## 184      0      0      0      0      0      0      0      0      0      0
    ## 185      0      0      0      0     13      0      0     35      0      0
    ## 186      0      0      0      0      0      8      0     26      0      0
    ## 187      0      0      0      0      0     67      0      0      0      0
    ## 188      0      0    287      0      0     15     65      0      0      0
    ## 189      0      0      0      0      0      0      0      0      0      0
    ## 190      0      0      0      0      0     17      0      0      0      0
    ## 191      0      0      0      0      0      0      0      0      0      0
    ## 192      0      0      0     19      0    219      0      0      0      0
    ## 193      0      0      0      0      0      0      0      0      0      0
    ## 194      0      0      0      0      0      0      0     14      0      0
    ## 195      0      0      0      0     18     11      0      0      0      0
    ## 196      0      9      0      0      0     14      0      0      0      0
    ## 197      0      0      0     19      0     78      0     49      0      0
    ## 198      0      0      0    747      0      0      0      0      0      0
    ## 199      0     13      0     23      0      0      0      8      0      0
    ## 200      0     28      0   1145      0     31      0      0      0      0
    ## 201      0    137      0      0      0     23      0      0      0      0
    ## 202      0     42      0    227    173      0      0      0      0      0
    ## 203      0    721      0    144      0      0      0      0      0      0
    ## 204      0    480      0   2108     42      0      0      0      0      0
    ## 205      0     17      0      0      0      0      0      0      0      0
    ## 206      0      0      0      7      0     52      0      0      0      0
    ## 207      0     44      0     22      0    116      0      0      0      0
    ## 208     10      0      0      0     16     18      0     97      0      0
    ## 209     87      0      0      0      0     22      0     60      0      0
    ## 210     88      0      0      0      0      0      0     10      0      0
    ## 211    146      0      0      0      0      5     32      0      0      0
    ## 212      0      0      0      0     31     26      0      0      0      0
    ## 213      0      0      0      0      0     10      0      0      0    353
    ## 214      0      0    181      0     20      0      8      0      0      0
    ## 215     34      0      0      0      0     10     79      0      0      0
    ## 216    221      0      0      0     28      0      0     67      0      0
    ## 217     68      0      0      0      0      0     13      0      0      0
    ## 218      0      0      0      0      0      0      0      4      0     33
    ## 219     42      0      0      0      0      0      0    193      0      0
    ## 220     19      0      0      0      0      9      0      0      0     72
    ## 221      0      0      0      0     11      0      0     24      0    337
    ## 222     86      0      0      0      7     60      0    201      0      0
    ## 223     15      0      0      0     24     12      0     10      0     45
    ## 224     46      0      0      0      8     23      0      0      0      0
    ## 225      0      0      0      0      0     30      0      0      0      0
    ##     otu_50 otu_51 otu_52 otu_53 otu_54 otu_55 otu_56 otu_57 otu_58 otu_59
    ## 1        0      0      0      0      0     20      0      0      0      0
    ## 2        0      0      0      0      3      0      0      0      0      0
    ## 3        0      0      0      0      0      0      0      0      6      0
    ## 4        0      0      0      0      0      0      0      0     41      0
    ## 5        0      0      0      0      0      0      0    202      0      0
    ## 6        0      0      0      0      0      0      0      0      0      0
    ## 7        0      0      0      0      0      0      0      0      0      0
    ## 8       11      0      0      0    217      0      0      0      0      0
    ## 9        0      0      0      0      7      0      0      0      0      0
    ## 10      78      0      0      0      5    288      0      0     18      0
    ## 11       0      0      0      0     10      0      0      0      0      0
    ## 12      47      6      0      0      0     16      0      0      0      0
    ## 13      10      0      0      0      0     12      0      0      0      0
    ## 14     306      0      0      0      0      0      0      0      0      0
    ## 15      79      0      0      0      4      0      0      0     40      0
    ## 16       0      0      0      0      0      0      0      0     19      0
    ## 17       0      0      0      0     88      0      0      0     73      0
    ## 18      25      0      0      0      0      0      0      0     36      0
    ## 19       0      0      0      0      0      0      0      0      0      0
    ## 20      60      0      0      0      0     54      0      0      0      0
    ## 21     105      0      0      0     11     52      0      0      0      0
    ## 22       7     82      0      0      0      0      0      0      0      0
    ## 23       0      0      0      0      0      7      0      0     50      0
    ## 24     271      0      0      0     64     17      0      0      0      0
    ## 25      32      0      0      0      0      0      0      0      0      0
    ## 26       0      0      0    102      0      0      0      0     14      0
    ## 27       0      0      0      0      0      0      0      0      0      0
    ## 28       0      0      0      0     29      0      0     70      0      0
    ## 29       0      0      0      0      0      0      0      0     25      0
    ## 30       0     40      0      0      0     15      0      0      0      0
    ## 31       0      0      0      0      0      0      0      0      0      0
    ## 32       0      0      0      0      8      0      0      0      0      0
    ## 33       0      0      0      0     32      0      0      0      0      0
    ## 34     160      0      0     10      8      0      0      0      0      0
    ## 35       0      0      0      0      0      0      0     64      0      0
    ## 36       0      0      0      0      0      0      0     68      0      0
    ## 37       0      0      0      0      0      0      0      0     52      0
    ## 38      97      0      0      0      7      0      0      0     50      0
    ## 39       0      0      0      0      0      0      0      0      0      0
    ## 40       0      0      0      0    127      0      0      0     28      0
    ## 41      19      0      0      0    178      0      0      0     31      0
    ## 42       0      0      0      0    298      0      0     35      0      0
    ## 43       0      0      0      0     59      0      0      0      0      0
    ## 44       0      0      0      0     47      0      0      0     60      0
    ## 45       0      0      0      0      0      0      0      0      0      0
    ## 46       0      0      0      0     65      0      0      0      0      0
    ## 47      26      0      0      0    140      0      0      0      5      0
    ## 48       0      0      0      0     13      0      0      0     40      0
    ## 49       0      0      0      0      0      0      0      0      0      0
    ## 50       0      0      0      0      0      0      0      0     76      0
    ## 51       0      0      0      0     80      0      0      0      0      0
    ## 52       0      0      0      0     25      0      0     53      0      0
    ## 53       0      0      0      0      6      0      0      0     47      0
    ## 54       0      0      0      0      3      0      0      0      0      0
    ## 55       0    415      0    276      0      0      0      0     11      0
    ## 56       0      0      0      0      0      0      0      0      0      0
    ## 57       0      0      0      0      0     44      0      0      0      0
    ## 58       0      0      0      0    252      0      0      0     47      0
    ## 59       0      0      0      0      0      0      0      0     13      0
    ## 60     142      0      0      0     21      0      0      0      0      0
    ## 61       0      0      0      0      0      0      0      0      0      0
    ## 62       0      0     34    405      0      0      0      0     94      0
    ## 63       0      0      0      0      0      0      0      0      0      0
    ## 64       0      0      0      0      0      0      0      0      0      0
    ## 65       0      0      0      0      0      0      0      0      0      0
    ## 66       0      0      0      0      0      0      0      0      0      0
    ## 67       0      0      0      0      0      0      0      0      0      0
    ## 68       0      0      0      0      0      0      0      0      0      3
    ## 69       0      0      0      0      0      0      0      0      0      0
    ## 70       0      0      0      0      0      0      0      0      0      0
    ## 71       0      0      0      0      0      0      0      0      0     48
    ## 72       0      0      0      0      0      0      0      0      0      0
    ## 73       0      0      0      0     76      0      0      0      0      0
    ## 74       0      0      0      0      0      0      0      0     13      0
    ## 75       0      0      0      0      0      0      0      0     55      0
    ## 76       0      0      0      0      0      0      0      0     95      0
    ## 77       0      0      0      0      0      0      0      0      0      0
    ## 78       0      0      0      0      0      0      0      0      0      0
    ## 79       0      0      0      0      0      0      0      0      0      0
    ## 80       0      0      0      0     38      0      0      0      0      0
    ## 81       0      0      0      0      9      0      0      0     28      0
    ## 82       0      0      0      0      0      0      0      0      0      0
    ## 83       0      0      0      0      0      0      0      0      0      0
    ## 84       0      0      0      0      0      0      0      0      0      0
    ## 85       0      0      0    169      0      0      0      0      0      0
    ## 86       0      0      0      0      0      0      0      0      0      0
    ## 87       0      0      0      0      0      0      0      0      0      0
    ## 88       0      0      0      0      0      0      0      0      0      0
    ## 89       0      0      0      0      0      0      0      0      0     60
    ## 90       0      0      0      0      0      0      0      0      0      0
    ## 91       0      0      0      0      0      0      0      0      0      0
    ## 92       0      0      0      0      0      0      0      0      0      0
    ## 93      12      0      0      0      0      0      0      0      0      0
    ## 94       0      0      0      0      0      0      0      0      0      0
    ## 95       0      6      0      0     19      0      0      0      0      0
    ## 96       0      0      0    110      0      0      0      0      0      0
    ## 97       0      0      0      0      0      0      0      0      0      0
    ## 98       0      0      0      6      0      0      0      0      0      0
    ## 99       0      0      0     37      0      0      0      0      0      0
    ## 100      0      0     18      0      0      0      0     13      0      0
    ## 101      0      0      0      0      0      0      0      0      0      0
    ## 102      0      0      7      0      0      0      0      0      0      0
    ## 103      0      0      0      0      0      0      0      0      0      0
    ## 104      0      0      0      0      0      0      0      0      0      0
    ## 105      0      0      0      0      0      0      0      0      0      0
    ## 106      0      0      0      0      0      0      0      0      0      0
    ## 107      0      0     12      0      0      0      0      0      0      0
    ## 108      0      0      0      0      0      0      0      0      0      0
    ## 109      0      0      0      0      0      0      0      0      0      0
    ## 110      5      0      0      0      0      0      0      0     15      0
    ## 111      0      0      0      0      0      0      0      0      0      0
    ## 112      0      0      0      0      0      0      0      0      0      0
    ## 113      0      0      0      0      0      0      0      0      0      0
    ## 114      0      0      0      0      5      0      0      0      0      0
    ## 115      0      0      0      0      0      0      0      0      0      0
    ## 116      0      0      0      0      0      0      0      0     24      0
    ## 117      0      0      0      0      0      0      0      0      0      0
    ## 118      0      0      0      0      0      0      0      0      0      0
    ## 119      0      0      0      0      0      0      0      0      0      0
    ## 120      0      0      0      0     39      0      0      0      0      0
    ## 121      0      0      0      0     10     21      0      0      0      0
    ## 122      0      0      0      0      0      0      0      0      0      0
    ## 123      0      0      0      0      0      0      0      0      0      0
    ## 124      0      0      0      0      0      0      0      0      0      0
    ## 125      0      0      0      0     25      0      0      0      0      0
    ## 126      0      0      0      0      0      0      0      0      0      0
    ## 127      0      0     16      0      0      0      0      0      0      0
    ## 128     10      0      0      0      0      0      0     14      0     63
    ## 129      0      0     94      0      0      0      0      0      0      0
    ## 130      0    457      8      0      0      0      0      0      0     19
    ## 131      0      0     90      0      0      0      0      0      0      0
    ## 132      0      0     10      0      0      0      0      0      0      0
    ## 133      0      0     24      0      0      0      0      0      0      0
    ## 134      0      0      0      0      0      0      0      0      0      0
    ## 135      0      0     42      0      0      0      0      0      0      9
    ## 136      0      0     28      0      0      0      0      0      0      0
    ## 137      0      0     14      0      0      0      0      0      0      0
    ## 138      0      0      0      0      0      0      0      0      0     22
    ## 139      0      0      0      0      0      0      0      0      0      0
    ## 140      0      0    193      0      0      0      0      0      0     12
    ## 141      0      0     86      0      0      0      0      0      0      0
    ## 142      0      0      0      0      0      0      0      0      0      0
    ## 143     27     15      0      0      0      0      0      0      0     18
    ## 144      0      0     91      0      0      0      0      0      0      0
    ## 145      0      0      0      0      0      0      0      0      0      0
    ## 146      0      0      0      0      0      0      0      0      0      0
    ## 147      0      0     30      0      0      0      0      0      0      0
    ## 148      0      0      0      0      0      0      0      0      0      0
    ## 149      0      0      0      0      0      0      0      0      0      0
    ## 150      0      0      0      0      0      0      0      0      0      0
    ## 151      0      0      0      0      0      0      0      0      0      0
    ## 152      0      0      0      0      0      0      0      0      0      0
    ## 153      0      0      0      0      0      0      0      0      0      0
    ## 154      0      0      0      0      0      0      0      0      0      0
    ## 155      0      0     31      0      0      0      0      0      0      0
    ## 156      0      0      0      0      0      0      0      0      0      0
    ## 157      0      0      0     79     20      0      0      0      0      0
    ## 158      0    437    143      0      0      0      0      0      0      0
    ## 159      0      0      0      0      0      0      0      0      0      0
    ## 160      0      0     23      0      0      0      0      0      0      0
    ## 161      0      0      0      0     33      0      0      0      0     13
    ## 162      0      0      0      0      0      0      0      0      0      0
    ## 163      0      0      0      0      0      0     21      0      0      0
    ## 164      0      0      0      0      0     37     30      0      0      0
    ## 165      0      0      0      0      0      0      0      0      0      0
    ## 166      0      0      0      0      0      0      0      0      0      0
    ## 167      0      0      0      0      0     37     54      0      0      0
    ## 168      0      0      0      0      0      0      0      0      0      0
    ## 169      0      0     26      0      0      0     20      0      0      0
    ## 170      0      0      0      0      0      0    170      0      0      0
    ## 171      0      0      0      0      0     40    689      0      0      0
    ## 172      0      0     31    351      0      0      0      0      0      0
    ## 173      0      0     18      0      0    229      0      0      0      0
    ## 174      0    657    239      0      0    314      0      0      0      0
    ## 175      0      0     89      0      0     11      0      0      0      0
    ## 176      0      0     54      0      0    268      0      0      0      0
    ## 177      0      0     62    943      0     10      0      0      0      0
    ## 178      0      0     11    601      0    213      0      0      0      0
    ## 179      0      0     90    135      0      0      0      0      0      0
    ## 180      0     66     15    161      0      0      0      0      0      0
    ## 181      0      0      0      0      0      0      0      0      0      0
    ## 182      0    790      0      0      0      0      0      0      0     74
    ## 183      0      0     37      0      0      0      0      0      0    260
    ## 184      0    551      0      0      0    175     12      0      0     28
    ## 185      0      0      0      0      0      0     26      0      0     43
    ## 186      0      0      0      0      0     13      0      0      0     15
    ## 187      0      0      0      0      0      0      0      0      0     87
    ## 188      0      0     36      0      0      0      6      0      0     84
    ## 189      0      0      0      0      0      0      0      0      0      0
    ## 190      0      0      0      0      0      0      0      0      0      0
    ## 191      0      0      0      0      0      0      0     30      0      0
    ## 192      0      0      0      0      0      0      0      0      0      0
    ## 193      0      0      0      0      0      0      0      0      0      0
    ## 194      0      0      0      0      0      0     73      0      0      0
    ## 195     91      0      0      0      0      0     33    134      0     67
    ## 196      0      0      0      0      0      0      0      0      0     23
    ## 197     43      0     13      0      0      0      0      0      0     29
    ## 198      0      0      0      0      0      0      0      0      0      0
    ## 199      0      0      0      0      0      8      0      0      0      0
    ## 200      0      0      0      0      0    119      0      0      0      0
    ## 201      0      0      0      0      0      0      0      0      0      0
    ## 202      0      0      0      0      0      0      0      0      0      0
    ## 203      0      0      0      0      0      0      0      0      0      0
    ## 204      0      0      0      0      0      0      0      0      0      0
    ## 205      0      0      0      0      0      0      0      0      0      0
    ## 206      0      0      0      0      0      0      0      0      0      0
    ## 207      0      0      0      0      0      0    160      0      0      0
    ## 208      7      0      0      0      0      0      0      0      0      0
    ## 209      0      0      6      0      0     61      0      0      0      0
    ## 210      8      0      0      0     16    617      0      6      0      0
    ## 211     37      0      0      0      9      0      0    372      0     93
    ## 212      0      0      0      0      0      0      0      0      0     47
    ## 213      6      0      0      0    198      0      0     56      0      0
    ## 214    149      0      0      0      0      0      0    270      0     33
    ## 215     24      0      0     29      0      0      0     34      0     39
    ## 216     11      0      0      0      0      0      0      0      0      0
    ## 217      7      0      0      0      0    176      0      0      0      0
    ## 218      0      0      0      0      0      0      0      0      0      0
    ## 219     10      0      0      0      0     51      0      0      0      0
    ## 220      0      0      0      0      8    450      0      0      0      0
    ## 221      0      0      0      0      0      0      0     61    121      0
    ## 222      0      0      0      0      0      0      0      0      0      0
    ## 223     58      0      0      0      0     48      0     13      0      0
    ## 224      0      0      0      0      0      9      0      0      0      0
    ## 225      0      0      0      0      0      0      0      0      0      0
    ##     otu_60 otu_61 otu_62 otu_63 otu_64 otu_65 otu_66 otu_67 otu_68 otu_69
    ## 1        0     67      0      0      0      0      0      0     16     18
    ## 2        0     14      0      0      0      0      0      0     61     10
    ## 3        0      0      0      0      0      0      0      0     10    101
    ## 4        0     97      0      0      0      0      0      0      3     64
    ## 5        0      0      0      0      0      0     14      0     20    220
    ## 6        0      0      0      0      0      0     12      0      0      0
    ## 7        0      0      0      0      0      0      0      0     38      9
    ## 8        0      0      0     93      0      0      0      0      0     34
    ## 9        0      0      0      0      0      0      0      0      0     10
    ## 10       0      0      0      0      0      0      0      0      0      0
    ## 11       0      0      0      0      0      0      0      0    113      0
    ## 12       0     40      0      0     13      0      0      0     22      0
    ## 13       0     57      0      0     18      0      0      0      0      0
    ## 14      35      0      0      0      0      0      0      0      0      0
    ## 15       0     87      0    272      0      0      0      0     30      0
    ## 16       0      0      0      0      0      0      0      0     80      0
    ## 17       0      0      0      0      0      0      0      0      9      0
    ## 18       0      0      0      0      0      0      0      0      0      0
    ## 19       0      0      0      0      0      0      0      0      0      0
    ## 20       0      0      0      0      0      0      0      0      0      0
    ## 21       0      0      0      0      0      0      0      0     17      0
    ## 22       0      0      0      0      0      0      0      0      0      0
    ## 23       0      0      0      0     17      0      0      0      0      0
    ## 24       0      0      0      0      0      0      0      0      0      0
    ## 25       0      0      0      0      0      0      0      0      0      0
    ## 26       0      0      0      0      0      0      0      0     10      0
    ## 27       0      0      0      0      0      0      0      0     10      0
    ## 28       0      0      0      0      0      0      0      0      0      0
    ## 29       0      0      0      0      0      0      0      0      0      0
    ## 30       0      0      0      0      0      0      0      0      0      0
    ## 31       0      0      0      0      0      9      0      0      0      0
    ## 32       0      0      0      0      0      0      0      0      0      0
    ## 33       0      0      0     80      0      0      0      0      0      0
    ## 34       0      0      0      0      0     35      0      0      0      0
    ## 35       0      0      0      0      0     75      0      0      0      0
    ## 36       0      0      0      0      0     22      0      0      0      0
    ## 37       0     11      0      0      0      0      0      0      0      0
    ## 38       0      0      0      0      0      0      0      0      0     16
    ## 39       0      0      0      0      0      0      0      0      0      0
    ## 40       0      0      0      0      0      0      0      0      0      0
    ## 41       0      0      0     37      0      0      0      0      0     41
    ## 42       0      0      0      0      0      0     22      0      0     37
    ## 43       0      6      0    170      0      0      0      0      0      8
    ## 44       0     19      0      0      0      0      0      0      0      8
    ## 45       0      0      0      0      0      0      0      0      0      0
    ## 46       0      0      0      0      0      0      0      0      0      0
    ## 47       0      0      0      0      0      0      0      0      0      9
    ## 48       0      0      0    169      0      0      0      0      0      0
    ## 49       0      0      0      0      0      0      0      0      0      0
    ## 50       0      0      0      0      0      0      0      0      0      0
    ## 51       0      0      0      0      0      0      0      0      0      0
    ## 52       0      0      0    190      0      0      0      0      0      0
    ## 53       0     71      0      0      0      0      0      0      0      0
    ## 54       0      0      0      0      0      0      0      0      0      0
    ## 55       0      0      0      0      0      0      0      0      0     85
    ## 56       0      0      0      0      0      0      0      0      0      0
    ## 57       0      0      0      0      0      0      0      0      0      0
    ## 58       0      0      0      0      0      0      0      0      0      0
    ## 59       0      0      0      0      0      0      0      0      0      0
    ## 60       0      0      0      0      0      0      0      0     34      0
    ## 61       0      0      0      0      0      0      0      0     48      0
    ## 62       0      0      0      0      0      0      0      0      0      0
    ## 63       0      0      0      0      0      0      0      0      0      0
    ## 64       0      0      0      0      0      0      0      0      0      0
    ## 65       0      0      0      0      0      0      0      0      0      0
    ## 66       0      0      0      0      0      0      0      0      0      0
    ## 67       0      0      0      0      0      0    199      0      0      0
    ## 68       0      0      0      0      0      0      0      0      0      0
    ## 69       0      0      0      0      0      0      0      0      0      0
    ## 70       0      0      0      0      0      0      0      0      0      0
    ## 71       0      0      0      0      0      0      0      0      0      0
    ## 72       0      0      0      0      0      0      0      0      0      0
    ## 73       0      0      0      0      0      0      0      0      0      0
    ## 74       0      0      0      0      0      0      0      0      0      0
    ## 75       0      0     26      0      0      0      0      0      0      0
    ## 76       0     56      0      0      0      0      0      0      0      0
    ## 77       0    244      0      0      0      0      0      0      0      0
    ## 78       0      0      0      0      0      0      0      0      0      0
    ## 79       0      0      0      0      0      0      0      0      0      0
    ## 80       0     10      0      0      0      0      0      0      0      0
    ## 81       0      0      0      0      0      0      0      0      0      0
    ## 82       0      0      0      0      0      0      0      0      0      0
    ## 83       0      0      0      0      0      0      0      0      0      0
    ## 84       0      0      0      0      0      0     11      0      0      0
    ## 85       0      0      0      0      0      0      0      0      0      0
    ## 86      58      0      0      0      0      0      0      0      0      0
    ## 87       0      0      0      0      0      0      0      0      0      0
    ## 88       0      0      0      0      0      0      0      0      0      0
    ## 89       0      0      0      0      0      0      0      0      0      0
    ## 90       0      0      0      0      0      0      0      0      0      0
    ## 91       0      0      0      0      0      0      0      0      0      0
    ## 92       0    241      0      0      0      0     20      0      0      0
    ## 93       0      0      0      0      0      0      0      0      0      0
    ## 94       0      0      0      0      0      0      0      0      0      0
    ## 95       0      0      0      0      0     85     67      0      0      0
    ## 96       0      0      0      0      0    161     24      0      0      0
    ## 97       0      0      0      0      0    201      0      0      0      0
    ## 98       0      0      0      0      0     39     11      0      0      0
    ## 99       0      0      0     84      0     64      0      0      0      0
    ## 100      0      0      0      0      0      0      0      0      0      0
    ## 101      0      0      0      0      0      0      0      0      0      0
    ## 102      0      0      0      0      0      0      0      0      0      0
    ## 103      0      0      0      0      0      0      0      0      0      0
    ## 104      0      0      0      0      0      0      0      0      0      0
    ## 105      0      0      0      0      0      0      0      0      0      0
    ## 106     15      0      0      0      0      0      0      0      0      0
    ## 107      0      0      0      0      0      0      0      0      0      0
    ## 108      0      0      0      0      0      0      0      0      0      0
    ## 109      0      0    316      0      0      0      0      0      0      0
    ## 110      0    265    566      0      0      0      0      0      0      0
    ## 111      0      0      0      0      0      0      0      0     14      0
    ## 112      0      0      0      0      0      0      0      0     14      0
    ## 113      0      0     34      0      0      0      0      0      0      0
    ## 114      0      0      0      0      0      0      0      0      0      0
    ## 115      0      0      0      0      0      0      0      0      0      0
    ## 116      0      0     12      0      0      0      0      0      0      0
    ## 117      0      0      0      0      0      0      0      0      0      0
    ## 118      0      0      0      0      0      0      0      0      0      0
    ## 119      0    105      0      0      0      0      0      0      0      0
    ## 120      0     14      3      0      0      0      0      0      0      0
    ## 121      0      0      0     68     23      0      0      0      0      0
    ## 122      0      0      0      0      0      0      0      0      0      0
    ## 123      0      0      0     16     52      0      0      0      0      0
    ## 124      0     56      0      0      0      0      0      0      0      0
    ## 125      0      0     11      0      0      0      0      0      0      0
    ## 126      0      0      0      0      0      0      0      0      0      0
    ## 127      0      0      0      0      0      0      7      0      0      0
    ## 128     34      0      0      0      0      0      0      0      0      0
    ## 129      0      0      0      0      0      0      0      0      0      0
    ## 130      0      0      0      0      0      0      0      0      0      0
    ## 131      0      0      0      0      0      0      0      0      0      0
    ## 132      0      0      0      0      0      0      0      0      0      0
    ## 133      0      0      0      0      0      0      0      0      0      0
    ## 134      0      0      0      0      0      0      0      0      0      0
    ## 135      0      0      0      0      0      0      0      0      0      0
    ## 136      0      0      0      0      0      0      0      0      0      0
    ## 137      0      0      0      0      0      0      0      0      0      0
    ## 138      0      0      0      0      0      0      0      0      0      0
    ## 139      0      0      0      0      0      0      0      0      0      0
    ## 140      0      0      0      0      0      0      0      0      0      0
    ## 141      0      0      0      0      0      0      0      0      0      0
    ## 142      0     24      0      0      0      0      0      0      0      0
    ## 143      0      0      0      0      0      0      0      0      0      0
    ## 144      0      0      0      0      0      0      0      0      0      0
    ## 145      0      0      0      0      0      0      0      0      0      0
    ## 146      0      0      0      0      0      0      0      0      0      0
    ## 147      0      0      0      0      0      0      0      0      0      0
    ## 148      0      0      0      0      0      0      0      0      0      0
    ## 149      0      0      0      0      0      0      0      0      0      0
    ## 150      0      0      0      0      0      0      0      0      0      0
    ## 151      0      0      0      0      0      0      0      0      0      0
    ## 152      0      0      0      0      0      0      0      0      0      0
    ## 153      0      0      0      0      0      0      0      0      0      0
    ## 154      0      0      0      0      0      0      0      0      0      0
    ## 155      0      0      0      0      0      0      0      0      0      0
    ## 156      0      0      0      0      0      0      0      0      0      0
    ## 157      0      0      0      0      0      0      0      0      0      0
    ## 158      0      0      0      0      0      0      0      0      0      0
    ## 159      0      0      0      0      0      0      0      0      0      0
    ## 160      0      0      0      0      0      0      0      0      0      0
    ## 161      0      0      0      0      0      0      0      0      0      0
    ## 162     71     36      0      0      0      0      0      0      0      0
    ## 163      0      0      0      0      0      0      0      0      0      0
    ## 164      0      0      0      0      0      0      0      0      0      0
    ## 165      0      0      0      0      0      0      0      0      0      0
    ## 166      0      0      0      0     17      0      0      0      0      0
    ## 167      0      0      0      0      0      0      0      0      0      0
    ## 168      0      0      0      0      0      0      0      0      0      0
    ## 169     27      0      0      0      0      0      0      0      0      0
    ## 170      0      0      0      0      0      0      0      0      0      0
    ## 171      0      0      0      0      0      0      0      0      0      0
    ## 172      0      0      0      0      0      0      0      0      0      0
    ## 173      0      0      0      0     48      0      0      0      0      0
    ## 174      0      0      0      0     76      0      0      0      0      0
    ## 175      0      0      0      0      0      0      0      0      0      0
    ## 176    157      0      0      0     50      0      0      0     67      0
    ## 177      0      0      0      0      0      0      0      0      0      0
    ## 178      0     20      0      0      0      0      0      0      0      0
    ## 179      0      0      0      0      0      0      0      0      0      0
    ## 180      0      0      0      0      0      0      0      0      0      0
    ## 181     44      0      0      0      0      0      0      0      0      0
    ## 182     74      0      0      0      0      0      0      0      0      0
    ## 183     39      0      0      0      0      0      0      0      0      0
    ## 184      9      0      0      0     30      0      0      0      0      0
    ## 185     70      0      0      0      0      0      0      0      0      0
    ## 186     65      0      0      0      0      0      0      0      0      0
    ## 187      8      0      0      0      0      0      0      0      0      0
    ## 188     65      0      0      0      0      0      0      0      0      0
    ## 189      0      0      0      0      0      0      0      0      0      0
    ## 190      0      0      0      0      0      0      0      0      0      0
    ## 191      0      0      0      0      0     17   1108      0      0      0
    ## 192      0      0      0     24      0     29     78      0      0      0
    ## 193      0      0      0      0      0      0      0      0      0      0
    ## 194      0      0      0      0      0      0     86      0      0      0
    ## 195     21      0      0      0      0      0     69      0      0      0
    ## 196      0      0      0      0      0      0      0      0      0      0
    ## 197     22      0      0      0      0      0     43      0      0      0
    ## 198      0      0      0      0      0      0     34    112      0      0
    ## 199      0      0      0      0      0      0      0      0      0      0
    ## 200      0      0      0      0      0      0      0    422      0      0
    ## 201      0      0      0      0      0      0      0      0      0      0
    ## 202      0      0      0      0      0      0      0     14      0      0
    ## 203      0      0      0      0      0      0      0    147      0      0
    ## 204      0      0      0      0      0      0      0      0      0      0
    ## 205      0      0      0      0      0      0      0      0      0      0
    ## 206      0      0      0      0      0      0      0      0      0      0
    ## 207      0      0      0      0      0      0      0      0      0      0
    ## 208      0      0      0      0      0      0      0      0      0      0
    ## 209      0      0      0      0     26      0      0      0      0      0
    ## 210      0      0      0      0    155      0      0      0      0      0
    ## 211      0      0      0      0      0      0      0      0      0      0
    ## 212      0      0      0      0      0      0      0      0      0      0
    ## 213      0      0      0      0      0      0      0      0     77      0
    ## 214      0      0      0      0      0      0      0      0      0      0
    ## 215      0      0      0      0      0      0      0      0      0      0
    ## 216      0      0      0      0      0      0      0      0      0      0
    ## 217      0      0      0      0    262      0      0      0      0      0
    ## 218      0      0      0      0      0      0      0      0      0      0
    ## 219     57     14      0      0     68      0      0      0      0      0
    ## 220      0      0      0      0    337      0      0      0      0      0
    ## 221      0      0      0      0      0      0      0      0     63      0
    ## 222    216      0      0      0      0      0      0      0      8      0
    ## 223      0      0      0      0     85      0      0      0      0      0
    ## 224      0      0      0      0      0      0      0      0      0      0
    ## 225      0      0      0      0      0      0      0      0      0      0
    ##     otu_70 otu_71 otu_72 otu_73 otu_74 otu_75 otu_76 otu_77 otu_78 otu_79
    ## 1        0      0     30      0     14      0     17      0      0      0
    ## 2        0      0      0      0      9      0      0      0      0      0
    ## 3        0      0      0      0      0      0      0      0      0      0
    ## 4        0      0      0      0     42      0     26      0      0      0
    ## 5        0      0      0      0      0      0      0      0      0      0
    ## 6        0      0      0      0      0      0      0      0      0      0
    ## 7        0      0      0      0      0      0      9      0      0      0
    ## 8      163      0      0      0      0      0      0      0     32      0
    ## 9       10      0      0      0      0      0      0      0      0      0
    ## 10       0      0      0      0      0      0      0      0      0      0
    ## 11       0      0      0      0      0      0      0      0      0      0
    ## 12       9      0      0      0      0     42     12      0      0      0
    ## 13     119      0      0      0     27    157     13      0      0     63
    ## 14       0      0      0      0      0      0      0      0      0      0
    ## 15      32      0     48      0     59      0     21      0      0      0
    ## 16       0      0      0      0      0      0      0      0      0      0
    ## 17      21      0     44      0      0     13      0      0     18      0
    ## 18      48     10     10      0      0     66      0      0      0      0
    ## 19       0      0      0      0      0      0      0      0      0      0
    ## 20      46      0      0      0      0      0      0      0      0      0
    ## 21       0      0    148      0      0      3      0      0      0      0
    ## 22       9      0     68      0      0      0      0      0      0    155
    ## 23      30      0     37      0      0      0      0      0      0      0
    ## 24       0      0      0      0      0      0      0      0     15      0
    ## 25       0      0    108      0      0      0      0      0      6      0
    ## 26       0      0     23      0      0      0      0      0      0      0
    ## 27       0      0      7      0      0      0      0      0      0      0
    ## 28       0      8      0      0      0      0      0      0     17      0
    ## 29       0      0      0      0      0      0      0      0      0      0
    ## 30       0    232      0      0      0      0      0      0      0      0
    ## 31       0     28      0      0      0      0      0      0      0      0
    ## 32       0    277      0      0      0      0      0      0      0      0
    ## 33       0     53      0      0      0      0      0      0      0      0
    ## 34       0      0      0      0      0      0      0      0      0      0
    ## 35       0    198      0      0      0      0      0      0      0      0
    ## 36       0      0      0      0      0      0      0      0      0      0
    ## 37       0      0      0      0      0      0      0      0      0      0
    ## 38       0      0      0      0      0      0      0      0      0      0
    ## 39       0      0      0      0      0      0      0      0      0      0
    ## 40       0      0      0      0      0      0      0      0     35      0
    ## 41       0      0      0      0      0      0      0      0     29      0
    ## 42       0      0      0      0      0      0      0      0     60      0
    ## 43       0      0      0      0      0      0      0      0     15      0
    ## 44       0      0      0      0      0      0     14      0     18      0
    ## 45       0      0      0      0      0      0      0      0      0      0
    ## 46       0      0      0      0      0      0      0      0      8      0
    ## 47       0     46      0      0      0      0      0      0     20      0
    ## 48       0      0      0      0      0      0      0      0      0      0
    ## 49       0      0      0      0      0      0      0      0      0      0
    ## 50       0      0      0      0      0      0      0      0      0      0
    ## 51       0      0     13      0      0      0      0      0      8      0
    ## 52       0    101      0      0      0      0      0      0      0      0
    ## 53       0      0      0      0     27      0     29      0      0      0
    ## 54       0    127    188      0      0      0      0      0      0      0
    ## 55       0      0      0      0      0      0      0      0      0      0
    ## 56       0      0      0      0      0      0      0      0      0      0
    ## 57       3      0     36      0      0      0      0      0      0      0
    ## 58       0      0      0      0      0      0      0      0     89      0
    ## 59       0      0      0      0      0      0      0      0      0      0
    ## 60       0      0      0      0      0      0      0      0      0      0
    ## 61       0      0      0      0      0      0      0      0      0      0
    ## 62      88      0      0      0      0      0      0      0      0      0
    ## 63      81      0     16      0      0      0      0      0      0      0
    ## 64       0      0      0      0      0      0      0      0      0      0
    ## 65       0      0      0      0      0      0      0      0      0      0
    ## 66      19      0      0      0      0      0      0      0      0      0
    ## 67       0      0      0      0      0      0      0      0      0      0
    ## 68       0      0      0      0      0      0      0      0      0      0
    ## 69       0      0      0      0      0      0      0      0      0      0
    ## 70       0      0      0      0      0      0      0      0      0      0
    ## 71       0      0      0      0      0      0      0      0      0      0
    ## 72       0      0      0      0      0      0      0      0      0      0
    ## 73       0      0      0      0      0      0      0      0     38      0
    ## 74       0      0      0      0      0      0      0      0      0     15
    ## 75       0      0      0      0      0      0      0      0      2      0
    ## 76       0      0      0      0     26      0     17      0      0      0
    ## 77       0      0      0      0     98      0     86      0      0      0
    ## 78       0      0     20      0      0      0      0      0      0      0
    ## 79       0      0      0      0      0      0      0      0      0      0
    ## 80       0      0      0      0     16      0      6      0      0      0
    ## 81       0      0      0      0      0      0      0      0      3      0
    ## 82       0      0      0      0      0      0      0      0      0      0
    ## 83       0      0     28      0      0      0      0      0      0      0
    ## 84       0      0      0      0      0      0      0      0      0      0
    ## 85       0      0      0      0      0      0      0      0      0      0
    ## 86       0      0      0      0      0      0      0      0      0      0
    ## 87       0      0      0      0      0      0      0      0      0      0
    ## 88       0      0      0      0      0      0      0      0      0      0
    ## 89       0      0      0      0      0      0      0      0      0      0
    ## 90       0      0      0      0      0      0      0      0      0      0
    ## 91       0      0      0      0      0      0      0      0      0      0
    ## 92       0      0      0      0     83      0     45      0      0      0
    ## 93       0      0      0      0      0      0      0      0      0      0
    ## 94       0      0      0      0      0      0      0      0      0      0
    ## 95       0      4      0      0      0      0      0      0      0      0
    ## 96       0      0      0      0      0      0      0      0      0      0
    ## 97       0      0      0      0      0      0      0      0      0      0
    ## 98       0      0      0      0      0      0      0      0      0      0
    ## 99       0      0      0      0      0      0      0      0      0      0
    ## 100      0      0      0      0      0      0      0      0      8      0
    ## 101      0      0      0      0      0      0      0      0      0      0
    ## 102      0      0      0      0      0      0      0      0      0      0
    ## 103      0      0      0      0      0      0      0      0      0      0
    ## 104      0      0      0      0      0      0      0      0      0      0
    ## 105      0      0      0      0      0      0      0      0      0      0
    ## 106      0      0      0      0      0      0      0      0      0      0
    ## 107      0      0      0      0      0      0      0      0      0      0
    ## 108      0      0      0      0      0      0      0      0      0      0
    ## 109      0      0      0      0      0      0      0    265      0      0
    ## 110    364      0      0      0     95      0    101     50      0      0
    ## 111      0      0      0      0      0      0      0     80      0      0
    ## 112      0      0      0      0      0      0      0     14      0      0
    ## 113      0      0      0      0      0      0      0      7      0      0
    ## 114      0      0      0      0      0      0      0      0      0      0
    ## 115      0      0      0      0      0      0      0      0      0      0
    ## 116      0      0      0      0      0      0      0      0      0      0
    ## 117      0      0      0      0      0      0      0      0     11      0
    ## 118      0      0      0      0      0      0      0      0      0      0
    ## 119      0      0      0      0     25      0     31      0      0      0
    ## 120      0      0      0      0      0      0      0      0     13      0
    ## 121      0      0      0      0      0      0      0      0      9      0
    ## 122      0      0      0      0      0      0      0      0      0      0
    ## 123      0      0      0      0      0      0      0      0      0      0
    ## 124      0      0      0      0     43      0     22      0      0      0
    ## 125      0      0      0      0      0      0      0      0      0      0
    ## 126      0      0      0      0      0      0      0      0      0      0
    ## 127      0      0      0      0      0      0      0      0      0      0
    ## 128      0      0     75      0      0     12      0      0      0      0
    ## 129      0      0      0      0      0      0      0      0      0      0
    ## 130      0      0      0      0      0      0      0      0      0      0
    ## 131      0      0      0      0      0      0      0      0      0      0
    ## 132      0      0      0      0      0      0      0      0      0      0
    ## 133      0      0      0      0      0      0      0      0      0      0
    ## 134      0      0      0      0      0      0      0      0      0      0
    ## 135      0      0      0      0      0      0      0      0      0      0
    ## 136      0      0      0      0      0      0      0      0      0      0
    ## 137      0      0      0      0      0      0      0      0      0      0
    ## 138      0      0      0      0      0      0      0      0      0      0
    ## 139      0      0      0      0      0      0      0      0      0      0
    ## 140      0      0      0      0      0      0      0      0      0      0
    ## 141      0      0      0      0      0      0      0      0      0      0
    ## 142      0      0      0      0     18      0      0      0      0      0
    ## 143      0      0      0      0      0      0      0      0      0      0
    ## 144      0      0      0      0      0      0      0      0      0      0
    ## 145      0      0      0      0      0      0      0      0      0      0
    ## 146      0      0      0      0      0      0      0      0      0      0
    ## 147      0      0      0      0      0      0      0      0      0      0
    ## 148      0      0      0      0      0      0      0      0      0      0
    ## 149      0      0      0      0      0      0      0      0      0      0
    ## 150      0      0      0      0      0      0      0      0      0      0
    ## 151      0      0      0      0      0      0      0      0      0      0
    ## 152      0      0      0      8      0      0      0      0      0      0
    ## 153      0      0      0      0      0      0      0      0      0      0
    ## 154      0      0      0      0      0      9      0      0      0      0
    ## 155      0      0      0      0      0      0      0      0      0      0
    ## 156      0      0      0      0      0      0      0      0     13      0
    ## 157      0      0      0      0      0      0      0      0      0      0
    ## 158      0      0      0      0      0      0      0      0      0      0
    ## 159      0      0      0      0      0      0      0      0      0      0
    ## 160      0      0      0      0      0      0      0      0      0      0
    ## 161      0      0      0      0      0      0      0      0     13     94
    ## 162      0      0      0      0     18      0      0      0      0     14
    ## 163      0      0      0     85      0      0      0      0      0      0
    ## 164      0      0      0     12      0      0      0      0      0      0
    ## 165      0      0      0     30      0      0      0      0      0      0
    ## 166      0      0      0      0      0      0      0      0      0      0
    ## 167      0      0      0    141      0      0      0      0      0      0
    ## 168      0      0      0      0      0      0      0      0      0      0
    ## 169      0      0      0    298      0      0      0      0      0      0
    ## 170      0      0      0     45      0      0      0      0      0      0
    ## 171      0      0      0      6      0      0      0      0      0      0
    ## 172      0     16      0      0      0     77      0      0      0      0
    ## 173      0      0      0      0      0      0      0      0      0      0
    ## 174      0      0      0      0      0      0      0      0      0      0
    ## 175      0      0      0      0      0      0      0      0      0      0
    ## 176      0      0      0      0      0      0      0      0      0      0
    ## 177      0      0      0      0      0      0      0      0      0      0
    ## 178      0      0      0      0      0      0     15      0      0     21
    ## 179      0      0      0      0      0      0      0      0      0     74
    ## 180      0     12      0      0      0      0      0      0      0      0
    ## 181      0      0      0      0      0      0      0      0      0      0
    ## 182      0      0      0      0      0      0      0      0      0      0
    ## 183      0      0      0      0      0      0      0      0      0      0
    ## 184      0      0      0      0      0      0      0      0      0      0
    ## 185      0      0      0      0      0      0      0      0      0      0
    ## 186      0      0      0      0      0      0      0      0      0      0
    ## 187      0      0      0     76      0      0      0      0      0      0
    ## 188      0      0      0      0      0      0      0      0      0      0
    ## 189      0      0      0     18      0      0      0      0      0      0
    ## 190      0      0      0      0      0      0      0      0      0     45
    ## 191      0      0      0      0      0      0      0      0      0      0
    ## 192      0      0      0      0      0      0      0      0      0      0
    ## 193      0      0      0      0      0      0      0      0      0      0
    ## 194      0      0      0      0      0      0      0      0      0      0
    ## 195      0      0      0      0      0      0      0      0      0      0
    ## 196      0      0      0      0      0      0      0      0      0      0
    ## 197      0      0      0      0      0      0      0      0      0      0
    ## 198      0      0      0      0      0      0      0      0      0      0
    ## 199      0      0      0      0      0      0      0      0      0      0
    ## 200      0      0      0      0      0      0      0      0      0      0
    ## 201      0      0      0      0      0      0      0      0     21      0
    ## 202      0      0      0      0      0      0      0      0      0      0
    ## 203      0      0      0      0      0      0      0      0      0      0
    ## 204      0      0      0      0      0      0      0      0      0      0
    ## 205      0      0      0      0      0      0      0      0      0      0
    ## 206      0      0      0      0      0      0      0      0      0      0
    ## 207      0      0      0      0      0      0      0      0      0      0
    ## 208      0      0      0      0      0      0      0      0      0      0
    ## 209    144      0      0      0      0      0      0      0      0      0
    ## 210      0      0      0      0      0      0      0      0      0      0
    ## 211      0      0      0      0      0      0      0      0      5      0
    ## 212      0      0      0      0      0      0      0      0      0      0
    ## 213    122      0      0      0      0      0      0      0     26      0
    ## 214     34      0      0      0      0      0      0      0      0      0
    ## 215      0      0      0      0      0      0      0      0      0      0
    ## 216      0      0      0      0      0      0      0      0      0      0
    ## 217     12      0      0      0      0    153      0      0      0      0
    ## 218      0      0      0      0      0      0      0      0      0     67
    ## 219    243      0     11      0      0    132     12      0      0      0
    ## 220      0      0      0      0      0      0      0      0      0      0
    ## 221      0      0     32      0      0    186      0      0      0      0
    ## 222     41      0     16      0      0      0      0      0      0      0
    ## 223      0      0      0      0      0    155      0      0      0      0
    ## 224     43      0      9      0      0      0      0      0      0      0
    ## 225      0      0      0      0      0      0      0      0      0      0
    ##     otu_80 otu_81 otu_82 otu_83 otu_84 otu_85 otu_86 otu_87 otu_88 otu_89
    ## 1        0      0     22      0      0      0      0      0      0      0
    ## 2        0      0      0      0      0     45      0      0      0      0
    ## 3        0      0      0      0      0      0      0      0      0      0
    ## 4        0      0      0      0      0      0      0      0      0      0
    ## 5        0      0     24      0      0     49      0      0      0      0
    ## 6        0      0      0      0      0      0      0      0      0      0
    ## 7        0      0      0      0      0     43      0      0      0      0
    ## 8        0      0      0      0      0     17      0      0      0      0
    ## 9        0      0      0      0      0      9      0      0      0      0
    ## 10       0     36      0      0      0      0      0      0      0      0
    ## 11       0      0      0      0      0      0      0      0      0      0
    ## 12       0      0      0      0      0      0      0      0      0      0
    ## 13       0      0      0      0      0      0      0      0      0     54
    ## 14       0      0      0      0      0      0      0      0      0      0
    ## 15       0      0      0      0      0      0      0      0      0      0
    ## 16       0      0      0      0      0      0      0      0      0      0
    ## 17       0      0      0      0      0      0      0      0      0      0
    ## 18       0      0      0      0      0      0      0      0      0      0
    ## 19       0      0      3      0      0      0      0      0      0      0
    ## 20       0     11      0      0      0      0      0      0      0      0
    ## 21       0      6      0      0      0      0      0      0      0      0
    ## 22       0      0      0      0      0      0      0      0      0    179
    ## 23       0      0      0      0      0      0      0      0      0      0
    ## 24       0      0      0      0      0      0      0      0      0      0
    ## 25       0      0      0      0      0      0      0      0      0      0
    ## 26       0      0      0      0      0      0      0      0      0      0
    ## 27       0      0      0      0      0      0      0      0      0      0
    ## 28       0      0      7      0      0      0      0      0      0      0
    ## 29       0      0      0      0      0      0      0      0      0      0
    ## 30       0      0      0      0      0      0      0      0      0      0
    ## 31       0      0      0      0      0      0      0      0      0      0
    ## 32       0      0      0      0     38      0      0      0      0      0
    ## 33       0      0      0      0      0      0      0      0      0      0
    ## 34       0      0      0      0      0      0      0      0      0      0
    ## 35       6      0      0      0      0      0      0     25      0      0
    ## 36       0      0      0      0      0      0      0      9      0      0
    ## 37       0      0      0      0      0      0      0      0      0      0
    ## 38       0      0      0      0      0      0      0      0      0      0
    ## 39       0      0      0      0      0      0      0      0      0      0
    ## 40       0      0      0      0      0      0      0      0      0      0
    ## 41       0      0      0      0      0      0      0      0      0      0
    ## 42       0      0      0      0      0      0      0      0      0      0
    ## 43       0      0      0      0      0      0      0      0      0      0
    ## 44       0      0      0      0      0      0      0      0      0      0
    ## 45       0      0      0      0      0      0      0      0      0      0
    ## 46       0      0      0      0      0      0      0      0      0      0
    ## 47       0      0      0      0      0      0      0      0      0      0
    ## 48       0      0      0      0      0      0      0      0      0      0
    ## 49       0      0      0      0      0      0      0      8      0      0
    ## 50       0      0     29      0      0      0      0      0      0      0
    ## 51       0      0      0      0      0      0      0      0      0      0
    ## 52       0      0      0      0      0      0      0      0      0      0
    ## 53       0      0      0      0      0      0      0      0      0      0
    ## 54       0      0      0      0      0      0      0      0      0      0
    ## 55       0      0      0      0      0      0      0      0      0      0
    ## 56       0      0      0      0      0      0      0      0      0      0
    ## 57       0     15      0     50      0      0      0      0      0      0
    ## 58       0      0      0      0      0      0      0     16      0      0
    ## 59       0      0      0      0      0      0      0      0      0      0
    ## 60       0      0      0      0      0      0      0      0      0      0
    ## 61       0      0      0      0      0      0      0      0      0      0
    ## 62       0      0      0     46      0      0      0      0      0      0
    ## 63       0      0      0    206      0      0      0      0      0      0
    ## 64       0      0      0      0      0      0      0      0      0      0
    ## 65       0      0      0      0      0      0      0      0      0      0
    ## 66       0      0      0      0      0      0      0      0      0      0
    ## 67       0      0      0      0      0      0      0      0      0      0
    ## 68       0      0      0      0      0      0      0      0      0      0
    ## 69       0      0      0      0      0      0      0      0      0      0
    ## 70       0      0      0      0      0      0      0      0      0      0
    ## 71       0      0      0      0      0      0      0      0      0      0
    ## 72       0      0      0      0      0      0      0      0      0      0
    ## 73       0      0      0      0      0      0      0      0      0      0
    ## 74       0      0      0      0      0      0      0      0      0      0
    ## 75       0      0      0      0      0      0      0      0      0      0
    ## 76       0      0      0      0      0      0      0      0      0      0
    ## 77       0      0      0      0      0      0      0      0      0      0
    ## 78       0      0      0      0      0      0      0      0      0      0
    ## 79       0      0      0      0      0      0      0      0      0      0
    ## 80       0      0      0      0      0      0      0      0      0      0
    ## 81       0      0      0      0      0      0      0      0      0      0
    ## 82       0      0      0      0     31      0      0      0      0      0
    ## 83       0      0      0      0      0      0      0      0      0      0
    ## 84       0      0      0      0      0      0      0      0      0      0
    ## 85       0      0      0      0      0      0      0      0      0      0
    ## 86       0      0      0      0      0      0      0      0      0      0
    ## 87       0      0      0      0      0      0      0      0      0      0
    ## 88       0      0      0      0      0      0      0      0      0      0
    ## 89       0      0      0      0      0      0      0      0      0      0
    ## 90       0      0      0      0      0      0      0      0      0      0
    ## 91       0      0      0      0      0      0      0      0      0      0
    ## 92       0      0      0      0      0      0      0      0      0      0
    ## 93       0      0      0      0      0      0      0     24      0      0
    ## 94       0      0      0      0      0      0      0      0      0      0
    ## 95       0      0      0      0      0      0      0     80      0      0
    ## 96       0      0      0      0      0      0      0      0      0      0
    ## 97       0      0      0      0      0      0      0      0      0      0
    ## 98       0      0      0      0      0      0      0      0      0      0
    ## 99       0      0      0      0      0      0      0    129      0      0
    ## 100      0      0    159      0      0      0      0      0      0      0
    ## 101      0      0     48      0      0      0      0      0      0      0
    ## 102      0      0      0      0      0      0      0      0      0      0
    ## 103      0      0     61      0     52      0      0      0      0      0
    ## 104      0      0      0      0     35      0      0      0      0      0
    ## 105    183      0      0      0      0      0      0      0      0      0
    ## 106      0      0      0      0      0      0      0      0      0      0
    ## 107    156      0      0      0      0      0      0      0      0      0
    ## 108      0      0      0      0      0      0      0      0      0      0
    ## 109      0      0      0      0      0      0      0      0      0      0
    ## 110      0      0      0      0     50      0      0      0      0      0
    ## 111      0      0      0      0      0     45      0      0      0      0
    ## 112      0      0      0      0      0     14    413      0      0      0
    ## 113      0      0      0      0      0      0    133      0      0      0
    ## 114      0      0      0      0      0     14      0      0      0      0
    ## 115      0      0      0      0      0     16      0      0      0      0
    ## 116      0      0      0      0      0     70      0      0      0      0
    ## 117      0      0      0      0      0      0     22      0      0      0
    ## 118      0      0      0      0      0      0      0      0      0      0
    ## 119      0      0      0      0      0      0      0      0      0      0
    ## 120      0      0      0      0      0      0      0      0      0      0
    ## 121      0      0      0      0    109      0      0      0      0      0
    ## 122      0      0      0      0      0      0      0      0      0      0
    ## 123      0      0      0      0    124      0      0      0      0      0
    ## 124      0      0      0      0     21      0      0      0      0      0
    ## 125      0      0      0      0      0      0      0      0      0      0
    ## 126      0      0      0      0      0      0      0      0      0      0
    ## 127      0      0      0      0      0      0      0      0      0      0
    ## 128      0      0      0      0      0      0      0      0      0      0
    ## 129      4      0      0      0      0      0      0      0     16      0
    ## 130      0      0      0      0      0      0      0      0      0      0
    ## 131      0      0      0      0      0      0      0      0     14      0
    ## 132      0      0      0      0      0      0      0      0    142      0
    ## 133     20      0      0      0      0      0      0     27      0      0
    ## 134      0      0      0      0      0      0      0      0      0      0
    ## 135      0      0      0      0      0      0      0      0      0      0
    ## 136      0      0      0      0     45      0      0      0      0      0
    ## 137      0      0      0      0      0      0      0      0     35      0
    ## 138      0      0      0      0      0      0      0      0      0      0
    ## 139      0      0      0      0      0      0      0      0      0      0
    ## 140      0      0      0      0      0      0      0      0      0      0
    ## 141      0      0      0      0      0      0      0      0      0      0
    ## 142      0      0      0      0     37      0      0      0      0      0
    ## 143      0      0      0      0      0      0      0      0      0      0
    ## 144      0      0      0      0     11      0      0      0      0      0
    ## 145      0      0      0      0      0      0      0      0      0      0
    ## 146      0      0      0      0      0      0      0      0      0      0
    ## 147      0      0      0      0      0      0      0      0      0      0
    ## 148      0      0      0      0      0      0      0      0      0      0
    ## 149      0      0      0      0      0      0      0      0      0      0
    ## 150      0      0      0      0      0      0      0      0      0      0
    ## 151      0      0      0      0     16      0      0      0      0      0
    ## 152      0      0      0      0      0      0      0      0     24      0
    ## 153      0      0      0      0      0      0      0      0      0      0
    ## 154      0      0      0      0     35      0      0      0      0      0
    ## 155      0      0     11      0      0      0      0      0      0      0
    ## 156      0      0      0      0    107      0      0      0      0      0
    ## 157      0      0      0      0     21      0      0      0      0      0
    ## 158      0      0      3      0      0      0      0      0      0      0
    ## 159      0      0      0      0    149      0      0      0      0      0
    ## 160      0      0      0      0      0      0      0      0      0      0
    ## 161      0      0      0      0      0      0      0      0      0      0
    ## 162      0      0      0      0     31      0      0      0      0      0
    ## 163      0      0      0      0      0      0      0      0      0      0
    ## 164      0      0      0      0      0      0      0      0      0      0
    ## 165      0      0      0      0      0      0      0      0      0      0
    ## 166      0      0      0      0      0      0      0      0      0      0
    ## 167      0     13      0      0      0      0      0      0      0      0
    ## 168      0      0      0      0      0      0      0      0      0      0
    ## 169      0      0      0      0      0      0      0      0      0      0
    ## 170      0      0      0      0      0      0      0      0      0      0
    ## 171      0      0      0      0      0      0      0      0      0      0
    ## 172      0      0      0      0      0      0      0      0      0      0
    ## 173      0    112      0     18      0      0      0      0      0      0
    ## 174      0     99      0      0      0      0      0      0      0      0
    ## 175      0      0      0      0      0      0      0      0      0      0
    ## 176      0    150      0      0      0      0      0      0      0      0
    ## 177      0      0      0      0      0      0      0      0      0      0
    ## 178      0     67      0      0      0      0      0      0      0      0
    ## 179      0      0      0      0      0      0      0      0      0      0
    ## 180      0      0      0      0      0      0      0      0      0      0
    ## 181      0      0      0      0      0      0      0      0      0      0
    ## 182      0      0      0      0      0      0      0      0      0      0
    ## 183      0      0      0      0      0      0      0      0      0      0
    ## 184      0     64      0      0      0      0      0      0      0      0
    ## 185      0      0      0      0      0      0      0      0      0      0
    ## 186      0      0      0      0      0      0      0      0      0      0
    ## 187      0      0      0      0      0      0      0      0      0      0
    ## 188     14      0      0      0      0      0      0      0      0      0
    ## 189      0      0      0      0      0      0      0      0      0      0
    ## 190      0      0      0      0      0      0      0      0      0      0
    ## 191      0      0      0      0      0      0      0      0      0      0
    ## 192      0      0      0      0      0      0      0      0      0      0
    ## 193      0      0      0      0      0      0      0      0      0      0
    ## 194      0      0      0      0     17      0      0      0      0      0
    ## 195      0      0      0      0      0      0      0      0      0      0
    ## 196      0      0      0      0      0      0      0      0      0      0
    ## 197      0      0      0      0      0      0      0      0      0      0
    ## 198      0      0      0      0      0      0      0      0      0      0
    ## 199      0      0      0      0      0      0      0      0      0      0
    ## 200      0      0      0      0      0      0      0      0      0      0
    ## 201      0      0      0      0      0      0      0      0      0      0
    ## 202      0      0      0      0      0      0      0      0      0      0
    ## 203      0      0      0      0      0      0      0      0      0      0
    ## 204      0      0      0      0      0      0      0      0      0      0
    ## 205      0      0      0      0      0      0      0      0      0      0
    ## 206      0      0      0      0      0      0      0      0      0      0
    ## 207      0      0      0      0      0      0      0      0      0      0
    ## 208      0      0      0      0      0      0      0      0      0      0
    ## 209      0     13      0      0      0      0      0      0      0      0
    ## 210      0    195      0      0      0      0      0      0      0      0
    ## 211      0      0      2      0      0      0      0      0      0      0
    ## 212      0      0      0      0      0      0      0      0      0      0
    ## 213      0      0      0      0      0      0      0      0      0      0
    ## 214      0      0      0      0      0      0      0      0      0      0
    ## 215     16      0      0      0      0      0      0      0      0      0
    ## 216      0      0      0      0      0      0      0      0      0    132
    ## 217      0      0      0      0      0      0      0      0      0      0
    ## 218      0      0      0      0      0      0      0      0      0    101
    ## 219      0      0      0      0      0      0      0      0      0      0
    ## 220      0      0      0      0      0      0      0      0      0      0
    ## 221      0      0      0      0      0      0      0      0      0    129
    ## 222      0      0      0      0      0      0      0      0      0      0
    ## 223      0      0      0      0      0      0      0      0      0      0
    ## 224      0      0      0      0      0      0      0      0      0      0
    ## 225      0      0      0      0      0      0      0      0      0      0
    ##     otu_90 otu_91 otu_92 otu_93 otu_94 otu_95 otu_96 otu_97 otu_98 otu_99
    ## 1        0      0      0      0      0      0      0      0      0      0
    ## 2        0      0      0     14      0      0      0      0      0      0
    ## 3        0      0      0      0      0      0      0      0      0      0
    ## 4        0      0      0      0      0      0      0      0      0      0
    ## 5        0      0      0      0      0      0      0      0      0      0
    ## 6        0      0      0      0      0      0      0      0      0      0
    ## 7        0      0      3      0      0      0      0      0      0      0
    ## 8        0      0      0      0      0      0      0      0      0      0
    ## 9        0      0      0      0      0      0      0      0      0      0
    ## 10       0      0      0      0      0      0      0     17      0      0
    ## 11       0      0      0      0      0      0      0      0      0      0
    ## 12       0      0      0      0      0      0      0      0      0      0
    ## 13       0      0      0      0      0      0      0      0      0      0
    ## 14       0      0      0      0      0      0      0      0      0      0
    ## 15       0      0      0      0      0      0      0     40      0      0
    ## 16       0      0      0      0      0      0      0      0      0      0
    ## 17       0      0      0      0      0      0      0      0      0      0
    ## 18       0      0      0      0      0      0      0      7      0      0
    ## 19       4      0      0      0      0      0      0      0      0      0
    ## 20       0      0      0      0      0      0      0      0      0      0
    ## 21       0      0      0      0      0      0      0      0      0      0
    ## 22       0      0      0      0      0      0      0      0      0      0
    ## 23       0      0      0      0      0      0      0      0      0      0
    ## 24       0      0      0      0      0      0      0      0      0      0
    ## 25       0      0      0      0      0      0      0      0      0      0
    ## 26      41      0      0     24      0      0      0      0     15      0
    ## 27       0      0      0      0      0      0      0      0      0      0
    ## 28       0      0      0      0      0      0      0      0     12      0
    ## 29       0      0      3      0      0      0      0      0      0      0
    ## 30       0      0      0      0      0      0      0      0      0      0
    ## 31       0      0      0      0      0      0      0      0      0      0
    ## 32       0      0      0      0      0      0      0      0      0      0
    ## 33       0      0      0      0      0      0      0     49      0      0
    ## 34       0      0      0      0      0      0      0      0      0      0
    ## 35       0      0      0      0      0      0      0      0      0      0
    ## 36       0      9      0      0      0      0      0      0      0      0
    ## 37       0      0      0      0      0      0      0      0      0      0
    ## 38       0      0      0      0      0      0      0      0     13      0
    ## 39       0      0      0      0      0      0      0      0     24      0
    ## 40     109      0      0      0      0      0      0      0    141      0
    ## 41      48      0      0      0      0      0      0      0     19      0
    ## 42       0      0      0      0      0      0      0     34      0      0
    ## 43       0      0      0      0      0      0      0      0      0      0
    ## 44       0      0      0      0      0      0      0      0      0      0
    ## 45       0      0      0      0      0      0      0     13      0      0
    ## 46       0     13      0      0      0      0      0      8      0      0
    ## 47       0      0      0      0      0      0      0      0      0      0
    ## 48       0     41      0      0      0      0      0      0      0      0
    ## 49       0      0      0     22      0      0      0      0      0      0
    ## 50       0      0      0      0      0      0      0      0      0      0
    ## 51       0      0      0     13      0      0      0      0      0      0
    ## 52       0      3      0      0      0      0      0     21      0      0
    ## 53       0     19      0      0      0      0      0      0      0      0
    ## 54       0      0      0      0      0      0      0      0      0      0
    ## 55       0      0      0      0      0      0      0      0      0      0
    ## 56       0      0      0      0      0      0      0      0      0      0
    ## 57       0      0      0      0      0      0      0      0      0      0
    ## 58       0      0      0      0      0      0      0      0      0      0
    ## 59       0      0      0      0      0      0      0      0      0      0
    ## 60       0      0      0      0      0      0      0      0      0      0
    ## 61       0      0      0      0      0      0      0      0      0      0
    ## 62       0      0      0      0      0      0      0      0      0      0
    ## 63       0      0      0      0      0      0      0      0      0      0
    ## 64       0      0      0      0      0      0      0      0      0      0
    ## 65       0      0      0      0      0      0      0      0      0      0
    ## 66       0      0      0      0      0      0      0      0      0      0
    ## 67       0      0      0      0      0      0      0      0      0      0
    ## 68       0      0      0      0      0      0      0      0      0      0
    ## 69       0      0      0      0      0      0      0      0      0      0
    ## 70       0      0      0      0      0      0      0      0      0      0
    ## 71       0      0    121      0      0      0      0      0      0      0
    ## 72       0      0      0      0      0      0      0      0      0      0
    ## 73       0    101      0      0      0      0      0      0      0      0
    ## 74      10      0      0      0      0      0      0      0      0      0
    ## 75       0     92      0      0      0      0      0      0      0      0
    ## 76       0      0      0      0      0      0      0      0      5      0
    ## 77       0      0      0      0      0      0      0      0      0      0
    ## 78       0      0      0      0      0     40      0      0      0      0
    ## 79       0      0      0      0      0      0      0      0     32      0
    ## 80       0      0      0      0      0      0      0      0      0      0
    ## 81       0      0      0      0      0      0      0      0     16      0
    ## 82       0      0      0      0      0      0      0      0      0      0
    ## 83       0      0      0      0      0      0      0      0      0      0
    ## 84       0      0      0      0      0      0      0      0      0      0
    ## 85       0      0      0    110      0      0      0      0      0      0
    ## 86       0      0      0      0      0      0      0      0      0      0
    ## 87       0      0      0      0      0      0      0      0      0      0
    ## 88       0      0      0      0      0      0      0      0      0      0
    ## 89       0      0      0      0      0      0      0      0      0      0
    ## 90     119      0      0      0      0      0      0      0      0      0
    ## 91       0      0      0      0    155      0      0      0      0      0
    ## 92       0      0      0      0     52      0      0      0      0      0
    ## 93       0      0      0      0      0      0      0      0      0      0
    ## 94       0      0      0      0      0      0      0      0      0      0
    ## 95       0      0      0      0      0      0      0      0      7      0
    ## 96       0      0     46    106      0      0      0      0      0      0
    ## 97       0      0      0      0      0      0      0      0      0      0
    ## 98       0      0      0      5      0      0      0      0      0      0
    ## 99       0      0      0     46      0      0      0      0      0      0
    ## 100     42      0      0      0      0      0      0      0      7      0
    ## 101      0      0      0      0      0      0      0      0     19      0
    ## 102      0      0     20      0      0      0      0      0      0      0
    ## 103      0      0      0      0      0      0      0      0      0      0
    ## 104      0      0      0      0      0      0      0      0      0      0
    ## 105      0      0      0      0      0      0      0      0      0     17
    ## 106     42      0      0      0      0      0      0      0     19      0
    ## 107      0      0      0      0      0      0      0      0      0      0
    ## 108     90      0      0      0      0      0      0      0      0      0
    ## 109      0      0      0      0      0      0      0      0      0      0
    ## 110      0      0      0      0      0      0      0      0      0      0
    ## 111      0      0      0      0      0      0      0      0      0      0
    ## 112      0      0      0      0      0      0      0      0      0      0
    ## 113      0      0      0      0      0      0      0      0      0      0
    ## 114      0      0      0      0      0      0      0      0      0      0
    ## 115      0      0      0      0      0      0      0      0      0      0
    ## 116      0      0      0      0      0      0      0      0      0      0
    ## 117      0      0      0      0      0      0      0      0      0      0
    ## 118      0      0      0      0      0      0      0      0      0      0
    ## 119      0      0      0      0      0      0      0      0      0      0
    ## 120      9      0     20      0      0      0      0      0      0      0
    ## 121      0      0      0      0      0      0      0      0      0      0
    ## 122      0      0      0      0      0      0      0      0      0      0
    ## 123      0      0      0      0      0      0      0      0      9      0
    ## 124      0      0      0      0      0      0      0      0      0      0
    ## 125      0      0      0      0      0      0      0      0      0      0
    ## 126      0      0      0      0      0      0      0      0      0      0
    ## 127      0      0      0      0      0      0      0      0      0     25
    ## 128      0      0      0      0      0      0      0      0      0      0
    ## 129      0      0      0      0      0      0      0      0      0      0
    ## 130      0      0      0      0      0      0      0      0      0      0
    ## 131      0      0      0      0      0      0      0      0      0      0
    ## 132      0      0      0      0      0      0      0      0      0      0
    ## 133      0      0      6      0      0      0      0      0      0      0
    ## 134      0      0      0      0      0      0      0      0      0      0
    ## 135      0      0      0      0      0      0      0      0      0      0
    ## 136      0      0      0      0      0      0      0      0      0      0
    ## 137      0      0      0      0      0      0      0      0      0      0
    ## 138      0      0      0      0      0      0      0      0      0      0
    ## 139      0      0      0      0      0      0      0      0      0      0
    ## 140      0      0      0      0      0      0      0      0      0      0
    ## 141      0      0      0      0      0      0      0      0      0     26
    ## 142      0      0      0      0      0      0      0      0      0      0
    ## 143      0      0     15      0      0      0      0      0      0      0
    ## 144      0      0      0      0      0      0      0      0      0     16
    ## 145      0      0      0      0      0      0      0      0      0      0
    ## 146      0      0      0      0      0      0      0      0      0     42
    ## 147      0      0      0      0      0      0      0      0      0      0
    ## 148      0      0      0      0      0      0      0      0      0      0
    ## 149      0      0      0      0      0      0      0      0      0      0
    ## 150      0      0      0      0      0      0    196      0      0      0
    ## 151      0      0      0      0      0      0      0      0      0     22
    ## 152      7      0      0      0      0      0      0      0      0      0
    ## 153      0      0      0      0      0      0      0      0      0      0
    ## 154      0      0     19      0      0      0      0      0      0      0
    ## 155      0      0      0      0      0      0      0      0      0      0
    ## 156      0      0      0      0      0      0      0      0      0      0
    ## 157      0      0     10     16      0      0      0      0      0      0
    ## 158      0      0      0      0      0      0      0      0      0      0
    ## 159      0      0      0      0      0      0      0      0      0      0
    ## 160      0      0      0      0      0      0      0      0      0      0
    ## 161      0      0      0      0      0      0      0      0      0      0
    ## 162      0      0    123      0      0      0      0      0      0      0
    ## 163      0      0      0      0      0      0      0      0      0      0
    ## 164      0      0      0      0      0      0      0      0      0      0
    ## 165      0      0      0      0      0      0      0      0      0      0
    ## 166      0      0      0      0      0     32      0      0      0      0
    ## 167      0      0      0      0      0     85      0      0      0      0
    ## 168      0      0      0      0      0    102      0      0      0      0
    ## 169      0      0      0      0      0     15      0      0      0      0
    ## 170      0      0      0      0      0      0      0      0      0      0
    ## 171      0      0      0      0      0      0      0      0      0      0
    ## 172      0      0      0      0      0    137      0      0      0      0
    ## 173      0      0      0      0      0     26      0      0      0      0
    ## 174      0      0      0      0      0      0      0      0      0      0
    ## 175      0      0      0      0      0      0      0      0      0      0
    ## 176      0      0      0      0      0      0      0      0      0      0
    ## 177      0      0      0      0      0      0      0      0      0      0
    ## 178      0      0      0      0      0      0      0      0      0      0
    ## 179      0      0      0      0      0      0      0      0      0      0
    ## 180      0      0      0      0      0      0      0      0      0      0
    ## 181      0      0      0      0      0      0      0      0      0      0
    ## 182      0      0      0      0      0      0      0      0      0      0
    ## 183      0      0      0      0      0      0      0      0      0      0
    ## 184      0      0      0      0      0      0      0      0      0      0
    ## 185      0      0      0      0      0      0      0      0      0      0
    ## 186      0      0      0      0      0      0      0      0      0      0
    ## 187      0      0      0      0      0      0      0      0      0      0
    ## 188      0      0      0      0      0      0      0      0      0      0
    ## 189      0      0      0      0      0      0      0      0      0      0
    ## 190      0      0      0      0      0      0      0      0      0      0
    ## 191      0      0      0      0      0      0      0      0      0      0
    ## 192      0      0      0      0      0      0      0      0      0      0
    ## 193      0      0      0      9      0      0      0      0      0      0
    ## 194      0      0      0      0      0      0      0      0      0      0
    ## 195      0      0      0      0      0      0      0      0      0      0
    ## 196      0      0      0      0      0      0      0      0      0      0
    ## 197      0      0     21      0      0      0      0      0      0      0
    ## 198      0      0      0      0      0      0      0      0      0      0
    ## 199      0      0      0      0      0      0      0      0      0      0
    ## 200      0      0      0      0      0      0      0      0      0      0
    ## 201      0      0      0      0      0      0      0      0      0      0
    ## 202      0      0      0      0      0      0      0      0      0      0
    ## 203      0      0      0      0      0      0      0      0      0      0
    ## 204      0      0      0      0      0      0      0      0      0      0
    ## 205      0      0      0      0      0      0      0      0      0      0
    ## 206      0      0      0      0      0      0      0      0      0      0
    ## 207      0      0      0      0      0      0      0      0      0      0
    ## 208      0     32      0      0      0      0      0      0      0      0
    ## 209      0      0      0      0      0      0      0      0      0      0
    ## 210      0      0      0      0      0      0      0      0      0      0
    ## 211      0      0      0      0      0      0      0      0      0      0
    ## 212      0      0      0      0      0      0      0      0      0      0
    ## 213      0      0      0      0      0      0      0      0      0      0
    ## 214      0      0      0      0      0      0      0      0      0      0
    ## 215      0      0      0     36      0      0      0      0      0      0
    ## 216      0      0      0      0      0      0      0      0      0      0
    ## 217      0      0      0      0      0      0      0      0      0      0
    ## 218      0      0      0      0      0      0      0      0      0      0
    ## 219      0      0      0      0      0      0      0      0      0      0
    ## 220      0      0      0      0      0      0      0      0      0      0
    ## 221      0      0      0      0      0      0      0      0      0      0
    ## 222      0      0      0      0      0      0      0      0      0      0
    ## 223      0      0      0      0      0      0      0      0      0      0
    ## 224      0      0      0      0      0      0      0      0      0      0
    ## 225      0      0      0      0      0      0      0      0      0      0
    ##     otu_100 otu_101 otu_102 otu_103 otu_104 otu_105 otu_106 otu_107 otu_108
    ## 1         0       0       0       0       0       0       0       0       0
    ## 2         0       0      13       0       0       0       0       0       0
    ## 3         0       0      26       0       0       0       0       0       0
    ## 4         0       0       0       0       0       0       0       0       0
    ## 5         0       0       0       0       0       0       0       0       0
    ## 6         0       0       4       0       0       0       0       0       0
    ## 7         0       0     186       0       0       0       0       0       0
    ## 8        43       0       0       0       0       0       0       0       0
    ## 9         0       0       0       0       0       0       0       0       0
    ## 10        0       0       0       0       0       0       0       0       0
    ## 11        0       0       0       0       0       0       0       0       0
    ## 12        0       0       0       0       0       0       0       0       0
    ## 13        0       0       0       0       0       0       0       0       0
    ## 14        0       0       0       0       0       0       0       0       0
    ## 15       12       0       0       0       0       0       0       0       0
    ## 16        0       0       0       0       0       0       0       0       0
    ## 17        0       0       0       0       0       0       0       0       0
    ## 18        0       0       0       0       0       0       0       0       0
    ## 19        0       0       0       0       0       0       0       0       0
    ## 20        0       0       0       0      34       0       0       0       0
    ## 21        0       0       0       0       0       0       0       0       0
    ## 22        0       0       0       0       0       0       0       0       0
    ## 23        0       0       0      13      88       0       0       0       0
    ## 24        0       0       0       0       0       0       0       0       0
    ## 25        0       0       0       0       0       0       0       0       0
    ## 26        0       0       0       8       0       0       0       0       0
    ## 27        0       0       0       0       0       0       0       0       0
    ## 28        0       0       0       0       0       0       0       0       0
    ## 29       39       0       0       0       0       0       0       0       0
    ## 30        0       0       0       0       0       0       0       0      29
    ## 31        0       0       0       0       0       0       0       0       0
    ## 32        0       0       0       0       0       0       0       0      36
    ## 33       18       0       0       0       0       0       0       0       0
    ## 34        0       0       0       0       0       0       0       0       0
    ## 35        0       0       0       0       0       0       0       0      10
    ## 36        0       0       0       0       0       0       0       0       0
    ## 37        0       0       0       0       0       0       0       0       0
    ## 38       22       0       0       0       0       0       0       0       0
    ## 39        0       0       0       0       0       0       0       0       0
    ## 40        0       0       0       0       0       0       0       0       0
    ## 41       10       0       0       0       0       0       0       0       0
    ## 42       15       5       0       0       0       0       0       0       0
    ## 43        0       0       0       0       0       0       0       0       0
    ## 44        0       0       0       0       0       0       0       0       0
    ## 45       22       0       0       0       0       0       0       0       0
    ## 46       16       0       0       0       0       0       0       0       0
    ## 47        0       7       0       0       0       0       0       0       0
    ## 48        0      73       0       0       0       0       0       0       0
    ## 49        0       0       0       0       0       0       0       0       0
    ## 50        0       0       0       0       0       0       0       0       0
    ## 51        0       0       0       0       0       0       0       0       0
    ## 52       20       0       0       0       0       0       0       0       0
    ## 53        0       0       0       0       0       0       0       0       0
    ## 54        0       0       0       0       0       0       0       0      21
    ## 55        3       0       0       0       0       0       0       0       0
    ## 56        0       0       0       0       0       0       0       0       0
    ## 57        0       0       0       0       0       0       0       0       0
    ## 58        0       0       0       0       0      18       0       0       0
    ## 59        0       0       0       0       0       0       0       0       0
    ## 60       14       0       0       0       0       0       0       0       0
    ## 61       49       0       0       0       0       0       0       0       0
    ## 62        0       0       0       0       0      20       0       0       0
    ## 63        0       0       0       0       0       0       0       0       0
    ## 64        0       0       0       0       0       0       0       0       0
    ## 65        0       0       0       0       0       0       0       0       0
    ## 66        0       0       0       0       0       0       0       0       0
    ## 67        0       0       0       0       0       0       0       0       0
    ## 68        0       0       0       0       0       0       0       0       0
    ## 69        0       0       0       0       0       0       0       0       0
    ## 70        0       0       0       0       0       0       0       0       0
    ## 71        0       0       0       0       0       0       0       0       0
    ## 72        0       0       0       0       0       0       0       0       0
    ## 73        0       0       0       0       0       0       0       0       0
    ## 74        0       0       0       0       0       0       0       0       0
    ## 75        0       0       0       0       0       0       0       0       0
    ## 76        0       0       0       0       0       0       0       0       0
    ## 77        0       0       0       0       0       0       0       0       0
    ## 78        0       0       0       0       0       0       0       0       0
    ## 79        0       0       0       0       0       0       0       0       0
    ## 80        0       0      49       0       0       0       0       0       0
    ## 81        0       0      36       0       0       0       0       0       0
    ## 82        0       0       0       0       0       0       0       0       0
    ## 83       19       0       0       0       0       0       0       0       0
    ## 84        0       0       0       0       0       0       0       0       0
    ## 85        0       0       0       0       0       0       0       0       0
    ## 86        0       0       0       0       0       0       0       0       0
    ## 87        0       0       0       0       0       0       0       0       0
    ## 88        0       0       0       0       0       0       0       0       0
    ## 89        0       0       0       0       0       0       0       0       0
    ## 90        0       0       0       0       0       0       0       0       0
    ## 91        0       0       0       0       0       0       0       0       0
    ## 92        0      11       0       0       0       0       0       0       0
    ## 93        0       0       0       0       0       0       0       0       0
    ## 94        0       0       0       0       0       0       0       0       0
    ## 95        0       0       0       0       0       0       0       0       0
    ## 96        0       0       0       0       0       0       0       0       0
    ## 97        0       0       0       0       0       0       0       0       0
    ## 98        0       0       0       0       0       0       0       0       0
    ## 99        0       0       0       0       0       0       0       0       0
    ## 100       0       0       0       0       0       0       0       0       0
    ## 101       0       0       0       0       0       0       0       0       0
    ## 102       0       0       0       0       0       0       0       0       0
    ## 103       0       0       0       0       0       0       0       0       0
    ## 104       0       0       0       0       0       0       0       0       0
    ## 105      30       0       0       0       0       0       0       0       0
    ## 106       0       0       0       0       0       0       0       0       0
    ## 107       0       0       0       0       0       0       0       0       0
    ## 108       0       0       0       0       0       0       0       0       0
    ## 109       0       7       0     155       0       0       0       0       0
    ## 110       0      19       0       0       0       0       0       0       0
    ## 111       0       0       0       0       0       0       0       0       0
    ## 112       0       0       0       0       0       0       0       0       0
    ## 113       0       0       0       0       0       9       0       0       0
    ## 114       0      17       0       0       0       0       0       0       0
    ## 115       0       0       0       0       0       0       0       0       0
    ## 116       0       0       0       0       0       0       0       0       0
    ## 117       0       0       0       0       0       0       0       0       0
    ## 118       0       0       0       0       0       0       0       0       0
    ## 119       0       0       0       0       0       0       0       0       0
    ## 120       0       0       0       0       0       0       0       0       0
    ## 121       0       0       0       0       0      12       0       0       0
    ## 122       0       0       0       0       0       0       0       0       0
    ## 123       0       0       0       0       0       0       0       0       0
    ## 124       0       0       0       0       0       0       0       0       0
    ## 125       0       0       0       0       0       0       0       0       0
    ## 126       0       0       0       0       0       0       0       0       0
    ## 127       0       0       0       0       0       0       0       0       0
    ## 128      20      19       0       0       0       0       0       0       0
    ## 129       0       0       0       0       0       0       0       0       0
    ## 130       0       0       0       0       0       0       0       0       0
    ## 131       0       0       0       0       0       0       0       0       0
    ## 132       0       0       0       0       0       0       0       0       0
    ## 133       0       0       0       0       0       0       0       0       0
    ## 134       0       0       0       0       0       0       0       0       0
    ## 135       0       0       0       0       0       0       0       0       0
    ## 136       0       0       0       0       0       0       0       0       0
    ## 137       0       0       0       0       0       0       0       0       0
    ## 138       0       0       0       0       0       0       0       0       0
    ## 139       0       0       0       0       0       0       0       0       0
    ## 140       0       0       0       0       0       0       0       0       0
    ## 141       0       0       0       0       0       0       0       0       0
    ## 142       0       0       0       0       0       0       0       0       0
    ## 143       0       0       0       0       0       0       0       0       0
    ## 144       0       0       0       0       0       0       0       0       0
    ## 145       0       0       0       0       0       0       0       0       0
    ## 146       0       0       0       0       0       0       0       0       0
    ## 147       0       0       0       0       0       0       0       0       0
    ## 148       0       0       0       0       0       0       0       0       0
    ## 149       0       0       0       0       0       0       0       0       0
    ## 150       0       0       0       0       0       0       0       0       0
    ## 151       0       0       0       0       0       0       0       0       0
    ## 152       0       0       0       0       0       0       0       0       0
    ## 153       0       0       0       0       0       0       0       0       0
    ## 154       0       0       0       0       0       0       0       0       0
    ## 155       0       0       0       0       0       0       0       0       0
    ## 156       0       0       0       0       0       0       0       0       0
    ## 157       0       0       0       0       0       0       0       0       0
    ## 158       0       0       0       0       0       0       0       0       0
    ## 159      13       0       0       0       0       0       0       0       0
    ## 160       0       0       0       0       0       0       0       0       0
    ## 161       0       0       0       0       0       0       0       0       0
    ## 162       0       0       0       0       0       0       0       0       0
    ## 163       0       0       0       0       0       0       0       0       0
    ## 164       0       0       0       0       0       0       0       0       0
    ## 165       0       0       0       0       0       0       0       0       0
    ## 166       0       0       0       0       0       0       0       0       0
    ## 167       0       0       0       0       0       0       0       0       0
    ## 168       0       0       0       0       0       0       0       0       0
    ## 169       0       0       0       0       0       0       0       0       0
    ## 170       0       0       0       0       0       0       0       0       0
    ## 171       0       0       0       0       0       0       0       0       0
    ## 172       0       0       0       0       0       0       0       0       0
    ## 173       0       0       0       0       0       0       0       0       0
    ## 174       0       0       0       0       0       0       0       0       0
    ## 175       0       0       0       0       0       0       0       0       0
    ## 176       0       0       0       0       0       0       0       0       0
    ## 177       0       0       0       0       0       0       0       0       0
    ## 178       0       4       0       0      29       0       0       0       0
    ## 179       0       9       0       0       0       0       0       0       0
    ## 180       0       0       0       0       0      31       0       0       0
    ## 181       0       0       0       0       0       0       0       0       0
    ## 182       0       0       0      16       0       0       0       0       0
    ## 183       0       0       0      15       0       0       0       0       0
    ## 184       0       0       0       0       0      79       0       0       0
    ## 185       0       0       0       0       0      38       0       0       0
    ## 186       0       0       0      59       0       0       0       0       0
    ## 187       8       0       0       0       0       0       0       0       0
    ## 188       0       0       0      33       0       0       0       0       0
    ## 189       0       0       0      16       0       0       0       0       0
    ## 190       0       0       0       0       0       0       0       0       0
    ## 191       0       0       0       0       0       0       0       0       0
    ## 192       0       0       0       0       0       0     116       0       0
    ## 193       0       0       0      15       0       0       0       0       0
    ## 194       0       0       0       0       0       0       0       0       0
    ## 195       0       0       0       0       0       0       0      41       0
    ## 196       0       0       0       0       0       0       0       0       0
    ## 197       0       0       0       0       0       0       0       0       0
    ## 198       0       0       0       0       0       0       0       0       0
    ## 199       0       0       0       0       0       0       0       0       0
    ## 200       0       0       0       0       0       0       0       0       0
    ## 201       0       0       0       0       0       0       0       0       0
    ## 202       0       0       0       0       0       0       0       0       0
    ## 203       0       0       0       0       0       0       0       0       0
    ## 204       0       0       0       0       0       0       0       0       0
    ## 205       0       0       0       0       0       0       0     129       0
    ## 206       0       0       0       0       0       0       0       0       0
    ## 207       0       0       0       0       0       0       0       0       0
    ## 208       0       0       0       0       0       0       0       0       0
    ## 209      13       0       0       0       0       0       0       0       0
    ## 210       0       0       0       0       0       0       0       0       0
    ## 211       0       0       0       0       0       0       0       0       0
    ## 212       0       0       0       0       0       0       0       0       0
    ## 213       0       0       0       0       0       0       0       0       0
    ## 214      64       0       0       0       0       0       0       0       0
    ## 215       0       0       0       0       0       0       0       0       0
    ## 216       0       0       0       0       0       0       0       0       0
    ## 217      46       0       0       0       0       0       0       0       0
    ## 218       0       0       0       0       0       0       0       0       0
    ## 219      14       0       0       0       0       0       0       0       0
    ## 220       0       0       0       0       0       0       0       0       0
    ## 221       0      27       0       0       0       0       0       0       0
    ## 222       0       0       0       0       0       0       0       0       0
    ## 223       0       0       0       0       0       0       0       0       0
    ## 224      12       0       0       0       0       0       0       0       0
    ## 225       0       0       0       0       0       0       0       0       0
    ##     otu_109 otu_110 otu_111 otu_112 otu_113 otu_114 otu_115 otu_116 otu_117
    ## 1         0       0       0       0       0       0       0       0       0
    ## 2         0       0       0       0       0       0       0       0       0
    ## 3         0       0       0       0       0       0       0       0       0
    ## 4         0       0       0       0       0       0       0       0       0
    ## 5         0       0       0       0       0       0       0       0       0
    ## 6         0       0       0       0       0       0       0       0       0
    ## 7         0       0       2       0       0       0       0       0       0
    ## 8         0       0       0       0       0       0       0       0       0
    ## 9         0       0       0       0       0       0       0       0       0
    ## 10        0       0       0       0       0       0       0       0       0
    ## 11        0       0       0       0       0       0       0       0       0
    ## 12        0       0       0       0       0       0       0       0       0
    ## 13        0       0       0       0       0       0       0       0       0
    ## 14        0       0       0       0       0       0       0       0       0
    ## 15        0       0       0       0       0       0       0       0       0
    ## 16        0       0       0       0       0       0       0       0       0
    ## 17        0       0       0       0       0       0       0       0       0
    ## 18        0       0       0       0       0       0       0       0       0
    ## 19        0       0       0       0       0       0       0       0       0
    ## 20        0       0       0       0       0       0       0       0       0
    ## 21        0       0       0       0       0       0       0       0       0
    ## 22        0       0       0       0       0       0       0       0       0
    ## 23        0       0       0       0       0      17       0       0       0
    ## 24        0       0       0       0       0       0       0       0       0
    ## 25        0       0       0       0       0       0       0       0       0
    ## 26        0       0       0       0       0       0       0       0       0
    ## 27        0       0       0       0       0       0       0       0       0
    ## 28        0       0       0       0       0       0      13       0       0
    ## 29        0       0       0       0       0       0       0       0       0
    ## 30        0       0       0       0       0       0      31       0       0
    ## 31        0       0       0       0       0       0       0       0       0
    ## 32        0       0       0       0       0       0       0       0       0
    ## 33        0       0       0       0       0       0       0       0       0
    ## 34        0       0       0       0       0       0       0       0      25
    ## 35        0       0       0       0       0       0       0       0      26
    ## 36        0       0       0       0       0       0       0       0       0
    ## 37        0       0       0       0       0       0       0       0       0
    ## 38        0       0       0       0       0       0       0       0       0
    ## 39        0       0       0       0       0       0       0       0       0
    ## 40        0       0       0       0       0       0       0       0       0
    ## 41        0       0       0       0       0       0       0       0       0
    ## 42        0       0       0       0       0       0       0       0       0
    ## 43        0       0       0       0       0       0       0       0       0
    ## 44        0       0      23       0       0       0       0       0       0
    ## 45        0       0      16       0       0       0       0       0       0
    ## 46        0       0       0       0       0       0       0       0       0
    ## 47        0       0       0       0       0       0       0       0       0
    ## 48        0       0       0       0       0       0       0       0       0
    ## 49        0       0       0       0       0       0      27       0       0
    ## 50        0       0       0       0       0       0       0       0       0
    ## 51        0       4       0       0       0       0      34       0       0
    ## 52        0       0       0       0       0       0      31       0       0
    ## 53        0       0       0       0       0       0       0       0       0
    ## 54        0       0       0       0       0       0       9       0       0
    ## 55        0       0       0       0       0       0       0       0       0
    ## 56        0       0       0       0       0       0       0       0       0
    ## 57        0       0       0       0       0       0       0       0       0
    ## 58        0       0       0       0       0       0       0       0       0
    ## 59        0       0       0       0       0       0       0       0       0
    ## 60        0       0       0       0       0       0      29       0       0
    ## 61        0       0       0       0       0       0       0       0       0
    ## 62        0       0       0       0       0       0       0       0       0
    ## 63        0       0       0       0       0       0       0       0       0
    ## 64        0       0       0       0       0       0       0       0       0
    ## 65        0       0       0       0       0       0       0       0       0
    ## 66        0       0       0       0       0       0       0       0       0
    ## 67        0       0       0       0       0       0       0       0       0
    ## 68        0       0       0       0       0       0       0       0       0
    ## 69        0       0       0       0       0       0       0       0       0
    ## 70        0       0       0       0       0       0       0       0       0
    ## 71        0       0       0       0       0       0       0       0       0
    ## 72        0       0       0       0       0       0       0       0       0
    ## 73        0       0       8       0       0       0       0       0       0
    ## 74        0       0       0       0       0       0       0       0       0
    ## 75       91       0       0       0       0       0       0       0       0
    ## 76        0       0       0       0       0       0       0       0       0
    ## 77        0       0       0       0       0       0       0       0       0
    ## 78        0       0       0       0       0       0       0       0       0
    ## 79        0       0       0       0       0       0       0       0       0
    ## 80        0       0       0       0       0       0       0       0       0
    ## 81        0       0       4       0       0       0       0       0       0
    ## 82        9       0       0       0       0       0       0       0       0
    ## 83        0       0       0       0       0       0       0       0       0
    ## 84        0       0       0       0       0       0       0       0       0
    ## 85        0       0       0       0       0       0       0       0       0
    ## 86        0       0       0       0       0       0       0       0       0
    ## 87        0      40       0       0       0       0       0       0       0
    ## 88        0       0       0       0       0       0       0       0       0
    ## 89        0       0       0       0       0       0       0       0       0
    ## 90        0      46       0       0       0       0       0       0       0
    ## 91        0       0       0       0       0       0       0       0       0
    ## 92        0       0       5       0       0       0       0       0       0
    ## 93        0       0       0       0       0       0       0       0       0
    ## 94        0       0       0       0       0       0       0       0       0
    ## 95        0       0       0       0       0       0      28       0       0
    ## 96        0       0      24       0       0       0       0       0       0
    ## 97        0       0       0       0       0       0       0       0       0
    ## 98        0       0      10       0       0       0       0       0       0
    ## 99        0       0      24       0       0       0       0       0       0
    ## 100       0       0      12       0       0       0       0       0       0
    ## 101       0       0       0       0       0       0       0       0       0
    ## 102       0       0       0       0       0       0       0       0       0
    ## 103       0       0       0       0       0       0       0       0       0
    ## 104       0       0       0       0       0       0       0       0       0
    ## 105       0       0       0       0       0       0      19       0       0
    ## 106       0       0       0       0       0       0       0       0       0
    ## 107       0       0       0       0       0       0       0       0       0
    ## 108       0       0       0       0       0       0       0       0       0
    ## 109       0       0       0       0       0       0      15       0       0
    ## 110       0       0       0       0       0       0       0       0       0
    ## 111       0       0       0       0       0       0       0       0       0
    ## 112       0       0       0       0       0       0       0       0       0
    ## 113       0       0       0       0       0       0       0       0       0
    ## 114       0       0       0       0       0       0       0       0       0
    ## 115       0       0       0       0       0       0       0       0       0
    ## 116       0       0       0       0       0       0       0       0       0
    ## 117       0       0       0       0       0       0       0       0       0
    ## 118       0       0       0       0       0       0       0       0       0
    ## 119       0       0       0       0       0       0       0       0       0
    ## 120       0      14       0       0       0       0       0       0       0
    ## 121       0       0       0       0       0       0       0       0       0
    ## 122       0       0       0       0       0       0       0       0       0
    ## 123       0      33       0       0       0       0       0       0       0
    ## 124       0       0       0       0       0       0       0       0       0
    ## 125       0       7       0       0       0       0       4       0       0
    ## 126       0       0       0       0       0       0       0       0       0
    ## 127       0       0       0       0       0       0       0       0       0
    ## 128       0       0       0       0       0       0       0       0       0
    ## 129       0       0       0       0       0       0       0       0       0
    ## 130       0       0       0       0       0       0       0       0       0
    ## 131       0       0       0       0       0       0       0       0       0
    ## 132       0       0       0       0       0       0       0       0       0
    ## 133       0       0       0       0       0       0       0      18       0
    ## 134       0       0       0       0       0       0       0       0       0
    ## 135       0       0       0       0       0       0       0       0       0
    ## 136       0       0       0       0       0       0       0       0       0
    ## 137       0       0       0       0       0       0       0      17       0
    ## 138       0       0       0       0      69       0       0       0       0
    ## 139       0       0       0       0       0       0       0       0       0
    ## 140       0       0       0       0       0       0       0       0       0
    ## 141       0       0       0       0       0       0       0       0       0
    ## 142       0       0      12       0       0       0       0       0       0
    ## 143       0       0       0       0       0       0      21       0       0
    ## 144       0       0       0       0       0       0       0       0       0
    ## 145       0       0       0       0       0       0       0       0       0
    ## 146       0       0       0       0       0       0       0       0       0
    ## 147       0       0       0       0       0       0       0       0       0
    ## 148       0       0       0       0       0       0       0       0       0
    ## 149       0       0       0       0       0       0       0       0       0
    ## 150       0       0       0       0       0       0       0       0       0
    ## 151       0       0       0       0       0       0       0       0       0
    ## 152       0       0       0       0       0       0       0       0       0
    ## 153       0       0       0       0       0       0       0       0       0
    ## 154       0       0       0       0       0       0       0       0       0
    ## 155       0       0       0       0       0       0       0       0       0
    ## 156       0       0       0       0       0       0       0       0       0
    ## 157       0       0       0       0       0       0       0       0       0
    ## 158       0       0       0       0       0       0       0       0       0
    ## 159       0       0       0       0       0       0       0       0       0
    ## 160       0       0       0       0       0       0       0       0       0
    ## 161       0       0       0       0       0       0       0       0       0
    ## 162       0       0       0       0       0       0       0       0       0
    ## 163       0       0       0       0       0       0       0       0       0
    ## 164       0       0       0       0       0       0       0       0       0
    ## 165       0       0       0       0       0       0       0       0       0
    ## 166       0       0       0       0       0       0       0       0       0
    ## 167       0       0       0       0       0       0       0       0       0
    ## 168       0       0       0       0       0       0       0       0       0
    ## 169       0       0       0      32       0       0       0       0       0
    ## 170       0       0       0       0       0       0       0       0       0
    ## 171       0       0       0       0       0       0       0       0       0
    ## 172       0       0       0       0       0       0       0       0       0
    ## 173       0       0       0       0       0       0       0       0       0
    ## 174       0       0       0       0       0       0       0       0       0
    ## 175       0       0       0       0       0       0       0       0       0
    ## 176       0       0       0       0       0       0       0       0       0
    ## 177       0       0       0       0       0       0       0       0       0
    ## 178       0       0       0       0       0       0       0       0       0
    ## 179       0       0       0       0       0       0       0       0       0
    ## 180       0       0       0       0       0       0       0       0       0
    ## 181       0       0       0       0       0       0       0       0       0
    ## 182       0       0       0       0       0       0       0       0       0
    ## 183       0       0       0      40       0       0       0       0       0
    ## 184       0       0       0       0       0       0       0       0       0
    ## 185       0       0       0       0       0       0       0       0       0
    ## 186       0       0       0       0      14       0       0       0       0
    ## 187       0       0       0       0       0       0       0       0       0
    ## 188       0       0       0       0       0       0       0      17       0
    ## 189       0       0       0       0       0       0       0       0       0
    ## 190       0       0       0       0       0       0       0       0       0
    ## 191       0       0       0       0       0       0       0       0       0
    ## 192       0       0       0       0       0       0       0       0       0
    ## 193       0       0       0       0       0       0       0       0       0
    ## 194       0       0       0       0       0       0       0       0       0
    ## 195       0       0       0       0       0       0       0       0       0
    ## 196       0       0       0       0       0       0       0       0       0
    ## 197       0       0       0       0       0       0      20       0       0
    ## 198       0       0       0       0       0       0       0       0       0
    ## 199       0       0       0       0       0       0       0       0       0
    ## 200       0       0       0       0       0       0       0       0       0
    ## 201       0       0       0       0       0       0       0       0       0
    ## 202       0       0       0       0       0       0       0       0       0
    ## 203       0       0       0       0       0       0       0       0       0
    ## 204       0       0       0       0       0       0       0       0       0
    ## 205       0       0       0       0       0       0       0       0       0
    ## 206       0       0       0       0       0       0       0       0       0
    ## 207       0       0       0       0       0       0       0       0       0
    ## 208       0       0       0       0       0       0       0       0       0
    ## 209       0       0       0       0       0       0       0       0       0
    ## 210       0       0       0       0       0       0       0       0       0
    ## 211       0       0       0       0       0       0       0       0       0
    ## 212       0       0       0       0       0       0       0       0       0
    ## 213       0       0       0       0       0       0       0       0       0
    ## 214       0       0       0       0       0       0       0       0       0
    ## 215       0       0       0       0       0       0       0       0       0
    ## 216       0       0       0       0       0       0       0       0       0
    ## 217       0       0       0       0       0       0       0       0       0
    ## 218       0       0       0       0       0       0       0       0       0
    ## 219       0       0       0       0       0      48       0       0       0
    ## 220       0       0       0       0       0       0       0       0       0
    ## 221       0       0       0       0       0       0       0       0       0
    ## 222       0       0       0       0       0       0       0       0       0
    ## 223       0       0       0       0       0       0       0       0       0
    ## 224       0       0       0       0       0       0       0       0       0
    ## 225       0       0       0       0       0       0       0       0       0
    ##     otu_118 otu_119 otu_120 otu_121 otu_122 otu_123 otu_124 otu_125 otu_126
    ## 1         0       0       0       0       0       0       0       0       0
    ## 2         0       0       0       0       0       0       0       0       0
    ## 3         0       0       0       0       0       0       0       0       0
    ## 4         0       0       0       0       0       0       0       0       0
    ## 5         0       0       0       0       0       0       0       0      26
    ## 6         0       0       0       0       0       0       0       0       0
    ## 7         0       0       0       0       0       0       0       0       0
    ## 8         0       0       0       0       0       0       0       0       0
    ## 9         0       0       0       0       0       0       0       0       0
    ## 10        0       0       0       0       0       0       0       0       0
    ## 11        0       0       0       0       0       0       0       0       0
    ## 12        0       0       0       0       0       0       0       0       0
    ## 13        0       0       0       0       0       0       0       0       0
    ## 14        0       0       0       0       0       0       0       0       0
    ## 15        0       0       0       0       0       0       0       0       0
    ## 16        0       0       0       0       0       0       0       0       0
    ## 17        0       0       0       0       0       0       0       0       0
    ## 18        0       0       0       0       0       0       0       0       0
    ## 19        0       0       0       0       0       0       0       0       0
    ## 20        0       0       0       0       0       0       0       0       0
    ## 21        0       0       0       0       0       0       0       0       0
    ## 22        0       0       0       0       0      18       0       0       0
    ## 23        0       0       0       0       0       0       0       0       0
    ## 24        0       0       0       0       0       0       0       0       0
    ## 25        0       0       0       0       0       0       0       0       0
    ## 26        0       0       0       0       0       0       0       0       0
    ## 27        0       0       0       0       0       0       0       0       0
    ## 28        0       0       0       0       0       0       0       0       0
    ## 29        0       0       0       0       0       0       0       0       0
    ## 30        0       0       0       0       0       0       0       0       0
    ## 31        0       0       0       0       0       0       0       0       0
    ## 32        0       0       0       0       0       0       0       0       0
    ## 33        0       0       0       0       0       0       0       0       0
    ## 34        0       0       0       0       0       0       0       0       0
    ## 35        0       0       0       0       0       0       0       0       0
    ## 36        0       0       0       0       0       0       0       0       0
    ## 37        0       0       0       0       0       0       0       0       0
    ## 38        0       0       0       0       0       0       0       0       0
    ## 39        0       0       0       0       0       0       0       0       0
    ## 40        0       0       0       0       0       0       0       0       0
    ## 41        0       0       0       0       0       0       0       0       0
    ## 42        0       0       0       0       0       0       0       0       0
    ## 43        0       0       0       0       0       0       0       0       0
    ## 44        0       0       0       0       0       0       0       0       0
    ## 45        0       0       0       0       0       0       0       0       0
    ## 46        0       0       0       0       0       0       0       0       0
    ## 47        0       0       0       0       0       0       0       0       0
    ## 48        0       0       0       0       0       0       0       0       0
    ## 49        0       0       0       0       0       0       0       0       0
    ## 50        0       0       0       0       0       0       0       0       0
    ## 51        0       0       0       0       0       0       0       0       0
    ## 52        0       0       0       0       0       0       0       0       0
    ## 53        0       0       0       0       0       0       0       0       0
    ## 54        0       0       0       0       0       0       0       0       0
    ## 55        0       0       0       0       0       0       0       0       0
    ## 56        0       0       0       0       0       0       0       0       0
    ## 57        0       0       0      25       0       0       0       0       0
    ## 58        0       0       0       0       0       0       0       0       0
    ## 59        0       0       0       0       0       0       0       0       0
    ## 60        0       0       0       0       0       0      27       0       0
    ## 61        0       0       0       0       0       0       0       0       0
    ## 62        0       0       0       0       0       0       0       0       0
    ## 63        0       0       0      33       0       0       0       0       0
    ## 64        0       0       0       0       0       0       0       0       0
    ## 65        0       0       0       0       0       0       0       0       0
    ## 66        0       0       0       0       0       0       0       0       0
    ## 67        0       0       0       0       0       0       0       0       0
    ## 68        0       0       0       0       0       0       0       0       0
    ## 69        0       0       0       0       0       0       0       0       0
    ## 70        0       0       0       0       0       0       0       0       0
    ## 71        0       0       0       0       0       0       0       0       0
    ## 72        0       0       0       0       0       0       0       0       0
    ## 73        0       0       0       0       0       0       0       0       0
    ## 74        0       0       0       0       0       0       0       0       0
    ## 75        0       0       0       0       0       0       0       0       0
    ## 76        0       0       0       0       0       0       0       0       0
    ## 77        0       0       0       0       0       0       0       0       0
    ## 78        0       0       0       0       0       0       0       0       0
    ## 79        0       0       0       0       0       0       0       0       0
    ## 80        0       0       0       0       0       0       0       0       0
    ## 81        0       0       0       0       0       0       0       0       0
    ## 82        0       0       0       0       0       0       0       0       0
    ## 83        0       0       0       0       0       0       0       0       0
    ## 84        0       0       0       0       0       0       0       0       0
    ## 85        0       0       0       0       0       0       0       0       0
    ## 86        0       0       0       0       0       0       0       0       0
    ## 87        0       0       0       0       0       0       0       0       0
    ## 88        0       0       0       0       0       0       0       0       0
    ## 89        0       0       0       0       0       0       0       0       0
    ## 90        0       0       0       0       0       0       0       0       0
    ## 91       32       0       0       0       0       0       0       0       0
    ## 92       16       0       0       0       0       0       0       0       0
    ## 93        0       0       0       0       0       0       0       0       0
    ## 94        0       0       0       0       0       0       0       0       0
    ## 95        0       0       0       0       0       0       0       0       0
    ## 96        0       0       0       0       0       0       0       0       0
    ## 97        0       0       0       0       0       0       0       0       0
    ## 98        0       0       0       0       0       0       0       0       0
    ## 99        0       0       0       0       0       0       0       0       0
    ## 100       0       0       0       0       0       0       0       0       0
    ## 101       0       0       0       0       0       0       0       0       0
    ## 102       0       0       0       0       0       0       0       0       0
    ## 103       0       0       0       0       0       0       0       0       0
    ## 104       0       0       0       0       0       0       0       0       0
    ## 105       0       0       0       0       0       0       0       0       0
    ## 106       0       0       0       0       0       0       0       0       0
    ## 107       0       0       0       0       0       0       0       0       0
    ## 108       0       0       0       0       0       0       0       0       0
    ## 109       0      45       0       0       0       0       0       0       0
    ## 110       0       0       0       0       0       0       0       0       0
    ## 111       0       0       0       0       0       0       0       0       0
    ## 112       0       0       0       0       0       0       0       0       0
    ## 113       0       0       0       0      31       0       0       0       0
    ## 114       0       0       0       0       0       0       0       0       0
    ## 115       0       0       0       0       0       0       0       0       0
    ## 116       0       0       0       0       0       0       0       0       0
    ## 117       0       0       0       0       0       0       0       0       0
    ## 118       0       0       0       0       0       0       0       0       0
    ## 119       0       0       0       0       0       0       0       0       0
    ## 120       0       0       0       0       0       0       0       0       0
    ## 121       0       0       0       0       0       0       0       0       0
    ## 122       0       0       0       0       0       0       0       0       0
    ## 123       0       0       0       0       0       0       0       0       0
    ## 124       0       0       0       0       0       0       0       0       0
    ## 125       0       0       0       0       0       0       0       0       0
    ## 126       0       0       0       0       0       0       0       0       0
    ## 127       0       0       0       0       0       0       0       0       0
    ## 128       0       0       0       0       0       0       0       0       0
    ## 129       0       0       0       0       0       0       0       0       0
    ## 130       0       0       0       0       0       0       0       0       0
    ## 131       0       0       0       0       0       0       0       0       0
    ## 132       0       0       0       0       0       0       0       0       0
    ## 133       0       0       0       0       0       0       0       0       0
    ## 134       0       0       0       0       0       0       0       0       0
    ## 135       0       0       0       0       0       0       0       0       0
    ## 136       0       0      34       0       0       0       0       0       0
    ## 137       0       0       0       0       0       0       0       0       0
    ## 138       0       0       0       0       0       0       0       0       0
    ## 139       0       0       0       0       0       0       0       0       0
    ## 140       0       0       0       0       0       0       0       0       0
    ## 141       0       0       0       0       0       0       0       0       0
    ## 142       0       0       0       0       0       0       0       0       0
    ## 143       0       0       0       0       0       0       0       0       0
    ## 144       0       0       0       0       0       0       0       0       0
    ## 145       0       0       0       0       0       0       0       0       0
    ## 146       0       0       0       0       0       0       0       0       0
    ## 147       0       0       0       0       0       0       0       0       0
    ## 148       0       0       0       0       0       0       0       0       0
    ## 149       0       0       0       0       0       0       0       0       0
    ## 150       0       0       0       0       0       0      27       0       0
    ## 151       0       0       0       0       0       0       0       0       0
    ## 152       0       0       0       0       0       0       0       0       0
    ## 153       0       0       0       0       0       0       0       0       0
    ## 154       0       0       0       0       0       0       0       0       0
    ## 155       0       0       0       0       0       0       0       0       0
    ## 156       0       0       0       0       0       0       0       0       0
    ## 157       0       0       0       0       0       0       0       0       0
    ## 158       0       0       0       0       0       0       0       0       0
    ## 159       0       0       0       0       0       0       0       0       0
    ## 160       0       0       0       0       0       0       0       0       0
    ## 161       0       0       0       0       0       0       0       0       0
    ## 162       0       0       0       0       0       0       0       0       0
    ## 163       0       0       0       0       0       0       0       0       0
    ## 164       0       0       0       0       0       0       0       0       0
    ## 165       0       0       0       0       0       0       0       0       0
    ## 166       0       0       0       0       0       0       0       0       0
    ## 167       0       0       0       0       0       0       0       0       0
    ## 168       0       0       0       0       0       0       0       0       0
    ## 169       0       0       0       0       0       0       0       0       0
    ## 170       0       0       0       0       0       0      18       0       0
    ## 171       0       0       0       0       0       0       0       0       0
    ## 172       0       0       0       0       0       0       0       4       0
    ## 173       0       0       0       0       0       0       0       0       0
    ## 174       0       0       0       0       0       0       0       0       0
    ## 175       0       0      12       0       0       0       0       0       0
    ## 176       0       0       0       0       0       0       0       0       0
    ## 177       0       0       0       0       0       0       0       0       0
    ## 178       0       0       0       0       0       0       0      23       0
    ## 179       0       0       0       0       0       0       0       0       0
    ## 180       0       0       0       0       0       0       0       0       0
    ## 181       0       0       0       0       0       0       0       0       0
    ## 182       0       0       0       0       0       0       0       0       0
    ## 183       0       0       0       0       0       0       0       0       0
    ## 184       0       0       0       0       0       0       0       0       0
    ## 185       0       0       0       0       0       0       0       0       0
    ## 186       0       0       0       0       0       0       0       0       0
    ## 187       0       0       0       0       0       0       0       0       0
    ## 188       0       0       0       0       0       0       0       0       0
    ## 189       0       0       0       0       0       0       0       0       0
    ## 190       0       0       0       0       0       0       0       0       0
    ## 191       0       0       0       0       0       0       0       0       0
    ## 192       0       0       0       0       0       0       0       0       0
    ## 193       0       0       0       0       0       0       0       0       0
    ## 194       0       0       0       0       0       0       0       0       0
    ## 195       0       0       0       0       0       0       0       0       0
    ## 196       0       0       0       0       0       0       0       0       0
    ## 197       0       0       0       0       0       0       0       0       0
    ## 198       0       0       0       0       0       0       0       0       0
    ## 199       0       0       0       0       0       0       0       0       0
    ## 200       0       0       0       0       0       0       0       0       0
    ## 201       0       0       0       0       0       0       0       0       0
    ## 202       0       0       0       0       0       0       0       0       0
    ## 203       0       0       0       0       0       0       0       0       0
    ## 204       0       0       0       0       0       0       0       0       0
    ## 205       0       0       0       0       0       0       0       0       0
    ## 206       0       0       0       0       0       0       0       0       0
    ## 207       0       0       0       0       0       0       0       0       0
    ## 208       0       0       0       0       0       0       0       0       0
    ## 209       0       0       0       0       0       0       0       0       0
    ## 210       0       0       0       0       0       0       0       0       0
    ## 211       0       0       0       0       0       0       0       0       0
    ## 212       0       0       0       0       0       0       0       0       0
    ## 213       0       0       0       0       0       0       0       0       0
    ## 214       0       0       0       0       0       0       0       0       0
    ## 215       0       0       0       0       0       0       0       0       0
    ## 216       0       0       0       0       0       0       0       0       0
    ## 217       0       0       0       0       0       0       0       0       0
    ## 218       0       0       0       0       0      11       0       0       0
    ## 219       0       0       0       0       0       0       0       0       0
    ## 220       0       0       0       0       0       0       0       0       0
    ## 221       0       0       0       0       0       0       0       0       0
    ## 222       0       0       0       0       0       0       0       0       0
    ## 223       0       0       0       0       0       0       0       0       0
    ## 224       0       0       0       0       0       0       0       0       0
    ## 225       0       0       0       0       0       0       0       0       0
    ##     otu_127 otu_128 otu_129 otu_130 otu_131 otu_132 otu_133 otu_134 otu_135
    ## 1         0       0       0       0       0       0       0       0       0
    ## 2         0       0       0       0       0       0       0       0       0
    ## 3         0       0       0       0       0       0       0       0       0
    ## 4         0       0       0       0       0       0       0       0       0
    ## 5         0       0      13       0       0       0       0       0       0
    ## 6         0       0       0       0       0       0       0       0       0
    ## 7         0       0       0       0       0       0       0       0       0
    ## 8         0       0       0       0       0       0       0       0       0
    ## 9         0       0       0       0       0       0       0       0       0
    ## 10        0       0       0       0       0       0       0       0       0
    ## 11        0       0       0       0       0       0       0       0       0
    ## 12        0       0       0       0       0       0       0       0       0
    ## 13        0       0       0       0       0       0       0       0       0
    ## 14        0       0       0       0       0       0       0       0       0
    ## 15        0       0       0       0       0       0       0       0       0
    ## 16        0       0       0       0       0       0       0       0       0
    ## 17        0       0       0       0       0       0       0       0       0
    ## 18        0       0       0       0       0       0       0       0       0
    ## 19        0       0       0       0       0       0       0       0       0
    ## 20        0       0       0       0       0       0       0       0       0
    ## 21        0       0       0       0       0       0       0       0       0
    ## 22        0       0       0       0       0       0       0       0       0
    ## 23        0       0       0       0       0       0       0       0       0
    ## 24        0       0       0       0       0       0       0       0       0
    ## 25        0       0       0       0       0       0       0       0       0
    ## 26        0       0       0       0       0       0       0       0       0
    ## 27        0       0       0       0       0       0       0       0       0
    ## 28        0       0       0       0       0       0       0       0       0
    ## 29        0       0       0       0       0       0       0       0       0
    ## 30        0       0       0       0       0       0       0       0       0
    ## 31        0       0       0       0       0       0       0       0       0
    ## 32        0       0       0       0       0       0       0       0       0
    ## 33        0       0       0       0       0       0       0       0       0
    ## 34       10       0       0       0       0       0       0       0       0
    ## 35        0       0       0       0       0       0       0       0       0
    ## 36        0       0       0       0       0       0       0       0       0
    ## 37        0       0       0       0       0       0       0       0       0
    ## 38        0       0       0       0       0       0       0       0       0
    ## 39        0       0       0       0       0       0       0       0       0
    ## 40        0       0       0       0       0       0       0       0       0
    ## 41        0       0       0       0       0       0       0       0       0
    ## 42        0       0       0       0       0       0       0       0       0
    ## 43        0       0       0       0       0       0       0       0       0
    ## 44        0       0       0       0       0       0       0       0       0
    ## 45        0       0       0       0       0       0       0       0       0
    ## 46        0       0       0       0       0       0       0       0       0
    ## 47       16       0       0       0       0       0       0       0       0
    ## 48        0       0       0       0       0       0       0       0       0
    ## 49        0       0       0       0       0       0       0       0       0
    ## 50        0       0       0       0       0       0       0       0       0
    ## 51        0       0       0       0       0       0       0       0       0
    ## 52        0       0       0       0       0       0       0       0       0
    ## 53        0       0       0       0       0       0       0       0       0
    ## 54        0       0       0       0       0       0       0       0       0
    ## 55        0       0       0       0       0       0       0       0       0
    ## 56        0       0       0       0       0       0       0       0       0
    ## 57        0       0       0       0       0       0       0       0       0
    ## 58        0       0       0       0       0       0       0       0       0
    ## 59        0       0       0       0       0       0       0       0       0
    ## 60        0       0       0       0       0       0       0       0       0
    ## 61        0       0       0       0       0       0       0       0       0
    ## 62        0       0       0       0       0       0       0       0       0
    ## 63        0       0       0       0       0       0       0       0       0
    ## 64        0       0       0       0       0       0       0       0       0
    ## 65        0       0       0       0       0       0       0       0       0
    ## 66        0       0       0       0       0       0       0       0       0
    ## 67        0       0       0       0       0       0       0       0       0
    ## 68        0       0       0       0       0       0       0       0       0
    ## 69        0       0       0       0       0       0       0       0       0
    ## 70        0       0       0       0       0       0       0       0       0
    ## 71        0       0       0       0       0       0       0       0       0
    ## 72        0       0       0       0       0       0       0       0       0
    ## 73        0       0       0       0       0       0       0       0       0
    ## 74        0       0       0       0       0       0       0       0       0
    ## 75        0       0       5       0       0       0       0       0       0
    ## 76        0       0       0       0       0       0       0       0       0
    ## 77        0       0       0       0       0       0       0       0       0
    ## 78        0       0       0       0       0       0       4       0       0
    ## 79        0       0       0       0       0       0       0       0       0
    ## 80        0       0       0       0       0       0       0       0       0
    ## 81        0       0       0       0       0       0       0       0       0
    ## 82        0       0       0       0       0       0       0       0       0
    ## 83        0       0       0       0       0       0       2       0       0
    ## 84        0       0       0       0       0       0       0       0       0
    ## 85        0       0       0       0       0       0       0       0       0
    ## 86        0       0       0       0       0       0       0       0       0
    ## 87        0       0       0       0       0       0       0       0       0
    ## 88        0       0       0       0       0       0       0       0       0
    ## 89        0       0       0       0       0       0       0       0       0
    ## 90        0       0       0       0       0       0       0       0       0
    ## 91        0       0       0       0       0       0       0       0       0
    ## 92        0       0       0       0       0       0       0       0       0
    ## 93        0       0       0       0       0       0       0       0       0
    ## 94        0       0       0       0       0       0       0       0       0
    ## 95        0       0       0       0       0       0       0       0       0
    ## 96        0       0       0      18       0       0       0       0       0
    ## 97        0       0       0       7       0       0       0       0       0
    ## 98        0       0       0       0       0       0       0       0       0
    ## 99        0       0       0       0       0       0       0       0       0
    ## 100       0       0       0       0       0       0       0       0       0
    ## 101       0       0       0       0       0       0       0       0       0
    ## 102       0       0       0       0       0       0       0       0       0
    ## 103       0       0       0       0       0       0       0       0       0
    ## 104       0       0       0       0       0       0       0       0       0
    ## 105       0       0       0       0       0       0       0       0       0
    ## 106       0       0       0       0       0       0       0       0       0
    ## 107       0       0       0       0       0       0       0       0       0
    ## 108       0       0       0       0       0       0       0       0       0
    ## 109       0       0       0       0       0       0       0       0       0
    ## 110       0       0       0       0       0       0       0       0       0
    ## 111       0       0       0       0       0       0       0       0       0
    ## 112       0       0       0       0       0       0       0       0       0
    ## 113       0       0       0       0       0       0       0       0       0
    ## 114       0       0       0       0       0       0       0       0       0
    ## 115       0       0       0       0       0       0       0       0       3
    ## 116       0       0       0       0       0       0       0       0       0
    ## 117       0       0       0       0       0       0       0       0       0
    ## 118       0       0       0       0       0       0       0       0       0
    ## 119       0       0       0       0       0       0       0       0       0
    ## 120       0       0       0       0       0       0       0       0       0
    ## 121       0       0       0       0       0       0       0       0      15
    ## 122       0       0       0       0       0       0       0       0       0
    ## 123       0       0       0       0       0       0       0       0       0
    ## 124       0       0       0       0       0       0       0       0       0
    ## 125       0       0       0       0       0       0       0       0       0
    ## 126       0       0       0       0       0       0       0       0       0
    ## 127       0       0       0       0       0       0       0       0       0
    ## 128       0       0       0       0       0       0       0       0       0
    ## 129       0       0       0       0       0       0       0       0       0
    ## 130       0       0       0       0       0       0       0       0       0
    ## 131       0       0       0       0       0       0       0       0       0
    ## 132       0       0       0       0       0       0       0       0       0
    ## 133       0       0       0       0       0       0       0       0       0
    ## 134       0       0       0       0       0       0       0       0       0
    ## 135       0       0       0       0       0       0       0       0       0
    ## 136       0       0       0       0      17      17       0       0       0
    ## 137       0       0       0       0       0       0       0       0       0
    ## 138       0       0       0       0       0       0       0       0       0
    ## 139       0       0       0       0       0       0       0       0       0
    ## 140       0       0       0       0       0       0       0       0       0
    ## 141       0       0       0       0       0       0       0       0       0
    ## 142       0       0       0       0       0       0       0       0       0
    ## 143       0       0       0       0       0       0       0      16       0
    ## 144       0       0       0       0       0       0       0       0       0
    ## 145       0       0       0       0       0       0       0       0       0
    ## 146       0       0       0       0       0       0       0       0       0
    ## 147       0       0       0       0       0       0       0       0       0
    ## 148       0       0       0       0       0       0       0       0       0
    ## 149       0       0       0       0       0       0       0       0       0
    ## 150       0       0       0       0       0       0       0       0       0
    ## 151       0       0       0       0       0       0       0       0       0
    ## 152       0       0       0       0       0       0       0       0       0
    ## 153       0       0       0       0       0       0       0       0       0
    ## 154       0       0       0       0       0       0       0       0       0
    ## 155       0       0       0       0       0       0       0       0       0
    ## 156       0       0       0       0       0       0       0       0       0
    ## 157       0       0       0       0       0       0       0       0       0
    ## 158       0       0       0       0       0       0       0       0       0
    ## 159       0       0       0       0       0       0       0       0       0
    ## 160       0       0       0       0       0       0       7       0       0
    ## 161       0       0       0       0       0       0       0       0       0
    ## 162       0       0       0       0       0       0       0       0       0
    ## 163       0       0       0       0       0       0       0       0       0
    ## 164       0       0       0       0       0       0       0       0       0
    ## 165       0       0       0       0       0       0       0       0       0
    ## 166       0       0       0       0       0       0       0       0       0
    ## 167       0       0       0       0       0       0       0       0       0
    ## 168       0       0       0       0       0       0       0       0       0
    ## 169       0       0       0       0       0       0       0       0       0
    ## 170       0       0       0       0       0       0       0       0       0
    ## 171       0       0       0       0       0       0       0       0       0
    ## 172       0       0       0       0       0       0       0       0       0
    ## 173       0       0       0       0       0       0       0       0       0
    ## 174       0       0       0       0       0       0       0       0       0
    ## 175       0       0       0       0       0       0       0       0       0
    ## 176       0       0       0       0       0       0       3       0       0
    ## 177       0       0       0       0       0       0       0       0       0
    ## 178       0       0       0       0       0       0       0       0       0
    ## 179       0       0       0       0       0       0       0       0       0
    ## 180       0       0       0       0       0       0       0       0       0
    ## 181       0       0       0       0       0       0       0       0       0
    ## 182       0       0       0       0       0       0       0       0       0
    ## 183       0       0       0       0       0       0       0       0       0
    ## 184       0       0       0       0       0       0       0       0       0
    ## 185       0       0       0       0       0       0       0       0       0
    ## 186       0       0       0       0       0       0       0       0       0
    ## 187       0       0       0       0       0       0       0       0       0
    ## 188       0       0       0       0       0       0       0       0       0
    ## 189       0       0       0       0       0       0       0       0       0
    ## 190       0       0       0       0       0       0       0       0       0
    ## 191       0       0       0       0       0       0       0       0       0
    ## 192       0      26       0       0       0       0       0       0       0
    ## 193       0       0       0       0       0       0       0       0       0
    ## 194       0       0       0       0       0       0       0       0       0
    ## 195       0       0       0       0       0       0       0       0       0
    ## 196       0       0       0       0       0       0       0       0       0
    ## 197       0       0       0       0       0       0       0       0       0
    ## 198       0       0       0       0       0       0       0       0       0
    ## 199       0       0       0       0       0       0       0       0       0
    ## 200       0       0       0       0       0       0       0       0       0
    ## 201       0       0       0       0       0       0       0       0       0
    ## 202       0       0       0       0       0       0       0       0       0
    ## 203       0       0       0       0       0       0       0       0       0
    ## 204       0       0       0       0       0       0       0       0       0
    ## 205       0       0       0       0       0       0       0       0       0
    ## 206       0       0       0       0       0       0       0       0       0
    ## 207       0       0       0       0       0       0       0       0       0
    ## 208       0       0       0       0       0       0       0       0       0
    ## 209       0       0       0       0       0       0       0       0       0
    ## 210       0       0       0       0       0       0       0       0       0
    ## 211       0       0       0       0       0       0       0       0       0
    ## 212       0       0       0       0       0       0       0       0       0
    ## 213       0       0       0       0       0       0       0       0       0
    ## 214       0       0       0       0       0       0       0       0       0
    ## 215       0       0       0       0       0       0       0       0       0
    ## 216       0       0       0       0       0       0       0       0       0
    ## 217       0       0       0       0       0       0       0       0       0
    ## 218       0       0       0       0       0       0       0       0       0
    ## 219       0       0       0       0       0       0       0       0       0
    ## 220       0       0       0       0       0       0       0       0       0
    ## 221       0       0       0       0       0       0       0       0       0
    ## 222       0       0       0       0       0       0       0       0       0
    ## 223       0       0       0       0       0       0       0       0       0
    ## 224       0       0       0       0       0       0       0       0       0
    ## 225       0       0       0       0       0       0       0       0       0
    ##     otu_136 otu_137 otu_138 otu_139 otu_140 otu_141 otu_142 otu_143 otu_144
    ## 1         0       0       0       0       0       0       0       0       0
    ## 2         0       0       0       0       0       0       0       0       0
    ## 3         0       0       0       0       0       0       0       0       0
    ## 4         0       0       0       0       0       0       0       0       0
    ## 5         0       0       0       0       0       0       0       0       0
    ## 6         0       0       0       0       0       0       0       0       0
    ## 7         0       0       0       0       0       0       0       0       0
    ## 8         0       0       0       0       0       0       0       7       0
    ## 9         0       0       0       0       0       0       0       0       0
    ## 10        0       0       0       0       0       0       0       0       0
    ## 11        0       0       0       0       0       0       0       0       0
    ## 12        0       0       0       0       0       0       0       0       0
    ## 13        0       0       0       0       0       0       0       0       0
    ## 14        0       0       0       0       0       0       0       0       0
    ## 15        0       0       0       0       0       0       0       0       0
    ## 16        0       0       0       0       0       0       0       0       0
    ## 17        0       0       0       0       0       0       0       0       0
    ## 18        0       0       0       0       0       0       0       0       0
    ## 19        0       0       0       0       0       0       0       0       0
    ## 20        0       0       0       0       0       0       0       0       0
    ## 21        0       0       0       0       0       0       0       0       0
    ## 22        0       0       0       0       0       0       0       0       0
    ## 23        0       0       0       0       0       0       0       0       0
    ## 24        0       0       0       0       0       0       0       0       0
    ## 25        0       0       0       0       0       0       0       0       0
    ## 26        0       0       0       0       0       0       0       0       0
    ## 27        0       0       0       0       0       0       0       0       0
    ## 28        0       0       0       0       0       0       0       0       0
    ## 29        0       0       0       0       0       0       0       0       0
    ## 30        0       0       0       0       0       0       0       0       0
    ## 31        0       0       0       0       0       0       0       0       0
    ## 32        0       0       0       0       0       0       0       0       0
    ## 33        0       0       0       0       0       0       0       0       0
    ## 34        0       0       0       0       0       0       0       0       0
    ## 35        0       0       0       0       0       0       0       0       0
    ## 36        0       0       0       0       0       0       0       0       0
    ## 37        0       0       0       0       0       0       0       0       0
    ## 38        0       0       0       0       0       0       0       0       0
    ## 39        0       0       0       0       0       0       0       0       0
    ## 40        0       0       0       0       0       0       0       0       0
    ## 41        0       0       0       0       0       0       0       0       0
    ## 42        0       0       0       0       0       0       0       0       0
    ## 43        0       0       0       0       0       0       0       0       0
    ## 44        0       0       0       0       0       0       0       0       0
    ## 45        0       0       0       0       0       0       0       0       0
    ## 46        0       0       0       0       0       0       0       0       0
    ## 47        0       0       0       0       0       0       0       0       0
    ## 48        0       0       0       0       0       0       0       0       0
    ## 49        0       0       0       0       0       0       0       0       0
    ## 50        0       0       0       0       0       0       0       0       0
    ## 51        0       0       0       0       0       0       0       0       0
    ## 52        0       0       0       0       0       0       0       0       0
    ## 53        0       0       0       0       0       0       0       0       0
    ## 54        0       0       0       0       0       0       0       0       0
    ## 55        0       0       0       0       0       0       0       0       0
    ## 56        0       0       0       0       0       0       0       0       0
    ## 57        0       0       0       0       0       0       0       0       0
    ## 58        0       0       0       0       0       0       0       0       0
    ## 59        0       0       0       0       0       0       0       0       0
    ## 60        0       0       0       0       0       0       0       0       0
    ## 61        0       0       0       0       0       0       0       0       0
    ## 62        0       0       0       0       0       0       0       0       0
    ## 63        0       0       0       0       0       0       0       0       0
    ## 64        0       0       0       0       0       0       0       0       0
    ## 65        0       0       0       0       0       0       0       0       0
    ## 66        0       0       0       0       0       0       0       0       0
    ## 67        0       0       0       0       0       0       0       0       0
    ## 68        0       0       0       0       0       0       0       0       0
    ## 69        0       0       0       0       0       0       0       0       0
    ## 70        0       0       0       0       0       0       0       0       0
    ## 71        0       0       0      12       0       0       0       0       0
    ## 72        0       0       0       0       0       0       0       0       0
    ## 73        0       0       0       0       0       0       0       0       0
    ## 74        0       0       0       0       0       0       0       0       0
    ## 75        0       0       0       0       0       0       0       0       0
    ## 76        0       0       0       0       0       0       0       0       0
    ## 77        0       0       0       0       0       0       0       0       0
    ## 78        0       0       0       0       0       0       0       0       0
    ## 79        0       0       0       0       0       0       0       0       0
    ## 80        0       0       0       0       0       0       0       0       0
    ## 81        0       0       0       0       0       0       0       0       0
    ## 82        0       0       0       0       0       0       0       0       0
    ## 83        0       0       0       0       0       0       0       0       0
    ## 84        0       0       0       0       0       0       0       0       0
    ## 85        0       0       0       0       0       0       0       0       0
    ## 86        0       0       0       0       0       0       0       0       0
    ## 87        0      13       0       0       0       0       0       0       0
    ## 88        0       0       0       0       0       0       0       0       0
    ## 89        0       0       0       0       0       0       0       0       0
    ## 90        0       0       0       0       0       0       0       0       0
    ## 91        0       0       0       0       0       0       0       0       0
    ## 92        0       0       0       0       0       0       0       0       0
    ## 93        0       0       0       0       0       0       0       0       0
    ## 94        0       0       0       0       0       0       0       0       0
    ## 95        0       0       0       0       0       0       0       0       0
    ## 96        0       0       0       0       0       0       0       0       0
    ## 97        0       0       0       0       0       0       0       0       0
    ## 98        0       0       0       0       0       0       0       0       0
    ## 99        0       0       0       0       0       0       0       0       0
    ## 100       0       0       0       0       0       0       0       0       0
    ## 101       0       0       0       0       0       0       0       0       0
    ## 102       0       0       0       0       0       0       0       0       0
    ## 103       0       0       0       0       0       0       0       0       0
    ## 104       0       0       0       0       0       0       0       0       0
    ## 105       0       0       0       0       0       0       0       0       0
    ## 106       0       0       0       0       0       0       0       0       0
    ## 107       0       0       0       0       0       0       0       0       0
    ## 108       0       0       0       0       0       0       0       0       0
    ## 109       0       0       0       0       0       0       0       0       0
    ## 110       0       0       0       0       0       0       0       0       0
    ## 111       0       0       0       0       0       0       0       0       0
    ## 112       0       0       0       0       0       0       8       0       0
    ## 113       0       0       0       0       0       0       0       0       0
    ## 114       0       0      13       0       0       0       0       0       0
    ## 115       0       0       0       0       0       0       0       0       0
    ## 116       0       0       0       0       0       0       0       0       0
    ## 117       0       0       0       0       0       0       0       0       0
    ## 118       0       0       0       0       0       0       0       0       0
    ## 119       0       0       0       0       0       0       0       0       0
    ## 120       0       6       0       0      11       0       0       0       0
    ## 121       0       0       0       0       0       0       0       0       0
    ## 122       0       0       0       0       0       0       0       0       0
    ## 123       0       0       0       0       0       0       0       0       0
    ## 124       0       0       0       0       0       0       0       0       0
    ## 125       0       0       0       0       0       0       0       0       0
    ## 126       0       0       0       0       0       0       0       0       0
    ## 127       5       0       0       0       0       0       0       0       0
    ## 128       0       0       0       0       0       0       0       0       0
    ## 129       0       0       0       0       0       0       0       0       0
    ## 130       0       0       0       0       0       0       0       0       0
    ## 131       0       0       0       0       0       0       0       0       0
    ## 132       0       0       0       0       0       0       0       0       0
    ## 133       0       0       0       0       0       0       0       0       0
    ## 134      10       0       0       0       0       0       0       0       0
    ## 135       0       0       0       0       0       0       0       0       0
    ## 136       0       0       0       0       0       0       0       0       0
    ## 137       0       0       0       0       0       0       0       0       0
    ## 138       0       0       0       0       0       0       0       0       0
    ## 139       0       0       0       0       0       0       0       0       0
    ## 140       0       0       0       0       0       9       0       0       0
    ## 141       0       0       0       0       0       0       0       0       0
    ## 142       0       0       0       0       0       0       0       0       0
    ## 143       0       0       0       0       0       0       0       0       0
    ## 144       0       0       0       0       0       0       0       0       0
    ## 145       0       0       0       0       0       0       0       0       0
    ## 146       0       0       0       0       0       0       0       0       0
    ## 147       0       0       0       0       0       0       0       0       0
    ## 148       0       0       0       0       0       0       0       0       0
    ## 149       0       0       0       0       0       0       0       0       0
    ## 150       0       0       0       0       0       0       0       0       0
    ## 151       0       0       0       0       0       0       0       0       0
    ## 152       0       0       0       0       0       0       0       0       0
    ## 153       0       0       0       0       0       0       0       0       0
    ## 154       0       0       0       0       0       0       0       0       0
    ## 155       0       0       0       0       0       0       0       0       0
    ## 156       0       0       0       0       0       0       0       0       0
    ## 157       0       0       0       0       0       0       0       0       0
    ## 158       0       0       0       0       0       0       0       0       0
    ## 159       0       0       0       0       0       0       0       0       0
    ## 160       0       0       0       0       0       0       0       0       0
    ## 161       0       0       0       0       0       0       0       0       0
    ## 162       0       0       0       0       0       0       0       0       0
    ## 163       0       0       0       0       0       0       0       0       0
    ## 164       0       0       0       0       0       0       0       0       0
    ## 165       0       0       0       0       0       0       0       0       0
    ## 166       0       0       0       0       0       0       0       0       0
    ## 167       0       0       0       0       0       0       0       0       0
    ## 168       0       0       0       0       0       0       0       0       0
    ## 169       0       0       0       0       0       0       0       0       0
    ## 170       0       0       0       0       0       0       0       0       0
    ## 171       0       0       0       0       0       0       0       0       0
    ## 172       0       0       0       0       0       0       0       0       0
    ## 173       0       0       0       0       0       0       0       0       0
    ## 174       0       0       0       0       0       0       0       0       0
    ## 175       0       0       0       0       0       0       0       0       0
    ## 176       0       0       0       0       0       0       0       0       0
    ## 177       0       0       0       0       0       0       0       0       0
    ## 178       0       0       0       0       0       0       0       0       0
    ## 179       0       0       0       0       0       0       0       0       0
    ## 180       0       0       0       0       0       0       0       0       0
    ## 181       0       0       0       0       0       0       0       0       0
    ## 182       0       0       0       0       0       0       0       0       0
    ## 183       0       0       0       0       0       0       0       0       0
    ## 184       0       0       0       0       0       0       0       0       0
    ## 185       0       0       0       0       0       0       0       0       0
    ## 186       0       0       0       0       0       0       0       0       0
    ## 187       0       0       0       0       0       0       0       0       0
    ## 188       0       0       0       0       0       0       0       0       0
    ## 189       0       0       0       0       0       0       0       0       0
    ## 190       0       0       0       0       0       0       0       0       0
    ## 191       0       0       0       0       0       0       0       0       0
    ## 192       0       0       0       0       0       0       0       0       0
    ## 193       0       0       0       0       0       0       0       0       0
    ## 194       0       0       0       0       0       0       0       0       0
    ## 195       0       0       0       0       0       0       0       0       0
    ## 196       0       0       0       0       0       0       0       0       6
    ## 197       0       0       0       0       0       0       0       0       0
    ## 198       0       0       0       0       0       0       0       0       0
    ## 199       0       0       0       0       0       0       0       0       0
    ## 200       0       0       0       0       0       0       0       0       0
    ## 201       0       0       0       0       0       0       0       0       0
    ## 202       0       0       0       0       0       0       0       0       0
    ## 203       0       0       0       0       0       0       0       0       0
    ## 204       0       0       0       0       0       0       0       0       0
    ## 205       0       0       0       0       0       0       0       0       0
    ## 206       0       0       0       0       0       0       0       0       0
    ## 207       0       0       0       0       0       0       0       0       0
    ## 208       0       0       0       0       0       0       0       0       0
    ## 209       0       0       0       0       0       0       0       0       0
    ## 210       0       0       0       0       0       0       0       0       0
    ## 211       0       0       0       0       0       0       0       0       0
    ## 212       0       0       0       0       0       0       0       0       0
    ## 213       0       0       0       0       0       0       0       0       0
    ## 214       0       0       0       0       0       0       0       0       0
    ## 215       0       0       0       0       0       0       0       0       0
    ## 216       0       0       0       0       0       0       0       0       0
    ## 217       0       0       0       0       0       0       0       0       0
    ## 218       0       0       0       0       0       0       0       0       0
    ## 219       0       0       0       0       0       0       0       0       0
    ## 220       0       0       0       0       0       0       0       0       0
    ## 221       0       0       0       0       0       0       0       0       0
    ## 222       0       0       0       0       0       0       0       0       0
    ## 223       0       0       0       0       0       0       0       0       0
    ## 224       0       0       0       0       0       0       0       0       0
    ## 225       0       0       0       0       0       0       0       0       0
    ##     otu_145 otu_146 otu_147 otu_148 otu_149 otu_150 otu_151 otu_152
    ## 1         0       0       0       0       0       0       0       0
    ## 2         0       0       0       0       0       0       0       0
    ## 3         0       0       0       0       0       0       0       0
    ## 4         0       0       0       0       0       0       0       0
    ## 5         0       0       0       0       0       0       0       0
    ## 6         0       0       0       0       0       0       0       0
    ## 7         0       0       0       0       0       0       0       0
    ## 8         0       0       0       0       0       0       0       0
    ## 9         0       0       0       0       0       0       0       0
    ## 10        0       0       0       0       0       0       0       0
    ## 11        0       0       0       0       0       0       0       0
    ## 12        0       0       0       0       0       0       0       0
    ## 13        0       0       0       0       0       0       0       0
    ## 14        0       0       0       0       0       0       0       0
    ## 15        0       0       0       0       0       0       0       0
    ## 16        0       0       0       0       0       0       0       0
    ## 17        0       0       0       0       0       0       0       0
    ## 18        0       0       0       0       0       0       0       0
    ## 19        0       0       0       0       0       0       0       0
    ## 20        0       0       0       0       0       0       0       0
    ## 21        0       0       0       0       0       0       0       0
    ## 22        0       0       0       0       0       0       0       0
    ## 23        0       0       0       0       0       0       0       0
    ## 24        0       0       0       0       0       0       0       0
    ## 25        0       0       0       0       0       0       0       0
    ## 26        0       0       0       0       0       0       0       0
    ## 27        0       0       0       0       0       0       0       0
    ## 28        0       0       0       0       0       0       0       0
    ## 29        0       0       0       0       0       0       0       0
    ## 30        0       0       0       0       0       0       0       0
    ## 31        0       0       0       0       0       0       0       0
    ## 32        0       0       0       0       0       0       0       0
    ## 33        0       0       0       0       0       0       0       0
    ## 34        0       0       0       0       0       0       0       0
    ## 35        0       0       0       0       0       0       0       0
    ## 36        0       0       0       0       2       0       0       0
    ## 37        0       0       0       0       0       0       0       0
    ## 38        0       0       0       0       0       0       0       0
    ## 39        0       0       0       0       0       0       0       0
    ## 40        0       0       0       0       0       0       0       0
    ## 41        0       0       0       0       0       0       0       0
    ## 42        0       0       0       0       0       0       0       0
    ## 43        0       0       0       0       0       0       0       0
    ## 44        0       0       0       0       0       0       0       0
    ## 45        0       0       0       0       0       0       0       0
    ## 46        0       0       0       0       0       0       0       0
    ## 47        0       0       0       0       0       0       0       0
    ## 48        0       0       0       0       0       0       0       0
    ## 49        0       0       0       0       0       0       0       0
    ## 50        0       0       0       0       0       0       0       0
    ## 51        0       0       0       0       0       0       0       0
    ## 52        0       0       0       0       0       0       0       0
    ## 53        0       0       0       0       0       0       0       0
    ## 54        0       0       0       0       0       0       0       0
    ## 55        0       0       0       0       0       0       0       0
    ## 56        0       0       0       0       0       0       0       0
    ## 57        0       0       0       0       0       0       0       0
    ## 58        0       0       0       0       0       0       0       0
    ## 59        0       0       0       0       0       0       0       0
    ## 60        0       0       0       0       0       0       0       0
    ## 61        0       0       0       0       0       0       0       0
    ## 62        0       0       0       0       0       0       0       0
    ## 63        0       0       0       0       0       0       0       0
    ## 64        0       0       0       0       0       0       0       0
    ## 65        0       0       0       0       0       0       0       0
    ## 66        0       0       0       0       0       0       0       0
    ## 67        0       0       0       0       0       0       0       0
    ## 68        0       0       0       0       0       0       0       0
    ## 69        0       0       0       0       0       0       0       0
    ## 70        0       0       0       0       0       0       0       0
    ## 71        0       0       0       0       0       0       0       0
    ## 72        0       0       0       0       0       0       0       0
    ## 73        0       0       0       0       0       0       0       0
    ## 74        0       0       0       0       0       0       0       0
    ## 75        0       0       0       0       0       0       0       0
    ## 76        0       0       0       0       0       0       0       0
    ## 77        0       0       0       0       0       0       0       0
    ## 78        0       0       0       0       0       2       0       0
    ## 79        0       0       0       0       0       0       0       0
    ## 80        0       0       0       0       0       0       0       0
    ## 81        0       0       0       0       0       0       0       0
    ## 82        0       0       0       0       0       0       0       0
    ## 83        0       0       0       0       0       0       0       0
    ## 84        0       0       0       0       0       0       0       0
    ## 85        0       0       0       0       0       0       0       0
    ## 86        0       0       0       0       0       0       0       0
    ## 87        0       0       0       0       0       0       0       0
    ## 88        0       0       0       0       0       0       0       0
    ## 89        0       0       0       0       0       0       0       0
    ## 90        0       0       0       0       0       0       0       0
    ## 91        0       0       0       0       0       0       0       0
    ## 92        0       0       0       0       0       0       0       0
    ## 93        0       0       0       0       0       0       0       0
    ## 94        0       0       0       0       0       0       0       0
    ## 95        0       0       0       0       0       0       0       0
    ## 96        0       0       0       0       0       0       0       0
    ## 97        0       0       0       0       0       0       0       0
    ## 98        0       0       0       0       0       0       0       0
    ## 99        0       0       4       0       0       0       0       0
    ## 100       0       0       0       0       0       0       0       0
    ## 101       0       0       0       0       0       0       0       0
    ## 102       0       0       0       0       0       0       0       0
    ## 103       0       0       0       0       0       0       0       0
    ## 104       0       0       0       0       0       0       0       0
    ## 105       0       0       0       0       0       0       0       0
    ## 106       0       0       0       0       0       0       0       0
    ## 107       0       0       0       0       0       0       0       0
    ## 108       0       0       0       0       0       0       0       0
    ## 109       0       0       0       0       0       0       0       0
    ## 110       0       0       0       0       0       0       0       0
    ## 111       0       0       0       0       0       0       0       0
    ## 112       5       0       0       0       0       0       0       0
    ## 113       0       0       0       0       0       0       0       0
    ## 114       0       0       0       0       0       0       0       0
    ## 115       0       0       0       0       0       0       0       0
    ## 116       0       0       0       0       0       0       0       0
    ## 117       0       0       0       0       0       0       0       0
    ## 118       0       0       0       0       0       0       0       0
    ## 119       0       0       0       0       0       0       0       0
    ## 120       0       0       0       0       0       0       0       0
    ## 121       0       0       0       0       0       0       0       0
    ## 122       0       0       0       0       0       0       0       0
    ## 123       0       0       0       0       0       0       0       0
    ## 124       0       0       0       0       0       0       0       0
    ## 125       0       0       0       0       0       0       0       0
    ## 126       0       0       0       0       0       0       0       0
    ## 127       0       0       0       0       0       0       0       0
    ## 128       0       0       0       0       0       0       0       0
    ## 129       0       0       0       0       0       0       0       0
    ## 130       0       0       0       0       0       0       0       0
    ## 131       0       0       0       0       0       0       0       0
    ## 132       0       0       0       0       0       0       0       0
    ## 133       0       0       0       0       0       0       0       0
    ## 134       0       0       0       0       0       0       0       0
    ## 135       0       0       0       0       0       0       0       0
    ## 136       0       0       0       0       0       0       0       0
    ## 137       0       0       0       0       0       0       0       0
    ## 138       0       0       0       0       0       0       0       0
    ## 139       0       0       0       0       0       0       0       0
    ## 140       0       0       0       0       0       0       0       0
    ## 141       0       0       0       0       0       0       0       0
    ## 142       0       0       0       0       0       0       0       0
    ## 143       0       0       0       0       0       0       0       0
    ## 144       0       0       0       0       0       0       0       0
    ## 145       0       0       0       0       0       0       0       0
    ## 146       0       0       0       0       0       0       0       0
    ## 147       0       0       0       0       0       0       0       0
    ## 148       0       0       0       0       0       0       0       0
    ## 149       0       0       0       0       0       0       0       0
    ## 150       0       0       0       0       0       0       0       0
    ## 151       0       0       0       0       0       0       0       0
    ## 152       0       0       0       0       0       0       0       0
    ## 153       0       0       0       0       0       0       0       0
    ## 154       0       0       0       0       0       0       0       0
    ## 155       0       0       0       0       0       0       0       0
    ## 156       0       0       0       0       0       0       0       0
    ## 157       0       0       0       0       0       0       0       0
    ## 158       0       0       0       0       0       0       0       0
    ## 159       0       5       0       0       0       0       0       0
    ## 160       0       0       0       0       0       0       0       0
    ## 161       0       0       0       0       0       0       0       0
    ## 162       0       0       0       0       0       0       0       0
    ## 163       0       0       0       0       0       0       0       0
    ## 164       0       0       0       0       0       0       0       0
    ## 165       0       0       0       0       0       0       0       0
    ## 166       0       0       0       0       0       0       0       0
    ## 167       0       0       0       0       0       0       0       0
    ## 168       0       0       0       0       0       0       0       0
    ## 169       0       0       0       0       0       0       0       0
    ## 170       0       0       0       0       0       0       0       0
    ## 171       0       0       0       0       0       0       0       0
    ## 172       0       0       0       0       0       0       0       0
    ## 173       0       0       0       0       0       0       0       0
    ## 174       0       0       0       0       0       0       0       0
    ## 175       0       0       0       0       0       0       0       0
    ## 176       0       0       0       0       0       0       2       0
    ## 177       0       0       0       0       0       0       0       0
    ## 178       0       0       0       0       0       0       0       0
    ## 179       0       0       0       0       0       0       0       0
    ## 180       0       0       0       0       0       0       0       0
    ## 181       0       0       0       0       0       0       0       0
    ## 182       0       0       0       0       0       0       0       0
    ## 183       0       0       0       0       0       0       0       0
    ## 184       0       0       0       3       0       0       0       0
    ## 185       0       0       0       0       0       0       0       0
    ## 186       0       0       0       0       0       0       0       0
    ## 187       0       0       0       0       0       0       0       0
    ## 188       0       0       0       0       0       0       0       0
    ## 189       0       0       0       0       0       0       0       0
    ## 190       0       0       0       0       0       0       0       0
    ## 191       0       0       0       0       0       0       0       0
    ## 192       0       0       0       0       0       0       0       0
    ## 193       0       0       0       0       0       0       0       0
    ## 194       0       0       0       0       0       0       0       0
    ## 195       0       0       0       0       0       0       0       0
    ## 196       0       0       0       0       0       0       0       0
    ## 197       0       0       0       0       0       0       0       0
    ## 198       0       0       0       0       0       0       0       0
    ## 199       0       0       0       0       0       0       0       0
    ## 200       0       0       0       0       0       0       0       0
    ## 201       0       0       0       0       0       0       0       0
    ## 202       0       0       0       0       0       0       0       0
    ## 203       0       0       0       0       0       0       0       0
    ## 204       0       0       0       0       0       0       0       0
    ## 205       0       0       0       0       0       0       0       0
    ## 206       0       0       0       0       0       0       0       0
    ## 207       0       0       0       0       0       0       0       0
    ## 208       0       0       0       0       0       0       0       0
    ## 209       0       0       0       0       0       0       0       0
    ## 210       0       0       0       0       0       0       0       0
    ## 211       0       0       0       0       0       0       0       2
    ## 212       0       0       0       0       0       0       0       0
    ## 213       0       0       0       0       0       0       0       0
    ## 214       0       0       0       0       0       0       0       0
    ## 215       0       0       0       0       0       0       0       0
    ## 216       0       0       0       0       0       0       0       0
    ## 217       0       0       0       0       0       0       0       0
    ## 218       0       0       0       0       0       0       0       0
    ## 219       0       0       0       0       0       0       0       0
    ## 220       0       0       0       0       0       0       0       0
    ## 221       0       0       0       0       0       0       0       0
    ## 222       0       0       0       0       0       0       0       0
    ## 223       0       0       0       0       0       0       0       0
    ## 224       0       0       0       0       0       0       0       0
    ## 225       0       0       0       0       0       0       0       0
    ## 
    ## $depth_spe_samps_rfy
    ## [1] 163
    ## 
    ## $spe_samps_rfy
    ## # A tibble: 225 × 139
    ##    field_key sample otu_1 otu_2 otu_3 otu_4 otu_5 otu_6 otu_7 otu_8 otu_9 otu_10
    ##        <dbl>  <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>  <dbl>
    ##  1         1      1    10    20     6     0    20    11     7    10     0      1
    ##  2         1      2     2     4    13    13    39     8     9     4     1      0
    ##  3         1      4     1    13     3     0    19     1    11    38     0      0
    ##  4         1      5     0     0    25     1    26    13    23     2     0      0
    ##  5         1      6     0     0    44     0    39     0     0     0     0      0
    ##  6         1      7     2    10    21     0    38    11     2    14     0      0
    ##  7         1      8     6    27     9     0     5     2     3    41     3      7
    ##  8         1      9    10     6    34     0     8     0    19    16     0      0
    ##  9         1     10     8    22     9     0    39     1     0     5     0      0
    ## 10         2      1    30     0    79     0     6     7     0     0     4      1
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
    ## # A tibble: 25 × 152
    ##    field_key otu_1 otu_2 otu_3 otu_4 otu_5 otu_6 otu_7 otu_8 otu_9 otu_10 otu_11
    ##        <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>  <dbl>  <dbl>
    ##  1         1   497  1453  1487   259  3162   685   957  2099    72    149   1297
    ##  2         2  2160    32  5002     0  2365   650   295     0   681    687      0
    ##  3         3    21   217    20  2541  1015   427  1154     0   393     13      0
    ##  4         4  2013  2306     0   639   676   701   369  1203  2833   1083    310
    ##  5         5   475  1397    16    21  2018   753   446     6  1197   1656     61
    ##  6         6   390    71  3415  3001   861  1088  1595     0   359     79      1
    ##  7         7    12     0  1221  1181  4056   874   136    70   110      0      0
    ##  8         8  1171  1506   239   560   947   913  1148   719   581   1117   1814
    ##  9         9   703  1653   129  1723   297   363   966  2733  1027    940   2864
    ## 10        10  1543  1106   119  1738   192   502  1272  3555   180    535   1365
    ## # ℹ 15 more rows
    ## # ℹ 140 more variables: otu_12 <dbl>, otu_13 <dbl>, otu_14 <dbl>, otu_15 <dbl>,
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
