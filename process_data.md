Database assembly: species data
================
Beau Larkin

Last updated: 24 February, 2023

- <a href="#description" id="toc-description">Description</a>
  - <a href="#its-data-all-fungi" id="toc-its-data-all-fungi">ITS data (all
    fungi)</a>
  - <a href="#18s-data-mycorrhizae" id="toc-18s-data-mycorrhizae">18S data
    (mycorrhizae)</a>
  - <a href="#desired-outcome" id="toc-desired-outcome">Desired outcome</a>
- <a href="#resources" id="toc-resources">Resources</a>
  - <a href="#packages-and-libraries"
    id="toc-packages-and-libraries">Packages and libraries</a>
  - <a href="#functions" id="toc-functions">Functions</a>
- <a href="#load-and-process-data" id="toc-load-and-process-data">Load and
  process data</a>
  - <a href="#import-files" id="toc-import-files">Import files</a>
  - <a href="#etl-using-etl" id="toc-etl-using-etl">ETL using
    <code>etl()</code></a>
  - <a href="#post-processing-18s-data"
    id="toc-post-processing-18s-data">Post-processing 18S data</a>

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
matrices. The six samples from each field with the greatest total
sequence abundance will be chosen, and sequences summed for each OTU
within fields. Rarefaction of sequencing depth to the minimum total
number of sequences will be applied.

For each taxonomy table, taxonomy strings must be parsed and unnecessary
characters removed. A function is used to streamline the pipeline and to
reduce errors. [Fungal
traits](https://link.springer.com/article/10.1007/s13225-020-00466-2)
data will be joined with the ITS taxonomy

For all tables, short and unique rownames must be created to allow for
easy joining of species and metadata tables.

For the 18S data, a second table is needed to produce a UNIFRAC distance
matrix. The table must have OTUs in rows with OTU ids.

# Resources

## Packages and libraries

``` r
packages_needed = c("tidyverse", "GUniFrac")
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

    spe_topn <-
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
        separate_wider_delim(cols = rowname, delim = "_", names = c("field_key", NA)) %>%
        select(field_key, everything()) %>%
        mutate(sum = rowSums(across(starts_with(cluster_type)))) %>%
        group_by(field_key) %>%
        slice_max(sum, n=samps) %>%
        ungroup() %>% 
        select(-sum) %>%
        mutate(field_key = as.numeric(field_key))
    zero_otu1 <- which(apply(spe_topn, 2, sum) == 0)
    spe_raw <- spe_topn[, -zero_otu1]
    
    spe_sum <-
        data.frame(
            spe_raw %>%
                group_by(field_key) %>%
                summarize(across(starts_with("otu"), ~ sum(.x)), .groups = "drop"),
            row.names = 1
        )
    depth <- min(rowSums(spe_sum))
    rfy <- Rarefy(spe_sum)
    zero_otu2 <- which(apply(rfy$otu.tab.rff, 2, sum) == 0)
    spe_rfy <- data.frame(rfy$otu.tab.rff[, -zero_otu2]) %>%
        rownames_to_column(var = "field_key") %>%
        mutate(field_key = as.numeric(field_key)) %>% 
        as_tibble()
    
    meta_raw <- meta %>% filter(!(otu_num %in% names(zero_otu1)))
    meta_rfy <- meta_raw %>% filter(!(otu_num %in% names(zero_otu2)))
    
    write_csv(meta_raw, paste0(getwd(), folder, "/spe_", gene, "_raw_taxonomy.csv"))
    write_csv(spe_raw, paste0(getwd(), folder, "/spe_", gene, "_raw.csv"))
    write_csv(meta_rfy, paste0(getwd(), folder, "/spe_", gene, "_rfy_taxonomy.csv"))
    write_csv(spe_rfy, paste0(getwd(), folder, "/spe_", gene, "_rfy.csv"))
    
    out <- list(
        spe_raw_meta = meta_raw,
        spe_raw      = spe_raw,
        spe_rfy_meta = meta_rfy,
        spe_rfy      = spe_rfy,
        depth_rfy    = depth
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
        samps = 6,
        traits = traits,
        varname = otu_num,
        gene = "ITS",
        cluster_type = "otu",
        colname_prefix = "ITS_TGP_",
        folder = "/clean_data"
    )
its
```

    ## $spe_raw_meta
    ## # A tibble: 2,795 × 9
    ##    otu_num otu_ID                phylum class order family genus species prima…¹
    ##    <chr>   <chr>                 <chr>  <chr> <chr> <chr>  <chr> <chr>   <chr>  
    ##  1 otu_1   352d386293a59777de3e… Ascom… Sord… Hypo… Nectr… Fusa… Fusari… plant_…
    ##  2 otu_2   a78342f18e3db5624f3d… Morti… Mort… Mort… Morti… Mort… Mortie… soil_s…
    ##  3 otu_3   dabfbac17ada76765b6d… Ascom… Sord… Glom… Plect… Gibe… <NA>    plant_…
    ##  4 otu_4   bdee5c3cbce33c88a5e8… Ascom… Euro… Chae… Herpo… unid… uniden… <NA>   
    ##  5 otu_5   73514f6e23fa430c5127… Ascom… Sord… Hypo… Nectr… <NA>  <NA>    <NA>   
    ##  6 otu_6   c92d481c3f8ce88e480e… Ascom… Doth… Pleo… Didym… <NA>  <NA>    <NA>   
    ##  7 otu_7   204f1bd97debac4d9c56… Ascom… Doth… Pleo… Peric… Peri… <NA>    plant_…
    ##  8 otu_8   0ab6be0adca17efdd24e… Ascom… Euro… Chae… Herpo… unid… uniden… <NA>   
    ##  9 otu_9   3c7865fa957956fc9c7f… Basid… Trem… Cyst… Mraki… Taus… Tauson… soil_s…
    ## 10 otu_10  caa87147e44034b05364… Ascom… <NA>  <NA>  <NA>   <NA>  <NA>    <NA>   
    ## # … with 2,785 more rows, and abbreviated variable name ¹​primary_lifestyle
    ## 
    ## $spe_raw
    ## # A tibble: 150 × 2,796
    ##    field_key otu_1 otu_2 otu_3 otu_4 otu_5 otu_6 otu_7 otu_8 otu_9 otu_10 otu_11
    ##        <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>  <dbl>  <dbl>
    ##  1         1   170   241     8     0    45   469  1502     0     0    108      0
    ##  2         1   153     0     4     0    27   203   498   274     0     98      0
    ##  3         1   117     0    13  2286    24   152   862     0     0     33      0
    ##  4         1   102   842     3     0    19    80   239     0     0     31      0
    ##  5         1    85   137    15    54    33   118   294     0     0    118      0
    ##  6         1   416    11    25     6    48   225   112     0     0    154      0
    ##  7        10   124   137    60     0    52   172   230     0     0    117      0
    ##  8        10   100   277    41     0    19   187    48     0     0      0      0
    ##  9        10   268    48    42     0    26   105   159  1023     0      0      0
    ## 10        10   243   186    40   230     0   166   262     0     0     91      0
    ## # … with 140 more rows, and 2,784 more variables: otu_12 <dbl>, otu_13 <dbl>,
    ## #   otu_14 <dbl>, otu_15 <dbl>, otu_16 <dbl>, otu_17 <dbl>, otu_18 <dbl>,
    ## #   otu_19 <dbl>, otu_20 <dbl>, otu_21 <dbl>, otu_22 <dbl>, otu_23 <dbl>,
    ## #   otu_24 <dbl>, otu_25 <dbl>, otu_26 <dbl>, otu_27 <dbl>, otu_28 <dbl>,
    ## #   otu_29 <dbl>, otu_30 <dbl>, otu_31 <dbl>, otu_32 <dbl>, otu_33 <dbl>,
    ## #   otu_34 <dbl>, otu_35 <dbl>, otu_36 <dbl>, otu_37 <dbl>, otu_38 <dbl>,
    ## #   otu_39 <dbl>, otu_40 <dbl>, otu_41 <dbl>, otu_42 <dbl>, otu_43 <dbl>, …
    ## 
    ## $spe_rfy_meta
    ## # A tibble: 2,789 × 9
    ##    otu_num otu_ID                phylum class order family genus species prima…¹
    ##    <chr>   <chr>                 <chr>  <chr> <chr> <chr>  <chr> <chr>   <chr>  
    ##  1 otu_1   352d386293a59777de3e… Ascom… Sord… Hypo… Nectr… Fusa… Fusari… plant_…
    ##  2 otu_2   a78342f18e3db5624f3d… Morti… Mort… Mort… Morti… Mort… Mortie… soil_s…
    ##  3 otu_3   dabfbac17ada76765b6d… Ascom… Sord… Glom… Plect… Gibe… <NA>    plant_…
    ##  4 otu_4   bdee5c3cbce33c88a5e8… Ascom… Euro… Chae… Herpo… unid… uniden… <NA>   
    ##  5 otu_5   73514f6e23fa430c5127… Ascom… Sord… Hypo… Nectr… <NA>  <NA>    <NA>   
    ##  6 otu_6   c92d481c3f8ce88e480e… Ascom… Doth… Pleo… Didym… <NA>  <NA>    <NA>   
    ##  7 otu_7   204f1bd97debac4d9c56… Ascom… Doth… Pleo… Peric… Peri… <NA>    plant_…
    ##  8 otu_8   0ab6be0adca17efdd24e… Ascom… Euro… Chae… Herpo… unid… uniden… <NA>   
    ##  9 otu_9   3c7865fa957956fc9c7f… Basid… Trem… Cyst… Mraki… Taus… Tauson… soil_s…
    ## 10 otu_10  caa87147e44034b05364… Ascom… <NA>  <NA>  <NA>   <NA>  <NA>    <NA>   
    ## # … with 2,779 more rows, and abbreviated variable name ¹​primary_lifestyle
    ## 
    ## $spe_rfy
    ## # A tibble: 25 × 2,790
    ##    field_key otu_1 otu_2 otu_3 otu_4 otu_5 otu_6 otu_7 otu_8 otu_9 otu_10 otu_11
    ##        <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>  <dbl>  <dbl>
    ##  1         1  1043  1231    68  2346   196  1247  3507   274     0    542      0
    ##  2         2  2889  1714   940    80   606    61   446     0     0   1790     81
    ##  3         3   351     7   309     0   540   271   161     0    27      0    427
    ##  4         4  1607     0     5     0  1395  2325  1491     0     0    820      0
    ##  5         5   441   102   422     0  1182  2772  4472     0    41   1011    120
    ##  6         6  1019   305   176     0   474   503   961     0  2815      0   1953
    ##  7         7  1036  1090  1030     8   599  2532    85     0  7377      0    264
    ##  8         8  1442   507   602  1134   510   912  1142     0    29   1187     90
    ##  9         9  1531  1029   790     0   152  1016   703    19     0    446     56
    ## 10        10   702   706   815   429   186   700   945  1420     0    210     28
    ## # … with 15 more rows, and 2,778 more variables: otu_12 <dbl>, otu_13 <dbl>,
    ## #   otu_14 <dbl>, otu_15 <dbl>, otu_16 <dbl>, otu_17 <dbl>, otu_18 <dbl>,
    ## #   otu_19 <dbl>, otu_20 <dbl>, otu_21 <dbl>, otu_22 <dbl>, otu_23 <dbl>,
    ## #   otu_24 <dbl>, otu_25 <dbl>, otu_26 <dbl>, otu_27 <dbl>, otu_28 <dbl>,
    ## #   otu_29 <dbl>, otu_30 <dbl>, otu_31 <dbl>, otu_32 <dbl>, otu_33 <dbl>,
    ## #   otu_34 <dbl>, otu_35 <dbl>, otu_36 <dbl>, otu_37 <dbl>, otu_38 <dbl>,
    ## #   otu_39 <dbl>, otu_40 <dbl>, otu_41 <dbl>, otu_42 <dbl>, otu_43 <dbl>, …
    ## 
    ## $depth_rfy
    ## [1] 45478

``` r
amf <- 
    etl(
        spe = amf_otu,
        taxa = amf_taxa,
        samps = 6,
        varname = otu_num,
        gene = "18S",
        cluster_type = "otu",
        colname_prefix = "X18S_TGP_",
        folder = "/clean_data"
    )
amf
```

    ## $spe_raw_meta
    ## # A tibble: 147 × 8
    ##    otu_num otu_ID                         class order family genus taxon acces…¹
    ##    <chr>   <chr>                          <chr> <chr> <chr>  <chr> <chr> <chr>  
    ##  1 otu_1   320f3edc7b48ba5691766ccc71b0d… Glom… Glom… Glome… Glom… Glom… VTX002…
    ##  2 otu_2   97cefc055a2fa1b8caf5b81635b3f… Glom… Glom… Glome… Glom… Glom… VTX001…
    ##  3 otu_3   2c47b9acb976f06611bddfe10f97e… Para… Para… Parag… Para… <NA>  <NA>   
    ##  4 otu_4   7487be824e0acea7ee8c22b94f137… Glom… Glom… Glome… Glom… Glom… VTX001…
    ##  5 otu_5   e8cbf85e2f78be8ff21c0c3e0f509… Glom… Glom… Glome… Glom… <NA>  <NA>   
    ##  6 otu_6   87ff74cc0f55f12dc4056a4a8a2c7… Glom… Glom… Claro… Clar… <NA>  <NA>   
    ##  7 otu_7   f10a50b1220da2d3e275f9772b3c2… Glom… Glom… Glome… Glom… Glom… VTX002…
    ##  8 otu_8   0d508e08bf3048aa561e7c9d96e3b… Glom… Glom… Glome… Glom… Glom… VTX003…
    ##  9 otu_9   4a4251fdd8c94240584c5d4ff6aeb… Glom… Glom… Glome… Glom… Glom… VTX002…
    ## 10 otu_10  210c71717f87c8e4912431ed98b55… Glom… Glom… Claro… Clar… Clar… VTX000…
    ## # … with 137 more rows, and abbreviated variable name ¹​accession
    ## 
    ## $spe_raw
    ## # A tibble: 150 × 148
    ##    field_key otu_1 otu_2 otu_3 otu_4 otu_5 otu_6 otu_7 otu_8 otu_9 otu_10 otu_11
    ##        <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>  <dbl>  <dbl>
    ##  1         1   102   195   378   267  1006   275   120   171    61      0      0
    ##  2         1    17   279    30     0   395    16   184   894     0      0   1021
    ##  3         1     0     0   486    30   566   223   550    33     0     16      0
    ##  4         1   130   505   129     0   113    26    55   758    26    110      0
    ##  5         1   144   359   193     0   648    30    16   151     0     28    361
    ##  6         1    55   111   275     0   516   110    21   148     0      0    103
    ##  7        10   745    21    94  1238    39   182   130     0    35    190    153
    ##  8        10    92   165    33    67    16    55  1037   625     0     57      0
    ##  9        10   235   195     0   101    14   125   118   948     0     48    887
    ## 10        10   399   256    38     0    10   190    59   599     0    321      0
    ## # … with 140 more rows, and 136 more variables: otu_12 <dbl>, otu_13 <dbl>,
    ## #   otu_14 <dbl>, otu_15 <dbl>, otu_16 <dbl>, otu_17 <dbl>, otu_18 <dbl>,
    ## #   otu_19 <dbl>, otu_20 <dbl>, otu_21 <dbl>, otu_22 <dbl>, otu_23 <dbl>,
    ## #   otu_24 <dbl>, otu_25 <dbl>, otu_26 <dbl>, otu_27 <dbl>, otu_28 <dbl>,
    ## #   otu_29 <dbl>, otu_30 <dbl>, otu_31 <dbl>, otu_32 <dbl>, otu_33 <dbl>,
    ## #   otu_34 <dbl>, otu_35 <dbl>, otu_36 <dbl>, otu_37 <dbl>, otu_38 <dbl>,
    ## #   otu_39 <dbl>, otu_40 <dbl>, otu_41 <dbl>, otu_42 <dbl>, otu_43 <dbl>, …
    ## 
    ## $spe_rfy_meta
    ## # A tibble: 146 × 8
    ##    otu_num otu_ID                         class order family genus taxon acces…¹
    ##    <chr>   <chr>                          <chr> <chr> <chr>  <chr> <chr> <chr>  
    ##  1 otu_1   320f3edc7b48ba5691766ccc71b0d… Glom… Glom… Glome… Glom… Glom… VTX002…
    ##  2 otu_2   97cefc055a2fa1b8caf5b81635b3f… Glom… Glom… Glome… Glom… Glom… VTX001…
    ##  3 otu_3   2c47b9acb976f06611bddfe10f97e… Para… Para… Parag… Para… <NA>  <NA>   
    ##  4 otu_4   7487be824e0acea7ee8c22b94f137… Glom… Glom… Glome… Glom… Glom… VTX001…
    ##  5 otu_5   e8cbf85e2f78be8ff21c0c3e0f509… Glom… Glom… Glome… Glom… <NA>  <NA>   
    ##  6 otu_6   87ff74cc0f55f12dc4056a4a8a2c7… Glom… Glom… Claro… Clar… <NA>  <NA>   
    ##  7 otu_7   f10a50b1220da2d3e275f9772b3c2… Glom… Glom… Glome… Glom… Glom… VTX002…
    ##  8 otu_8   0d508e08bf3048aa561e7c9d96e3b… Glom… Glom… Glome… Glom… Glom… VTX003…
    ##  9 otu_9   4a4251fdd8c94240584c5d4ff6aeb… Glom… Glom… Glome… Glom… Glom… VTX002…
    ## 10 otu_10  210c71717f87c8e4912431ed98b55… Glom… Glom… Claro… Clar… Clar… VTX000…
    ## # … with 136 more rows, and abbreviated variable name ¹​accession
    ## 
    ## $spe_rfy
    ## # A tibble: 25 × 147
    ##    field_key otu_1 otu_2 otu_3 otu_4 otu_5 otu_6 otu_7 otu_8 otu_9 otu_10 otu_11
    ##        <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>  <dbl>  <dbl>
    ##  1         1   366  1233  1255   241  2732   562   800  1814    74    129   1265
    ##  2         2  1969    12  4298     0  1675   607   225     0   700    592      0
    ##  3         3    30   218    23  2208   856   360   941     0   382     18      0
    ##  4         4  1422  1949     0   666   626   687   238   930  2525    766    313
    ##  5         5   348  1362    16    32  1886   711   322    11  1047   1275     62
    ##  6         6   399    86  3329  2465   832   523  1473     2   376     62      4
    ##  7         7    11     0  1034   861  3716   193    37    83    47      0      0
    ##  8         8  1121  1056   201   316   836   753  1105   461   492   1192   1557
    ##  9         9   363  1501   130  1473   287   279   933  2203  1089    796   2570
    ## 10        10  1353   890   115  1184    74   484  1207  3196   165    522   1146
    ## # … with 15 more rows, and 135 more variables: otu_12 <dbl>, otu_13 <dbl>,
    ## #   otu_14 <dbl>, otu_15 <dbl>, otu_16 <dbl>, otu_17 <dbl>, otu_18 <dbl>,
    ## #   otu_19 <dbl>, otu_20 <dbl>, otu_21 <dbl>, otu_22 <dbl>, otu_23 <dbl>,
    ## #   otu_24 <dbl>, otu_25 <dbl>, otu_26 <dbl>, otu_27 <dbl>, otu_28 <dbl>,
    ## #   otu_29 <dbl>, otu_30 <dbl>, otu_31 <dbl>, otu_32 <dbl>, otu_33 <dbl>,
    ## #   otu_34 <dbl>, otu_35 <dbl>, otu_36 <dbl>, otu_37 <dbl>, otu_38 <dbl>,
    ## #   otu_39 <dbl>, otu_40 <dbl>, otu_41 <dbl>, otu_42 <dbl>, otu_43 <dbl>, …
    ## 
    ## $depth_rfy
    ## [1] 16611

## Post-processing 18S data

To produce a UNIFRAC distance table, the trimmed table `amf$spe_rfy`
must be imported back into QIIME, where sequence data can be used with
abundances to create a phylogeny and distance matrix. The data frame
must be transposed and use legal column names (i.e., non-numeric). Site
metadata is used to create better column names.

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
    left_join(amf$spe_rfy_meta %>% select(otu_num, otu_ID), by = join_by(otu_num)) %>% 
    select(otu_ID, everything(), -otu_num)
write_tsv(amf_export, paste0(getwd(), "/otu_tables/18S/spe_18S_rfy_export.tsv"))
```
