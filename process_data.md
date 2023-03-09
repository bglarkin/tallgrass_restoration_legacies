Database assembly: species data
================
Beau Larkin

Last updated: 09 March, 2023

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
        select(-sum) %>%
        summarize(across(starts_with("otu"), ~ sum(.x)), .groups = "drop") %>%
        mutate(field_key = as.numeric(field_key)) %>% 
        arrange(field_key)
    
    zero_otu1 <- which(apply(spe_topn, 2, sum) == 0)
    if(length(zero_otu1) == 0) {
        spe_raw <- spe_topn
    } else {
        spe_raw <- spe_topn[, -zero_otu1]
    }
    
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
    if(length(zero_otu2) == 0) {
        spe_rfy <- data.frame(rfy$otu.tab.rff) %>%
            rownames_to_column(var = "field_key") %>%
            mutate(field_key = as.numeric(field_key)) %>% 
            arrange(field_key) %>% 
            as_tibble()
    } else {
        spe_rfy <- data.frame(rfy$otu.tab.rff[, -zero_otu2]) %>%
            rownames_to_column(var = "field_key") %>%
            mutate(field_key = as.numeric(field_key)) %>% 
            arrange(field_key) %>% 
            as_tibble()
    }
    
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
    ## # A tibble: 25 × 2,796
    ##    field_key otu_1 otu_2 otu_3 otu_4 otu_5 otu_6 otu_7 otu_8 otu_9 otu_10 otu_11
    ##        <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>  <dbl>  <dbl>
    ##  1         1  1043  1231    68  2346   196  1247  3507   274     0    542      0
    ##  2         2  3256  1921  1062    92   695    65   509     0     0   2039     90
    ##  3         3   404     7   361     0   628   312   193     0    30      0    497
    ##  4         4  1695     0     5     0  1471  2439  1562     0     0    860      0
    ##  5         5   542   123   527     0  1435  3396  5575     0    50   1248    147
    ##  6         6  1369   409   235     0   635   689  1312     0  3778      0   2675
    ##  7         7  1235  1324  1244     9   697  2970    98     0  8845      0    312
    ##  8         8  1760   601   744  1395   615  1118  1398     0    38   1469    108
    ##  9         9  2199  1453  1181     0   220  1503  1054    27     0    611     79
    ## 10        10   842   848   962   521   214   834  1114  1721     0    256     34
    ## # … with 15 more rows, and 2,784 more variables: otu_12 <dbl>, otu_13 <dbl>,
    ## #   otu_14 <dbl>, otu_15 <dbl>, otu_16 <dbl>, otu_17 <dbl>, otu_18 <dbl>,
    ## #   otu_19 <dbl>, otu_20 <dbl>, otu_21 <dbl>, otu_22 <dbl>, otu_23 <dbl>,
    ## #   otu_24 <dbl>, otu_25 <dbl>, otu_26 <dbl>, otu_27 <dbl>, otu_28 <dbl>,
    ## #   otu_29 <dbl>, otu_30 <dbl>, otu_31 <dbl>, otu_32 <dbl>, otu_33 <dbl>,
    ## #   otu_34 <dbl>, otu_35 <dbl>, otu_36 <dbl>, otu_37 <dbl>, otu_38 <dbl>,
    ## #   otu_39 <dbl>, otu_40 <dbl>, otu_41 <dbl>, otu_42 <dbl>, otu_43 <dbl>, …
    ## 
    ## $spe_rfy_meta
    ## # A tibble: 2,785 × 9
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
    ## # … with 2,775 more rows, and abbreviated variable name ¹​primary_lifestyle
    ## 
    ## $spe_rfy
    ## # A tibble: 25 × 2,786
    ##    field_key otu_1 otu_2 otu_3 otu_4 otu_5 otu_6 otu_7 otu_8 otu_9 otu_10 otu_11
    ##        <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>  <dbl>  <dbl>
    ##  1         1  1043  1231    68  2346   196  1247  3507   274     0    542      0
    ##  2         2  2898  1690   930    78   632    56   449     0     0   1810     82
    ##  3         3   356     7   305     0   541   276   165     0    27      0    432
    ##  4         4  1633     0     5     0  1412  2311  1471     0     0    826      0
    ##  5         5   452   102   415     0  1163  2765  4511     0    41   1014    125
    ##  6         6  1012   289   178     0   484   519   990     0  2768      0   2000
    ##  7         7  1032  1092  1069     8   596  2493    87     0  7476      0    252
    ##  8         8  1450   487   618  1134   496   903  1141     0    33   1192    100
    ##  9         9  1505  1016   812     0   151  1035   713    23     0    435     57
    ## 10        10   680   723   811   450   180   692   931  1471     0    217     31
    ## # … with 15 more rows, and 2,774 more variables: otu_12 <dbl>, otu_13 <dbl>,
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
    ## # A tibble: 25 × 148
    ##    field_key otu_1 otu_2 otu_3 otu_4 otu_5 otu_6 otu_7 otu_8 otu_9 otu_10 otu_11
    ##        <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>  <dbl>  <dbl>
    ##  1         1   448  1449  1491   297  3244   680   946  2155    87    154   1485
    ##  2         2  4071    25  8546     0  3449  1276   462     0  1400   1187      0
    ##  3         3    44   345    31  3448  1341   555  1489     0   570     25      0
    ##  4         4  2493  3480     0  1223  1109  1222   416  1616  4457   1313    565
    ##  5         5   726  2861    27    58  3854  1447   664    18  2196   2621    123
    ##  6         6   725   125  6271  4679  1467   965  2782     5   725    149      6
    ##  7         7    16     0  1547  1294  5483   283    67   112    67      0      0
    ##  8         8  2082  2001   425   624  1657  1503  2228   919   941   2343   3033
    ##  9         9   629  2585   211  2490   487   452  1593  3743  1811   1380   4372
    ## 10        10  1829  1192   165  1596    96   649  1627  4344   228    730   1551
    ## # … with 15 more rows, and 136 more variables: otu_12 <dbl>, otu_13 <dbl>,
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
    ##  1         1   375  1237  1228   254  2725   572   812  1824    71    129   1241
    ##  2         2  1979    16  4251     0  1673   605   232     0   672    599      0
    ##  3         3    29   227    19  2222   887   369   924     0   367     15      0
    ##  4         4  1416  2004     0   703   639   713   223   894  2569    735    319
    ##  5         5   367  1374    14    32  1886   693   315     7  1070   1298     68
    ##  6         6   411    72  3349  2513   804   501  1500     2   376     74      2
    ##  7         7     8     0  1056   873  3655   195    51    76    44      0      0
    ##  8         8  1097  1048   201   303   824   797  1119   466   477   1217   1526
    ##  9         9   355  1501   126  1455   264   280   909  2159  1107    837   2555
    ## 10        10  1333   891   113  1200    70   474  1197  3242   169    542   1182
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
