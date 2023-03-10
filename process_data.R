#' ---
#' title: "Database assembly: species data"
#' author: "Beau Larkin\n"
#' date: "Last updated: `r format(Sys.time(), '%d %B, %Y')`"
#' output:
#'   github_document:
#'     toc: true
#'     toc_depth: 3
#'     df_print: paged
#' ---
#'
#' # Description
#' Microbial sequence abundances were produced by Lorinda Bullington in QIIME2. ETL must
#' be performed on output text files to allow downstream analysis. 
#' 
#' ## ITS data (all fungi)
#' Sequence abundances in 97% similar OTUs in individual samples form the base data. 
#' The abundances are raw (not rarefied).  
#' 
#' - ITS taxonomy are included in a separate file. 
#' 
#' ## 18S data (mycorrhizae)
#' Sequence abundance in 97% similar OTUs in individual samples. The abundances are raw
#' (not rarefied). 
#' 
#' - 18S taxonomy are included in a separate file.
#' - A unifrac distance matrix will be created and included after sample selection and 
#' sequence depth rarefaction.
#' 
#' ## Desired outcome
#' For each raw table, species OTU codes must be aligned with short, unique keys, and then species
#' tables must be transposed into sites-species matrices. The six samples from each field with the 
#' greatest total sequence abundance will be chosen, and sequences summed for each OTU within fields.
#' Rarefaction of sequencing depth to the minimum total number of sequences will be applied.
#' 
#' For each taxonomy table, taxonomy strings must be parsed and unnecessary characters
#' removed. A function is used to streamline the pipeline and to reduce errors. 
#' [Fungal traits](https://link.springer.com/article/10.1007/s13225-020-00466-2) data will be joined
#' with the ITS taxonomy
#' 
#' For all tables, short and unique rownames must be created to allow for easy joining of species
#' and metadata tables. 
#' 
#' For all sequences, zero-abundance and singleton OTUs must be removed after rarefying.
#' 
#' For the 18S data, a second table is needed to produce a UNIFRAC distance matrix. The table
#' must have OTUs in rows with OTU ids. 
#' 
#' # Resources
#' ## Packages and libraries
packages_needed = c("tidyverse", "GUniFrac")
packages_installed = packages_needed %in% rownames(installed.packages())
#+ packages,message=FALSE
if (any(!packages_installed)) {
    install.packages(packages_needed[!packages_installed])
}
#+ libraries,message=FALSE
for (i in 1:length(packages_needed)) {
    library(packages_needed[i], character.only = T)
}
#' 
#' ## Functions
#' *NOTE:* all `write_csv()` steps have been commented out as of 2023-03-13 to prevent overwriting
#' existing files. This is because function `Rarefy()` produces inconsistent results. Due to rounding,
#' a very few OTUs are retained or lost (<1%) when this function is rerun. These different outcomes 
#' change nothing about how results would be interpreted, but they do change axis limits and other trivial
#' parameters that cause headaches later. 
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
    
    spe_topn <- 
        spe_t %>%
        mutate(sum = rowSums(across(starts_with(cluster_type)))) %>%
        group_by(field_key) %>% 
        slice_max(sum, n=samps) %>%
        select(-sum)
    
    spe_topn_sum <- 
        spe_topn %>% 
        group_by(field_key) %>% 
        summarize(across(starts_with(cluster_type), ~ sum(.x)), .groups = "drop") %>%
        mutate(field_key = as.numeric(field_key)) %>% 
        arrange(field_key)
    
    zero_otu <- which(apply(spe_topn_sum, 2, sum) == 0)
    if(length(zero_otu) == 0) {
        spe_raw <- spe_topn_sum
    } else {
        spe_raw <- spe_topn_sum[, -zero_otu]
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
    
    single_zero_otus <- which(apply(rfy$otu.tab.rff, 2, sum) <= 1)
    if(length(single_zero_otus) == 0) {
        spe_rfy <- data.frame(rfy$otu.tab.rff) %>%
            rownames_to_column(var = "field_key") %>%
            mutate(field_key = as.numeric(field_key)) %>% 
            arrange(field_key) %>% 
            as_tibble()
    } else {
        spe_rfy <- data.frame(rfy$otu.tab.rff[, -single_zero_otus]) %>%
            rownames_to_column(var = "field_key") %>%
            mutate(field_key = as.numeric(field_key)) %>% 
            arrange(field_key) %>% 
            as_tibble()
    }
    
    meta_raw <- meta %>% filter(!(otu_num %in% names(zero_otu)))
    meta_rfy <- meta_raw %>% filter(!(otu_num %in% names(single_zero_otus)))
    
    # Commented out 2023-03-13, see note above.
    # write_csv(spe_t, paste0(getwd(), folder, "/spe_", gene, "_raw_samps_all.csv"))
    # write_csv(spe_topn, paste0(getwd(), folder, "/spe_", gene, "_raw_samps_topn.csv"))
    # write_csv(meta_raw, paste0(getwd(), folder, "/spe_", gene, "_raw_taxonomy.csv"))
    # write_csv(spe_raw, paste0(getwd(), folder, "/spe_", gene, "_raw.csv"))
    # write_csv(meta_rfy, paste0(getwd(), folder, "/spe_", gene, "_rfy_taxonomy.csv"))
    # write_csv(spe_rfy, paste0(getwd(), folder, "/spe_", gene, "_rfy.csv"))
    
    out <- list(
        spe_raw_meta = meta_raw,
        spe_raw      = spe_raw,
        spe_rfy_meta = meta_rfy,
        spe_rfy      = spe_rfy,
        depth_rfy    = depth
    )
    
    return(out)
    
}
#'
#' # Load and process data
#' ## Import files
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
#' 
#' ## ETL using `etl()`
#' Schema: `process_qiime(spe, taxa, samps, traits=NULL, varname, gene, cluster_type, colname_prefix, folder)`
#+ otu_its,message=FALSE
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
#+ otu_18S,message=FALSE
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
#' 
#' ## Post-processing 18S data
#' To produce a UNIFRAC distance table, the trimmed table `amf$spe_rfy` must be 
#' imported back into QIIME, where sequence data can be used with abundances to 
#' create a phylogeny and distance matrix. The data frame must be transposed and 
#' use legal column names (i.e., non-numeric). 
#' Site metadata is used to create better column names.
#' 
#' The resultant distance matrix will be imported again when needed for multivariate
#' analysis and ordination. 
#+ import_sites,message=FALSE
sites <- read_csv(paste0(getwd(), "/clean_data/sites.csv"), show_col_types = FALSE)
#+ wrangle_amf_spe
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
