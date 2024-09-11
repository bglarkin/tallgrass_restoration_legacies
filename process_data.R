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
#' be performed on output text files to allow downstream analysis. Iterative steps are needed
#' to find out the optimal number of samples to keep from each field to ensure equal 
#' sampling effort and adequate representation of diversity. 
#' 
#' ## Workflow
#' 1. The script `process_data.R` is run first. A few samples failed to amplify, resulting in
#' some fields characterized by 9 samples and others by 10. Also, one ITS sample was ambiguously assigned
#' to site and was removed. To balance sampling effort
#' across fields, the top 9 samples by sequence abundance are chosen from each field. 
#' **Assign "pre" to the argument `process_step` in the etl function so that files are created in the /clean_data/pre/... directory.**
#' 1. Next, `microbial_diagnostics_pre.R` is run to investigate sequencing depth in samples 
#' and species accumulation in fields. A few samples are known to have low sequence abundance 
#' (an order of magnitude lower than the maximum), and the consequence of rarefying to this 
#' small depth must be known. A new cutoff for sequence depth, and definition of further 
#' samples which must be cut, is recommended. 
#' 1. Then, `process_data.R` is run again, this time with the number of samples retained
#' per field set to the levels recommended in `microbial_diagnostics_pre.R`. 
#' As of 2023-10-11, the recommended number of samples to keep from all fields is **8 from the ITS dataset** and 
#' **7 from the 18S dataset.**
#'    - If downstream analyses have been completed, then it's likely that the `process_data.R`
#'    script has been left at this step. 
#'    - **Assign "post" to the argument `process_step` in the etl function so that files are created in the /clean_data/... directory.**
#' 1. Finally, `microbial_diagnostics_post.R` is run. It is very similar to the "_pre" script,
#' but a different file is used so that the two may be compared. 
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
#' tables must be transposed into sites-species matrices. Some fields only retained nine samples. 
#' To correct for survey effort, the nine samples from each field with the 
#' greatest total sequence abundance will be chosen, and sequences summed for each OTU within fields.
#' Rarefaction of sequencing depth to the minimum total number of sequences will be applied to summed OTUs.
#' 
#' For each taxonomy table, taxonomy strings must be parsed and unnecessary characters
#' removed. A function is used to streamline the pipeline and to reduce errors. 
#' [Fungal traits](https://link.springer.com/article/10.1007/s13225-020-00466-2) data will be joined
#' with the ITS taxonomy
#' 
#' For all tables, short and unique rownames must be created to allow for easy joining of species
#' and metadata tables. 
#' 
#' For all sequences, zero-abundance and singleton OTUs must be removed after OTUs have been 
#' rarefied and summed within fields.
#' 
#' For the 18S data, a second table is needed to produce a UNIFRAC distance matrix. The table
#' must have OTUs in rows with OTU ids. 
#' 
#' # Resources
#' ## Packages and libraries
packages_needed = c("tidyverse", "vegan", "knitr")
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
#' *NOTE:* Due to rounding,
#' a very few OTUs are retained or lost (<1%) when function `Rarefy()` is rerun. These different outcomes 
#' change nothing about how results would be interpreted, but they do change axis limits and other trivial
#' parameters that cause headaches later. Be advised that downstream changes will be needed if `Rarefy()` 
#' is rerun and new files are created by `write_csv()` in the steps at the end of this function.
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
# Site metadata
#+ import_sites,message=FALSE
sites    <- read_csv(paste0(getwd(), "/clean_data/sites.csv"), show_col_types = FALSE)
#' 
#' ## ETL using `etl()`
#' Schema: `process_qiime(spe, taxa, samps, traits=NULL, varname, gene, cluster_type, colname_prefix, folder)`
#' **Note:** If doing the first step of this workflow, retain 9 samples per field from each dataset and proceed
#' to further diagnostics in `microbial_diagnostics_pre.R`.
#' **Note:** If doing the third step of this workflow, Retain 8 samples from ITS and 7 samples from 18S
#' per field based on results from `microbial_diagnostics_pre.R`. Proceed to `microbial_diagnostics_post.R` for 
#' final exploration of the datasets. 
#+ otu_its,message=FALSE
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
        folder = "/clean_data",
        process_step = "pre"
    )
its
#+ otu_18S,message=FALSE,warning=FALSE
amf <-
    etl(
        spe = amf_otu,
        taxa = amf_taxa,
        samps = 9,
        varname = otu_num,
        gene = "18S",
        cluster_type = "otu",
        colname_prefix = "X18S_TGP_",
        folder = "/clean_data",
        process_step = "pre"
    )
amf
#' 
#' ## Post-processing 18S data
#' To produce a UNIFRAC distance table, the trimmed table `amf$spe_rfy` must be 
#' imported back into QIIME, where sequence data can be used with abundances to 
#' create a phylogeny and distance matrix. The data frame must be transposed and 
#' use legal column names (i.e., non-numeric). 
#' 
#' Site metadata is used to create better column names.
#' 
#' The resultant distance matrix will be imported again when needed for multivariate
#' analysis and ordination. 
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
    left_join(amf$spe_meta %>% select(otu_num, otu_ID), by = join_by(otu_num)) %>%
    select(otu_ID, everything(), -otu_num)
write_tsv(amf_export, paste0(getwd(), "/otu_tables/18S/spe_18S_rfy_export.tsv"))

