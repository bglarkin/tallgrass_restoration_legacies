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
#' The microbial species data included here were produced by Lorinda Bullington in 2022
#' using QIIME II. See methods for details. 
#' 
#' ## ITS data (all fungi)
#' Lorinda included two datasets as of 2022-01-06. The first is a table of sequence
#' variants (SVs), assigned based on 100% similarity of ITS sequences in each cluster. 
#' The second is a table of operational taxonomic units (OTUs) based on 97% sequence 
#' similarity. 
#' Each table also includes various other data and metadata, including taxonomy,
#' trophic guilds, and references. 
#' 
#' ## 18S data (mycorrhizae)
#' Created by Lorinda on 2022-02-14. As with the ITS data, files with 97% similar OTUs 
#' and 100% similar SVs were created. Distance matrices for each were created, both 
#' weighted and unweighted. UNIFRAC distance was used. UNIFRAC is different from BC and others in that it accounts for phylogenetic distance, 
#' which can be informative for 18s and 16S, but not so much for ITS. 
#' Weighted UNIFRAC considers abundances where as non-weighted is based on presence/absence. 
#' 
#' ## Desired outcome
#' For each raw table, the metadata must be separated from the species abundances.
#' Species OTU codes must be aligned with short, unique keys, and then species
#' tables must be transposed into sites-species matrices. Rownames must be cleaned 
#' to align with site metadata files. Taxonomy strings must be parsed and unnecessary characters
#' removed. A function is used to streamline the pipeline and to reduce errors.
#' 
#' # Resources
#' ## Packages and libraries
packages_needed = c("tidyverse", "readxl")
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
process_qiime <- function(data, varname, gene, cluster_type, colname_prefix, folder) {
    
    # Variable definitions
    # data            = Dataframe or tibble with Qiime and FunGuild output
    # varname         = An unique key will be created to replace the cumbersome OTU key 
    #                   which is produced by Qiime. Varname is the desired column
    #                   name for this new OTU key. Unquoted. 
    # gene            = Gene region targeted in sequencing, e.g.: "ITS", "18S", "16S". 
    #                   Quoted. Used to select() column names, so must match text in 
    #                   column names. Also used to create distinct file names.
    # cluster_type    = Clustering algorithm output, e.g.: "otu", "sv". Quoted.
    #                   Used to create filenames in output files. 
    # colname_prefix  = Prefix to text of OTU column names. The function removes
    #                   this prefix to make OTU names more concise. 
    # folder          = The function creates output files in the working directory by
    #                   default. To use a subfolder, use this variable. Quoted. 
    #                   Include the "/" before the folder name. If no folder 
    #                   name is desired, use "".
    
    varname <- enquo(varname)
    
    if(gene == "ITS") {
        meta <-
            data %>%
            mutate(!!varname := paste0(cluster_type, "_", row_number())) %>% 
            select(!starts_with(gene)) %>%
            rename(
                otu_ID = `#OTU ID`,
                taxon = Taxon,
                taxon_level = `Taxon Level`,
                trophic_mode = `Trophic Mode`,
                guild = Guild,
                growth_morphology = `Growth Morphology`,
                trait = Trait,
                confidence = `Confidence Ranking`,
                notes = Notes,
                citation = `Citation/Source`
            ) %>% 
            select(!!varname, everything()) %>% 
            separate(
                taxonomy,
                into = c(
                    "kingdom",
                    "phylum",
                    "class",
                    "order",
                    "family",
                    "genus",
                    "species"
                ),
                sep = "; ",
                remove = TRUE,
                fill = "right"
            ) %>% 
            mutate(kingdom = str_sub(kingdom, 4, nchar(kingdom)),
                   phylum  = str_sub(phylum,  4, nchar(phylum)),
                   class   = str_sub(class,   4, nchar(class)),
                   order   = str_sub(order,   4, nchar(order)),
                   family  = str_sub(family,  4, nchar(family)),
                   genus   = str_sub(genus,   4, nchar(genus)),
                   species = str_sub(species, 4, nchar(species)))
        write_csv(meta, 
                  paste0(getwd(), folder, "/spe_", gene, "_", cluster_type, "_funGuild.csv"))
    } else {
        meta <-
            data %>%
            mutate(!!varname := paste0(cluster_type, "_", row_number())) %>% 
            select(!starts_with(gene)) %>%
            rename(otu_ID = `#OTU ID`) %>% 
            select(!!varname, everything()) %>% 
            separate(taxonomy, 
                     c("class", "order", "family", "genus", "taxon", "accession"), 
                     sep = ";", remove = TRUE, fill = "right")
        write_csv(meta, 
                  paste0(getwd(), folder, "/spe_", gene, "_", cluster_type, "_taxonomy.csv"))
    }
    
    spe_all <- 
        data.frame(
            data %>% 
                mutate(!!varname := paste0(cluster_type, "_", row_number())) %>% 
                select(!!varname, starts_with(gene)),
            row.names = 1
        ) %>%
        t() %>% 
        as.data.frame() %>% 
        rownames_to_column() %>% 
        mutate(site = str_remove(rowname, colname_prefix)) %>% 
        separate(col = site, into = c("site_key", "sample"), sep = "_") %>% 
        select(site_key, sample, everything(), -rowname) %>% 
        arrange(as.numeric(site_key), as.numeric(sample))
    write_csv(spe_all, 
              paste0(getwd(), folder, "/spe_", gene, "_", cluster_type, "_siteSpeMatrix_allReps.csv"))
    
    spe_count <- 
        spe_all %>% 
        group_by(site_key) %>% 
        summarize(samples = n()) %>% 
        arrange(as.numeric(site_key))
    write_csv(spe_count, 
              paste0(getwd(), folder, "/spe_", gene, "_", cluster_type, "_samples.csv"))
}
#'
#' # Load and process data
#' ## Import files
otu_its <- read_excel(paste0(getwd(), "/otu_tables/ITS/OTU_table_rrfd_3200_w_taxa.guilds.xlsx"), na = "-")
sv_its  <- read_excel(paste0(getwd(), "/otu_tables/ITS/SV_table_rrfd_3200.guilds.xlsx"), na = "-")
otu_18s <- read_delim(paste0(getwd(), "/otu_tables/TGP_18S_tables_021722/OTUs_18S_TGP_table_rarefied.txt"), 
                      delim = "\t", show_col_types = FALSE)
sv_18s  <- read_delim(paste0(getwd(), "/otu_tables/TGP_18S_tables_021722/SVs_18S_TGP_rarefied_table.txt"), 
                      delim = "\t", show_col_types = FALSE)
#' 
#' ## ETL using `process_qiime()`
#' Schema: `process_qiime(data, varname, "gene", "cluster_type", "colname_prefix", "folder")`
#+ otu_its
process_qiime(otu_its, otu_num, "ITS", "otu", "ITS_TGP_",  "/clean_data")
#+ sv_its
process_qiime(sv_its,  sv_num,  "ITS", "sv",  "ITS_TGP_",  "/clean_data")
#+ otu_18S
process_qiime(otu_18s, otu_num, "18S", "otu", "X18S_TGP_", "/clean_data")
#+ sv_18S
process_qiime(sv_18s,  sv_num,  "18S", "sv",  "X18S_TGP_", "/clean_data")

