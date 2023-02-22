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
#' similarity. Based on previous work, we've decided to use OTS-based data exclusively 
#' in this project.
#' Each table also includes various other data and metadata, including taxonomy,
#' trophic guilds, and references. 
#' 
#' For later examinations of sequence abundances in guilds, we will need to rarefy abundances
#' within guild subsets. The raw (un-rarefied) sequence data are loaded here for later use. 
#' 
#' ## 18S data (mycorrhizae)
#' Created by Lorinda on 2022-02-14. As with the ITS data, files with 97% similar OTUs 
#' and 100% similar SVs were created. Based on previous work, we've decided to use OTS-based data exclusively 
#' in this project. 
#' 
#' Weighted and unweighted UNIFRAC distance matrices were created. 
#' UNIFRAC is different from Bray-Curtis 
#' and others in that it accounts for phylogenetic distance, 
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
#' UNIFRAC distances must be coerced to distance objects.
#' 
#' The Fungal Traits data needs basic ETL for ease of later use.
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
process_qiime <- function(data, traits=NULL, varname, gene, cluster_type, append="", colname_prefix, folder) {
    
    # Variable definitions
    # data            = Dataframe or tibble with Qiime and FunGuild output
    # traits          = Additional dataframe of traits or guilds.
    # varname         = An unique key will be created to replace the cumbersome OTU key 
    #                   which is produced by Qiime. Varname is the desired column
    #                   name for this new OTU key. Unquoted. 
    # gene            = Gene region targeted in sequencing, e.g.: "ITS", "18S", "16S". 
    #                   Quoted. Used to select() column names, so must match text in 
    #                   column names. Also used to create distinct file names.
    # cluster_type    = Clustering algorithm output, e.g.: "otu", "sv". Quoted.
    #                   Used to create simple cluster IDs. 
    # append          = Additional text desired to differentiate among output files.
    #                   Defaults to NULL (""). Lead the string with "_".
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
                   species = str_sub(species, 4, nchar(species))) %>% 
            left_join(traits, by = join_by(phylum, class, order, family, genus)) %>% 
            select(-otu_ID, -kingdom, -growth_morphology, -trait, -notes, -citation)
        write_csv(meta, 
                  paste0(getwd(), folder, "/spe_", gene, append, "_guilds.csv"))
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
                  paste0(getwd(), folder, "/spe_", gene, append, "_taxonomy.csv"))
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
              paste0(getwd(), folder, "/spe_", gene, append, "_abund.csv"))
    
}
#'
#' # Load and process data
#' ## Import files
otu_its <- read_excel(paste0(getwd(), "/otu_tables/ITS/OTU_table_rrfd_3200_w_taxa.guilds.xlsx"), na = "-")
otu_18S <- read_delim(paste0(getwd(), "/otu_tables/TGP_18S_tables_021722/OTUs_18S_TGP_table_rarefied.txt"), 
                      delim = "\t", show_col_types = FALSE)
traits  <- read_excel(paste0(getwd(), "/otu_tables/13225_2020_466_MOESM4_ESM.xlsx"), sheet = "export_ver") %>% 
    select(phylum:primary_lifestyle)
#' 
#' ## ETL using `process_qiime()`
#' Schema: `process_qiime(data, varname, "gene", "cluster_type", "colname_prefix", "folder")`
#+ otu_its
process_qiime(
    data = otu_its,
    traits = traits,
    varname = otu_num,
    gene = "ITS",
    cluster_type = "otu",
    colname_prefix = "ITS_TGP_",
    folder = "/clean_data"
)
#+ otu_18S
process_qiime(
    data = otu_18S,
    varname = otu_num,
    gene = "18S",
    cluster_type = "otu",
    colname_prefix = "X18S_TGP_",
    folder = "/clean_data"
)
#' 



#' ## Resample and produce field averages
#' We examine diversity at the field level, so diversity obtained at samples should be averaged 
#' for each field. We collected ten samples from each field, but processing failed for some samples
#' at one step or another in the pipeline. Samples must be randomly resampled to the smallest 
#' number obtained in a series to produce comparable diversity metrics. Some OTUs or SVs may be "lost"
#' as a result, these were rare. Resultant zero sum columns will be removed. The following function 
#' will resample the site-species data to the correct number of samples and remove zero sum columns.
#' 
#' ### Import sites-species tables
its_otu_all <- read_csv(paste0(getwd(), "/clean_data/spe_ITS_otu_siteSpeMatrix_allReps.csv"), 
                        show_col_types = FALSE)
its_sv_all  <- read_csv(paste0(getwd(), "/clean_data/spe_ITS_sv_siteSpeMatrix_allReps.csv"), 
                        show_col_types = FALSE)
amf_otu_all <- read_csv(paste0(getwd(), "/clean_data/spe_18S_otu_siteSpeMatrix_allReps.csv"), 
                        show_col_types = FALSE)
amf_sv_all  <- read_csv(paste0(getwd(), "/clean_data/spe_18S_sv_siteSpeMatrix_allReps.csv"), 
                        show_col_types = FALSE)
#+ resample_fields_function 
resample_fields <- function(data, min, cluster_type) {
    set.seed(482)
    avg_df <-
        data %>%
        group_by(site_key) %>%
        slice_sample(n = min) %>%
        summarize(across(starts_with(cluster_type), ~ mean(.x, na.rm = TRUE)))
    null_spe <- which(apply(avg_df, 2, sum) == 0)
    out <- avg_df[,-null_spe]
    return(out)
}
#' 
#' The minimum number of samples in a field for each gene is:
#' 
#' - ITS = 8 samples
#' - 18S = 7 samples
#' 
#' With this, we can run the function for each dataset:
its_otu_avg <- resample_fields(its_otu_all, 8, "otu") %>% 
    write_csv(paste0(getwd(), "/clean_data/spe_ITS_otu_siteSpeMatrix_avg.csv"))
its_sv_avg  <- resample_fields(its_sv_all,  8, "sv") %>% 
    write_csv(paste0(getwd(), "/clean_data/spe_ITS_sv_siteSpeMatrix_avg.csv"))
amf_otu_avg <- resample_fields(amf_otu_all, 7, "otu") %>% 
    write_csv(paste0(getwd(), "/clean_data/spe_18S_otu_siteSpeMatrix_avg.csv"))
amf_sv_avg  <- resample_fields(amf_sv_all,  7, "sv") %>% 
    write_csv(paste0(getwd(), "/clean_data/spe_18S_sv_siteSpeMatrix_avg.csv"))





# Unifrac