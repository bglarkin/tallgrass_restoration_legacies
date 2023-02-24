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
#' # Resources
#' ## Packages and libraries
packages_needed = c("tidyverse")
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
process_qiime <- function(spe, taxa, samps, traits=NULL, varname, gene, cluster_type, colname_prefix, folder) {
    
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
            select(-otu_ID, -kingdom, -Confidence)
        write_csv(meta,
                  paste0(getwd(), folder, "/spe_", gene, "_taxaguild.csv"))
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
            select(-otu_ID, -Confidence)
        write_csv(meta,
                  paste0(getwd(), folder, "/spe_", gene, "_taxonomy.csv"))
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
        mutate(field_key = as.numeric(field_key))
    null_spe <- which(apply(spe_topn, 2, sum) == 0)
    spe_all <- spe_topn[,-null_spe]
    write_csv(spe_all, paste0(getwd(), folder, "/spe_", gene, "_raw.csv"))
    out <- list(
        spe_meta = meta %>% filter(!(otu_num %in% names(null_spe))),
        spe_all  = spe_all
    )
    return(out)
    
}


#'
#' # Load and process data
#' ## Import files
its_otu  <- read_delim(paste0(getwd(), "/otu_tables/ITS/ITS_otu_raw.txt"), show_col_types = FALSE)
its_taxa <- read_delim(paste0(getwd(), "/otu_tables/ITS/ITS_otu_taxonomy.txt"), show_col_types = FALSE)
# The 18S OTU file contains an unknown site label in the last column; remove it
amf_otu  <- read_delim(paste0(getwd(), "/otu_tables/18S/18S_otu_raw.txt"), show_col_types = FALSE) %>% select(-last_col())
amf_taxa <- read_delim(paste0(getwd(), "/otu_tables/18S/18S_otu_taxonomy.txt"), show_col_types = FALSE)
traits   <- read_csv(paste0(getwd(),   "/otu_tables/2023-02-23_fungal_traits.csv"), show_col_types = FALSE) %>% 
    select(phylum:primary_lifestyle)
#' 
#' ## ETL using `process_qiime()`
#' Schema: `process_qiime(data, varname, "gene", "cluster_type", "colname_prefix", "folder")`
#+ otu_its,message=FALSE
its <- 
    process_qiime(
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
#+ otu_18S,message=FALSE
amf <- 
    process_qiime(
        spe = amf_otu,
        taxa = amf_taxa,
        samps = 6,
        varname = otu_num,
        gene = "18S",
        cluster_type = "otu",
        colname_prefix = "X18S_TGP_",
        folder = "/clean_data"
    )
#' 



## Sum in fields, rarefy, and go



#' ## Resample and produce field averages
#' We examine diversity at the field level, so diversity obtained at samples should be averaged 
#' for each field. We collected ten samples from each field, but processing failed for some samples
#' at one step or another in the pipeline. Samples must be randomly resampled to the smallest 
#' number obtained in a series to produce comparable diversity metrics. Some OTUs or SVs may be "lost"
#' as a result, these were rare. Resultant zero sum columns will be removed. The following function 
#' will resample the site-species data to the correct number of samples and remove zero sum columns.



#+ resample_fields_function 
# resample_fields <- function(data, min, cluster_type) {
#     set.seed(482)
#     avg_df <-
#         data %>%
#         group_by(site_key) %>%
#         slice_sample(n = min) %>%
#         summarize(across(starts_with(cluster_type), ~ mean(.x, na.rm = TRUE)))
#     null_spe <- which(apply(avg_df, 2, sum) == 0)
#     out <- avg_df[,-null_spe]
#     return(out)
# }
#' 
#' The minimum number of samples in a field for each gene is:
#' 
#' - ITS = 8 samples
#' - 18S = 7 samples
#' 
#' With this, we can run the function for each dataset:
# its_otu_avg <- resample_fields(its_otu_all, 8, "otu") %>% 
#     write_csv(paste0(getwd(), "/clean_data/spe_ITS_otu_siteSpeMatrix_avg.csv"))
# its_sv_avg  <- resample_fields(its_sv_all,  8, "sv") %>% 
#     write_csv(paste0(getwd(), "/clean_data/spe_ITS_sv_siteSpeMatrix_avg.csv"))
# amf_otu_avg <- resample_fields(amf_otu_all, 7, "otu") %>% 
#     write_csv(paste0(getwd(), "/clean_data/spe_18S_otu_siteSpeMatrix_avg.csv"))
# amf_sv_avg  <- resample_fields(amf_sv_all,  7, "sv") %>% 
#     write_csv(paste0(getwd(), "/clean_data/spe_18S_sv_siteSpeMatrix_avg.csv"))





# Unifrac


# but before unifrac, you need to fix the site_key and sample columns; must be one column with a continuous numeric
# that can join back to the sites table



