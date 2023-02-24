#' ---
#' title: "Microbial data: microbial guilds and taxonomy"
#' author: "Beau Larkin\n"
#' date: "Last updated: `r format(Sys.time(), '%d %B, %Y')`"
#' output:
#'   github_document:
#'     toc: true
#'     toc_depth: 3
#'     df_print: paged
#'     fig_width: 5.5
#'     fig_height: 5
#' ---
#'
#' # Description
#' Sequence clusters identified in QIIME2 are annotated with taxonomic information and
#' metadata from [Fungal traits](https://link.springer.com/article/10.1007/s13225-020-00466-2). 
#' In this report, sequence abundances in taxonomic groups or fungal guilds are compared 
#' across field types and with time since restoration. 
#' 
#' The full sequence abundance tables were rarefied to make sequencing depth equvalent
#' across fields. This can result in lower-abundance OTUs dropping to zero, making
#' comparisons within guilds less informative. Here, the raw abundances are used, filtered
#' into particular guilds, and then rarefied within those guilds to better reveal 
#' trends.
#' 
#' # Packages and libraries
packages_needed = c("tidyverse",
                    "knitr",
                    "conflicted",
                    "ggbeeswarm",
                    "colorspace",
                    "rsq",
                    "lme4",
                    "multcomp")
packages_installed = packages_needed %in% rownames(installed.packages())
#+ packages,message=FALSE
if (any(!packages_installed)) {
    install.packages(packages_needed[!packages_installed])
}
#+ libraries,message=FALSE
for (i in 1:length(packages_needed)) {
    library(packages_needed[i], character.only = T)
}
#+ conflicts,message=FALSE
conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")
#'
#' # Data
#' ## Sites-species tables
#' CSV files were produced in `process_data.R`
spe <- list(
    its_raw = read_csv(
        paste0(getwd(), "/clean_data/spe_ITS_raw.csv"),
        show_col_types = FALSE
    ),
    its_rfy = read_csv(
        paste0(getwd(), "/clean_data/spe_ITS_rfy.csv"),
        show_col_types = FALSE
    ),
    amf_raw = read_csv(
        paste0(getwd(), "/clean_data/spe_18S_raw.csv"),
        show_col_types = FALSE
    ),
    amf_rfy = read_csv(
        paste0(getwd(), "/clean_data/spe_18S_rfy.csv"),
        show_col_types = FALSE
    )
)
#' ## Species metadata
#' Load taxonomy for all and guilds (`primary_lifestyle` in Fungal Traits)
#' for ITS OTUs.
meta <- list(
    its_raw = 
        read_csv(
            paste0(getwd(), "/clean_data/spe_ITS_raw_taxonomy.csv"),
            show_col_types = FALSE
        ),
    its_rfy = 
        read_csv(
            paste0(getwd(), "/clean_data/spe_ITS_rfy_taxonomy.csv"),
            show_col_types = FALSE
        ),
    amf_raw = 
        read_csv(
            paste0(getwd(), "/clean_data/spe_18S_raw_taxonomy.csv"),
            show_col_types = FALSE
        ),
    amf_rfy = 
        read_csv(
            paste0(getwd(), "/clean_data/spe_18S_rfy_taxonomy.csv"),
            show_col_types = FALSE
        )
)
#' ## Site metadata and design
#' Set remnants to 50 years old as a placeholder. This number will not be used in
#' a quantitative sense, for example in models.
rem_age <- 50
sites   <-
    read_csv(paste0(getwd(), "/clean_data/sites.csv"), show_col_types = FALSE) %>%
    mutate(
        field_type = factor(
            field_type,
            ordered = TRUE,
            levels = c("corn", "restored", "remnant")
        ),
        yr_since = replace(yr_since, which(field_type == "remnant"), rem_age)
    ) %>%
    select(-lat, -long, -yr_restore, -yr_rank)
#' 
#' ## Joined species, metadata, and design tables
#' Functions streamline this process
#+ join_spe_meta_fun
join_spe_meta <-
    function(spe, meta) {
        spe %>%
            pivot_longer(starts_with("otu"),
                         names_to = "otu_num",
                         values_to = "seq_abund") %>%
            filter(seq_abund != 0) %>%
            left_join(meta, by = join_by(otu_num)) %>%
            left_join(sites, by = join_by(field_key))
    }
#+ spe_meta_list
spe_meta <- list(
    its_raw = 
        join_spe_meta(spe$its_raw, meta$its_raw) %>%
        write_csv(paste0(getwd(), "/clean_data/speTaxa_ITS_raw.csv")),
    its_rfy = 
        join_spe_meta(spe$its_rfy, meta$its_rfy) %>%
        write_csv(paste0(getwd(), "/clean_data/speTaxa_ITS_rfy.csv")),
    amf_raw = 
        join_spe_meta(spe$amf_raw, meta$amf_raw) %>%
        write_csv(paste0(getwd(), "/clean_data/speTaxa_18S_raw.csv")),
    amf_rfy = 
        join_spe_meta(spe$amf_rfy, meta$amf_rfy) %>%
        write_csv(paste0( getwd(), "/clean_data/speTaxa_18S_rfy.csv" ))
)
#' 
#' # Analysis and Results
#' ## ITS-based data
#' The first pass uses the entire rarefied tables for ITS or AMF. Interesting
#' patterns can then be explored with the raw sequence abundance data, rarefied
#' in particular guilds or taxa. 
#' 
#' Functions streamline data processing, model fitting, and results output.
#' ### Function: ITS taxonomy
#+ its_tax_trophic
its_taxaGuild <- function(data) {
    # What is the distribution among site types at the class level?
    taxonomy_df <-
        data %>%
        group_by(phylum, class, field_type, field_name) %>%
        summarize(abund = sum(seq_abund), .groups = "drop") %>%
        group_by(phylum, class, field_type) %>%
        summarize(mean = mean(abund) %>% round(., 2),
                  .groups = "drop") %>%
        pivot_wider(
            names_from = field_type,
            values_from = mean,
            values_fill = 0
        ) %>%
        select(phylum, class, corn, restored, remnant) %>%
        arrange(-remnant)
    print(kable(
        taxonomy_df,
        format = "pandoc",
        caption = "Distribution of ITS OTUs in classes; mean sequence abundance by field type"
    ))
    # What is the distribution of `primary_lifestyles` among site types?
    guild_df <-
        data %>%
        group_by(primary_lifestyle, field_type, field_name) %>%
        summarize(abund = sum(seq_abund), .groups = "drop") %>%
        group_by(primary_lifestyle, field_type) %>%
        summarize(mean = round(mean(abund), 1), .groups = "drop") %>%
        pivot_wider(
            names_from = field_type,
            values_from = mean,
            values_fill = 0
        ) %>%
        select(primary_lifestyle, corn, restored, remnant) %>%
        arrange(-remnant)
    print(kable(
        guild_df,
        format = "pandoc",
        caption = "Distribution of ITS OTUs by Fungal Trait 'primary_lifestyle'; mean sequence abundance by field type"
    ))
}
#' 






#' ### Function: ITS and guilds
#' Sequence abundances are analyzed in the guilds "Arbuscular Mycorrhizal", "Plant Pathogen",
#' and "Undefined Saprotroph". Other guilds are inappropriate or contain few sequences. 
#' 
#' FunGuild identifies the level of taxonomy assignment, the trophic mode, and the confidence
#' of assignment. View the [README](https://github.com/UMNFuN/FUNGuild) on GitHub.
#' Note on confidence ranking: I don't know what this refers to. Is it the taxonomic assignment,
#' guild/trophic mode, or all data? This is not explained.
#' 
#' To conduct summaries of FunGuild metadata, it would seem appropriate to choose OTUs with
#' higher confidence rankings and more specific taxonomic assignments.
#+ its_guilds
its_test_taxaGuild <- function(data) {
    pl <- c("soil_saprotroph", "wood_saprotroph", "litter_saprotroph", "plant_pathogen")
    df1 <- data.frame()
    for (i in 1:length(pl)) {
        cat("---------------------------------\n")
        print(pl[i])
        cat("---------------------------------\n")
        mod_data <- data %>%
            filter(primary_lifestyle == pl[i]) %>%
            group_by(field_type, region, field_name, yr_since, primary_lifestyle) %>%
            summarize(seq_sum = sum(seq_abund), .groups = "drop")
        print(kable(mod_data %>% arrange(-seq_sum), format = "pandoc"))
        cat("----------------------------------------------------\n\n")
        mmod <-
            lmer(seq_sum ~ field_type + (1 | region),
                 data = mod_data,
                 REML = FALSE)
        print(mmod)
        cat("----------------------------------------------------\n\n")
        mmod_null <-
            lmer(seq_sum ~ 1 + (1 | region),
                 data = mod_data,
                 REML = FALSE)
        print(mmod_null)
        cat("----------------------------------------------------\n\n")
        print(anova(mmod, mmod_null))
        cat("----------------------------------------------------\n\n")
        mod_tuk <-
            glht(mmod,
                 linfct = mcp(field_type = "Tukey"),
                 test = adjusted("holm"))
        print(summary(mod_tuk))
        print(cld(mod_tuk))
        cat("----------------------------------------------------\n\n")
        print(paste(
            "Years since restoration and",
            pl[i],
            "sequence abundance in Blue Mounds Area"
        ))
        mod_data2 <- mod_data %>%
            filter(region == "BM", field_type == "restored")
        print(summary(lm(seq_sum ~ yr_since,
                         data = mod_data2)))
        cat("\n\n\n")
        df1 <- rbind(df1, mod_data)
    }
    
    return(df1)

}




#' 
#' ### ITS sequences in OTU clusters
#' Function outputs are verbose, but details may be necessary later so they are displayed here.
#+ its_tax_trophic_otu,message=FALSE
its_taxaGuild(spe_meta$its_rfy)
#' 
#+ its_guilds_otu,message=FALSE
its_rfy_guilds <- its_test_taxaGuild(spe_meta$its_rfy)
#' 
#' Soil saprotroph increases with years since
#' Wood saprotroph differs among field types and decreases with years since
#' Plant pathogens decrease with years since
#' Re-rarefy in just these groups, test, and make figures. 





#' 
#' ## 18S-based data (AMF)
#' A function streamlines analysis and results output.
#+ amf_taxa_function
amf_tax <- function(data, cluster_type) {
    cat("---------------------------------\n")
    print(paste("AMF", cluster_type))
    cat("---------------------------------\n")
    amf_df <-
        data %>%
        group_by(family, field_type, region, site_name, yr_since) %>%
        summarize(seq_sum = sum(seq_abund) %>% round(., 1),
                  .groups = "drop")
    amf_df_summary <-
        amf_df %>%
        group_by(family, field_type) %>%
        summarize(seq_avg = mean(seq_sum) %>% round(., 1),
                  .groups = "drop") %>%
        pivot_wider(
            names_from = field_type,
            values_from = seq_avg,
            names_sort = TRUE,
            values_fill = 0
        ) %>%
        arrange(-remnant)
    print(kable(amf_df_summary, format = "pandoc"))
    write_csv(
        amf_df_summary,
        paste0(
            getwd(),
            "/microbial_guild_taxonomy_files/amf_",
            cluster_type,
            "_taxonomy.csv"
        )
    )
    cat("\n---------------------------------\n")
    print("Compare abundances across field types with mixed model")
    cat("---------------------------------\n")
    test_families <-
        amf_df %>% 
        count(region, family, field_type) %>% 
        count(region, family) %>% 
        filter(n == 3) %>% 
        pull(family) %>% 
        unique()
    for (i in 1:length(test_families)) {
        cat("\n---------------------------------\n")
        print(test_families[i])
        cat("---------------------------------\n")
        mmod <-
            lmer(
                seq_sum ~ field_type + (1 | region),
                data = amf_df %>% filter(family == test_families[i]),
                REML = FALSE
            )
        print(mmod)
        cat("----------------------------------------------------\n\n")
        mmod_null <-
            lmer(
                seq_sum ~ 1 + (1 | region),
                data = amf_df %>% filter(family == test_families[i]),
                REML = FALSE
            )
        print(mmod_null)
        cat("----------------------------------------------------\n\n")
        print(anova(mmod, mmod_null))
        cat("----------------------------------------------------\n\n")
        mod_tuk <-
            glht(mmod,
                 linfct = mcp(field_type = "Tukey"),
                 test = adjusted("holm"))
        print(summary(mod_tuk))
        print(cld(mod_tuk))
        cat("\n")
    }
    cat("\n---------------------------------\n")
    print("Test abundances with years since restoration")
    cat("---------------------------------\n")
    all7 <-
        amf_df %>%
        filter(field_type == "restored", region == "BM") %>%
        count(family) %>%
        filter(n == 7) %>%
        pull(family)
    mod_data <-
        amf_df %>%
        filter(field_type == "restored", region == "BM", family %in% all7)
        for (i in 1:length(all7)) {
        print(all7[i])
        print(summary(lm(
            seq_sum ~ yr_since, data = mod_data %>% filter(family == all7[i])
        )))
    }
    return(amf_df)
}
#' 
#' ## AMF OTUs
#' Function output is verbose but retained as explained previously.
#+ amf_otu_summary,message=FALSE
amf_otu_summary <- amf_tax(spe_meta$amf_otu, "otu")
#' Claroideoglomeraceae differs across field types with a likelihood ratio test result p<0.01. 
#' Tukey's post-hoc test with Holm correction performed, letters on the figure show differences.
#+ claroideoglomeraceae_otu_fields_fig,message=FALSE,fig.align='center'
amf_otu_summary %>% 
    filter(family == "Claroideoglomeraceae") %>% 
    ggplot(aes(x = field_type, y = seq_sum)) +
    geom_boxplot(varwidth = TRUE, fill = "gray90", outlier.shape = NA) +
    geom_beeswarm(aes(fill = region), shape = 21, size = 2, dodge.width = 0.2) +
    annotate("text", label = c("a", "b", "ab"), x = c(1,2,3), y = rep(350, 3)) +
    labs(x = "", y = "Sequence abundance", title = "AMF variation in field types, 97% OTU",
         caption = "Likelihood ratio test p<0.01, Tukey's post-hoc with Holm correction at p<0.05") +
    scale_fill_discrete_qualitative(palette = "Dark3") +
    theme_classic()
#' Gigasporaceae increased with time since restoration by a simple linear regression, 
#' $R^2_{adj}$ = 0.81, p < 0.01
#+ gigasporaceae_otu_time_fig,message=FALSE,fig.align='center'
amf_otu_summary %>% 
    filter(field_type == "restored", region == "BM", family == "Gigasporaceae") %>% 
    ggplot(aes(x = yr_since, y = seq_sum)) +
    geom_smooth(method = "lm", linewidth = 0.4, se = FALSE) +
    geom_point(size = 2, shape = 21, fill = "gray60") +
    labs(x = "Years since restoration", 
         y = "Sequence abundance", 
         title = "Gigasporaceae abundance since restoration, 97% OTU",
         caption = "R2Adj = 0.81, p<0.01") +
    theme_classic()






#' 
#' # Conclusions: taxa and guilds
#' Little variation exists here for ITS or AMF sequences among field types, although 
#' classes of fungi identified through ITS sequences remain to be closely examined. 
#' It's striking that plant pathogens decline as restorations age while 
#' the AMF family *Gigasporaceae* increases, but this contrast was not found in any 
#' other group of AMF and the *Gigasporaceae* aren't particularly abundant to begin with.
#' 
#' 
#' 
#' 
#' # New section: INDVAL; let's look for species with stories...