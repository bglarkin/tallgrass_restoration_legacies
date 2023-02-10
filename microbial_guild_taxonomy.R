#' ---
#' title: "Microbial data: microbial guilds and taxonomy"
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
#' Sequence clusters identified in QIIME are annotated with taxonomic information and
#' metadata from [FunGuild](https://github.com/UMNFuN/FUNGuild). In this report, sequence abundances
#' in taxonomic groups or fungal guilds are compared across field types and with time
#' since restoration.
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
    its_otu = read_csv(
        paste0(getwd(), "/clean_data/spe_ITS_otu_siteSpeMatrix_avg.csv"),
        show_col_types = FALSE
    ),
    its_sv  = read_csv(
        paste0(getwd(), "/clean_data/spe_ITS_sv_siteSpeMatrix_avg.csv"),
        show_col_types = FALSE
    ),
    amf_otu = read_csv(
        paste0(getwd(), "/clean_data/spe_18S_otu_siteSpeMatrix_avg.csv"),
        show_col_types = FALSE
    ),
    amf_sv  = read_csv(
        paste0(getwd(), "/clean_data/spe_18S_sv_siteSpeMatrix_avg.csv"),
        show_col_types = FALSE
    )
)
#' ## Species metadata
#' Load guild and taxonomy data (ITS sequences)
guild <- list(
    its_otu =
        read_csv(
            paste0(getwd(), "/clean_data/spe_ITS_otu_funGuild.csv"),
            show_col_types = FALSE
        ),
    its_sv  =
        read_csv(
            paste0(getwd(), "/clean_data/spe_ITS_sv_funGuild.csv"),
            show_col_types = FALSE
        )
)
#' Load taxonomy data (AMF/18S sequences)
taxonomy <- list(
    amf_otu =
        read_csv(
            paste0(getwd(), "/clean_data/spe_18S_otu_taxonomy.csv"),
            show_col_types = FALSE
        ),
    amf_sv  =
        read_csv(
            paste0(getwd(), "/clean_data/spe_18S_sv_taxonomy.csv"),
            show_col_types = FALSE
        )
)
#' ## Site metadata and design
#' Oldfields are filtered out because they could not be replicated in regions.
#' Set remnants to 50 years old as a placeholder. This number will not be used in
#' a quantitative sense, for example in models.
rem_age <- 50
sites   <-
    read_csv(paste0(getwd(), "/clean_data/site.csv"), show_col_types = FALSE) %>%
    mutate(
        field_type = factor(
            site_type,
            ordered = TRUE,
            levels = c("corn", "restored", "remnant")
        ),
        yr_since = replace(yr_since, which(site_type == "remnant"), rem_age)
    ) %>%
    filter(field_type != "oldfield") %>% 
    select(-lat, -long, -yr_restore, -site_type)
#' 
#' ## Joined species, metadata, and design tables
#' Functions streamline this process
#+ join_spe_meta
join_spe_meta <-
    function(spe,
             meta,
             filter_char = "otu",
             clust_name = "otu_num",
             abund_name = "seq_abund") {
        spe %>%
            pivot_longer(starts_with(filter_char),
                         names_to = clust_name,
                         values_to = abund_name) %>%
            filter(seq_abund != 0) %>%
            left_join(meta, by = clust_name) %>%
            left_join(sites, by = "site_key")
    }
#+ spe_meta_list
spe_meta <- list(
    its_otu =
        join_spe_meta(
            spe$its_otu,
            guild$its_otu,
            filter_char = "otu",
            clust_name = "otu_num"
        ) %>%
        select(-otu_ID,-trait,-notes,-citation) %>%
        write_csv(paste0(
            getwd(), "/clean_data/speGuild_ITS_otu.csv"
        )),
    its_sv =
        join_spe_meta(
            spe$its_sv,
            guild$its_sv,
            filter_char = "sv",
            clust_name = "sv_num"
        ) %>%
        select(-otu_ID,-trait,-notes,-citation) %>%
        write_csv(paste0(
            getwd(), "/clean_data/speGuild_ITS_sv.csv"
        )),
    amf_otu =
        join_spe_meta(
            spe$amf_otu,
            taxonomy$amf_otu,
            filter_char = "otu",
            clust_name = "otu_num"
        ) %>%
        select(-otu_ID,-accession) %>%
        write_csv(paste0(
            getwd(), "/clean_data/speTaxa_18S_otu.csv"
        )),
    amf_sv =
        join_spe_meta(
            spe$amf_sv,
            taxonomy$amf_sv,
            filter_char = "sv",
            clust_name = "sv_num"
        ) %>%
        select(-otu_ID,-accession) %>%
        write_csv(paste0(
            getwd(), "/clean_data/speTaxa_18S_sv.csv"
        ))
)
#' 
#' # Analysis and Results
#' ## ITS-based data
#' Functions streamline data processing, model fitting, and results output.
#' ### Function: ITS taxonomy
#+ its_tax_trophic
its_tax_trophic <- function(data, taxon_level = 9, cluster_type) {
    # What is the distribution among site types at the class level?
    taxonomy_df <-
        data %>%
        filter(taxon_level >= taxon_level &
                   confidence %in% c("Highly Probable", "Probable")) %>%
        group_by(phylum, class, field_type, site_name) %>%
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
        caption = paste(
            "Distribution of",
            cluster_type,
            "clusters in classes, mean sequence abundance by field type"
        )
    ))
    write_csv(taxonomy_df, 
              paste0(getwd(), "/microbial_diversity_files/its_", cluster_type, "_taxonomy.csv"))
    # What is the distribution of trophic modes among site types?
    trophic_df <-
        data %>%
        filter(taxon_level >= taxon_level &
                   confidence %in% c("Highly Probable", "Probable")) %>%
        group_by(trophic_mode, field_type, site_name) %>%
        summarize(abund = sum(seq_abund), .groups = "drop") %>%
        group_by(trophic_mode, field_type) %>%
        summarize(mean = round(mean(abund), 1), .groups = "drop") %>%
        pivot_wider(
            names_from = field_type,
            values_from = mean,
            values_fill = 0
        ) %>%
        select(trophic_mode, corn, restored, remnant) %>%
        arrange(-remnant)
    print(kable(
        trophic_df,
        format = "pandoc",
        caption = paste(
            "Distribution of",
            cluster_type,
            "clusters in classes, mean sequence abundance by field type"
        )
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
its_guilds <- function(data) {
    guilds <-
        c("Arbuscular Mycorrhizal",
          "Plant Pathogen",
          "Undefined Saprotroph")
    df1 <- data.frame()
    for (i in 1:length(guilds)) {
        cat("\n---------------------------------\n")
        print(guilds[i])
        cat("---------------------------------\n")
        mod_data <- data %>%
            filter(
                taxon_level >= 9 &
                    confidence %in% c("Highly Probable", "Probable") &
                    guild == guilds[i]
            ) %>%
            group_by(field_type, region, site_name, yr_since, guild) %>%
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
            guilds[i],
            "sequence abundance"
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
its_tax_trophic(spe_meta$its_otu, cluster_type = "OTU")
#' 
#' All comparisions across field types are non-significant
#+ its_guilds_otu,message=FALSE
its_otu_guilds <- its_guilds(spe_meta$its_otu)
#' 
#' Plant pathogens correlate with restoration age in Blue Mounds area. 
#+ its_otu_pathogen_correlation,message=FALSE
its_otu_guilds %>% 
    filter(field_type == "restored", guild == "Plant Pathogen", region == "BM") %>% 
    ggplot(aes(x = yr_since, y = seq_sum)) +
    geom_smooth(method = "lm", se = TRUE) +
    geom_label(aes(label = site_name)) +
    labs(x = "Years since restoration", y = "Sum of ITS sequences (97% OTUs)", caption = "R2adj=0.59, p<0.05", title = "Plant pathogen sequence abundance in restored fields") +
    theme_classic()
#' 
#' ### ITS sequences in SV clusters
#' Function outputs are verbose, but details may be necessary later so they are displayed here.
#+ its_tax_trophic_sv,message=FALSE
its_tax_trophic(spe_meta$its_sv, cluster_type = "SV")
#' 
#' All comparisions across field types are non-significant
#+ its_guilds_sv,message=FALSE
its_sv_guilds <- its_guilds(spe_meta$its_sv)
#' 
#' Plant pathogens correlate with restoration age in Blue Mounds area. 
#+ its_sv_pathogen_correlation,message=FALSE
its_sv_guilds %>% 
    filter(field_type == "restored", guild == "Plant Pathogen", region == "BM") %>% 
    ggplot(aes(x = yr_since, y = seq_sum)) +
    geom_smooth(method = "lm", se = TRUE) +
    geom_label(aes(label = site_name)) +
    labs(x = "Years since restoration", y = "Sum of ITS sequences (100% SVs)", caption = "R2adj=0.58, p<0.05", title = "Plant pathogen sequence abundance in restored fields") +
    theme_classic()
#' 
#' ## 18S-based data (AMF)
#' A function streamlines analysis and results output.
#+ amf_taxa_function
amf_tax <- function(data, cluster_type) {
    cat("\n---------------------------------\n")
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
            "/microbial_diversity_files/amf_",
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
#+ claroideoglomeraceae_otu_fields_fig,message=FALSE
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
#' Gigasporaceae increased with time since restoration by a simple linear regression, $R^2_{adj}$ = 0.81, p < 0.01
#+ gigasporaceae_otu_time_fig,message=FALSE
amf_otu_summary %>% 
    filter(field_type == "restored", region == "BM", family == "Gigasporaceae") %>% 
    ggplot(aes(x = yr_since, y = seq_sum)) +
    geom_smooth(method = "lm", linewidth = 0.4, se = FALSE) +
    geom_point(size = 2, shape = 21, fill = "gray60") +
    labs(x = "Years since restoration", y = "Sequence abundance", title = "Gigasporaceae abundance since restoration, 97% OTU",
         caption = "R2Adj = 0.81, p<0.01") +
    theme_classic()





amf_sv_summary  <- amf_tax(spe_meta$amf_sv,  "sv")
# Claroideoglomeraceae different in fields a, b, ab, test p<0.01
# Gigasporaceae different over time R2adj 0.65, p<0.05
amf_sv_summary %>% 
    filter(family == "Claroideoglomeraceae") %>% 
    ggplot(aes(x = field_type, y = seq_sum)) +
    geom_boxplot(varwidth = TRUE, fill = "gray90", outlier.shape = NA) +
    geom_beeswarm(aes(fill = region), shape = 21, size = 2, dodge.width = 0.2) +
    annotate("text", label = c("a", "b", "ab"), x = c(1,2,3), y = rep(350, 3)) +
    labs(x = "", y = "Sequence abundance", title = "AMF variation in field types, 100% SV",
         caption = "Likelihood ratio test p<0.01, Tukey's post-hoc with Holm correction at p<0.05") +
    scale_fill_discrete_qualitative(palette = "Dark3") +
    theme_classic()
amf_sv_summary %>% 
    filter(field_type == "restored", region == "BM", family == "Gigasporaceae") %>% 
    ggplot(aes(x = yr_since, y = seq_sum)) +
    geom_smooth(method = "lm", linewidth = 0.4, se = FALSE) +
    geom_point(size = 2, shape = 21, fill = "gray60") +
    labs(x = "Years since restoration", y = "Sequence abundance", title = "Gigasporaceae abundance since restoration, 100% SV",
         caption = "R2Adj = 0.65, p<0.05") +
    theme_classic()
