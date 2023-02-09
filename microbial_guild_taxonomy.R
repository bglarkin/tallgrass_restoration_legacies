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
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#' # Packages and libraries
packages_needed = c("tidyverse",
                    "knitr",
                    "conflicted",
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
#' Oldfields are filtered out because they could not be replicated in regions.
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


# Need a unified OTU and OTU-based metadata table to analyze taxonomy and guilds
# These unified tables will also be necessary to filter species and produce ordinations
# Note: function process_taxa() embedded in this function
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

# we are ready to summarize and explore the data.
# FunGuild identifies the level of taxonomy assignment, the trophic mode, and the confidence
# of assignment. View the [README.md](https://github.com/UMNFuN/FUNGuild) on GitHub.
# Note on confidence ranking: I don't know what this refers to. Is it the taxonomic assignment,
# guild/trophic mode, or all data? This is not explained.
#
# To conduct summaries of FunGuild metadata, it would seem appropriate to choose OTUs with
# higher confidence rankings and more specific taxonomic assignments.

# ____Summary Analysis of ITS OTU data ----------------------


its_tax_trophic <- function(data, taxon_level = 9, cluster_type) {
    # What is the distribution among site types at the class level?
    taxonomy_df <-
        data %>%
        filter(taxon_level >= taxon_level &
                   confidence %in% c("Highly Probable", "Probable")) %>%
        group_by(phylum, class, field_type, site_name) %>%
        summarize(abund = sum(seq_abund), .groups = "drop") %>%
        group_by(phylum, class, field_type) %>%
        summarize(median = median(abund) %>% round(., 1),
                  .groups = "drop") %>%
        pivot_wider(
            names_from = field_type,
            values_from = median,
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
            "clusters in classes, median sequence abundance by field type"
        )
    ))
    # What is the distribution of trophic modes among site types?
    trophic_df <-
        data %>%
        filter(taxon_level >= taxon_level &
                   confidence %in% c("Highly Probable", "Probable")) %>%
        group_by(trophic_mode, field_type, site_name) %>%
        summarize(abund = sum(seq_abund), .groups = "drop") %>%
        group_by(trophic_mode, field_type) %>%
        summarize(median = round(median(abund), 1), .groups = "drop") %>%
        pivot_wider(
            names_from = field_type,
            values_from = median,
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
            "clusters in classes, median sequence abundance by field type"
        )
    ))
    
}


# NEEDS TO BE SEQUENCE ABNDANCE INSTEAD OF RICHNESS; CONSIDER A FUNCTION HERE...
# Guilds are more specific than trophic modes. Filter to plant pathogens and AMF.


its_guilds <- function(data) {
    guilds <-
        c("Arbuscular Mycorrhizal",
          "Plant Pathogen",
          "Undefined Saprotroph")
    for (i in 1:length(guilds)) {
        cat("---------------------------------\n")
        print(guilds[i])
        cat("---------------------------------\n")
        mod_data <- data %>%
            filter(
                taxon_level >= 9 &
                    confidence %in% c("Highly Probable", "Probable") &
                    guild == guilds[i]
            ) %>%
            group_by(field_type, region, site_name, guild) %>%
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
            lmer(seq_sum ~ 1 + (1 |
                                    region),
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
        cat("\n\n\n")
    }
}


# ITS OTU
its_tax_trophic(spe_meta$its_otu, cluster_type = "OTU")
its_guilds(spe_meta$its_otu)
# all NS
# CHECK THE USE OF SUMS...DIFFERENT NUMBER OF SITES FOR FIELD TYPES!
# LOOK AT YEARS SINCE RESTORATION IN RESTORED FIELDS

# ITS SV
its_tax_trophic(spe_meta$its_sv, cluster_type = "SV")
its_guilds(spe_meta$its_sv)
# all NS

# LOOK AT YEARS SINCE RESTORATION IN RESTORED FIELDS





# Significance labels needed for plot
sig_labs_otu <- data.frame(
    guild = rep("Arbuscular Mycorrhizal", 3),
    lab = c("a", "b", "ab"),
    xpos = c(1,2,3),
    ypos = rep(46, 3)
)
its_otu_spe_meta %>% 
    filter(taxon_level >= taxon_level & 
               confidence %in% c("Highly Probable", "Probable") &
               guild %in% c("Plant Pathogen", "Arbuscular Mycorrhizal")) %>% 
    count(field_type, region, site_name, guild) %>% 
    ggplot(aes(x = field_type, y = n)) +
    facet_wrap(vars(guild), scales = "free_y") +
    geom_boxplot(fill = "gray90", varwidth = TRUE, outlier.shape = NA) +
    geom_beeswarm(aes(color = region), dodge.width = 0.2) +
    geom_text(data = sig_labs_otu, aes(x = xpos, y = ypos, label = lab)) +
    labs(x = "", y = "Richness", caption = "OTUs at 97% similarity; mixed linear model with region as a random effect.") +
    theme_bw()







# 18S
spe_meta$amf_otu %>% glimpse()
spe_meta$amf_otu %>% group_by(site_key) %>% summarize(sum = sum(seq_abund))
apply(spe_meta$amf_otu, 2, unique)

# Don't sum, different number of sites!
spe_meta$amf_otu %>% 
    group_by(genus, taxon, field_type) %>% 
    summarize(seq_sum = sum(seq_abund) %>% round(., 1), .groups = "drop") %>% 
    pivot_wider(names_from = field_type, values_from = seq_sum, names_sort = TRUE, values_fill = 0) %>% 
    arrange(corn) %>%
    kable(format = "pandoc")
# 
