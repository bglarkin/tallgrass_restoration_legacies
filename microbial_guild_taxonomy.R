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
sites   <- read_csv(paste0(getwd(), "/clean_data/site.csv"), show_col_types = FALSE) %>% 
    mutate(field_type = factor(site_type, ordered = TRUE, levels = c("corn", "restored", "remnant")),
           yr_since = replace(yr_since, which(site_type == "remnant"), rem_age)) %>% 
    filter(site_type != "oldfield")
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
            left_join(sites %>% select(-lat,-long,-yr_restore), by = "site_key")
    }



spe_meta <- list(
    its_otu =
        join_spe_meta(
            spe$its_otu,
            guild$its_otu,
            filter_char = "otu",
            clust_name = "otu_num"
        ) %>%
        select(-otu_ID, -trait, -notes, -citation) %>%
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
        select(-otu_ID, -trait, -notes, -citation) %>%
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
        select(-otu_ID, -accession) %>%
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
        select(-otu_ID, -accession) %>%
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


# CONSIDER A FUNCTION FOR THESE...
# Filter for quality of taxa assignment
taxon_level <- 9
# What is the distribution among site types at the class level?
spe_meta$its_otu %>%
    filter(taxon_level >= taxon_level & confidence %in% c("Highly Probable", "Probable")) %>%
    group_by(phylum, class, field_type, site_name) %>% 
    summarize(abund = sum(seq_abund), .groups = "drop") %>% 
    group_by(phylum, class, field_type) %>% 
    summarize(median = median(abund) %>% round(., 1), .groups = "drop") %>% 
    pivot_wider(names_from = field_type, values_from = median, values_fill = 0) %>% 
    select(phylum, class, corn, restored, remnant) %>%
    arrange(-remnant) %>% 
    kable(format = "pandoc", caption = "Distribution of OTUs in classes, median sequence abundance by field type")
# What is the distribution of trophic modes among site types?
spe_meta$its_otu %>% 
    filter(taxon_level >= taxon_level & confidence %in% c("Highly Probable", "Probable")) %>% 
    group_by(trophic_mode, field_type, site_name) %>% 
    summarize(abund = sum(seq_abund), .groups = "drop") %>% 
    group_by(trophic_mode, field_type) %>% 
    summarize(median = round(median(abund), 1), .groups = "drop") %>% 
    pivot_wider(names_from = field_type, values_from = median, values_fill = 0) %>% 
    select(trophic_mode, corn, restored, remnant) %>%
    arrange(-remnant) %>% 
    kable(format = "pandoc", caption = "Distribution of OTUs by trophic traits, median sequence abundance by field type")


# NEEDS TO BE SEQUENCE ABNDANCE INSTEAD OF RICHNESS; CONSIDER A FUNCTION HERE...
# Guilds are more specific than trophic modes. Filter to plant pathogens and AMF.
guilds <- c("Arbuscular Mycorrhizal", "Plant Pathogen")
for(i in 1:length(guilds)) {
    print(guilds[i])
    mod_data <- its_otu_spe_meta %>% 
        filter(taxon_level >= taxon_level & 
                   confidence %in% c("Highly Probable", "Probable") &
                   guild == guilds[i]) %>% 
        count(field_type, region, site_name, guild)
    mmod <- lmer(n ~ field_type + (1 | region), data = mod_data)
    mmod_null <- lmer(n ~ 1 + (1 | region), data = mod_data)
    print(anova(mmod, mmod_null))
    mod_tuk <- glht(mmod, linfct = mcp(field_type = "Tukey"), test = adjusted("holm"))
    print(summary(mod_tuk))
    print(cld(mod_tuk))
}
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
# It's possible that species richness can't be compared in these guilds if there are 
# too few species identified. Sites with low richness may just not have enough
# species detections. Check the rarefaction in these guilds and consider plotting
# the number of sequences per OTU rather than the count of OTUs.
amf_otu <- 
    its_otu_spe_meta %>% 
    filter(taxon_level >= taxon_level & 
               confidence %in% c("Highly Probable", "Probable") &
               guild == "Arbuscular Mycorrhizal") %>% 
    left_join(read_csv(paste0(getwd(), "/clean_data/spe_ITS_otu_samples.csv"))) %>% 
    mutate(seqs = as.integer(seq_abund * samples)) %>% 
    select(site_name, otu_num, seqs) %>% 
    pivot_wider(names_from = otu_num, values_from = seqs, values_fill = 0)
rarecurve(data.frame(amf_otu, row.names = 1), step = 1, label = TRUE, col = "blue", xlab = "Sequences", ylab = "AMF OTUs")
# Number of sequences is indeed very uneven, meaning that comparing species richness isn't appropriate
# Examine boxplots of sequence abundance for differences from richness plots
data.frame(
    site_name = amf_otu[, 1],
    seqs = apply(amf_otu[, -1], MARGIN = 1, sum)
) %>% 
    left_join(sites) %>% 
    ggplot(aes(x = field_type, y = seqs)) +
    geom_boxplot(fill = "gray90", varwidth = TRUE, outlier.shape = NA) +
    geom_beeswarm(aes(color = region), dodge.width = 0.2) +
    theme_bw()
# Some rank order changes, honestly I don't think the story would change at all

