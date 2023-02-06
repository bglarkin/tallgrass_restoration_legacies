# guilds and taxonomy.....


#' ## Species metadata
#' Load guild data
its_otu_guild <- read_csv(paste0(getwd(), "/clean_data/spe_ITS_otu_funGuild.csv"), show_col_types = FALSE)
its_sv_guild  <- read_csv(paste0(getwd(), "/clean_data/spe_ITS_sv_funGuild.csv"), show_col_types = FALSE)
#' Load taxonomy data
amf_otu_tax <- read_csv(paste0(getwd(), "/clean_data/spe_18S_otu_taxonomy.csv"), show_col_types = FALSE)
amf_sv_tax  <- read_csv(paste0(getwd(), "/clean_data/spe_18S_sv_taxonomy.csv"), show_col_types = FALSE)
#'




# __Trophic Modes, Guilds and Taxonomy ----------------------
# Do any species register 0 abundance across all sites?
which(apply(its_otu, MARGIN = 2, FUN = sum) == 0) # No
which(apply(its_sv,  MARGIN = 2, FUN = sum) == 0) # No
# Function to clean up Qiime taxomomy from guilds files
process_taxa <- function(data) {
    data %>% 
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
}
# Need a unified OTU and OTU-based metadata table to analyze taxonomy and guilds
# These unified tables will also be necessary to filter species and produce ordinations
# Note: function process_taxa() embedded in this function
join_spe_meta <- function(spe, meta, filter_char = "otu", clust_name = "otu_num", abund_name = "seq_abund") {
    spe %>%
        na_if(., 0) %>% 
        pivot_longer(starts_with(filter_char), names_to = clust_name, values_to = abund_name) %>% 
        drop_na() %>% 
        left_join(process_taxa(meta) %>% select(-otu_ID, -trait, -notes, -citation), by = clust_name) %>% 
        left_join(sites %>% select(starts_with("site"), region, yr_rank, yr_since), by = "site_key")
}

# With these functions in hand, we are ready to summarize and explore the data. 
# FunGuild identifies the level of taxonomy assignment, the trophic mode, and the confidence
# of assignment. View the [README.md](https://github.com/UMNFuN/FUNGuild) on GitHub.
# Note on confidence ranking: I don't know what this refers to. Is it the taxonomic assignment,
# guild/trophic mode, or all data? This is not explained. 
# 
# To conduct summaries of FunGuild metadata, it would seem appropriate to choose OTUs with
# higher confidence rankings and more specific taxonomic assignments. 

# ____Summary Analysis of ITS OTU data ----------------------
# Create the summary data using functions; write csv for use in other scripts
its_otu_spe_meta <- join_spe_meta(its_otu, its_otu_guild) %>% glimpse()
write_csv(its_otu_spe_meta, paste0(getwd(), "/clean_data/speMeta_ITS_otu.csv"))
# Filter for quality of taxa assignment and investigate NAs
taxon_level <- 9
its_otu_spe_meta %>% 
    filter(taxon_level >= taxon_level & confidence %in% c("Highly Probable", "Probable")) %>% 
    filter(is.na(species))
# What is the distribution among site types at the class level?
its_otu_spe_meta %>%
    filter(taxon_level >= taxon_level & confidence %in% c("Highly Probable", "Probable")) %>%
    group_by(phylum, class, field_type, site_name) %>% 
    summarize(abund = sum(seq_abund), .groups = "drop") %>% 
    group_by(phylum, class, field_type) %>% 
    summarize(median = median(abund) %>% round(., 1), .groups = "drop") %>% 
    pivot_wider(names_from = field_type, values_from = median, values_fill = 0) %>% 
    rownames_to_column(var = "number") %>% 
    select(number, phylum, class, corn, restored, remnant) %>%
    arrange(-remnant) %>% 
    kable(format = "pandoc", caption = "Distribution of OTUs in classes, median by site")
# What is the distribution of trophic modes among site types?
its_otu_spe_meta %>% 
    filter(taxon_level >= taxon_level & confidence %in% c("Highly Probable", "Probable")) %>% 
    group_by(trophic_mode, field_type, site_name) %>% 
    summarize(abund = sum(seq_abund), .groups = "drop") %>% 
    group_by(trophic_mode, field_type) %>% 
    summarize(median = round(median(abund), 1), .groups = "drop") %>% 
    pivot_wider(names_from = field_type, values_from = median, values_fill = 0) %>% 
    rownames_to_column(var = "number") %>% 
    select(number, trophic_mode, corn, restored, remnant) %>%
    kable(format = "pandoc", caption = "Distribution of OTUs by trophic traits, median by site")
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


# ____Summary Analysis of ITS SV data ----------------------
# Use functions to produce table; write csv to wd for use in other scripts
its_sv_spe_meta <- join_spe_meta(its_sv, its_sv_guild, filter_char = "sv", clust_name = "sv_num") %>% glimpse()
write_csv(its_sv_spe_meta, paste0(getwd(), "/clean_data/speMeta_ITS_sv.csv"))
# Highest quality data are with specific taxon level (13), and good confidence
# All NAs in this filtered quality set are from no specific identification. All are identified to genus
taxon_level <- 9
its_sv_spe_meta %>% 
    filter(taxon_level >= taxon_level & confidence %in% c("Highly Probable", "Probable")) %>% 
    filter(is.na(species))
its_sv_spe_meta %>% 
    filter(taxon_level >= taxon_level & confidence %in% c("Highly Probable", "Probable")) %>% 
    filter(is.na(genus))
# What is the distribution among site types at the class level?
its_sv_spe_meta %>% 
    filter(taxon_level >= taxon_level & confidence %in% c("Highly Probable", "Probable")) %>% 
    group_by(phylum, class, field_type, site_name) %>% 
    summarize(abund = sum(seq_abund), .groups = "drop") %>% 
    group_by(phylum, class, field_type) %>% 
    summarize(median = round(median(abund), 1), .groups = "drop") %>% 
    pivot_wider(names_from = field_type, values_from = median, values_fill = 0) %>% 
    rownames_to_column(var = "number") %>% 
    select(number, phylum, class, corn, restored, remnant) %>%
    kable(format = "pandoc", caption = "Distribution of SVs in classes, median by site")
# What is the distribution of guilds among site types?
its_sv_spe_meta %>% 
    filter(taxon_level >= taxon_level & confidence %in% c("Highly Probable", "Probable")) %>% 
    group_by(trophic_mode, field_type, site_name) %>% 
    summarize(abund = sum(seq_abund), .groups = "drop") %>% 
    group_by(trophic_mode, field_type) %>% 
    summarize(median = round(median(abund), 1), .groups = "drop") %>% 
    pivot_wider(names_from = field_type, values_from = median, values_fill = 0) %>% 
    rownames_to_column(var = "number") %>% 
    select(number, trophic_mode, corn, restored, remnant) %>%
    kable(format = "pandoc", caption = "Distribution of SVs by trophic traits, median by site")
# Guilds are more specific than trophic modes. Filter to plant pathogens and AMF.
guilds <- c("Arbuscular Mycorrhizal", "Plant Pathogen")
for(i in 1:length(guilds)) {
    print(guilds[i])
    mod_data <- its_sv_spe_meta %>% 
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
its_sv_spe_meta %>% 
    filter(taxon_level >= taxon_level & 
               confidence %in% c("Highly Probable", "Probable") &
               guild %in% c("Plant Pathogen", "Arbuscular Mycorrhizal")) %>% 
    count(field_type, region, site_name, guild) %>% 
    ggplot(aes(x = field_type, y = n)) +
    facet_wrap(vars(guild), scales = "free_y") +
    geom_boxplot(fill = "gray90", varwidth = TRUE) +
    geom_beeswarm(aes(color = region), dodge.width = 0.2) +
    labs(x = "", y = "Richness", caption = "SVs at 100% similarity") +
    theme_bw()
# Look at change over time in guilds...
its_sv_spe_meta %>% 
    filter(taxon_level >= taxon_level & 
               confidence %in% c("Highly Probable", "Probable") &
               guild %in% c("Plant Pathogen", "Arbuscular Mycorrhizal")) %>% 
    count(field_type, region, site_name, guild) %>% 
    left_join(its_sv_spe_meta %>% select(site_name, yr_since), by = "site_name") %>% 
    ggplot(aes(x = yr_since, y = n)) +
    facet_wrap(vars(guild), scales = "free_y") +
    geom_point(aes(color = region, shape = field_type)) +
    theme_bw()
