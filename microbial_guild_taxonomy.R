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
                    "multcomp",
                    "indicspecies",
                    "GUniFrac")
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
sites   <-
    read_csv(paste0(getwd(), "/clean_data/sites.csv"), show_col_types = FALSE) %>%
    mutate(
        field_type = factor(
            field_type,
            ordered = TRUE,
            levels = c("corn", "restored", "remnant")
        )) %>%
    select(-lat, -long, -yr_restore, -yr_rank) %>% 
    arrange(field_key)
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
#' # Functions 
#' Functions streamline data processing, model fitting, and results output.
#' ### Function: ITS taxonomy
#' This function simplifies and displays the sequence distribution among taxa and
#' across primary lifestyles.
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
#' Several primary lifestyles have been chosen from Fungal Traits for further
#' examination. These lifestyles were the largest by sequence abundance and thought to be
#' most informative given the habitat and questions applied. 
#' 
#' This function filters those groups and tests them among field types. 
#' These tests aren't technically valid due to pseudoreplication, but this analysis can
#' help us find trends worthy of further study. 
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
#' ## 18S-based data (AMF)
#' This function simplifies and displays taxonomic information about the AMF OTUs. 
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
#' ## Re-rarefy in guilds (or groups)
#' To examine trends or differences within subgroups of OTUs, the raw sequence
#' data should be re-rarefied within those groups. Otherwise, especially with low-abundance
#' groups, data and sites may have been lost when the entire species matrix was 
#' rarefied. This function automates the process.
#' 
#' Outputs are:
#' 
#' 1. Sequencing depth used for the subset of OTUs
#' 1. Number of OTUs excluded by rarefying
#' 1. The re-rarefied samples-species matrix
#' 1. The OTU list in long form, with abundances, species, and site metadata
#' 
#+ rerare_function
rerare <- function(spe, meta, grp_var, grp, site) {
    # spe       = species matrix with raw abundances
    # meta      = species metadata matching the OTU list with raw abundances
    # grp_var   = variable name from `meta` desired for grouping and filtering
    #             the OTUs (e.g., `primary_lifestyle`, `family`)
    # grp       = string or factor level name of the group desired from `grp_var`
    # site      = site metadata to combine with sequence abundance long-form
    #             output table
    
    grp_var <- enquo(grp_var)
    
    data <- 
        spe %>% 
        column_to_rownames(var = "field_key") %>% 
        t() %>% 
        as.data.frame() %>% 
        rownames_to_column(var = "otu_num") %>% 
        as_tibble() %>% 
        left_join(meta, by = join_by(otu_num)) %>% 
        filter(!!grp_var == grp) %>% 
        column_to_rownames(var = "otu_num") %>% 
        select(-colnames(meta)[-1]) %>% 
        t() %>% 
        as.data.frame()
    
    depth <- min(rowSums(data))
    rfy <- Rarefy(data)
    zero_otu <- which(apply(rfy$otu.tab.rff, 2, sum) == 0)
    rrfd <- data.frame(rfy$otu.tab.rff[, -zero_otu]) %>%
        rownames_to_column(var = "field_key") %>%
        mutate(field_key = as.numeric(field_key)) %>% 
        arrange(field_key) %>% 
        as_tibble()
    
    rrfd_speTaxa <- 
        rrfd %>% 
        pivot_longer(cols = starts_with("otu"), 
                     names_to = "otu_num", 
                     values_to = "seq_abund") %>% 
        filter(seq_abund > 0) %>% 
        left_join(meta, by = join_by(otu_num)) %>% 
        left_join(site, by = join_by(field_key)) %>% 
        select(-otu_ID)
    
    return(list(
        seq_depth = depth,
        zero_otu_num = length(zero_otu),
        rrfd = rrfd,
        rrfd_speTaxa = rrfd_speTaxa
    ))
    
}
#' 
#' ## Perform Indicator Species Analysis
#' Function `inspan()` takes a combined species and sites data frame and 
#' wrangles it through the analysis to filter OTUs for indicators of field types. 
#' The output is top candidate OTUs joined with species metadata for further analysis. 
#+ inspan_function
inspan <- function(data, np, meta) {
    # data is the samples-species matrix joined with the sites data frame
    # the join aligns the grouping vector with field numbers
    # np is the desired number of permutations
    # meta is the appropriate species metadata table for the original data
    spe <- data.frame(
        data %>% select(field_key, starts_with("otu")),
        row.names = 1
    )
    grp = data$field_type
    mp <- multipatt(
        spe, 
        grp, 
        max.order = 1, 
        control = how(nperm = np))
    si <- mp$sign %>% 
        select(index, stat, p.value) %>% 
        mutate(field_type = case_when(index == 1 ~ "corn", 
                                      index == 2 ~ "restored", 
                                      index == 3 ~ "remnant")) %>% 
        filter(p.value < 0.05) %>% 
        rownames_to_column(var = "otu_num") %>%
        select(-index) %>% 
        as_tibble()
    A  <- mp$A %>% 
        as.data.frame() %>% 
        rownames_to_column(var = "otu_num") %>% 
        pivot_longer(cols = corn:remnant, 
                     names_to = "field_type", 
                     values_to = "A")
    B <- mp$B %>% 
        as.data.frame() %>% 
        rownames_to_column(var = "otu_num") %>% 
        pivot_longer(cols = corn:remnant, 
                     names_to = "field_type", 
                     values_to = "B")
    out <- 
        si %>% 
        left_join(A, by = join_by(otu_num, field_type)) %>% 
        left_join(B, by = join_by(otu_num, field_type)) %>% 
        left_join(meta %>% select(-otu_ID), by = join_by(otu_num)) %>% 
        select(otu_num, A, B, stat, p.value, 
               field_type, primary_lifestyle, everything()) %>% 
        arrange(field_type, -stat)
    
    return(out)
    
}
#' 
#' # Analysis and Results
#' ## ITS sequences in OTU clusters
#' Function outputs are verbose, but details may be necessary later so they are displayed here.
#+ its_tax_trophic_otu,message=FALSE
its_taxaGuild(spe_meta$its_rfy)
#' 
#+ its_guilds_otu,message=FALSE
its_rfy_guilds <- its_test_taxaGuild(spe_meta$its_rfy)
#' 
#' Model tests on `field_type` are technically invalid due to pseudoreplication, but are included here
#' to point out trends that we may be able to present in some other valid way. Trends 
#' with restoration age in Blue Mounds are clearly justified. These trends are:
#' 
#' - Soil saprotroph increases with years since
#' - Wood saprotroph differs among field types and decreases with years since
#' - Plant pathogens decrease with years since
#' 
#' ## ITS-based indicators
#' An indicator species analysis is warranted, identifying which species correlate strongly with `field_type`. 
#' Performing this with all ITS data may identify 
#' particular species to further examine. The analysis should also be done with data re-rarefied into 
#' the guilds identified here, again to showcase particular species which seem to drive differences among
#' field types. It's also of value because this approach avoids the problem we have with pseudoreplication.
#' 
#' With indicator species analysis performed using package [indicspecies](http://sites.google.com/site/miqueldecaceres/),
#' the index values A and B show the specificity and fidelity components of the IndVal combined index. The 
#' combined index value is noted as 'stat' in the output table below.  

# whole dataset doesn't need to be rarefied
its_inspan <- 
    spe$its_rfy %>% 
    left_join(sites, by = join_by(field_key)) %>% 
    inspan(., 1999, meta$its_rfy)


# soil saprotrophs, re-rarefy
ssap <- rerare(spe$its_raw, meta$its_raw, primary_lifestyle, "soil_saprotroph", sites)
ssap$rrfd_speTaxa %>% 
    group_by(field_type, order, field_key) %>% 
    summarize(seq_sum = sum(seq_abund)) %>% 
    group_by(field_type, order) %>% 
    summarize(seq_avg = mean(seq_sum)) %>% 
    mutate(seq_max = max(seq_avg),
           seq_prop = seq_avg / seq_max) %>% group_by(field_type) %>% summarize(sum = sum(seq_prop))
# How to make the columns all the same height? 
    ungroup() %>% 
    filter(seq_prop >= 0.01) %>% 
    ggplot(aes(x = field_type, y = seq_prop)) +
    geom_col(aes(fill = order), color = "black") +
    labs(x = "", y = "Proportion of sequence abundance") +
    scale_fill_discrete_sequential(palette = "Viridis") +
    theme_classic()
    
# recode the low abundance species to "other" so that the bars line up. 
# finally, consider the other basic outputs you might want here

ssap_inspan <- 
    ssap$rrfd %>% 
    left_join(sites, by = join_by(field_key)) %>% 
    inspan(., 1999, meta$its_raw)










wsap <- rerare(spe$its_raw, meta$its_raw, primary_lifestyle, "wood_saprotroph", sites)
wsap_inspan <- 
    wsap$rrfd %>% 
    left_join(sites, by = join_by(field_key)) %>% 
    inspan(., 1999, meta$its_raw)

ppat <- rerare(spe$its_raw, meta$its_raw, primary_lifestyle, "plant_pathogen", sites)
ppat_inspan <- 
    ppat$rrfd %>% 
    left_join(sites, by = join_by(field_key)) %>% 
    inspan(., 1999, meta$its_raw)






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