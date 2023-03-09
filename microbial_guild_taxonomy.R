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
#' The full sequence abundance tables were rarefied to make sequencing depth equivalent
#' across fields. This can result in lower-abundance OTUs dropping to zero. Within guilds, loss 
#' of OTUs could change or bias interpretations of richness, diversity, and composition. We 
#' tried using raw sequence data and rarefying within guilds to address this problem, but 
#' in each case the sequence depth was so small that additional OTUs were lost and abundances were
#' significantly lowered. 
#' 
#' We may try a different approach which is described in [Semchenko et al. 2018](https://www.science.org/doi/10.1126/sciadv.aau4578),
#' but for now, the analysis uses data from the entire rarefied tables for ITS and 18S sequences.
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
                    "GUniFrac",
                    "vegan")
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
#' across primary lifestyles. Use the argument `other_threshold` to 
#' choose a small (e.g., 2, the default) cutoff, below which orders are relabeled as "other".
#+ its_tax_trophic
its_taxaGuild <- function(data, other_threshold=2) {
    # What is the distribution among site types at the class level?
    taxonomy_df <-
        data %>%
        group_by(phylum, order, field_type, field_name) %>%
        summarize(abund = sum(seq_abund), .groups = "drop") %>%
        group_by(phylum, order, field_type) %>%
        summarize(mean = mean(abund) %>% round(., 2),
                  .groups = "drop") %>%
        pivot_wider(
            names_from = field_type,
            values_from = mean,
            values_fill = 0
        ) %>%
        select(phylum, order, corn, restored, remnant) %>%
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
    # Create table
    table <- kable(guild_df, format = "pandoc",
        caption = "Distribution of ITS OTUs by Fungal Trait 'primary_lifestyle'; mean sequence abundance by field type")
    # Plot the most abundant orders across field types
    plot_orders <- 
        data %>% 
        filter(order != is.na(order), order != "unidentified") %>% 
        group_by(field_type, order, field_key) %>% 
        summarize(seq_sum = sum(seq_abund), .groups = "drop_last") %>% 
        summarize(seq_avg = mean(seq_sum), .groups = "drop_last") %>% 
        mutate(seq_comp = (seq_avg / sum(seq_avg)) * 100,
               order = replace(order, which(seq_comp < 2), paste0("Other (OTU<", other_threshold, "%)"))) %>% 
        group_by(field_type, order) %>% 
        summarize(seq_comp = sum(seq_comp), .groups = "drop") %>% 
        ggplot(., aes(x = field_type, y = seq_comp)) +
        geom_col(aes(fill = order), color = "black") +
        labs(x = "", y = "Proportion of sequence abundance",
             title = "Composition of fungi by order") +
        scale_fill_discrete_sequential(name = "Order", palette = "Plasma") +
        theme_classic()
    # Plot the composition of primary lifestyles
    plot_guilds <- 
        data %>% 
        filter(primary_lifestyle != is.na(primary_lifestyle)) %>% 
        group_by(field_type, primary_lifestyle, field_key) %>% 
        summarize(seq_sum = sum(seq_abund), .groups = "drop_last") %>% 
        summarize(seq_avg = mean(seq_sum), .groups = "drop_last") %>% 
        mutate(seq_comp = (seq_avg / sum(seq_avg)) * 100,
               primary_lifestyle = replace(primary_lifestyle, which(seq_comp < 2), paste0("Other (OTU<", other_threshold, "%)"))) %>% 
        group_by(field_type, primary_lifestyle) %>% 
        summarize(seq_comp = sum(seq_comp), .groups = "drop") %>% 
        ggplot(., aes(x = field_type, y = seq_comp)) +
        geom_col(aes(fill = primary_lifestyle), color = "black") +
        labs(x = "", y = "Proportion of sequence abundance",
             title = "Composition of fungi by primary lifestyle") +
        scale_fill_discrete_sequential(name = "Primary lifestyle", palette = "Inferno") +
        theme_classic()
    
    print(list(table,
               plot_orders,
               plot_guilds))
    
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
    pl <- c("soil_saprotroph", "plant_pathogen", "ectomycorrhizal", "wood_saprotroph", "litter_saprotroph")
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
amf_tax <- function(data) {
    cat("---------------------------------\n")
    print(paste("AMF"))
    cat("---------------------------------\n")
    amf_df <-
        data %>%
        group_by(family, field_type, region, field_name, yr_since) %>%
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
#' ## Examine change over time in guilds
#' Function `guiltime()` filters ITS data to a user-specified guild and 
#' produces linear models and plots of change in sequence abundance over time since restoration
#' in Blue Mounds and Fermilab. 
guiltime <- function(pl) {
    d <- spe_meta$its_rfy %>%
        filter(
            primary_lifestyle == pl,
            region %in% c("BM", "FL"),
            field_type == "restored"
        ) %>% 
        group_by(field_key, field_name, region, yr_since) %>% 
        summarize(seq_sum = sum(seq_abund), .groups = "drop")
    
    bm <- summary(
        lm(seq_sum ~ yr_since, data = d %>% filter(region == "BM"))
    )
    fl <- summary(
        lm(seq_sum ~ yr_since, data = d %>% filter(region == "FL"))
    )
    
    fits <- data.frame(
        rbind(BM = c(coef(bm)[1,1], coef(bm)[2,1], coef(bm)[2,4]),
              FL = c(coef(fl)[1,1], coef(fl)[2,1], coef(fl)[2,4]))) %>% 
        mutate(lty = case_when(X3 < 0.05 ~ "a", TRUE ~ NA_character_)) %>% 
        rownames_to_column(var = "region")
    plot <- 
        ggplot(d, aes(x = yr_since, y = seq_sum)) +
        facet_wrap(vars(region), scales = "free_y") +
        geom_point() +
        geom_abline(data = fits, aes(slope = X2, intercept = X1, linetype = lty), color = "blue") +
        labs(x = "Years since restoration", 
             y = "Sum of ITS sequences",
             caption = "Solid line, if present, shows linear relationship at p<0.05") +
        theme_bw() +
        theme(legend.position = "none")
    
    out <- list(
        bm_summary = bm,
        fl_summary = fl,
        plot = plot
    )
    
    print(out)
    
}
#' 
#' ## Re-rarefy in guilds (or groups)
#' To examine richness and composition within subgroups of OTUs, the raw sequence
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
#' ## Filter to guilds or taxonomic groups
#' To examine richness and composition within subgroups of OTUs, the rarefied table
#' must be transposed, filtered, and transposed back. The function `filgu()` or "filter 
#' guilds" automates this process. 
#' 
#' Outputs are:
#' 
#' 1. The resulting samples-species matrix
#' 1. Sequence abundances in long-form, with site and species metadata
#' 
#+ filgu_function
filgu <- function(spe, meta, grp_var, grp, site) {
    # spe       = species matrix with raw abundances
    # meta      = species metadata matching the OTU list with raw abundances
    # grp_var   = variable name from `meta` desired for grouping and filtering
    #             the OTUs (e.g., `primary_lifestyle`, `family`)
    # grp       = string or factor level name of the group desired from `grp_var`
    # site      = site metadata to combine with sequence abundance long-form
    #             output table
    
    grp_var <- enquo(grp_var)
    
    filspe <- 
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
        as.data.frame() %>%
        rownames_to_column(var = "field_key") %>%
        mutate(field_key = as.numeric(field_key)) %>% 
        arrange(field_key) %>% 
        as_tibble()
    
    cs <- colSums(filspe %>% select(-field_key))
    rs <- rowSums(filspe %>% select(-field_key))
    
    hist(cs,
         breaks = length(cs),
         main = "Histogram of OTU sequence sums",
         xlab = "Number of sequences")
    
    hist(rs,
         breaks = length(rs),
         main = "Histogram of sequence abundance in samples",
         xlab = "Number of sequences")
    
    filspeTaxa <- 
        filspe %>% 
        pivot_longer(cols = starts_with("otu"), 
                     names_to = "otu_num", 
                     values_to = "seq_abund") %>% 
        filter(seq_abund > 0) %>% 
        left_join(meta, by = join_by(otu_num)) %>% 
        left_join(site, by = join_by(field_key)) %>% 
        select(-otu_ID)
    
    print(list(
        OTUs_n = length(cs),
        Sites_n = length(rs)
    ))
    
    return(list(
        filspe = filspe,
        filspeTaxa = filspeTaxa
    ))
    
}
#' 
#' ### Calculate Hill's series on a samples-species matrix
#' The objects `$rrfd` from **rerare()** or `$filspe` from **filgu()** can be passed to this function
#+ calc_diversity_function
calc_diversity <- function(spe) {
    spe_mat <- data.frame(spe, row.names = 1)
    
    N0  <- apply(spe_mat > 0, MARGIN = 1, FUN = sum)
    N1  <- exp(diversity(spe_mat))
    N2  <- diversity(spe_mat, "inv")
    E10 <- N1 / N0
    E20 <- N2 / N0
    
    return(
        data.frame(N0, N1, N2, E10, E20) %>%
            rownames_to_column(var = "field_key") %>%
            mutate(field_key = as.integer(field_key)) %>%
            left_join(sites, by = join_by(field_key)) %>%
            pivot_longer(
                cols = N0:E20,
                names_to = "hill_index",
                values_to = "value"
            ) %>%
            mutate(hill_index = factor(
                hill_index,
                ordered = TRUE,
                levels = c("N0", "N1", "N2", "E10", "E20")
            ))
    )
}
#' 
#' ### Results from re-rarefied data
#' After re-rarefying into a guild (or taxonomic group), produce diversity statistics 
#' and calculate percent composition; display results. For plotting, it's convenient 
#' to limit the number of taxonomic orders displayed. Use the argument `other_threshold` to 
#' choose a small (e.g., 2, the default) cutoff, below which orders are relabeled as "other".
gudicom <- function(div, rrfd, grp_var, gene="its", other_threshold=2) {
    hillfield <-     
        ggplot(div, aes(x = field_type, y = value)) +
        facet_wrap(vars(hill_index), scales = "free_y") +
        geom_boxplot(varwidth = TRUE, fill = "gray90", outlier.shape = NA) +
        geom_beeswarm(aes(fill = region), shape = 21, size = 2, dodge.width = 0.2) +
        labs(x = "", y = "Index value", title = paste("Microbial diversity (Hill's):", grp_var),
             caption = "Re-rarefied in the group; N0-richness, N1-e^Shannon, N2-Simpson, E10=N1/N0, E20=N2/N0, width=n") +
        scale_fill_discrete_qualitative(palette = "Dark3") +
        theme_bw()
    hilltime <- 
        div %>% 
        filter(field_type == "restored", region %in% c("BM", "FL")) %>% 
        ggplot(aes(x = yr_since, y = value)) +
        facet_grid(rows = vars(hill_index), cols = vars(region), scales = "free") +
        geom_smooth(method = "lm") +
        geom_point() +
        labs(x = "Years since restoration", y = "Index value", title = paste("Microbial diversity (Hill's) in restored fields:", grp_var),
             caption = "Re-rarefied in the group; N0-richness, N1-e^Shannon, N2-Simpson, E10=N1/N0, E20=N2/N0, width=n") +
        theme_bw()
    if(gene == "its") {
        comp <- 
            rrfd %>% 
            filter(order != is.na(order), order != "unidentified") %>% 
            group_by(field_type, order, field_key) %>% 
            summarize(seq_sum = sum(seq_abund), .groups = "drop_last") %>% 
            summarize(seq_avg = mean(seq_sum), .groups = "drop_last") %>% 
            mutate(seq_comp = (seq_avg / sum(seq_avg)) * 100,
                   order = replace(order, which(seq_comp < other_threshold), paste0("Other (OTU<", other_threshold, "%)"))) %>% 
            group_by(field_type, order) %>% 
            summarize(seq_comp = sum(seq_comp), .groups = "drop")
        comp_plot <-
            ggplot(comp, aes(x = field_type, y = seq_comp)) +
            geom_col(aes(fill = order), color = "black") +
            labs(x = "", y = "Proportion of sequence abundance",
                 title = paste("Composition of", grp_var)) +
            scale_fill_discrete_sequential(name = "Order", palette = "Plasma") +
            theme_classic()
        
        print(list(
            Hills_field_type = hillfield,
            Hills_yrs_since_restoration = hilltime,
            Composition = comp_plot
        ))
        
        return(comp)
        
    } else {
        print(list(
            Hills_field_type = hillfield,
            Hills_yrs_since_restoration = hilltime
        ))
    }
    
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
#' ## ITS sequences
#' Recall the number of OTUs recovered in each dataset. The effect of rarefing did not change
#' richness or diversity very much. 
Map(function(x) ncol(x)-1, spe)
#' 
#' ### Composition in field types
#' Function outputs are verbose, but details may be necessary later so they are displayed here.
#+ its_tax_trophic_otu,message=FALSE,fig.height=7,fig.align='center'
its_taxaGuild(spe_meta$its_rfy)
#' 
#+ its_guilds_otu,message=FALSE
its_rfy_guilds <- its_test_taxaGuild(spe_meta$its_rfy)
#' 
#' Model tests on `field_type` are technically invalid due to pseudoreplication, but are included here
#' to point out trends that we may be able to present in some other valid way. Trends 
#' with restoration age in Blue Mounds are clearly justified. Results are shown in descending 
#' order based on sequence abundance in remnants:
#' 
#' - Soil saprotroph increases with years since
#' - Plant pathogens decrease with years since
#' - Ectomycorrhizal abundance is very low in corn/restored and with 
#' little replication; nothing can be said except that it's relatively abundant in remnants.
#' - Wood saprotroph differs among field types (corn vs. remnant; restored intermediate) and decreases with years since
#' - Litter saprotroph is abundant everywhere, but differences over time or field type are weak.
#' 
#' #### ITS-based indicators
#' An indicator species analysis is warranted, identifying which species correlate strongly with `field_type`. 
#' Performing this with all ITS data may identify particular species to further examine, although it remains
#' a problem that we cannot distinguish field type from an individual field due to pseudoreplication. 
#' 
#' Following the indicator species analysis, richness and composition of selected guilds is 
#' calculated. These calculations are done with data re-rarefied into 
#' the guilds identified here, again to showcase particular species which seem to drive differences among
#' field types. It's also of value because this approach avoids the problem we have with pseudoreplication.
#' 
#' With indicator species analysis performed using package [indicspecies](http://sites.google.com/site/miqueldecaceres/),
#' the index values A and B show the specificity and fidelity components of the IndVal combined index. The 
#' combined index value is noted as 'stat' in the output table below.  
#+ its_inspan_all
its_inspan <- 
    spe$its_rfy %>% 
    left_join(sites, by = join_by(field_key)) %>% 
    inspan(., 1999, meta$its_rfy)
#+ its_inspan_stats
its_inspan %>%
    mutate(field_type = factor(
        field_type,
        ordered = TRUE,
        levels = c("corn", "restored", "remnant")
    )) %>%
    group_by(field_type) %>%
    summarize(
        n_otu = n(),
        stat_avg = mean(stat),
        stat_sd = sd(stat)
    ) %>% 
    kable(format = "pandoc", caption = "Indicator species stats of entire rarefied ITS table")
#' Potential indicators were filtered to p.value<0.05 before this summary was produced. 
#' Cornfields are a restrictive habitat for soil microbes, and that is reflected in the results here.
#' More species have higher specificity and fidelity to cornfields than the other field types. The top ten
#' indicators for each field type are printed here; the entire table is available for further use.
#+ its_inspan_top10
its_inspan %>% 
    mutate(field_type = factor(
    field_type,
    ordered = TRUE,
    levels = c("corn", "restored", "remnant")
)) %>%
    group_by(field_type) %>% 
    slice_max(order_by = stat, n = 10) %>% 
    arrange(field_type, -stat) %>% 
    kable(format = "pandoc", caption = "Indicator species of ITS OTUs (top 10 per field type)")
#' 
#' ### Soil saprotrophs
#' #### Trends over time
#+ ssap_guiltime,fig.width=8,fig.height=4.5,fig.align='center'
guiltime("soil_saprotroph")
#' Sequence abundance of soil saprotrophs increases over time in the Blue Mounds 
#' area ($R^2_{Adj}=0.58, p<0.05$), but
#' this appears to be leveraged by Karla Ott's property, though. With all
#' that big bluestem...maybe there is more litter and soil carbon? It will be good 
#' to look at trends in soil chemistry. 
#' 
#' #### Diversity
#+ ssap_filgu
ssap <- filgu(spe$its_rfy, meta$its_rfy, primary_lifestyle, "soil_saprotroph", sites)
#' Most OTUs contain few sequences, but several range from hundreds to 25,000 sequences.
#' The 25 samples are all retained, and vary from 4000 to 14000 sequences. None are so small that 
#' results would be biased by poor representation bias from being rarefied.
#+ ssap_div
ssap_div <- calc_diversity(ssap$filspe)
#' Diversity measures are stored in this data frame for further use...
#+ ssap_composition,message=FALSE,fig.width=7,fig.height=7,fig.align='center'
ssap_comp <- gudicom(ssap_div, ssap$filspeTaxa, "soil_saprotroph")
#' Richness increases from corn to remnant, but within-group variability is high. Diversity 
#' indices look muddy. Diversity indices increase with years since restoration, but the 
#' significance of this remains to be seen. 
#' 
#' Composition of soil saprotrophs by order can be modified somewhat by choosing the 
#' threshold for lumping rare orders into an "other" category. Leaving this at the default
#' of <2%, nine named orders are left. *Agarics* increase strongly from corn to remnant; *Cystofilobasidiales*
#' and *Filobasidiales* aren't found outside of cornfields. Generally, cornfield composition looks 
#' different than the other two, but remnants do appear somewhat intermediate. *Mortierellales* appear less 
#' in remnants than corn or former corn fields.
#' 
#' #### Indicators
#+ ssap_inspan
ssap_inspan <- 
    ssap$filspe %>% 
    left_join(sites, by = join_by(field_key)) %>% 
    inspan(., 1999, meta$its_raw)
#+ ssap_inspan_stats
ssap_inspan %>%
    mutate(field_type = factor(
        field_type,
        ordered = TRUE,
        levels = c("corn", "restored", "remnant")
    )) %>%
    group_by(field_type) %>%
    summarize(
        n_otu = n(),
        stat_avg = mean(stat),
        stat_sd = sd(stat)
    ) %>% 
    kable(format = "pandoc", caption = "Indicator species stats: soil saprotrophs")
#' We see the same trend as before, where more indicators are found in 
#' cornfields, and their indicator stats are stronger. 
#+ ssap_inspan_table
ssap_inspan %>% 
    mutate(field_type = factor(
        field_type,
        ordered = TRUE,
        levels = c("corn", "restored", "remnant")
    )) %>%
    arrange(field_type, -stat) %>% 
    kable(format = "pandoc", caption = "Indicator species of soil saprotrophs")
#' A later task will be to comb these tables for species with good stories...
#' 
#' ### Plant pathogens 
#' #### Trends over time
#+ ppat_guiltime,fig.width=8,fig.height=4.5,fig.align='center'
guiltime("plant_pathogen")
#' A strong decline in pathogens is seen in Blue Mounds' restored fields 
#' ($R^2_{Adj}=0.75, p<0.01$), and although two distinct groups are apparent, no single
#' site displays undue leverage. It's possible that a signal like this will be found
#' in soil chemistry or plant data and can help explain what we are seeing here. Recall also
#' that AMF were previously found to increase along this same sequence...maybe that will
#' still hold up. 
#' 
#' #### Diversity
#+ ppat_filgu
ppat <- filgu(spe$its_rfy, meta$its_rfy, primary_lifestyle, "plant_pathogen", sites)
#' All samples are retained and contain 2000-12000 sequences, so none are so limited as to bias 
#' results. 
#+ ppat_div
ppat_div <- calc_diversity(ppat$filspe)
#+ ppat_composition,message=FALSE,fig.width=7,fig.height=7,fig.align='center'
ppat_comp <- gudicom(ppat_div, ppat$filspeTaxa, "plant_pathogen", other_threshold = 1)
#' Richness and diversity look flat or declining from corn to remnants and evenness takes a 
#' hit in restored and remnant fields. It looks like we have fewer pathogens, but more dominant 
#' individual taxa become established. Pathogen diversity decreases with years since restoration 
#' in Blue Mounds, but if the dumbbell plots can be believed, the opposite appears true in Fermi.
#' 
#' Many pathogen orders are rare, so the argument `other_threshold` was adjusted to 
#' show more diversity. Shifts don't appear pronounced. *Diaporthales* decreases in composition
#' from corn to remnant while *Hypocreales* pathogens increase. *Cantharellales* appear a
#' small component but are possibly "late successional" pathogens, possibly associated with some
#' native plant in a plant-soil feedback. 
#' 
#' #### Indicators
#+ ppat_inspan
ppat_inspan <- 
    ppat$filspe %>% 
    left_join(sites, by = join_by(field_key)) %>% 
    inspan(., 1999, meta$its_raw)
#+ ppat_inspan_stats
ppat_inspan %>%
    mutate(field_type = factor(
        field_type,
        ordered = TRUE,
        levels = c("corn", "restored", "remnant")
    )) %>%
    group_by(field_type) %>%
    summarize(
        n_otu = n(),
        stat_avg = mean(stat),
        stat_sd = sd(stat)
    ) %>% 
    kable(format = "pandoc", caption = "Indicator species stats: plant pathogens")
#' We see the same trend as before, where more indicators are found in 
#' cornfields, and their indicator stats are stronger. Composition at the level of taxonomic order
#' isn't telling the whole story.
#' 
#' Plant pathogen indicators are nearly all in *Ascomycota.*
#+ ppat_inspan_table
ppat_inspan %>% 
    mutate(field_type = factor(
        field_type,
        ordered = TRUE,
        levels = c("corn", "restored", "remnant")
    )) %>%
    arrange(field_type, -stat) %>% 
    kable(format = "pandoc", caption = "Indicator species of plant pathogens")
#' 
#' ### Wood saprotrophs
#' #### Trends over time
#+ wsap_guiltime,fig.width=8,fig.height=4.5,fig.align='center'
guiltime("wood_saprotroph") 
#' Interestingly a strong negative relationship over time since restoration ($R^2_{Adj}=0.73, p<0.01$)
#' in sharp contrast to the increasing relationship found with soil saprotrophs. Apparently many wood
#' saprotrophs live in cornfield soil...let's see:
#' 
#' #### Diversity
#+ wsap_filgu
wsap <- filgu(spe$its_rfy, meta$its_rfy, primary_lifestyle, "wood_saprotroph", sites)
#' Samples contain 800-4400 sequences. Sequence depth is low; these aren't abundant or numerous taxa.
#' Only 123 OTUs comprise this group. 
#+ wsap_div
wsap_div <- calc_diversity(wsap$filspe)
#+ wsap_composition,message=FALSE,fig.width=7,fig.height=7,fig.align='center'
wasp_comp <- gudicom(wsap_div, wsap$filspeTaxa, "wood_saprotroph")
#' With diversity, not much jumps out. 
#' 
#' Diversity appears high across fields and years compared with other guilds.
#' While *Agaric* soil saprotrophs increased strongly from corn to remnants, 
#' they declined when characterized as wood saprotrophs.
#' 
#' #### Indicators
#+ wsap_inspan
wsap_inspan <- 
    wsap$filspe %>% 
    left_join(sites, by = join_by(field_key)) %>% 
    inspan(., 1999, meta$its_raw)
#+ wsap_inspan_stats
wsap_inspan %>%
    mutate(field_type = factor(
        field_type,
        ordered = TRUE,
        levels = c("corn", "restored", "remnant")
    )) %>%
    group_by(field_type) %>%
    summarize(
        n_otu = n(),
        stat_avg = mean(stat),
        stat_sd = sd(stat)
    ) %>% 
    kable(format = "pandoc", caption = "Indicator species stats: wood saprotrophs")
#' Few species show specificity or fidelity. Corn fields have a few unusual taxa, though. 
#' Less so with remnants, and none with restored fields. 
#+ wsap_inspan_table
wsap_inspan %>% 
    mutate(field_type = factor(
        field_type,
        ordered = TRUE,
        levels = c("corn", "restored", "remnant")
    )) %>%
    arrange(field_type, -stat) %>% 
    kable(format = "pandoc", caption = "Indicator species of wood saprotrophs")
#' 
#' ### Litter saprotrophs
#' #### Trends over time
#+ lsap_guiltime,fig.width=8,fig.height=4.5,fig.align='center'
guiltime("litter_saprotroph") 
#' 
#' #### Diversity
#+ lsap_filgu
lsap <- filgu(spe$its_rfy, meta$its_rfy, primary_lifestyle, "litter_saprotroph", sites)
#' Slightly more numerous than the wood saprotrophs, but similarly not abundant or numerous. Recall that 
#' when this group was rarefied in the guild, sampling depth was 297, or an order of magnitude less 
#' than what we have here. Several OTUs were lost. 
#+ lsap_div
lsap_div <- calc_diversity(lsap$filspe)
#+ lsap_composition,message=FALSE,fig.width=7,fig.height=7,fig.align='center'
lsap_comp <- gudicom(lsap_div, lsap$filspeTaxa, "litter_saprotroph")
#' With no litter in cornfields, it's perhaps not surprising to see increasing trends across field types
#' with this guild. Trends over time aren't convincing, except possibly in Fermi.
#' 
#' #### Indicators
#+ lsap_inspan
lsap_inspan <- 
    lsap$filspe %>% 
    left_join(sites, by = join_by(field_key)) %>% 
    inspan(., 1999, meta$its_raw)
#+ lsap_inspan_stats
lsap_inspan %>%
    mutate(field_type = factor(
        field_type,
        ordered = TRUE,
        levels = c("corn", "restored", "remnant")
    )) %>%
    group_by(field_type) %>%
    summarize(
        n_otu = n(),
        stat_avg = mean(stat),
        stat_sd = sd(stat)
    ) %>% 
    kable(format = "pandoc", caption = "Indicator species stats: litter saprotrophs")
#+ lsap_inspan_table
lsap_inspan %>% 
    mutate(field_type = factor(
        field_type,
        ordered = TRUE,
        levels = c("corn", "restored", "remnant")
    )) %>%
    arrange(field_type, -stat) %>% 
    kable(format = "pandoc", caption = "Indicator species of litter saprotrophs")
#' 





# 2023-03-08 this is where I left off





#' ## AMF
#' Function output is verbose but retained as explained previously.
#+ amf_otu_summary,message=FALSE
amf_summary <- amf_tax(spe_meta$amf_rfy)
#' 
#' Let's look at abundances across field types in four families:
#+ amf_boxplot,fig.width=7,fig.height=6,fig.align='center'
amf_summary %>% 
    filter(family %in% c("Claroideoglomeraceae", "Paraglomeraceae", "Diversisporaceae", "Gigasporaceae")) %>% 
    ggplot(aes(x = field_type, y = seq_sum)) +
    facet_wrap(vars(family), scales = "free_y") +
    geom_boxplot(varwidth = TRUE, fill = "gray90", outlier.shape = NA) +
    geom_beeswarm(aes(fill = region), shape = 21, size = 2, dodge.width = 0.2) +
    labs(x = "", y = "Sum of sequence abundance", title = "AMF abundance in families and field types") +
    scale_fill_discrete_qualitative(palette = "Dark3") +
    theme_bw()
#' 
#' Let's also look at change over time in Fermi and Blue Mounds:
#+ amf_change_time_all,message=FALSE,fig.height=8,fig.width=5,fig.align='center'
amf_summary %>% 
    filter(family %in% c("Claroideoglomeraceae", "Paraglomeraceae", "Diversisporaceae", "Gigasporaceae"),
           region %in% c("BM", "FL"),
           field_type == "restored") %>% 
    ggplot(aes(x = yr_since, y = seq_sum)) +
    facet_grid(rows = vars(family), cols = vars(region), scales = "free") +
    geom_smooth(method = "lm") +
    geom_point() +
    labs(x = "", y = "Sum of sequences abundances") +
    theme_bw()
#' 
#' And finally composition at the family level. *Glomeraceae* truly dominates.
#+ amf_composition,fig.height=7,fig.align='center'
amf_summary %>% 
    group_by(field_type, family) %>% 
    summarize(seq_avg = mean(seq_sum), .groups = "drop_last") %>% 
    mutate(seq_comp = (seq_avg / sum(seq_avg)) * 100,
           order = replace(family, which(seq_comp < 0), "Other")) %>% 
    group_by(field_type, order) %>% 
    summarize(seq_comp = sum(seq_comp), .groups = "drop") %>% 
    ggplot(aes(x = field_type, y = seq_comp)) +
    geom_col(aes(fill = order), color = "black") +
    labs(x = "", y = "Proportion of sequence abundance",
         title = "Composition of AMF by order") +
    scale_fill_discrete_sequential(name = "Order", palette = "Plasma") +
    theme_classic()
    
    

    
    
    

#' From the mean sequence abundances in field types and trends over time, the following families look interesting:
#' 
#' - *Claroideoglomeraceae:* low in corn; significantly by likelihood ratio test
#' - *Paraglomeraceae:* highest in corn, declines through restoration and remnant, declines in BM and FL but 
#' likely not a significant trend
#' - *Diversisporaceae:* highest in corn, declines through restoration and remnant
#' - *Gigasporaceae:* low in corn, and also the only one with a significant change with years
#' since restoration, and this only in Blue Mounds. These increase over time (recall that pathogens decline over time). 
#' 
#' In the next section, we will examine these families more closely by first re-rarefying abundances within
#' families.
#' 
#' ### Claroideoglomeraceae
#' 
#' #### Diversity

# The re-rarefy thing isn't working. Pause to figure out why. 
# You can still do diversity and so forth below.

#' Sequencing depth of 290, perhaps too rare to justify examination.
#+ claroid_div
# claroid_div <- calc_diversity(claroid$rrfd)


#+ claroid_divplot,message=FALSE,results=FALSE,fig.width=7,fig.height=7,fig.align='center'
# gudicom(claroid_div, claroid$rrfd_speTaxa, "Claroideoglomeraceae", gene = "amf")


#' With no litter in cornfields, it's perhaps not surprising to see increasing trends across field types
#' with this guild. Trends over time aren't convincing, except possibly in Fermi.
#' 



#' 
#' # Conclusions: taxa and guilds
#' Little variation exists here for ITS or AMF sequences among field types, although 
#' classes of fungi identified through ITS sequences remain to be closely examined. 
#' It's striking that plant pathogens decline as restorations age while 
#' the AMF family *Gigasporaceae* increases, but this contrast was not found in any 
#' other group of AMF and the *Gigasporaceae* aren't particularly abundant to begin with.
#' 

