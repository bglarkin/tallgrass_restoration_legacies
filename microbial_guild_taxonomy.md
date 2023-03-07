Microbial data: microbial guilds and taxonomy
================
Beau Larkin

Last updated: 06 March, 2023

- <a href="#description" id="toc-description">Description</a>
- <a href="#packages-and-libraries"
  id="toc-packages-and-libraries">Packages and libraries</a>
- <a href="#data" id="toc-data">Data</a>
  - <a href="#sites-species-tables"
    id="toc-sites-species-tables">Sites-species tables</a>
  - <a href="#species-metadata" id="toc-species-metadata">Species
    metadata</a>
  - <a href="#site-metadata-and-design"
    id="toc-site-metadata-and-design">Site metadata and design</a>
  - <a href="#joined-species-metadata-and-design-tables"
    id="toc-joined-species-metadata-and-design-tables">Joined species,
    metadata, and design tables</a>
- <a href="#functions" id="toc-functions">Functions</a>
  - <a href="#function-its-and-guilds"
    id="toc-function-its-and-guilds">Function: ITS and guilds</a>
  - <a href="#18s-based-data-amf" id="toc-18s-based-data-amf">18S-based data
    (AMF)</a>
  - <a href="#examine-change-over-time-in-guilds"
    id="toc-examine-change-over-time-in-guilds">Examine change over time in
    guilds</a>
  - <a href="#re-rarefy-in-guilds-or-groups"
    id="toc-re-rarefy-in-guilds-or-groups">Re-rarefy in guilds (or
    groups)</a>
    - <a href="#calculate-hills-series-on-a-samples-species-matrix"
      id="toc-calculate-hills-series-on-a-samples-species-matrix">Calculate
      Hill’s series on a samples-species matrix</a>
    - <a href="#results-from-re-rarefied-data"
      id="toc-results-from-re-rarefied-data">Results from re-rarefied data</a>
  - <a href="#perform-indicator-species-analysis"
    id="toc-perform-indicator-species-analysis">Perform Indicator Species
    Analysis</a>
- <a href="#analysis-and-results" id="toc-analysis-and-results">Analysis
  and Results</a>
  - <a href="#its-sequences" id="toc-its-sequences">ITS sequences</a>
    - <a href="#soil-saprotrophs" id="toc-soil-saprotrophs">Soil
      saprotrophs</a>
    - <a href="#plant-pathogens" id="toc-plant-pathogens">Plant pathogens</a>
    - <a href="#wood-saprotrophs" id="toc-wood-saprotrophs">Wood
      saprotrophs</a>
    - <a href="#litter-saprotrophs" id="toc-litter-saprotrophs">Litter
      saprotrophs</a>
  - <a href="#amf-otus" id="toc-amf-otus">AMF OTUs</a>
- <a href="#conclusions-taxa-and-guilds"
  id="toc-conclusions-taxa-and-guilds">Conclusions: taxa and guilds</a>

# Description

Sequence clusters identified in QIIME2 are annotated with taxonomic
information and metadata from [Fungal
traits](https://link.springer.com/article/10.1007/s13225-020-00466-2).
In this report, sequence abundances in taxonomic groups or fungal guilds
are compared across field types and with time since restoration.

The full sequence abundance tables were rarefied to make sequencing
depth equivalent across fields. This can result in lower-abundance OTUs
dropping to zero. Here, when richness and composition within guilds is
calculated, the raw abundances are used, filtered into particular
guilds, and then rarefied within those guilds to produce accurate
results.

# Packages and libraries

``` r
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
```

``` r
if (any(!packages_installed)) {
    install.packages(packages_needed[!packages_installed])
}
```

``` r
for (i in 1:length(packages_needed)) {
    library(packages_needed[i], character.only = T)
}
```

``` r
conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")
```

# Data

## Sites-species tables

CSV files were produced in `process_data.R`

``` r
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
```

## Species metadata

Load taxonomy for all and guilds (`primary_lifestyle` in Fungal Traits)
for ITS OTUs.

``` r
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
```

## Site metadata and design

``` r
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
```

## Joined species, metadata, and design tables

Functions streamline this process

``` r
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
```

``` r
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
```

# Functions

Functions streamline data processing, model fitting, and results output.
\### Function: ITS taxonomy This function simplifies and displays the
sequence distribution among taxa and across primary lifestyles. Use the
argument `other_threshold` to choose a small (e.g., 2, the default)
cutoff, below which orders are relabeled as “other”.

``` r
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
    plot <- 
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
             title = "Composition of fungi") +
        scale_fill_discrete_sequential(name = "Order", palette = "Plasma") +
        theme_classic()
    
    print(list(table, plot))
    
}
```

### Function: ITS and guilds

Several primary lifestyles have been chosen from Fungal Traits for
further examination. These lifestyles were the largest by sequence
abundance and thought to be most informative given the habitat and
questions applied.

This function filters those groups and tests them among field types.
These tests aren’t technically valid due to pseudoreplication, but this
analysis can help us find trends worthy of further study.

``` r
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
```

## 18S-based data (AMF)

This function simplifies and displays taxonomic information about the
AMF OTUs.

``` r
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
```

## Examine change over time in guilds

Function `guiltime()` filters ITS data to a user-specified guild and
produces linear models and plots of change in sequence abundance over
time since restoration in Blue Mounds and Fermilab.

``` r
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
```

## Re-rarefy in guilds (or groups)

To examine richness and composition within subgroups of OTUs, the raw
sequence data should be re-rarefied within those groups. Otherwise,
especially with low-abundance groups, data and sites may have been lost
when the entire species matrix was rarefied. This function automates the
process.

Outputs are:

1.  Sequencing depth used for the subset of OTUs
2.  Number of OTUs excluded by rarefying
3.  The re-rarefied samples-species matrix
4.  The OTU list in long form, with abundances, species, and site
    metadata

``` r
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
```

### Calculate Hill’s series on a samples-species matrix

The object `$rrfd` from rerare() can be passed to this function

``` r
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
```

### Results from re-rarefied data

After re-rarefying into a guild (or taxonomic group), produce diversity
statistics and calculate percent composition; display results. For
plotting, it’s convenient to limit the number of taxonomic orders
displayed. Use the argument `other_threshold` to choose a small (e.g.,
2, the default) cutoff, below which orders are relabeled as “other”.

``` r
gudicom <- function(div, rrfd, grp_var, other_threshold=2) {
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
    
}
```

## Perform Indicator Species Analysis

Function `inspan()` takes a combined species and sites data frame and
wrangles it through the analysis to filter OTUs for indicators of field
types. The output is top candidate OTUs joined with species metadata for
further analysis.

``` r
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
```

# Analysis and Results

## ITS sequences

Function outputs are verbose, but details may be necessary later so they
are displayed here.

``` r
its_taxaGuild(spe_meta$its_rfy)
```

    ## 
    ## 
    ## Table: Distribution of ITS OTUs in classes; mean sequence abundance by field type
    ## 
    ## phylum                             order                                        corn   restored   remnant
    ## ---------------------------------  --------------------------------------  ---------  ---------  --------
    ## Ascomycota                         Hypocreales                               7295.00    6209.00   6917.25
    ## Ascomycota                         Chaetothyriales                            624.00    5573.44   5824.50
    ## Basidiomycota                      Agaricales                                2766.40    2082.00   5000.75
    ## Ascomycota                         Pleosporales                              6398.40    6405.56   4929.00
    ## Ascomycota                         Sordariales                              10384.80    4295.44   3317.25
    ## Ascomycota                         Helotiales                                2631.60    3266.38   3221.25
    ## Ascomycota                         NA                                        1267.80    2366.94   2874.25
    ## NA                                 NA                                        1215.00    1911.19   2318.00
    ## Ascomycota                         Onygenales                                  69.40    1267.56   2063.75
    ## Mortierellomycota                  Mortierellales                            3310.80    2572.88   1641.50
    ## Basidiomycota                      Thelephorales                                2.00      26.43   1041.75
    ## Ascomycota                         Geoglossales                                 5.67    1461.58    912.00
    ## Ascomycota                         GS34                                         0.00      72.00    692.00
    ## Ascomycota                         Glomerellales                             1740.00    1117.06    660.50
    ## Ascomycota                         unidentified                                 0.00     828.80    608.75
    ## Basidiomycota                      Cantharellales                             264.80     754.75    551.00
    ## Glomeromycota                      Glomerales                                  78.80     414.44    421.00
    ## Ascomycota                         Coniochaetales                             709.80     162.12    300.25
    ## Ascomycota                         Sordariomycetes_ord_Incertae_sedis          38.80      51.44    239.75
    ## Basidiomycota                      Russulales                                   4.00       3.50    182.33
    ## Ascomycota                         Xylariales                                  63.20     285.88    174.00
    ## Ascomycota                         Pezizales                                  927.60     311.31    168.00
    ## Basidiomycota                      NA                                          75.80     548.12    165.25
    ## Ascomycota                         Magnaporthales                              76.80     118.88    155.75
    ## Basidiomycota                      Sebacinales                                 27.20     621.31    155.50
    ## Ascomycota                         Capnodiales                                542.80     505.50    147.75
    ## Basidiomycota                      Boletales                                    5.00       8.00    143.00
    ## Glomeromycota                      NA                                           5.00     144.23    128.00
    ## Ascomycota                         Minutisphaerales                             0.00      58.00    106.50
    ## Ascomycota                         Chaetosphaeriales                          259.50     257.40    101.67
    ## Ascomycota                         Branch06                                    10.00     132.57     88.67
    ## Mucoromycota                       NA                                           0.00      24.00     83.00
    ## Ascomycota                         Mytilinidales                                0.00       0.00     82.00
    ## Basidiomycota                      Trichosporonales                           115.00      39.44     76.33
    ## Basidiomycota                      Auriculariales                             105.00     223.25     76.00
    ## Basidiomycota                      Tremellales                                 15.00      90.81     75.25
    ## Basidiomycota                      Filobasidiales                             809.50     292.21     71.67
    ## Chytridiomycota                    Spizellomycetales                          164.40      87.27     68.50
    ## Mucoromycota                       Umbelopsidales                               0.00       2.00     64.00
    ## Basidiomycota                      Cystofilobasidiales                       2273.80      84.15     62.50
    ## Ascomycota                         Mytilinidiales                               0.00      14.00     62.00
    ## Ascomycota                         Thelebolales                                90.75      35.69     62.00
    ## Basidiomycota                      Hymenochaetales                             11.80     179.43     60.50
    ## Ascomycota                         Dothideomycetes_ord_Incertae_sedis           0.00       0.00     52.00
    ## Chytridiomycota                    Rhizophlyctidales                          211.80     107.47     51.75
    ## Ascomycota                         Venturiales                                 31.00      72.92     45.67
    ## Basidiomycota                      Ustilaginales                                2.00      84.75     45.33
    ## Ascomycota                         Verrucariales                                0.00       0.00     45.00
    ## Ascomycota                         Orbiliales                                  22.67      81.29     31.33
    ## Basidiomycota                      Atheliales                                   0.00     148.50     30.00
    ## Basidiomycota                      Tremellodendropsidales                       9.00      28.08     30.00
    ## Basidiomycota                      Polyporales                                 19.80      28.00     28.25
    ## Basidiomycota                      Geminibasidiales                            63.00      76.75     28.00
    ## Ascomycota                         Tubeufiales                                 64.60     166.12     27.00
    ## Ascomycota                         Saccharomycetales                          263.50      38.10     26.00
    ## Basidiomycota                      Trechisporales                             118.40     364.69     24.67
    ## Ascomycota                         Myrmecridiales                               0.00      71.00     24.00
    ## Ascomycota                         Microascales                               108.00      72.45     23.00
    ## Chytridiomycota                    Chytridiales                                 5.00      85.25     21.00
    ## Rozellomycota                      GS11                                         0.00       2.00     20.00
    ## Basidiomycota                      Phallales                                  213.25      41.38     19.33
    ## Ascomycota                         GS32                                         0.00       0.00     19.00
    ## Basidiomycota                      Erythrobasidiales                            0.00       4.00     19.00
    ## Ascomycota                         Diaporthales                               124.00      10.43     14.50
    ## Chytridiomycota                    Rhizophydiales                              19.75      16.00     14.50
    ## Ascomycota                         Eurotiales                                 105.00      35.40     14.25
    ## Ascomycota                         Ostropales                                   0.00      52.00     13.00
    ## Basidiomycota                      Leucosporidiales                            33.33       8.57     12.50
    ## Ascomycota                         Savoryellales                               14.50      18.67     12.00
    ## Ascomycota                         Archaeorhizomycetales                        0.00      33.00     11.00
    ## Ascomycota                         Acrospermales                                0.00       8.00     10.00
    ## Glomeromycota                      Archaeosporales                              3.00       8.50      9.00
    ## Chlorophyta                        Chaetopeltidales                            11.50       5.50      7.50
    ## Mortierellomycota                  NA                                           0.00       4.50      7.00
    ## Ascomycota                         Rhytismatales                                0.00       1.00      6.00
    ## Basidiobolomycota                  Basidiobolales                               6.00      15.67      6.00
    ## Basidiomycota                      Agaricomycetes_ord_Incertae_sedis            0.00       0.00      6.00
    ## Ascomycota                         Dothideales                                  0.00      29.38      5.00
    ## Basidiomycota                      Geastrales                                  63.00      50.00      5.00
    ## Chlorophyta                        NA                                           9.00       7.50      5.00
    ## Anthophyta                         Poales                                       3.00       3.33      4.00
    ## Basidiomycota                      Microbotryomycetes_ord_Incertae_sedis        7.67      11.00      4.00
    ## Basidiomycota                      unidentified                                11.00     119.75      4.00
    ## Basidiomycota                      Atractiellales                              13.50       0.00      3.00
    ## Glomeromycota                      Diversisporales                              0.00       7.75      2.50
    ## Ichthyosporia_phy_Incertae_sedis   unidentified                                 0.00       0.00      2.00
    ## Mucoromycota                       Mucorales                                    0.00       5.83      2.00
    ## Ascomycota                         Candelariales                                0.00       4.00      1.00
    ## Entorrhizomycota                   Entorrhizales                                0.00       0.00      1.00
    ## Anthophyta                         Asterales                                    0.00       5.00      0.00
    ## Anthophyta                         Brassicales                                  8.00       5.00      0.00
    ## Anthophyta                         Commelinales                                 5.00     217.00      0.00
    ## Anthophyta                         Fabales                                      0.00      15.00      0.00
    ## Ascomycota                         Boliniales                                  38.50      44.50      0.00
    ## Ascomycota                         Botryosphaeriales                           13.50      19.80      0.00
    ## Ascomycota                         Jahnulales                                  11.00       0.00      0.00
    ## Ascomycota                         Pezizomycotina_ord_Incertae_sedis           32.50     489.00      0.00
    ## Ascomycota                         Phacidiales                                  0.00       5.50      0.00
    ## Ascomycota                         Phomatosporales                            665.00      10.67      0.00
    ## Ascomycota                         Trichosphaeriales                            5.00      19.50      0.00
    ## Basidiomycota                      Agaricostilbales                             2.00       2.00      0.00
    ## Basidiomycota                      Corticiales                                  0.00      43.43      0.00
    ## Basidiomycota                      Cystobasidiales                             32.50      11.67      0.00
    ## Basidiomycota                      Entylomatales                                0.00      11.20      0.00
    ## Basidiomycota                      Holtermanniales                              4.00      14.00      0.00
    ## Basidiomycota                      Kriegeriales                                 0.00       4.00      0.00
    ## Basidiomycota                      Platygloeales                                0.00      48.25      0.00
    ## Basidiomycota                      Pucciniales                                  0.00       9.00      0.00
    ## Basidiomycota                      Sporidiobolales                             44.00      16.25      0.00
    ## Basidiomycota                      Tilletiales                                  0.00      14.50      0.00
    ## Basidiomycota                      Urocystidales                               58.00       4.25      0.00
    ## Calcarisporiellomycota             Calcarisporiellales                          0.00       6.00      0.00
    ## Cercozoa                           unidentified                                 7.00       2.00      0.00
    ## Chlorophyta                        Chaetophorales                               0.00      23.00      0.00
    ## Chlorophyta                        Chlamydomonadales                            0.00       2.00      0.00
    ## Chlorophyta                        Chlorellales                                 0.00       6.00      0.00
    ## Chlorophyta                        Sphaeropleales                               0.00       5.00      0.00
    ## Chytridiomycota                    unidentified                                 8.00       2.00      0.00
    ## Chytridiomycota                    NA                                          14.00       8.00      0.00
    ## Glomeromycota                      Paraglomerales                               0.00      15.00      0.00
    ## Glomeromycota                      unidentified                                 2.00      20.40      0.00
    ## Haplosporidia                      Haplosporidia_ord_Incertae_sedis             3.00      10.00      0.00
    ## Mucoromycota                       GS22                                         0.00       4.00      0.00
    ## [[1]]
    ## 
    ## 
    ## Table: Distribution of ITS OTUs by Fungal Trait 'primary_lifestyle'; mean sequence abundance by field type
    ## 
    ## primary_lifestyle            corn   restored   remnant
    ## -----------------------  --------  ---------  --------
    ## NA                        20151.0    24785.1   27184.5
    ## soil_saprotroph            7123.4     5736.2    6067.2
    ## plant_pathogen             7068.4     7165.0    5622.2
    ## ectomycorrhizal              11.0      142.8    1979.2
    ## wood_saprotroph            2841.6     2290.1    1127.2
    ## litter_saprotroph          2063.2     1434.7    1124.2
    ## dung_saprotroph            2735.0     1559.6     936.0
    ## animal_parasite             660.6     1082.9     649.2
    ## mycoparasite               1820.2      598.2     213.2
    ## unspecified_saprotroph      773.0      159.0     175.8
    ## arbuscular_mycorrhizal       54.0      188.4     153.5
    ## root_endophyte               16.0      329.7     134.8
    ## pollen_saprotroph           152.4       75.6      52.5
    ## nectar/tap_saprotroph        37.5       36.0      26.0
    ## lichen_parasite              12.3       55.2      23.5
    ## epiphyte                      0.0        6.0      23.0
    ## lichenized                    0.0      100.0      13.0
    ## unspecified_pathotroph        0.0       13.2      10.5
    ## foliar_endophyte              3.0       31.5       4.0
    ## algal_parasite                2.0        4.2       0.0
    ## 
    ## [[2]]

<img src="microbial_guild_taxonomy_files/figure-gfm/its_tax_trophic_otu-1.png" style="display: block; margin: auto;" />

``` r
its_rfy_guilds <- its_test_taxaGuild(spe_meta$its_rfy)
```

    ## ---------------------------------
    ## [1] "soil_saprotroph"
    ## ---------------------------------
    ## 
    ## 
    ## field_type   region   field_name    yr_since  primary_lifestyle    seq_sum
    ## -----------  -------  -----------  ---------  ------------------  --------
    ## corn         FL       FLC2                 0  soil_saprotroph        12046
    ## restored     BM       KORP1               28  soil_saprotroph        10616
    ## remnant      LP       LPREM1              NA  soil_saprotroph        10360
    ## restored     FL       FLRSP1              10  soil_saprotroph         9609
    ## restored     FL       FLRP5               35  soil_saprotroph         8813
    ## corn         FL       FLC1                 0  soil_saprotroph         7845
    ## corn         LP       LPC1                 0  soil_saprotroph         7329
    ## restored     BM       PHRP1               11  soil_saprotroph         7306
    ## restored     FL       FLRSP3              10  soil_saprotroph         7074
    ## remnant      BM       MBREM1              NA  soil_saprotroph         6558
    ## restored     BM       BBRP1               16  soil_saprotroph         6020
    ## corn         BM       PHC1                 0  soil_saprotroph         5167
    ## restored     FL       FLRP1               40  soil_saprotroph         5093
    ## restored     FL       FLRSP2              10  soil_saprotroph         4837
    ## restored     FL       FLRP4               36  soil_saprotroph         4550
    ## restored     BM       MBRP1               18  soil_saprotroph         4508
    ## remnant      FG       FGREM1              NA  soil_saprotroph         4427
    ## restored     LP       LPRP1                4  soil_saprotroph         4302
    ## restored     BM       ERRP1                3  soil_saprotroph         4171
    ## restored     BM       MHRP2                2  soil_saprotroph         4073
    ## restored     LP       LPRP2                4  soil_saprotroph         3956
    ## restored     BM       MHRP1                7  soil_saprotroph         3835
    ## corn         FG       FGC1                 0  soil_saprotroph         3230
    ## restored     FG       FGRP1               15  soil_saprotroph         3017
    ## remnant      FL       FLREM1              NA  soil_saprotroph         2924
    ## ----------------------------------------------------
    ## 
    ## Linear mixed model fit by maximum likelihood  ['lmerMod']
    ## Formula: seq_sum ~ field_type + (1 | region)
    ##    Data: mod_data
    ##       AIC       BIC    logLik  deviance  df.resid 
    ##  471.6598  477.7542 -230.8299  461.6598        20 
    ## Random effects:
    ##  Groups   Name        Std.Dev.
    ##  region   (Intercept)    0    
    ##  Residual             2476    
    ## Number of obs: 25, groups:  region, 4
    ## Fixed Effects:
    ##  (Intercept)  field_type.L  field_type.Q  
    ##       6309.0        -746.8         701.4  
    ## optimizer (nloptwrap) convergence code: 0 (OK) ; 0 optimizer warnings; 1 lme4 warnings 
    ## ----------------------------------------------------
    ## 
    ## Linear mixed model fit by maximum likelihood  ['lmerMod']
    ## Formula: seq_sum ~ 1 + (1 | region)
    ##    Data: mod_data
    ##       AIC       BIC    logLik  deviance  df.resid 
    ##  468.8281  472.4847 -231.4140  462.8281        22 
    ## Random effects:
    ##  Groups   Name        Std.Dev.
    ##  region   (Intercept)    0    
    ##  Residual             2534    
    ## Number of obs: 25, groups:  region, 4
    ## Fixed Effects:
    ## (Intercept)  
    ##        6067  
    ## optimizer (nloptwrap) convergence code: 0 (OK) ; 0 optimizer warnings; 1 lme4 warnings 
    ## ----------------------------------------------------
    ## 
    ## Data: mod_data
    ## Models:
    ## mmod_null: seq_sum ~ 1 + (1 | region)
    ## mmod: seq_sum ~ field_type + (1 | region)
    ##           npar    AIC    BIC  logLik deviance  Chisq Df Pr(>Chisq)
    ## mmod_null    3 468.83 472.48 -231.41   462.83                     
    ## mmod         5 471.66 477.75 -230.83   461.66 1.1683  2     0.5576
    ## ----------------------------------------------------
    ## 
    ## 
    ##   Simultaneous Tests for General Linear Hypotheses
    ## 
    ## Multiple Comparisons of Means: Tukey Contrasts
    ## 
    ## 
    ## Fit: lmer(formula = seq_sum ~ field_type + (1 | region), data = mod_data, 
    ##     REML = FALSE)
    ## 
    ## Linear Hypotheses:
    ##                         Estimate Std. Error z value Pr(>|z|)
    ## restored - corn == 0       -1387       1268  -1.094    0.513
    ## remnant - corn == 0        -1056       1661  -0.636    0.797
    ## remnant - restored == 0      331       1384   0.239    0.968
    ## (Adjusted p values reported -- single-step method)
    ## 
    ##     corn restored  remnant 
    ##      "a"      "a"      "a" 
    ## ----------------------------------------------------
    ## 
    ## [1] "Years since restoration and soil_saprotroph sequence abundance in Blue Mounds Area"
    ## 
    ## Call:
    ## lm(formula = seq_sum ~ yr_since, data = mod_data2)
    ## 
    ## Residuals:
    ##       1       2       3       4       5       6       7 
    ##  -597.3   342.4  1424.6 -2538.3  -851.6   458.9  1761.3 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)  
    ## (Intercept)  3185.04    1055.79   3.017   0.0295 *
    ## yr_since      214.51      71.02   3.020   0.0294 *
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 1611 on 5 degrees of freedom
    ## Multiple R-squared:  0.646,  Adjusted R-squared:  0.5752 
    ## F-statistic: 9.123 on 1 and 5 DF,  p-value: 0.0294
    ## 
    ## 
    ## 
    ## 
    ## ---------------------------------
    ## [1] "plant_pathogen"
    ## ---------------------------------
    ## 
    ## 
    ## field_type   region   field_name    yr_since  primary_lifestyle    seq_sum
    ## -----------  -------  -----------  ---------  ------------------  --------
    ## restored     BM       MHRP2                2  plant_pathogen         12342
    ## restored     BM       MHRP1                7  plant_pathogen         12155
    ## restored     LP       LPRP1                4  plant_pathogen         11571
    ## restored     BM       PHRP1               11  plant_pathogen         11314
    ## restored     FG       FGRP1               15  plant_pathogen         10631
    ## corn         LP       LPC1                 0  plant_pathogen         10322
    ## restored     BM       ERRP1                3  plant_pathogen          9998
    ## remnant      LP       LPREM1              NA  plant_pathogen          8111
    ## corn         FG       FGC1                 0  plant_pathogen          7598
    ## corn         FL       FLC2                 0  plant_pathogen          6751
    ## remnant      FG       FGREM1              NA  plant_pathogen          6293
    ## restored     BM       BBRP1               16  plant_pathogen          6215
    ## restored     LP       LPRP2                4  plant_pathogen          6111
    ## remnant      FL       FLREM1              NA  plant_pathogen          5713
    ## corn         FL       FLC1                 0  plant_pathogen          5469
    ## restored     FL       FLRP1               40  plant_pathogen          5444
    ## corn         BM       PHC1                 0  plant_pathogen          5202
    ## restored     FL       FLRSP2              10  plant_pathogen          5168
    ## restored     FL       FLRSP1              10  plant_pathogen          4884
    ## restored     BM       MBRP1               18  plant_pathogen          4364
    ## restored     FL       FLRP5               35  plant_pathogen          3994
    ## restored     FL       FLRP4               36  plant_pathogen          3797
    ## restored     BM       KORP1               28  plant_pathogen          3552
    ## restored     FL       FLRSP3              10  plant_pathogen          3100
    ## remnant      BM       MBREM1              NA  plant_pathogen          2372
    ## ----------------------------------------------------
    ## 
    ## Linear mixed model fit by maximum likelihood  ['lmerMod']
    ## Formula: seq_sum ~ field_type + (1 | region)
    ##    Data: mod_data
    ##       AIC       BIC    logLik  deviance  df.resid 
    ##  478.2810  484.3754 -234.1405  468.2810        20 
    ## Random effects:
    ##  Groups   Name        Std.Dev.
    ##  region   (Intercept) 1334    
    ##  Residual             2625    
    ## Number of obs: 25, groups:  region, 4
    ## Fixed Effects:
    ##  (Intercept)  field_type.L  field_type.Q  
    ##       6871.4       -1260.3        -876.9  
    ## ----------------------------------------------------
    ## 
    ## Linear mixed model fit by maximum likelihood  ['lmerMod']
    ## Formula: seq_sum ~ 1 + (1 | region)
    ##    Data: mod_data
    ##       AIC       BIC    logLik  deviance  df.resid 
    ##  475.8970  479.5537 -234.9485  469.8970        22 
    ## Random effects:
    ##  Groups   Name        Std.Dev.
    ##  region   (Intercept) 1146    
    ##  Residual             2759    
    ## Number of obs: 25, groups:  region, 4
    ## Fixed Effects:
    ## (Intercept)  
    ##        7138  
    ## ----------------------------------------------------
    ## 
    ## Data: mod_data
    ## Models:
    ## mmod_null: seq_sum ~ 1 + (1 | region)
    ## mmod: seq_sum ~ field_type + (1 | region)
    ##           npar    AIC    BIC  logLik deviance  Chisq Df Pr(>Chisq)
    ## mmod_null    3 475.90 479.55 -234.95   469.90                     
    ## mmod         5 478.28 484.38 -234.14   468.28 1.6161  2     0.4457
    ## ----------------------------------------------------
    ## 
    ## 
    ##   Simultaneous Tests for General Linear Hypotheses
    ## 
    ## Multiple Comparisons of Means: Tukey Contrasts
    ## 
    ## 
    ## Fit: lmer(formula = seq_sum ~ field_type + (1 | region), data = mod_data, 
    ##     REML = FALSE)
    ## 
    ## Linear Hypotheses:
    ##                         Estimate Std. Error z value Pr(>|z|)
    ## restored - corn == 0       182.8     1367.0   0.134    0.990
    ## remnant - corn == 0      -1782.3     1766.4  -1.009    0.567
    ## remnant - restored == 0  -1965.1     1496.1  -1.313    0.383
    ## (Adjusted p values reported -- single-step method)
    ## 
    ##     corn restored  remnant 
    ##      "a"      "a"      "a" 
    ## ----------------------------------------------------
    ## 
    ## [1] "Years since restoration and plant_pathogen sequence abundance in Blue Mounds Area"
    ## 
    ## Call:
    ## lm(formula = seq_sum ~ yr_since, data = mod_data2)
    ## 
    ## Residuals:
    ##       1       2       3       4       5       6       7 
    ##  -953.1 -1871.0   723.3 -2080.9  1732.4   111.4  2337.9 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept) 12953.84    1234.13  10.496 0.000135 ***
    ## yr_since     -361.61      83.02  -4.356 0.007319 ** 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 1884 on 5 degrees of freedom
    ## Multiple R-squared:  0.7914, Adjusted R-squared:  0.7497 
    ## F-statistic: 18.97 on 1 and 5 DF,  p-value: 0.007319
    ## 
    ## 
    ## 
    ## 
    ## ---------------------------------
    ## [1] "ectomycorrhizal"
    ## ---------------------------------
    ## 
    ## 
    ## field_type   region   field_name    yr_since  primary_lifestyle    seq_sum
    ## -----------  -------  -----------  ---------  ------------------  --------
    ## remnant      FG       FGREM1              NA  ectomycorrhizal         3537
    ## remnant      FL       FLREM1              NA  ectomycorrhizal         2727
    ## remnant      BM       MBREM1              NA  ectomycorrhizal         1120
    ## restored     BM       MBRP1               18  ectomycorrhizal          543
    ## remnant      LP       LPREM1              NA  ectomycorrhizal          533
    ## restored     FL       FLRP1               40  ectomycorrhizal           14
    ## corn         LP       LPC1                 0  ectomycorrhizal           11
    ## restored     BM       MHRP2                2  ectomycorrhizal            8
    ## restored     FG       FGRP1               15  ectomycorrhizal            6
    ## ----------------------------------------------------
    ## 
    ## Linear mixed model fit by maximum likelihood  ['lmerMod']
    ## Formula: seq_sum ~ field_type + (1 | region)
    ##    Data: mod_data
    ##      AIC      BIC   logLik deviance df.resid 
    ## 156.2751 157.2613 -73.1376 146.2751        4 
    ## Random effects:
    ##  Groups   Name        Std.Dev.
    ##  region   (Intercept)   0.0   
    ##  Residual             818.5   
    ## Number of obs: 9, groups:  region, 4
    ## Fixed Effects:
    ##  (Intercept)  field_type.L  field_type.Q  
    ##          711          1392           696  
    ## optimizer (nloptwrap) convergence code: 0 (OK) ; 0 optimizer warnings; 1 lme4 warnings 
    ## ----------------------------------------------------
    ## 
    ## Linear mixed model fit by maximum likelihood  ['lmerMod']
    ## Formula: seq_sum ~ 1 + (1 | region)
    ##    Data: mod_data
    ##      AIC      BIC   logLik deviance df.resid 
    ## 159.6979 160.2896 -76.8490 153.6979        6 
    ## Random effects:
    ##  Groups   Name        Std.Dev.
    ##  region   (Intercept)    0    
    ##  Residual             1236    
    ## Number of obs: 9, groups:  region, 4
    ## Fixed Effects:
    ## (Intercept)  
    ##       944.3  
    ## optimizer (nloptwrap) convergence code: 0 (OK) ; 0 optimizer warnings; 1 lme4 warnings 
    ## ----------------------------------------------------
    ## 
    ## Data: mod_data
    ## Models:
    ## mmod_null: seq_sum ~ 1 + (1 | region)
    ## mmod: seq_sum ~ field_type + (1 | region)
    ##           npar    AIC    BIC  logLik deviance  Chisq Df Pr(>Chisq)  
    ## mmod_null    3 159.70 160.29 -76.849   153.70                       
    ## mmod         5 156.28 157.26 -73.138   146.28 7.4228  2    0.02444 *
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## ----------------------------------------------------
    ## 
    ## 
    ##   Simultaneous Tests for General Linear Hypotheses
    ## 
    ## Multiple Comparisons of Means: Tukey Contrasts
    ## 
    ## 
    ## Fit: lmer(formula = seq_sum ~ field_type + (1 | region), data = mod_data, 
    ##     REML = FALSE)
    ## 
    ## Linear Hypotheses:
    ##                         Estimate Std. Error z value Pr(>|z|)   
    ## restored - corn == 0       131.7      915.1   0.144  0.98836   
    ## remnant - corn == 0       1968.2      915.1   2.151  0.07677 . 
    ## remnant - restored == 0   1836.5      578.8   3.173  0.00411 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## (Adjusted p values reported -- single-step method)
    ## 
    ##     corn restored  remnant 
    ##     "ab"      "a"      "b" 
    ## ----------------------------------------------------
    ## 
    ## [1] "Years since restoration and ectomycorrhizal sequence abundance in Blue Mounds Area"
    ## 
    ## Call:
    ## lm(formula = seq_sum ~ yr_since, data = mod_data2)
    ## 
    ## Residuals:
    ## ALL 2 residuals are 0: no residual degrees of freedom!
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)
    ## (Intercept)   -58.88        NaN     NaN      NaN
    ## yr_since       33.44        NaN     NaN      NaN
    ## 
    ## Residual standard error: NaN on 0 degrees of freedom
    ## Multiple R-squared:      1,  Adjusted R-squared:    NaN 
    ## F-statistic:   NaN on 1 and 0 DF,  p-value: NA
    ## 
    ## 
    ## 
    ## 
    ## ---------------------------------
    ## [1] "wood_saprotroph"
    ## ---------------------------------
    ## 
    ## 
    ## field_type   region   field_name    yr_since  primary_lifestyle    seq_sum
    ## -----------  -------  -----------  ---------  ------------------  --------
    ## restored     LP       LPRP2                4  wood_saprotroph         4706
    ## restored     FL       FLRSP2              10  wood_saprotroph         4056
    ## corn         LP       LPC1                 0  wood_saprotroph         3819
    ## restored     FG       FGRP1               15  wood_saprotroph         3783
    ## corn         FG       FGC1                 0  wood_saprotroph         3186
    ## corn         BM       PHC1                 0  wood_saprotroph         3065
    ## restored     BM       MHRP2                2  wood_saprotroph         2802
    ## corn         FL       FLC1                 0  wood_saprotroph         2789
    ## restored     BM       ERRP1                3  wood_saprotroph         2717
    ## restored     LP       LPRP1                4  wood_saprotroph         2532
    ## restored     BM       MHRP1                7  wood_saprotroph         2326
    ## restored     BM       PHRP1               11  wood_saprotroph         2233
    ## restored     FL       FLRP4               36  wood_saprotroph         2154
    ## remnant      FL       FLREM1              NA  wood_saprotroph         1992
    ## restored     FL       FLRSP1              10  wood_saprotroph         1909
    ## restored     BM       MBRP1               18  wood_saprotroph         1857
    ## restored     FL       FLRP5               35  wood_saprotroph         1528
    ## restored     FL       FLRP1               40  wood_saprotroph         1490
    ## corn         FL       FLC2                 0  wood_saprotroph         1349
    ## remnant      LP       LPREM1              NA  wood_saprotroph          965
    ## restored     FL       FLRSP3              10  wood_saprotroph          882
    ## remnant      BM       MBREM1              NA  wood_saprotroph          882
    ## restored     BM       KORP1               28  wood_saprotroph          855
    ## restored     BM       BBRP1               16  wood_saprotroph          812
    ## remnant      FG       FGREM1              NA  wood_saprotroph          670
    ## ----------------------------------------------------
    ## 
    ## Linear mixed model fit by maximum likelihood  ['lmerMod']
    ## Formula: seq_sum ~ field_type + (1 | region)
    ##    Data: mod_data
    ##       AIC       BIC    logLik  deviance  df.resid 
    ##  425.1543  431.2487 -207.5772  415.1543        20 
    ## Random effects:
    ##  Groups   Name        Std.Dev.
    ##  region   (Intercept) 206.4   
    ##  Residual             957.3   
    ## Number of obs: 25, groups:  region, 4
    ## Fixed Effects:
    ##  (Intercept)  field_type.L  field_type.Q  
    ##       2110.8       -1225.1        -287.3  
    ## ----------------------------------------------------
    ## 
    ## Linear mixed model fit by maximum likelihood  ['lmerMod']
    ## Formula: seq_sum ~ 1 + (1 | region)
    ##    Data: mod_data
    ##       AIC       BIC    logLik  deviance  df.resid 
    ##  427.4485  431.1052 -210.7243  421.4485        22 
    ## Random effects:
    ##  Groups   Name        Std.Dev.
    ##  region   (Intercept)    0    
    ##  Residual             1108    
    ## Number of obs: 25, groups:  region, 4
    ## Fixed Effects:
    ## (Intercept)  
    ##        2214  
    ## optimizer (nloptwrap) convergence code: 0 (OK) ; 0 optimizer warnings; 1 lme4 warnings 
    ## ----------------------------------------------------
    ## 
    ## Data: mod_data
    ## Models:
    ## mmod_null: seq_sum ~ 1 + (1 | region)
    ## mmod: seq_sum ~ field_type + (1 | region)
    ##           npar    AIC    BIC  logLik deviance  Chisq Df Pr(>Chisq)  
    ## mmod_null    3 427.45 431.11 -210.72   421.45                       
    ## mmod         5 425.15 431.25 -207.58   415.15 6.2942  2    0.04298 *
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## ----------------------------------------------------
    ## 
    ## 
    ##   Simultaneous Tests for General Linear Hypotheses
    ## 
    ## Multiple Comparisons of Means: Tukey Contrasts
    ## 
    ## 
    ## Fit: lmer(formula = seq_sum ~ field_type + (1 | region), data = mod_data, 
    ##     REML = FALSE)
    ## 
    ## Linear Hypotheses:
    ##                         Estimate Std. Error z value Pr(>|z|)  
    ## restored - corn == 0      -514.4      493.3  -1.043   0.5451  
    ## remnant - corn == 0      -1732.5      642.9  -2.695   0.0188 *
    ## remnant - restored == 0  -1218.2      538.5  -2.262   0.0598 .
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## (Adjusted p values reported -- single-step method)
    ## 
    ##     corn restored  remnant 
    ##      "a"     "ab"      "b" 
    ## ----------------------------------------------------
    ## 
    ## [1] "Years since restoration and wood_saprotroph sequence abundance in Blue Mounds Area"
    ## 
    ## Call:
    ## lm(formula = seq_sum ~ yr_since, data = mod_data2)
    ## 
    ## Residuals:
    ##       1       2       3       4       5       6       7 
    ## -829.91   59.82  150.26  371.29  -18.79   66.72  200.60 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)  2891.47     277.63  10.415 0.000141 ***
    ## yr_since      -78.10      18.68  -4.182 0.008639 ** 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 423.7 on 5 degrees of freedom
    ## Multiple R-squared:  0.7777, Adjusted R-squared:  0.7332 
    ## F-statistic: 17.49 on 1 and 5 DF,  p-value: 0.008639
    ## 
    ## 
    ## 
    ## 
    ## ---------------------------------
    ## [1] "litter_saprotroph"
    ## ---------------------------------
    ## 
    ## 
    ## field_type   region   field_name    yr_since  primary_lifestyle    seq_sum
    ## -----------  -------  -----------  ---------  ------------------  --------
    ## corn         FG       FGC1                 0  litter_saprotroph       4932
    ## restored     BM       ERRP1                3  litter_saprotroph       4833
    ## remnant      FL       FLREM1              NA  litter_saprotroph       2526
    ## corn         BM       PHC1                 0  litter_saprotroph       2349
    ## restored     BM       MBRP1               18  litter_saprotroph       2221
    ## restored     BM       BBRP1               16  litter_saprotroph       2076
    ## restored     LP       LPRP2                4  litter_saprotroph       2070
    ## restored     BM       MHRP1                7  litter_saprotroph       2017
    ## restored     LP       LPRP1                4  litter_saprotroph       1975
    ## restored     BM       MHRP2                2  litter_saprotroph       1749
    ## restored     FL       FLRSP2              10  litter_saprotroph       1497
    ## restored     BM       PHRP1               11  litter_saprotroph       1118
    ## corn         FL       FLC1                 0  litter_saprotroph       1116
    ## corn         FL       FLC2                 0  litter_saprotroph       1023
    ## corn         LP       LPC1                 0  litter_saprotroph        896
    ## restored     FL       FLRSP3              10  litter_saprotroph        828
    ## remnant      LP       LPREM1              NA  litter_saprotroph        820
    ## restored     BM       KORP1               28  litter_saprotroph        754
    ## remnant      FG       FGREM1              NA  litter_saprotroph        708
    ## remnant      BM       MBREM1              NA  litter_saprotroph        443
    ## restored     FL       FLRP4               36  litter_saprotroph        425
    ## restored     FL       FLRP5               35  litter_saprotroph        421
    ## restored     FL       FLRSP1              10  litter_saprotroph        396
    ## restored     FG       FGRP1               15  litter_saprotroph        369
    ## restored     FL       FLRP1               40  litter_saprotroph        206
    ## ----------------------------------------------------
    ## 
    ## Linear mixed model fit by maximum likelihood  ['lmerMod']
    ## Formula: seq_sum ~ field_type + (1 | region)
    ##    Data: mod_data
    ##       AIC       BIC    logLik  deviance  df.resid 
    ##  434.3499  440.4442 -212.1749  424.3499        20 
    ## Random effects:
    ##  Groups   Name        Std.Dev.
    ##  region   (Intercept)  263    
    ##  Residual             1148    
    ## Number of obs: 25, groups:  region, 4
    ## Fixed Effects:
    ##  (Intercept)  field_type.L  field_type.Q  
    ##       1555.8        -691.8         141.2  
    ## ----------------------------------------------------
    ## 
    ## Linear mixed model fit by maximum likelihood  ['lmerMod']
    ## Formula: seq_sum ~ 1 + (1 | region)
    ##    Data: mod_data
    ##       AIC       BIC    logLik  deviance  df.resid 
    ##  432.0696  435.7262 -213.0348  426.0696        22 
    ## Random effects:
    ##  Groups   Name        Std.Dev.
    ##  region   (Intercept)  185.5  
    ##  Residual             1201.7  
    ## Number of obs: 25, groups:  region, 4
    ## Fixed Effects:
    ## (Intercept)  
    ##        1517  
    ## ----------------------------------------------------
    ## 
    ## Data: mod_data
    ## Models:
    ## mmod_null: seq_sum ~ 1 + (1 | region)
    ## mmod: seq_sum ~ field_type + (1 | region)
    ##           npar    AIC    BIC  logLik deviance  Chisq Df Pr(>Chisq)
    ## mmod_null    3 432.07 435.73 -213.03   426.07                     
    ## mmod         5 434.35 440.44 -212.18   424.35 1.7197  2     0.4232
    ## ----------------------------------------------------
    ## 
    ## 
    ##   Simultaneous Tests for General Linear Hypotheses
    ## 
    ## Multiple Comparisons of Means: Tukey Contrasts
    ## 
    ## 
    ## Fit: lmer(formula = seq_sum ~ field_type + (1 | region), data = mod_data, 
    ##     REML = FALSE)
    ## 
    ## Linear Hypotheses:
    ##                         Estimate Std. Error z value Pr(>|z|)
    ## restored - corn == 0      -662.2      591.9  -1.119    0.498
    ## remnant - corn == 0       -978.4      771.1  -1.269    0.408
    ## remnant - restored == 0   -316.2      646.2  -0.489    0.875
    ## (Adjusted p values reported -- single-step method)
    ## 
    ##     corn restored  remnant 
    ##      "a"      "a"      "a" 
    ## ----------------------------------------------------
    ## 
    ## [1] "Years since restoration and litter_saprotroph sequence abundance in Blue Mounds Area"
    ## 
    ## Call:
    ## lm(formula = seq_sum ~ yr_since, data = mod_data2)
    ## 
    ## Residuals:
    ##       1       2       3       4       5       6       7 
    ##   270.5  2002.2  -105.1   573.2  -498.3 -1160.7 -1081.8 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)  
    ## (Intercept)  3067.40     785.52   3.905   0.0114 *
    ## yr_since      -78.87      52.84  -1.493   0.1958  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 1199 on 5 degrees of freedom
    ## Multiple R-squared:  0.3082, Adjusted R-squared:  0.1699 
    ## F-statistic: 2.228 on 1 and 5 DF,  p-value: 0.1958

Model tests on `field_type` are technically invalid due to
pseudoreplication, but are included here to point out trends that we may
be able to present in some other valid way. Trends with restoration age
in Blue Mounds are clearly justified. Results are shown in descending
order based on sequence abundance in remnants:

- Soil saprotroph increases with years since
- Plant pathogens decrease with years since
- Ectomycorrhizal abundance is very low in corn/restored and with little
  replication; nothing can be said except that it’s relatively abundant
  in remnants.
- Wood saprotroph differs among field types (corn vs. remnant; restored
  intermediate) and decreases with years since
- Litter saprotroph is abundant everywhere, but differences over time or
  field type are weak.

#### ITS-based indicators

An indicator species analysis is warranted, identifying which species
correlate strongly with `field_type`. Performing this with all ITS data
may identify particular species to further examine, although it remains
a problem that we cannot distinguish field type from an individual field
due to pseudoreplication.

Following the indicator species analysis, richness and composition of
selected guilds is calculated. These calculations are done with data
re-rarefied into the guilds identified here, again to showcase
particular species which seem to drive differences among field types.
It’s also of value because this approach avoids the problem we have with
pseudoreplication.

With indicator species analysis performed using package
[indicspecies](http://sites.google.com/site/miqueldecaceres/), the index
values A and B show the specificity and fidelity components of the
IndVal combined index. The combined index value is noted as ‘stat’ in
the output table below.

``` r
its_inspan <- 
    spe$its_rfy %>% 
    left_join(sites, by = join_by(field_key)) %>% 
    inspan(., 1999, meta$its_rfy)
```

``` r
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
```

| field_type | n_otu |  stat_avg |   stat_sd |
|:-----------|------:|----------:|----------:|
| corn       |    88 | 0.8299548 | 0.0965605 |
| restored   |    10 | 0.8150214 | 0.0363025 |
| remnant    |    47 | 0.7462101 | 0.0822420 |

Indicator species stats of entire rarefied ITS table

Potential indicators were filtered to p.value\<0.05 before this summary
was produced. Cornfields are a restrictive habitat for soil microbes,
and that is reflected in the results here. More species have higher
specificity and fidelity to cornfields than the other field types. The
top ten indicators for each field type are printed here; the entire
table is available for further use.

``` r
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
```

| otu_num  |         A |      B |      stat | p.value | field_type | primary_lifestyle | phylum            | class              | order                              | family                             | genus            | species                 |
|:---------|----------:|-------:|----------:|--------:|:-----------|:------------------|:------------------|:-------------------|:-----------------------------------|:-----------------------------------|:-----------------|:------------------------|
| otu_537  | 1.0000000 | 1.0000 | 1.0000000 |  0.0005 | corn       | soil_saprotroph   | Basidiomycota     | Agaricomycetes     | Agaricales                         | Bolbitiaceae                       | Conocybe         | Conocybe_apala          |
| otu_204  | 0.9937578 | 1.0000 | 0.9968740 |  0.0005 | corn       | NA                | Mortierellomycota | Mortierellomycetes | Mortierellales                     | Mortierellaceae                    | NA               | NA                      |
| otu_172  | 0.9772048 | 1.0000 | 0.9885367 |  0.0005 | corn       | plant_pathogen    | Ascomycota        | Dothideomycetes    | Pleosporales                       | Corynesporascaceae                 | Corynespora      | Corynespora_cassiicola  |
| otu_188  | 0.9759492 | 1.0000 | 0.9879014 |  0.0005 | corn       | NA                | NA                | NA                 | NA                                 | NA                                 | NA               | NA                      |
| otu_9    | 0.9753667 | 1.0000 | 0.9876066 |  0.0025 | corn       | soil_saprotroph   | Basidiomycota     | Tremellomycetes    | Cystofilobasidiales                | Mrakiaceae                         | Tausonia         | Tausonia_pullulans      |
| otu_200  | 0.9724757 | 1.0000 | 0.9861418 |  0.0005 | corn       | plant_pathogen    | Ascomycota        | Dothideomycetes    | Pleosporales                       | Phaeosphaeriaceae                  | Ophiosphaerella  | unidentified            |
| otu_59   | 0.9602305 | 1.0000 | 0.9799135 |  0.0005 | corn       | soil_saprotroph   | Mortierellomycota | Mortierellomycetes | Mortierellales                     | Mortierellaceae                    | Mortierella      | NA                      |
| otu_694  | 0.9400850 | 1.0000 | 0.9695798 |  0.0010 | corn       | NA                | NA                | NA                 | NA                                 | NA                                 | NA               | NA                      |
| otu_553  | 0.9378783 | 1.0000 | 0.9684412 |  0.0015 | corn       | plant_pathogen    | Ascomycota        | Sordariomycetes    | Magnaporthales                     | Magnaporthaceae                    | Gaeumannomyces   | NA                      |
| otu_364  | 0.9318632 | 1.0000 | 0.9653306 |  0.0005 | corn       | NA                | Ascomycota        | Sordariomycetes    | Sordariales                        | Lasiosphaeriaceae                  | Cladorrhinum     | NA                      |
| otu_332  | 0.9219288 | 0.8125 | 0.8654867 |  0.0370 | restored   | plant_pathogen    | Ascomycota        | Sordariomycetes    | Glomerellales                      | Plectosphaerellaceae               | Plectosphaerella | NA                      |
| otu_177  | 0.9809886 | 0.7500 | 0.8577537 |  0.0280 | restored   | NA                | Ascomycota        | Dothideomycetes    | Pleosporales                       | NA                                 | NA               | NA                      |
| otu_817  | 1.0000000 | 0.6875 | 0.8291562 |  0.0190 | restored   | NA                | Ascomycota        | NA                 | NA                                 | NA                                 | NA               | NA                      |
| otu_461  | 0.8351648 | 0.8125 | 0.8237545 |  0.0180 | restored   | NA                | Ascomycota        | Dothideomycetes    | Pleosporales                       | Phaeosphaeriaceae                  | NA               | NA                      |
| otu_35   | 0.7234228 | 0.9375 | 0.8235344 |  0.0440 | restored   | animal_parasite   | Ascomycota        | Sordariomycetes    | Hypocreales                        | Clavicipitaceae                    | Metarhizium      | NA                      |
| otu_193  | 0.8307978 | 0.8125 | 0.8215979 |  0.0480 | restored   | NA                | Basidiomycota     | Agaricomycetes     | Sebacinales                        | unidentified                       | unidentified     | unidentified            |
| otu_107  | 0.8061297 | 0.8125 | 0.8093086 |  0.0310 | restored   | NA                | Ascomycota        | Dothideomycetes    | Pleosporales                       | NA                                 | NA               | NA                      |
| otu_114  | 0.6963432 | 0.9375 | 0.8079739 |  0.0025 | restored   | soil_saprotroph   | Mortierellomycota | Mortierellomycetes | Mortierellales                     | Mortierellaceae                    | Mortierella      | unidentified            |
| otu_238  | 0.9179982 | 0.6250 | 0.7574621 |  0.0455 | restored   | NA                | Ascomycota        | Leotiomycetes      | NA                                 | NA                                 | NA               | NA                      |
| otu_10   | 0.5687968 | 1.0000 | 0.7541862 |  0.0150 | restored   | NA                | Ascomycota        | NA                 | NA                                 | NA                                 | NA               | NA                      |
| otu_772  | 0.9272420 | 1.0000 | 0.9629340 |  0.0010 | remnant    | NA                | Ascomycota        | Sordariomycetes    | NA                                 | NA                                 | NA               | NA                      |
| otu_629  | 0.9159892 | 1.0000 | 0.9570732 |  0.0015 | remnant    | NA                | Ascomycota        | Leotiomycetes      | Helotiales                         | Hyaloscyphaceae                    | Microscypha      | unidentified            |
| otu_159  | 0.8185686 | 1.0000 | 0.9047478 |  0.0020 | remnant    | NA                | Ascomycota        | Sordariomycetes    | Sordariomycetes_ord_Incertae_sedis | Sordariomycetes_fam_Incertae_sedis | Pleurophragmium  | unidentified            |
| otu_135  | 0.7768230 | 1.0000 | 0.8813757 |  0.0070 | remnant    | plant_pathogen    | Ascomycota        | Sordariomycetes    | Hypocreales                        | Nectriaceae                        | Ilyonectria      | NA                      |
| otu_854  | 1.0000000 | 0.7500 | 0.8660254 |  0.0030 | remnant    | NA                | Ascomycota        | NA                 | NA                                 | NA                                 | NA               | NA                      |
| otu_1740 | 1.0000000 | 0.7500 | 0.8660254 |  0.0040 | remnant    | NA                | Glomeromycota     | Glomeromycetes     | Glomerales                         | Glomeraceae                        | NA               | NA                      |
| otu_1098 | 0.9716841 | 0.7500 | 0.8536762 |  0.0075 | remnant    | NA                | NA                | NA                 | NA                                 | NA                                 | NA               | NA                      |
| otu_235  | 0.7275292 | 1.0000 | 0.8529532 |  0.0480 | remnant    | NA                | Ascomycota        | Leotiomycetes      | Helotiales                         | Hyaloscyphaceae                    | NA               | NA                      |
| otu_1468 | 0.9332261 | 0.7500 | 0.8366119 |  0.0055 | remnant    | NA                | Ascomycota        | Sordariomycetes    | NA                                 | NA                                 | NA               | NA                      |
| otu_140  | 0.9276552 | 0.7500 | 0.8341111 |  0.0310 | remnant    | soil_saprotroph   | Ascomycota        | Sordariomycetes    | Hypocreales                        | Stachybotryaceae                   | Striaticonidium  | Striaticonidium_cinctum |

Indicator species of ITS OTUs (top 10 per field type)

### Soil saprotrophs

#### Trends over time

``` r
guiltime("soil_saprotroph")
```

    ## $bm_summary
    ## 
    ## Call:
    ## lm(formula = seq_sum ~ yr_since, data = d %>% filter(region == 
    ##     "BM"))
    ## 
    ## Residuals:
    ##       1       2       3       4       5       6       7 
    ##  -597.3   342.4  1424.6 -2538.3  -851.6   458.9  1761.3 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)  
    ## (Intercept)  3185.04    1055.79   3.017   0.0295 *
    ## yr_since      214.51      71.02   3.020   0.0294 *
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 1611 on 5 degrees of freedom
    ## Multiple R-squared:  0.646,  Adjusted R-squared:  0.5752 
    ## F-statistic: 9.123 on 1 and 5 DF,  p-value: 0.0294
    ## 
    ## 
    ## $fl_summary
    ## 
    ## Call:
    ## lm(formula = seq_sum ~ yr_since, data = d %>% filter(region == 
    ##     "FL"))
    ## 
    ## Residuals:
    ##       1       2       3       4       5       6 
    ##  -850.7 -1568.0  2651.5  2358.1 -2413.9  -176.9 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)  
    ## (Intercept)  7686.71    1896.81   4.052   0.0154 *
    ## yr_since      -43.58      69.88  -0.624   0.5667  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 2325 on 4 degrees of freedom
    ## Multiple R-squared:  0.08861,    Adjusted R-squared:  -0.1392 
    ## F-statistic: 0.3889 on 1 and 4 DF,  p-value: 0.5667
    ## 
    ## 
    ## $plot

<img src="microbial_guild_taxonomy_files/figure-gfm/ssap_guiltime-1.png" style="display: block; margin: auto;" />

Sequence abundance of soil saprotrophs increases over time in the Blue
Mounds area ($R^2_{Adj}=0.58, p<0.05$), but this appears to be leveraged
by Karla Ott’s property, though. With all that big bluestem…maybe there
is more litter and soil carbon? It will be good to look at trends in
soil chemistry.

#### Diversity

``` r
(ssap <- rerare(spe$its_raw, meta$its_raw, primary_lifestyle, "soil_saprotroph", sites))
```

    ## $seq_depth
    ## [1] 3586
    ## 
    ## $zero_otu_num
    ## [1] 6
    ## 
    ## $rrfd
    ## # A tibble: 25 × 245
    ##    field_key otu_2 otu_9 otu_14 otu_27 otu_37 otu_41 otu_47 otu_49 otu_55 otu_59
    ##        <dbl> <dbl> <dbl>  <dbl>  <dbl>  <dbl>  <dbl>  <dbl>  <dbl>  <dbl>  <dbl>
    ##  1         1   716     0     97      2    181     18    841      0    613      0
    ##  2         2  1464     0    452      2    142     66      0      0      0      0
    ##  3         3     7    29    413      0    246    471      0      0      0    280
    ##  4         4     0     0      0      0    281    220     20      0     20      0
    ##  5         5   118    47      0      0    122     35      0      0      0    220
    ##  6         6   143  1289    360    143     22    261      0      0      0     74
    ##  7         7   326  2261     44     46    100     50      0      2      0    207
    ##  8         8   601    38    373     71    219     62     19      0      0     25
    ##  9         9   708     0    295      2     96     48    978      0    743      9
    ## 10        10   570     0    307      0     45     68    383      0   1225      0
    ## # … with 15 more rows, and 234 more variables: otu_61 <dbl>, otu_66 <dbl>,
    ## #   otu_70 <dbl>, otu_75 <dbl>, otu_79 <dbl>, otu_88 <dbl>, otu_89 <dbl>,
    ## #   otu_100 <dbl>, otu_102 <dbl>, otu_106 <dbl>, otu_114 <dbl>, otu_132 <dbl>,
    ## #   otu_134 <dbl>, otu_140 <dbl>, otu_144 <dbl>, otu_154 <dbl>, otu_168 <dbl>,
    ## #   otu_170 <dbl>, otu_180 <dbl>, otu_182 <dbl>, otu_192 <dbl>, otu_195 <dbl>,
    ## #   otu_209 <dbl>, otu_216 <dbl>, otu_221 <dbl>, otu_223 <dbl>, otu_224 <dbl>,
    ## #   otu_234 <dbl>, otu_276 <dbl>, otu_280 <dbl>, otu_283 <dbl>, …
    ## 
    ## $rrfd_speTaxa
    ## # A tibble: 935 × 14
    ##    field_key otu_num seq_abund phylum   class order family genus species prima…¹
    ##        <dbl> <chr>       <dbl> <chr>    <chr> <chr> <chr>  <chr> <chr>   <chr>  
    ##  1         1 otu_2         716 Mortier… Mort… Mort… Morti… Mort… Mortie… soil_s…
    ##  2         1 otu_14         97 Mortier… Mort… Mort… Morti… Mort… <NA>    soil_s…
    ##  3         1 otu_27          2 Basidio… Trem… Filo… Pisku… Soli… <NA>    soil_s…
    ##  4         1 otu_37        181 Mortier… Mort… Mort… Morti… Mort… <NA>    soil_s…
    ##  5         1 otu_41         18 Mortier… Mort… Mort… Morti… Mort… Mortie… soil_s…
    ##  6         1 otu_47        841 Ascomyc… Geog… Geog… Geogl… Geog… uniden… soil_s…
    ##  7         1 otu_55        613 Basidio… Agar… Agar… Hygro… Hygr… <NA>    soil_s…
    ##  8         1 otu_75          4 Basidio… Trem… Filo… Pisku… Soli… Solico… soil_s…
    ##  9         1 otu_79         55 Mortier… Mort… Mort… Morti… Mort… <NA>    soil_s…
    ## 10         1 otu_114        65 Mortier… Mort… Mort… Morti… Mort… uniden… soil_s…
    ## # … with 925 more rows, 4 more variables: field_name <chr>, region <chr>,
    ## #   field_type <ord>, yr_since <dbl>, and abbreviated variable name
    ## #   ¹​primary_lifestyle

``` r
(ssap_div <- calc_diversity(ssap$rrfd))
```

<div data-pagedtable="false">

<script data-pagedtable-source type="application/json">
{"columns":[{"label":["field_key"],"name":[1],"type":["dbl"],"align":["right"]},{"label":["field_name"],"name":[2],"type":["chr"],"align":["left"]},{"label":["region"],"name":[3],"type":["chr"],"align":["left"]},{"label":["field_type"],"name":[4],"type":["ord"],"align":["right"]},{"label":["yr_since"],"name":[5],"type":["dbl"],"align":["right"]},{"label":["hill_index"],"name":[6],"type":["ord"],"align":["right"]},{"label":["value"],"name":[7],"type":["dbl"],"align":["right"]}],"data":[{"1":"1","2":"BBRP1","3":"BM","4":"restored","5":"16","6":"N0","7":"31.00000000"},{"1":"1","2":"BBRP1","3":"BM","4":"restored","5":"16","6":"N1","7":"10.86224023"},{"1":"1","2":"BBRP1","3":"BM","4":"restored","5":"16","6":"N2","7":"7.23303504"},{"1":"1","2":"BBRP1","3":"BM","4":"restored","5":"16","6":"E10","7":"0.35039485"},{"1":"1","2":"BBRP1","3":"BM","4":"restored","5":"16","6":"E20","7":"0.23332371"},{"1":"2","2":"ERRP1","3":"BM","4":"restored","5":"3","6":"N0","7":"36.00000000"},{"1":"2","2":"ERRP1","3":"BM","4":"restored","5":"3","6":"N1","7":"7.39856050"},{"1":"2","2":"ERRP1","3":"BM","4":"restored","5":"3","6":"N2","7":"4.18034042"},{"1":"2","2":"ERRP1","3":"BM","4":"restored","5":"3","6":"E10","7":"0.20551557"},{"1":"2","2":"ERRP1","3":"BM","4":"restored","5":"3","6":"E20","7":"0.11612057"},{"1":"3","2":"FGC1","3":"FG","4":"corn","5":"0","6":"N0","7":"34.00000000"},{"1":"3","2":"FGC1","3":"FG","4":"corn","5":"0","6":"N1","7":"16.74643781"},{"1":"3","2":"FGC1","3":"FG","4":"corn","5":"0","6":"N2","7":"12.61583470"},{"1":"3","2":"FGC1","3":"FG","4":"corn","5":"0","6":"E10","7":"0.49254229"},{"1":"3","2":"FGC1","3":"FG","4":"corn","5":"0","6":"E20","7":"0.37105396"},{"1":"4","2":"FGREM1","3":"FG","4":"remnant","5":"NA","6":"N0","7":"32.00000000"},{"1":"4","2":"FGREM1","3":"FG","4":"remnant","5":"NA","6":"N1","7":"6.41412570"},{"1":"4","2":"FGREM1","3":"FG","4":"remnant","5":"NA","6":"N2","7":"3.14150276"},{"1":"4","2":"FGREM1","3":"FG","4":"remnant","5":"NA","6":"E10","7":"0.20044143"},{"1":"4","2":"FGREM1","3":"FG","4":"remnant","5":"NA","6":"E20","7":"0.09817196"},{"1":"5","2":"FGRP1","3":"FG","4":"restored","5":"15","6":"N0","7":"29.00000000"},{"1":"5","2":"FGRP1","3":"FG","4":"restored","5":"15","6":"N1","7":"10.96471698"},{"1":"5","2":"FGRP1","3":"FG","4":"restored","5":"15","6":"N2","7":"6.34555593"},{"1":"5","2":"FGRP1","3":"FG","4":"restored","5":"15","6":"E10","7":"0.37809369"},{"1":"5","2":"FGRP1","3":"FG","4":"restored","5":"15","6":"E20","7":"0.21881227"},{"1":"6","2":"FLC1","3":"FL","4":"corn","5":"0","6":"N0","7":"30.00000000"},{"1":"6","2":"FLC1","3":"FL","4":"corn","5":"0","6":"N1","7":"9.42357851"},{"1":"6","2":"FLC1","3":"FL","4":"corn","5":"0","6":"N2","7":"5.55324586"},{"1":"6","2":"FLC1","3":"FL","4":"corn","5":"0","6":"E10","7":"0.31411928"},{"1":"6","2":"FLC1","3":"FL","4":"corn","5":"0","6":"E20","7":"0.18510820"},{"1":"7","2":"FLC2","3":"FL","4":"corn","5":"0","6":"N0","7":"29.00000000"},{"1":"7","2":"FLC2","3":"FL","4":"corn","5":"0","6":"N1","7":"4.55609794"},{"1":"7","2":"FLC2","3":"FL","4":"corn","5":"0","6":"N2","7":"2.39787143"},{"1":"7","2":"FLC2","3":"FL","4":"corn","5":"0","6":"E10","7":"0.15710683"},{"1":"7","2":"FLC2","3":"FL","4":"corn","5":"0","6":"E20","7":"0.08268522"},{"1":"8","2":"FLREM1","3":"FL","4":"remnant","5":"NA","6":"N0","7":"44.00000000"},{"1":"8","2":"FLREM1","3":"FL","4":"remnant","5":"NA","6":"N1","7":"17.87530791"},{"1":"8","2":"FLREM1","3":"FL","4":"remnant","5":"NA","6":"N2","7":"11.61712984"},{"1":"8","2":"FLREM1","3":"FL","4":"remnant","5":"NA","6":"E10","7":"0.40625700"},{"1":"8","2":"FLREM1","3":"FL","4":"remnant","5":"NA","6":"E20","7":"0.26402568"},{"1":"9","2":"FLRP1","3":"FL","4":"restored","5":"40","6":"N0","7":"37.00000000"},{"1":"9","2":"FLRP1","3":"FL","4":"restored","5":"40","6":"N1","7":"9.40699910"},{"1":"9","2":"FLRP1","3":"FL","4":"restored","5":"40","6":"N2","7":"5.97985170"},{"1":"9","2":"FLRP1","3":"FL","4":"restored","5":"40","6":"E10","7":"0.25424322"},{"1":"9","2":"FLRP1","3":"FL","4":"restored","5":"40","6":"E20","7":"0.16161761"},{"1":"10","2":"FLRP4","3":"FL","4":"restored","5":"36","6":"N0","7":"38.00000000"},{"1":"10","2":"FLRP4","3":"FL","4":"restored","5":"36","6":"N1","7":"10.62099112"},{"1":"10","2":"FLRP4","3":"FL","4":"restored","5":"36","6":"N2","7":"5.94897886"},{"1":"10","2":"FLRP4","3":"FL","4":"restored","5":"36","6":"E10","7":"0.27949977"},{"1":"10","2":"FLRP4","3":"FL","4":"restored","5":"36","6":"E20","7":"0.15655208"},{"1":"11","2":"FLRP5","3":"FL","4":"restored","5":"35","6":"N0","7":"36.00000000"},{"1":"11","2":"FLRP5","3":"FL","4":"restored","5":"35","6":"N1","7":"10.43799584"},{"1":"11","2":"FLRP5","3":"FL","4":"restored","5":"35","6":"N2","7":"7.25301697"},{"1":"11","2":"FLRP5","3":"FL","4":"restored","5":"35","6":"E10","7":"0.28994433"},{"1":"11","2":"FLRP5","3":"FL","4":"restored","5":"35","6":"E20","7":"0.20147269"},{"1":"12","2":"FLRSP1","3":"FL","4":"restored","5":"10","6":"N0","7":"36.00000000"},{"1":"12","2":"FLRSP1","3":"FL","4":"restored","5":"10","6":"N1","7":"8.65822702"},{"1":"12","2":"FLRSP1","3":"FL","4":"restored","5":"10","6":"N2","7":"5.70549395"},{"1":"12","2":"FLRSP1","3":"FL","4":"restored","5":"10","6":"E10","7":"0.24050631"},{"1":"12","2":"FLRSP1","3":"FL","4":"restored","5":"10","6":"E20","7":"0.15848594"},{"1":"13","2":"FLRSP2","3":"FL","4":"restored","5":"10","6":"N0","7":"43.00000000"},{"1":"13","2":"FLRSP2","3":"FL","4":"restored","5":"10","6":"N1","7":"11.64522116"},{"1":"13","2":"FLRSP2","3":"FL","4":"restored","5":"10","6":"N2","7":"5.87998416"},{"1":"13","2":"FLRSP2","3":"FL","4":"restored","5":"10","6":"E10","7":"0.27081910"},{"1":"13","2":"FLRSP2","3":"FL","4":"restored","5":"10","6":"E20","7":"0.13674382"},{"1":"14","2":"FLRSP3","3":"FL","4":"restored","5":"10","6":"N0","7":"32.00000000"},{"1":"14","2":"FLRSP3","3":"FL","4":"restored","5":"10","6":"N1","7":"4.17268922"},{"1":"14","2":"FLRSP3","3":"FL","4":"restored","5":"10","6":"N2","7":"2.44287674"},{"1":"14","2":"FLRSP3","3":"FL","4":"restored","5":"10","6":"E10","7":"0.13039654"},{"1":"14","2":"FLRSP3","3":"FL","4":"restored","5":"10","6":"E20","7":"0.07633990"},{"1":"15","2":"KORP1","3":"BM","4":"restored","5":"28","6":"N0","7":"35.00000000"},{"1":"15","2":"KORP1","3":"BM","4":"restored","5":"28","6":"N1","7":"13.19041481"},{"1":"15","2":"KORP1","3":"BM","4":"restored","5":"28","6":"N2","7":"9.63317023"},{"1":"15","2":"KORP1","3":"BM","4":"restored","5":"28","6":"E10","7":"0.37686899"},{"1":"15","2":"KORP1","3":"BM","4":"restored","5":"28","6":"E20","7":"0.27523344"},{"1":"16","2":"LPC1","3":"LP","4":"corn","5":"0","6":"N0","7":"40.00000000"},{"1":"16","2":"LPC1","3":"LP","4":"corn","5":"0","6":"N1","7":"11.59630336"},{"1":"16","2":"LPC1","3":"LP","4":"corn","5":"0","6":"N2","7":"6.69893489"},{"1":"16","2":"LPC1","3":"LP","4":"corn","5":"0","6":"E10","7":"0.28990758"},{"1":"16","2":"LPC1","3":"LP","4":"corn","5":"0","6":"E20","7":"0.16747337"},{"1":"17","2":"LPREM1","3":"LP","4":"remnant","5":"NA","6":"N0","7":"48.00000000"},{"1":"17","2":"LPREM1","3":"LP","4":"remnant","5":"NA","6":"N1","7":"11.79998712"},{"1":"17","2":"LPREM1","3":"LP","4":"remnant","5":"NA","6":"N2","7":"6.63429296"},{"1":"17","2":"LPREM1","3":"LP","4":"remnant","5":"NA","6":"E10","7":"0.24583307"},{"1":"17","2":"LPREM1","3":"LP","4":"remnant","5":"NA","6":"E20","7":"0.13821444"},{"1":"18","2":"LPRP1","3":"LP","4":"restored","5":"4","6":"N0","7":"53.00000000"},{"1":"18","2":"LPRP1","3":"LP","4":"restored","5":"4","6":"N1","7":"15.95453677"},{"1":"18","2":"LPRP1","3":"LP","4":"restored","5":"4","6":"N2","7":"9.79506781"},{"1":"18","2":"LPRP1","3":"LP","4":"restored","5":"4","6":"E10","7":"0.30102900"},{"1":"18","2":"LPRP1","3":"LP","4":"restored","5":"4","6":"E20","7":"0.18481260"},{"1":"19","2":"LPRP2","3":"LP","4":"restored","5":"4","6":"N0","7":"35.00000000"},{"1":"19","2":"LPRP2","3":"LP","4":"restored","5":"4","6":"N1","7":"13.42663395"},{"1":"19","2":"LPRP2","3":"LP","4":"restored","5":"4","6":"N2","7":"9.24891540"},{"1":"19","2":"LPRP2","3":"LP","4":"restored","5":"4","6":"E10","7":"0.38361811"},{"1":"19","2":"LPRP2","3":"LP","4":"restored","5":"4","6":"E20","7":"0.26425473"},{"1":"20","2":"MBREM1","3":"BM","4":"remnant","5":"NA","6":"N0","7":"36.00000000"},{"1":"20","2":"MBREM1","3":"BM","4":"remnant","5":"NA","6":"N1","7":"6.60123947"},{"1":"20","2":"MBREM1","3":"BM","4":"remnant","5":"NA","6":"N2","7":"4.22371778"},{"1":"20","2":"MBREM1","3":"BM","4":"remnant","5":"NA","6":"E10","7":"0.18336776"},{"1":"20","2":"MBREM1","3":"BM","4":"remnant","5":"NA","6":"E20","7":"0.11732549"},{"1":"21","2":"MBRP1","3":"BM","4":"restored","5":"18","6":"N0","7":"47.00000000"},{"1":"21","2":"MBRP1","3":"BM","4":"restored","5":"18","6":"N1","7":"14.87007531"},{"1":"21","2":"MBRP1","3":"BM","4":"restored","5":"18","6":"N2","7":"7.99806196"},{"1":"21","2":"MBRP1","3":"BM","4":"restored","5":"18","6":"E10","7":"0.31638458"},{"1":"21","2":"MBRP1","3":"BM","4":"restored","5":"18","6":"E20","7":"0.17017153"},{"1":"22","2":"MHRP1","3":"BM","4":"restored","5":"7","6":"N0","7":"38.00000000"},{"1":"22","2":"MHRP1","3":"BM","4":"restored","5":"7","6":"N1","7":"13.81476801"},{"1":"22","2":"MHRP1","3":"BM","4":"restored","5":"7","6":"N2","7":"9.40389658"},{"1":"22","2":"MHRP1","3":"BM","4":"restored","5":"7","6":"E10","7":"0.36354653"},{"1":"22","2":"MHRP1","3":"BM","4":"restored","5":"7","6":"E20","7":"0.24747096"},{"1":"23","2":"MHRP2","3":"BM","4":"restored","5":"2","6":"N0","7":"37.00000000"},{"1":"23","2":"MHRP2","3":"BM","4":"restored","5":"2","6":"N1","7":"10.14585020"},{"1":"23","2":"MHRP2","3":"BM","4":"restored","5":"2","6":"N2","7":"5.67423090"},{"1":"23","2":"MHRP2","3":"BM","4":"restored","5":"2","6":"E10","7":"0.27421217"},{"1":"23","2":"MHRP2","3":"BM","4":"restored","5":"2","6":"E20","7":"0.15335759"},{"1":"24","2":"PHC1","3":"BM","4":"corn","5":"0","6":"N0","7":"38.00000000"},{"1":"24","2":"PHC1","3":"BM","4":"corn","5":"0","6":"N1","7":"8.28354349"},{"1":"24","2":"PHC1","3":"BM","4":"corn","5":"0","6":"N2","7":"4.91962802"},{"1":"24","2":"PHC1","3":"BM","4":"corn","5":"0","6":"E10","7":"0.21798799"},{"1":"24","2":"PHC1","3":"BM","4":"corn","5":"0","6":"E20","7":"0.12946390"},{"1":"25","2":"PHRP1","3":"BM","4":"restored","5":"11","6":"N0","7":"41.00000000"},{"1":"25","2":"PHRP1","3":"BM","4":"restored","5":"11","6":"N1","7":"9.60375740"},{"1":"25","2":"PHRP1","3":"BM","4":"restored","5":"11","6":"N2","7":"4.99606668"},{"1":"25","2":"PHRP1","3":"BM","4":"restored","5":"11","6":"E10","7":"0.23423799"},{"1":"25","2":"PHRP1","3":"BM","4":"restored","5":"11","6":"E20","7":"0.12185528"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>

</div>

Diversity measures are stored in this data frame for further use…

``` r
(ssap_comp <- gudicom(ssap_div, ssap$rrfd_speTaxa, "soil_saprotroph"))
```

    ## $Hills_field_type

<img src="microbial_guild_taxonomy_files/figure-gfm/ssap_composition-1.png" style="display: block; margin: auto;" />

    ## 
    ## $Hills_yrs_since_restoration

<img src="microbial_guild_taxonomy_files/figure-gfm/ssap_composition-2.png" style="display: block; margin: auto;" />

    ## 
    ## $Composition

<img src="microbial_guild_taxonomy_files/figure-gfm/ssap_composition-3.png" style="display: block; margin: auto;" />

<div data-pagedtable="false">

<script data-pagedtable-source type="application/json">
{"columns":[{"label":["field_type"],"name":[1],"type":["ord"],"align":["right"]},{"label":["order"],"name":[2],"type":["chr"],"align":["left"]},{"label":["seq_comp"],"name":[3],"type":["dbl"],"align":["right"]}],"data":[{"1":"corn","2":"Agaricales","3":"11.085075"},{"1":"corn","2":"Cystofilobasidiales","3":"20.510270"},{"1":"corn","2":"Filobasidiales","3":"10.197536"},{"1":"corn","2":"Mortierellales","3":"46.455734"},{"1":"corn","2":"Other (OTU<2%)","3":"4.338644"},{"1":"corn","2":"Pezizales","3":"2.455787"},{"1":"corn","2":"Phallales","3":"4.956954"},{"1":"restored","2":"Agaricales","3":"21.061653"},{"1":"restored","2":"Filobasidiales","3":"6.077823"},{"1":"restored","2":"Geoglossales","3":"16.333898"},{"1":"restored","2":"Hypocreales","3":"4.470648"},{"1":"restored","2":"Mortierellales","3":"38.646304"},{"1":"restored","2":"Other (OTU<2%)","3":"10.106596"},{"1":"restored","2":"Pezizales","3":"3.303077"},{"1":"remnant","2":"Agaricales","3":"40.701440"},{"1":"remnant","2":"Geoglossales","3":"9.057892"},{"1":"remnant","2":"Helotiales","3":"4.227016"},{"1":"remnant","2":"Hypocreales","3":"10.692250"},{"1":"remnant","2":"Mortierellales","3":"28.473286"},{"1":"remnant","2":"Other (OTU<2%)","3":"6.848116"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>

</div>

Richness increases from corn to remnant, but within-group variability is
high. Diversity indices look muddy. Diversity indices increase with
years since restoration, but the significance of this remains to be
seen.

Composition of soil saprotrophs by order can be modified somewhat by
choosing the threshold for lumping rare orders into an “other” category.
Leaving this at the default of \<2%, nine named orders are left. Agarics
increase strongly from corn to remnant; Cystofilobasidiales and
Filobasidiales aren’t found outside of cornfields. Generally, cornfield
composition looks different than the other two, but remnants do appear
somewhat intermediate.

#### Indicators

``` r
ssap_inspan <- 
    ssap$rrfd %>% 
    left_join(sites, by = join_by(field_key)) %>% 
    inspan(., 1999, meta$its_raw)
```

``` r
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
```

| field_type | n_otu |  stat_avg |   stat_sd |
|:-----------|------:|----------:|----------:|
| corn       |     5 | 0.9342112 | 0.0705360 |
| restored   |     2 | 0.7794787 | 0.0024830 |
| remnant    |     2 | 0.7515430 | 0.0628423 |

Indicator species stats: soil saprotrophs

We see the same trend as before, where more indicators are found in
cornfields, and their indicator stats are stronger.

``` r
ssap_inspan %>% 
    mutate(field_type = factor(
        field_type,
        ordered = TRUE,
        levels = c("corn", "restored", "remnant")
    )) %>%
    arrange(field_type, -stat) %>% 
    kable(format = "pandoc", caption = "Indicator species of soil saprotrophs")
```

| otu_num  |         A |      B |      stat | p.value | field_type | primary_lifestyle | phylum            | class              | order               | family           | genus              | species                 |
|:---------|----------:|-------:|----------:|--------:|:-----------|:------------------|:------------------|:-------------------|:--------------------|:-----------------|:-------------------|:------------------------|
| otu_537  | 1.0000000 | 1.0000 | 1.0000000 |  0.0005 | corn       | soil_saprotroph   | Basidiomycota     | Agaricomycetes     | Agaricales          | Bolbitiaceae     | Conocybe           | Conocybe_apala          |
| otu_9    | 0.9414168 | 1.0000 | 0.9702664 |  0.0015 | corn       | soil_saprotroph   | Basidiomycota     | Tremellomycetes    | Cystofilobasidiales | Mrakiaceae       | Tausonia           | Tausonia_pullulans      |
| otu_59   | 0.9303537 | 1.0000 | 0.9645485 |  0.0005 | corn       | soil_saprotroph   | Mortierellomycota | Mortierellomycetes | Mortierellales      | Mortierellaceae  | Mortierella        | NA                      |
| otu_134  | 0.8393619 | 1.0000 | 0.9161670 |  0.0020 | corn       | soil_saprotroph   | Mortierellomycota | Mortierellomycetes | Mortierellales      | Mortierellaceae  | Mortierella        | NA                      |
| otu_41   | 0.6725220 | 1.0000 | 0.8200744 |  0.0165 | corn       | soil_saprotroph   | Mortierellomycota | Mortierellomycetes | Mortierellales      | Mortierellaceae  | Mortierella        | Mortierella_minutissima |
| otu_114  | 0.6510157 | 0.9375 | 0.7812344 |  0.0155 | restored   | soil_saprotroph   | Mortierellomycota | Mortierellomycetes | Mortierellales      | Mortierellaceae  | Mortierella        | unidentified            |
| otu_2    | 0.6048530 | 1.0000 | 0.7777230 |  0.0295 | restored   | soil_saprotroph   | Mortierellomycota | Mortierellomycetes | Mortierellales      | Mortierellaceae  | Mortierella        | Mortierella_exigua      |
| otu_140  | 0.8447774 | 0.7500 | 0.7959793 |  0.0435 | remnant    | soil_saprotroph   | Ascomycota        | Sordariomycetes    | Hypocreales         | Stachybotryaceae | Striaticonidium    | Striaticonidium_cinctum |
| otu_2138 | 1.0000000 | 0.5000 | 0.7071068 |  0.0175 | remnant    | soil_saprotroph   | Ascomycota        | Leotiomycetes      | Thelebolales        | Pseudeurotiaceae | Gymnostellatospora | NA                      |

Indicator species of soil saprotrophs

A later task will be to comb these tables for species with good stories…

### Plant pathogens

#### Trends over time

``` r
guiltime("plant_pathogen")
```

    ## $bm_summary
    ## 
    ## Call:
    ## lm(formula = seq_sum ~ yr_since, data = d %>% filter(region == 
    ##     "BM"))
    ## 
    ## Residuals:
    ##       1       2       3       4       5       6       7 
    ##  -953.1 -1871.0   723.3 -2080.9  1732.4   111.4  2337.9 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept) 12953.84    1234.13  10.496 0.000135 ***
    ## yr_since     -361.61      83.02  -4.356 0.007319 ** 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 1884 on 5 degrees of freedom
    ## Multiple R-squared:  0.7914, Adjusted R-squared:  0.7497 
    ## F-statistic: 18.97 on 1 and 5 DF,  p-value: 0.007319
    ## 
    ## 
    ## $fl_summary
    ## 
    ## Call:
    ## lm(formula = seq_sum ~ yr_since, data = d %>% filter(region == 
    ##     "FL"))
    ## 
    ## Residuals:
    ##       1       2       3       4       5       6 
    ##   961.7  -664.8  -462.7   555.3   839.3 -1228.7 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)   
    ## (Intercept) 4277.575    826.145   5.178  0.00662 **
    ## yr_since       5.117     30.435   0.168  0.87463   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 1013 on 4 degrees of freedom
    ## Multiple R-squared:  0.007018,   Adjusted R-squared:  -0.2412 
    ## F-statistic: 0.02827 on 1 and 4 DF,  p-value: 0.8746
    ## 
    ## 
    ## $plot

<img src="microbial_guild_taxonomy_files/figure-gfm/ppat_guiltime-1.png" style="display: block; margin: auto;" />

A strong decline in pathogens is seen in Blue Mounds’ restored fields
($R^2_{Adj}=0.75, p<0.01$), and although two distinct groups are
apparent, no single site displays undue leverage. It’s possible that a
signal like this will be found in soil chemistry or plant data and can
help explain what we are seeing here. Recall also that AMF were
previously found to increase along this same sequence…maybe that will
still hold up.

#### Diversity

``` r
(ppat <- rerare(spe$its_raw, meta$its_raw, primary_lifestyle, "plant_pathogen", sites))
```

    ## $seq_depth
    ## [1] 2786
    ## 
    ## $zero_otu_num
    ## [1] 7
    ## 
    ## $rrfd
    ## # A tibble: 25 × 154
    ##    field_key otu_1 otu_3 otu_7 otu_13 otu_16 otu_21 otu_23 otu_28 otu_33 otu_43
    ##        <dbl> <dbl> <dbl> <dbl>  <dbl>  <dbl>  <dbl>  <dbl>  <dbl>  <dbl>  <dbl>
    ##  1         1   477    32  1583      8     56      0    164      0     29     12
    ##  2         2   819   255   139    200     64      6    100      9    226     29
    ##  3         3   115   121    61    749    230    755     71     55     40     48
    ##  4         4   728     2   691      0      0      0     60      0     37      0
    ##  5         5   120   100  1179     46    117     36     27     29     45     71
    ##  6         6   521    79   519      5    258    514    127     56     16    103
    ##  7         7   417   449    39    232    312    212     88    212     83     68
    ##  8         8   679   289   572     70     46     50     29     41    170     43
    ##  9         9   722   426   389    100     17      0    331      0    110      5
    ## 10        10   497   600   699     60     11      1     36      0    142      0
    ## # … with 15 more rows, and 143 more variables: otu_53 <dbl>, otu_58 <dbl>,
    ## #   otu_65 <dbl>, otu_68 <dbl>, otu_87 <dbl>, otu_99 <dbl>, otu_135 <dbl>,
    ## #   otu_137 <dbl>, otu_153 <dbl>, otu_172 <dbl>, otu_179 <dbl>, otu_200 <dbl>,
    ## #   otu_212 <dbl>, otu_279 <dbl>, otu_285 <dbl>, otu_289 <dbl>, otu_294 <dbl>,
    ## #   otu_304 <dbl>, otu_315 <dbl>, otu_319 <dbl>, otu_325 <dbl>, otu_332 <dbl>,
    ## #   otu_383 <dbl>, otu_391 <dbl>, otu_408 <dbl>, otu_416 <dbl>, otu_425 <dbl>,
    ## #   otu_432 <dbl>, otu_504 <dbl>, otu_511 <dbl>, otu_521 <dbl>, …
    ## 
    ## $rrfd_speTaxa
    ## # A tibble: 848 × 14
    ##    field_key otu_num seq_abund phylum   class order family genus species prima…¹
    ##        <dbl> <chr>       <dbl> <chr>    <chr> <chr> <chr>  <chr> <chr>   <chr>  
    ##  1         1 otu_1         477 Ascomyc… Sord… Hypo… Nectr… Fusa… Fusari… plant_…
    ##  2         1 otu_3          32 Ascomyc… Sord… Glom… Plect… Gibe… <NA>    plant_…
    ##  3         1 otu_7        1583 Ascomyc… Doth… Pleo… Peric… Peri… <NA>    plant_…
    ##  4         1 otu_13          8 Ascomyc… Sord… Glom… Plect… Plec… Plecto… plant_…
    ##  5         1 otu_16         56 Ascomyc… Sord… Hypo… Nectr… Nect… Nectri… plant_…
    ##  6         1 otu_23        164 Ascomyc… Doth… Pleo… Pleos… Alte… <NA>    plant_…
    ##  7         1 otu_33         29 Ascomyc… Sord… Hypo… Nectr… Fusa… <NA>    plant_…
    ##  8         1 otu_43         12 Ascomyc… Sord… Hypo… Nectr… Fusa… Fusari… plant_…
    ##  9         1 otu_58         26 Ascomyc… Doth… Pleo… Phaeo… Para… <NA>    plant_…
    ## 10         1 otu_65          4 Ascomyc… Sord… Hypo… Nectr… Gibb… Gibber… plant_…
    ## # … with 838 more rows, 4 more variables: field_name <chr>, region <chr>,
    ## #   field_type <ord>, yr_since <dbl>, and abbreviated variable name
    ## #   ¹​primary_lifestyle

``` r
(ppat_div <- calc_diversity(ppat$rrfd))
```

<div data-pagedtable="false">

<script data-pagedtable-source type="application/json">
{"columns":[{"label":["field_key"],"name":[1],"type":["dbl"],"align":["right"]},{"label":["field_name"],"name":[2],"type":["chr"],"align":["left"]},{"label":["region"],"name":[3],"type":["chr"],"align":["left"]},{"label":["field_type"],"name":[4],"type":["ord"],"align":["right"]},{"label":["yr_since"],"name":[5],"type":["dbl"],"align":["right"]},{"label":["hill_index"],"name":[6],"type":["ord"],"align":["right"]},{"label":["value"],"name":[7],"type":["dbl"],"align":["right"]}],"data":[{"1":"1","2":"BBRP1","3":"BM","4":"restored","5":"16","6":"N0","7":"24.00000000"},{"1":"1","2":"BBRP1","3":"BM","4":"restored","5":"16","6":"N1","7":"4.70238168"},{"1":"1","2":"BBRP1","3":"BM","4":"restored","5":"16","6":"N2","7":"2.75564635"},{"1":"1","2":"BBRP1","3":"BM","4":"restored","5":"16","6":"E10","7":"0.19593257"},{"1":"1","2":"BBRP1","3":"BM","4":"restored","5":"16","6":"E20","7":"0.11481860"},{"1":"2","2":"ERRP1","3":"BM","4":"restored","5":"3","6":"N0","7":"42.00000000"},{"1":"2","2":"ERRP1","3":"BM","4":"restored","5":"3","6":"N1","7":"13.13196914"},{"1":"2","2":"ERRP1","3":"BM","4":"restored","5":"3","6":"N2","7":"7.80282324"},{"1":"2","2":"ERRP1","3":"BM","4":"restored","5":"3","6":"E10","7":"0.31266593"},{"1":"2","2":"ERRP1","3":"BM","4":"restored","5":"3","6":"E20","7":"0.18578151"},{"1":"3","2":"FGC1","3":"FG","4":"corn","5":"0","6":"N0","7":"38.00000000"},{"1":"3","2":"FGC1","3":"FG","4":"corn","5":"0","6":"N1","7":"10.84358940"},{"1":"3","2":"FGC1","3":"FG","4":"corn","5":"0","6":"N2","7":"6.14939899"},{"1":"3","2":"FGC1","3":"FG","4":"corn","5":"0","6":"E10","7":"0.28535762"},{"1":"3","2":"FGC1","3":"FG","4":"corn","5":"0","6":"E20","7":"0.16182629"},{"1":"4","2":"FGREM1","3":"FG","4":"remnant","5":"NA","6":"N0","7":"27.00000000"},{"1":"4","2":"FGREM1","3":"FG","4":"remnant","5":"NA","6":"N1","7":"8.95621863"},{"1":"4","2":"FGREM1","3":"FG","4":"remnant","5":"NA","6":"N2","7":"6.07629307"},{"1":"4","2":"FGREM1","3":"FG","4":"remnant","5":"NA","6":"E10","7":"0.33171180"},{"1":"4","2":"FGREM1","3":"FG","4":"remnant","5":"NA","6":"E20","7":"0.22504789"},{"1":"5","2":"FGRP1","3":"FG","4":"restored","5":"15","6":"N0","7":"33.00000000"},{"1":"5","2":"FGRP1","3":"FG","4":"restored","5":"15","6":"N1","7":"7.22763545"},{"1":"5","2":"FGRP1","3":"FG","4":"restored","5":"15","6":"N2","7":"4.23173877"},{"1":"5","2":"FGRP1","3":"FG","4":"restored","5":"15","6":"E10","7":"0.21901926"},{"1":"5","2":"FGRP1","3":"FG","4":"restored","5":"15","6":"E20","7":"0.12823451"},{"1":"6","2":"FLC1","3":"FL","4":"corn","5":"0","6":"N0","7":"28.00000000"},{"1":"6","2":"FLC1","3":"FL","4":"corn","5":"0","6":"N1","7":"9.66205863"},{"1":"6","2":"FLC1","3":"FL","4":"corn","5":"0","6":"N2","7":"7.45879011"},{"1":"6","2":"FLC1","3":"FL","4":"corn","5":"0","6":"E10","7":"0.34507352"},{"1":"6","2":"FLC1","3":"FL","4":"corn","5":"0","6":"E20","7":"0.26638536"},{"1":"7","2":"FLC2","3":"FL","4":"corn","5":"0","6":"N0","7":"36.00000000"},{"1":"7","2":"FLC2","3":"FL","4":"corn","5":"0","6":"N1","7":"15.53008412"},{"1":"7","2":"FLC2","3":"FL","4":"corn","5":"0","6":"N2","7":"11.24826243"},{"1":"7","2":"FLC2","3":"FL","4":"corn","5":"0","6":"E10","7":"0.43139123"},{"1":"7","2":"FLC2","3":"FL","4":"corn","5":"0","6":"E20","7":"0.31245173"},{"1":"8","2":"FLREM1","3":"FL","4":"remnant","5":"NA","6":"N0","7":"38.00000000"},{"1":"8","2":"FLREM1","3":"FL","4":"remnant","5":"NA","6":"N1","7":"13.78493004"},{"1":"8","2":"FLREM1","3":"FL","4":"remnant","5":"NA","6":"N2","7":"8.00953079"},{"1":"8","2":"FLREM1","3":"FL","4":"remnant","5":"NA","6":"E10","7":"0.36276132"},{"1":"8","2":"FLREM1","3":"FL","4":"remnant","5":"NA","6":"E20","7":"0.21077713"},{"1":"9","2":"FLRP1","3":"FL","4":"restored","5":"40","6":"N0","7":"31.00000000"},{"1":"9","2":"FLRP1","3":"FL","4":"restored","5":"40","6":"N1","7":"11.09263807"},{"1":"9","2":"FLRP1","3":"FL","4":"restored","5":"40","6":"N2","7":"7.44109949"},{"1":"9","2":"FLRP1","3":"FL","4":"restored","5":"40","6":"E10","7":"0.35782703"},{"1":"9","2":"FLRP1","3":"FL","4":"restored","5":"40","6":"E20","7":"0.24003547"},{"1":"10","2":"FLRP4","3":"FL","4":"restored","5":"36","6":"N0","7":"31.00000000"},{"1":"10","2":"FLRP4","3":"FL","4":"restored","5":"36","6":"N1","7":"10.47022731"},{"1":"10","2":"FLRP4","3":"FL","4":"restored","5":"36","6":"N2","7":"6.62503307"},{"1":"10","2":"FLRP4","3":"FL","4":"restored","5":"36","6":"E10","7":"0.33774927"},{"1":"10","2":"FLRP4","3":"FL","4":"restored","5":"36","6":"E20","7":"0.21371074"},{"1":"11","2":"FLRP5","3":"FL","4":"restored","5":"35","6":"N0","7":"35.00000000"},{"1":"11","2":"FLRP5","3":"FL","4":"restored","5":"35","6":"N1","7":"16.09347537"},{"1":"11","2":"FLRP5","3":"FL","4":"restored","5":"35","6":"N2","7":"10.63341297"},{"1":"11","2":"FLRP5","3":"FL","4":"restored","5":"35","6":"E10","7":"0.45981358"},{"1":"11","2":"FLRP5","3":"FL","4":"restored","5":"35","6":"E20","7":"0.30381180"},{"1":"12","2":"FLRSP1","3":"FL","4":"restored","5":"10","6":"N0","7":"33.00000000"},{"1":"12","2":"FLRSP1","3":"FL","4":"restored","5":"10","6":"N1","7":"8.82851684"},{"1":"12","2":"FLRSP1","3":"FL","4":"restored","5":"10","6":"N2","7":"6.01566501"},{"1":"12","2":"FLRSP1","3":"FL","4":"restored","5":"10","6":"E10","7":"0.26753081"},{"1":"12","2":"FLRSP1","3":"FL","4":"restored","5":"10","6":"E20","7":"0.18229288"},{"1":"13","2":"FLRSP2","3":"FL","4":"restored","5":"10","6":"N0","7":"32.00000000"},{"1":"13","2":"FLRSP2","3":"FL","4":"restored","5":"10","6":"N1","7":"9.44485637"},{"1":"13","2":"FLRSP2","3":"FL","4":"restored","5":"10","6":"N2","7":"6.28784879"},{"1":"13","2":"FLRSP2","3":"FL","4":"restored","5":"10","6":"E10","7":"0.29515176"},{"1":"13","2":"FLRSP2","3":"FL","4":"restored","5":"10","6":"E20","7":"0.19649527"},{"1":"14","2":"FLRSP3","3":"FL","4":"restored","5":"10","6":"N0","7":"33.00000000"},{"1":"14","2":"FLRSP3","3":"FL","4":"restored","5":"10","6":"N1","7":"8.74171331"},{"1":"14","2":"FLRSP3","3":"FL","4":"restored","5":"10","6":"N2","7":"3.91821508"},{"1":"14","2":"FLRSP3","3":"FL","4":"restored","5":"10","6":"E10","7":"0.26490040"},{"1":"14","2":"FLRSP3","3":"FL","4":"restored","5":"10","6":"E20","7":"0.11873379"},{"1":"15","2":"KORP1","3":"BM","4":"restored","5":"28","6":"N0","7":"23.00000000"},{"1":"15","2":"KORP1","3":"BM","4":"restored","5":"28","6":"N1","7":"7.10127100"},{"1":"15","2":"KORP1","3":"BM","4":"restored","5":"28","6":"N2","7":"4.66777240"},{"1":"15","2":"KORP1","3":"BM","4":"restored","5":"28","6":"E10","7":"0.30875091"},{"1":"15","2":"KORP1","3":"BM","4":"restored","5":"28","6":"E20","7":"0.20294663"},{"1":"16","2":"LPC1","3":"LP","4":"corn","5":"0","6":"N0","7":"37.00000000"},{"1":"16","2":"LPC1","3":"LP","4":"corn","5":"0","6":"N1","7":"10.36015814"},{"1":"16","2":"LPC1","3":"LP","4":"corn","5":"0","6":"N2","7":"7.04256682"},{"1":"16","2":"LPC1","3":"LP","4":"corn","5":"0","6":"E10","7":"0.28000427"},{"1":"16","2":"LPC1","3":"LP","4":"corn","5":"0","6":"E20","7":"0.19033964"},{"1":"17","2":"LPREM1","3":"LP","4":"remnant","5":"NA","6":"N0","7":"29.00000000"},{"1":"17","2":"LPREM1","3":"LP","4":"remnant","5":"NA","6":"N1","7":"8.98197294"},{"1":"17","2":"LPREM1","3":"LP","4":"remnant","5":"NA","6":"N2","7":"5.90425754"},{"1":"17","2":"LPREM1","3":"LP","4":"remnant","5":"NA","6":"E10","7":"0.30972320"},{"1":"17","2":"LPREM1","3":"LP","4":"remnant","5":"NA","6":"E20","7":"0.20359509"},{"1":"18","2":"LPRP1","3":"LP","4":"restored","5":"4","6":"N0","7":"35.00000000"},{"1":"18","2":"LPRP1","3":"LP","4":"restored","5":"4","6":"N1","7":"10.62316038"},{"1":"18","2":"LPRP1","3":"LP","4":"restored","5":"4","6":"N2","7":"5.92056983"},{"1":"18","2":"LPRP1","3":"LP","4":"restored","5":"4","6":"E10","7":"0.30351887"},{"1":"18","2":"LPRP1","3":"LP","4":"restored","5":"4","6":"E20","7":"0.16915914"},{"1":"19","2":"LPRP2","3":"LP","4":"restored","5":"4","6":"N0","7":"42.00000000"},{"1":"19","2":"LPRP2","3":"LP","4":"restored","5":"4","6":"N1","7":"12.36279068"},{"1":"19","2":"LPRP2","3":"LP","4":"restored","5":"4","6":"N2","7":"6.65897630"},{"1":"19","2":"LPRP2","3":"LP","4":"restored","5":"4","6":"E10","7":"0.29435216"},{"1":"19","2":"LPRP2","3":"LP","4":"restored","5":"4","6":"E20","7":"0.15854705"},{"1":"20","2":"MBREM1","3":"BM","4":"remnant","5":"NA","6":"N0","7":"32.00000000"},{"1":"20","2":"MBREM1","3":"BM","4":"remnant","5":"NA","6":"N1","7":"6.66992069"},{"1":"20","2":"MBREM1","3":"BM","4":"remnant","5":"NA","6":"N2","7":"3.12776880"},{"1":"20","2":"MBREM1","3":"BM","4":"remnant","5":"NA","6":"E10","7":"0.20843502"},{"1":"20","2":"MBREM1","3":"BM","4":"remnant","5":"NA","6":"E20","7":"0.09774278"},{"1":"21","2":"MBRP1","3":"BM","4":"restored","5":"18","6":"N0","7":"39.00000000"},{"1":"21","2":"MBRP1","3":"BM","4":"restored","5":"18","6":"N1","7":"7.30842183"},{"1":"21","2":"MBRP1","3":"BM","4":"restored","5":"18","6":"N2","7":"3.39703756"},{"1":"21","2":"MBRP1","3":"BM","4":"restored","5":"18","6":"E10","7":"0.18739543"},{"1":"21","2":"MBRP1","3":"BM","4":"restored","5":"18","6":"E20","7":"0.08710353"},{"1":"22","2":"MHRP1","3":"BM","4":"restored","5":"7","6":"N0","7":"41.00000000"},{"1":"22","2":"MHRP1","3":"BM","4":"restored","5":"7","6":"N1","7":"8.96811439"},{"1":"22","2":"MHRP1","3":"BM","4":"restored","5":"7","6":"N2","7":"4.96135132"},{"1":"22","2":"MHRP1","3":"BM","4":"restored","5":"7","6":"E10","7":"0.21873450"},{"1":"22","2":"MHRP1","3":"BM","4":"restored","5":"7","6":"E20","7":"0.12100857"},{"1":"23","2":"MHRP2","3":"BM","4":"restored","5":"2","6":"N0","7":"44.00000000"},{"1":"23","2":"MHRP2","3":"BM","4":"restored","5":"2","6":"N1","7":"13.54135534"},{"1":"23","2":"MHRP2","3":"BM","4":"restored","5":"2","6":"N2","7":"7.83722142"},{"1":"23","2":"MHRP2","3":"BM","4":"restored","5":"2","6":"E10","7":"0.30775808"},{"1":"23","2":"MHRP2","3":"BM","4":"restored","5":"2","6":"E20","7":"0.17811867"},{"1":"24","2":"PHC1","3":"BM","4":"corn","5":"0","6":"N0","7":"27.00000000"},{"1":"24","2":"PHC1","3":"BM","4":"corn","5":"0","6":"N1","7":"11.79504627"},{"1":"24","2":"PHC1","3":"BM","4":"corn","5":"0","6":"N2","7":"9.15705673"},{"1":"24","2":"PHC1","3":"BM","4":"corn","5":"0","6":"E10","7":"0.43685357"},{"1":"24","2":"PHC1","3":"BM","4":"corn","5":"0","6":"E20","7":"0.33915025"},{"1":"25","2":"PHRP1","3":"BM","4":"restored","5":"11","6":"N0","7":"38.00000000"},{"1":"25","2":"PHRP1","3":"BM","4":"restored","5":"11","6":"N1","7":"13.27340226"},{"1":"25","2":"PHRP1","3":"BM","4":"restored","5":"11","6":"N2","7":"9.26432057"},{"1":"25","2":"PHRP1","3":"BM","4":"restored","5":"11","6":"E10","7":"0.34930006"},{"1":"25","2":"PHRP1","3":"BM","4":"restored","5":"11","6":"E20","7":"0.24379791"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>

</div>

``` r
(ppat_comp <- gudicom(ppat_div, ppat$rrfd_speTaxa, "plant_pathogen", other_threshold = 1))
```

    ## $Hills_field_type

<img src="microbial_guild_taxonomy_files/figure-gfm/ppat_composition-1.png" style="display: block; margin: auto;" />

    ## 
    ## $Hills_yrs_since_restoration

<img src="microbial_guild_taxonomy_files/figure-gfm/ppat_composition-2.png" style="display: block; margin: auto;" />

    ## 
    ## $Composition

<img src="microbial_guild_taxonomy_files/figure-gfm/ppat_composition-3.png" style="display: block; margin: auto;" />

<div data-pagedtable="false">

<script data-pagedtable-source type="application/json">
{"columns":[{"label":["field_type"],"name":[1],"type":["ord"],"align":["right"]},{"label":["order"],"name":[2],"type":["chr"],"align":["left"]},{"label":["seq_comp"],"name":[3],"type":["dbl"],"align":["right"]}],"data":[{"1":"corn","2":"Diaporthales","3":"1.292612"},{"1":"corn","2":"Glomerellales","3":"22.999095"},{"1":"corn","2":"Hypocreales","3":"36.867648"},{"1":"corn","2":"Magnaporthales","3":"1.233857"},{"1":"corn","2":"Other (OTU<1%)","3":"2.367830"},{"1":"corn","2":"Pleosporales","3":"35.238957"},{"1":"restored","2":"Cantharellales","3":"2.569449"},{"1":"restored","2":"Capnodiales","3":"1.831766"},{"1":"restored","2":"Glomerellales","3":"13.487499"},{"1":"restored","2":"Hypocreales","3":"42.904617"},{"1":"restored","2":"Other (OTU<1%)","3":"3.501145"},{"1":"restored","2":"Pleosporales","3":"32.927534"},{"1":"restored","2":"Ustilaginales","3":"1.146258"},{"1":"restored","2":"Xylariales","3":"1.631732"},{"1":"remnant","2":"Cantharellales","3":"3.668935"},{"1":"remnant","2":"Glomerellales","3":"9.318749"},{"1":"remnant","2":"Hypocreales","3":"54.017742"},{"1":"remnant","2":"Other (OTU<1%)","3":"2.067005"},{"1":"remnant","2":"Pezizales","3":"1.722505"},{"1":"remnant","2":"Pleosporales","3":"29.205064"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>

</div>

Richness and diversity look flat or declining from corn to remnants and
evenness takes a hit in restored and remnant fields. It looks like we
have fewer pathogens, but more dominant individual taxa become
established. Pathogen diversity decreases with years since restoration
in Blue Mounds, but if the dumbbell plots can be believed, the opposite
appears true in Fermi.

Many pathogen orders are rare, so the argument `other_threshold` was
adjusted to show more diversity. Shifts don’t appear pronounced.
*Diaporthales* decreases in composition from corn to remnant while
*Hypocreales* pathogens increase. *Cantharellales* appear a small
component but are possibly “late successional” pathogens, possibly
associated with some native plant in a plant-soil feedback.

#### Indicators

``` r
ppat_inspan <- 
    ppat$rrfd %>% 
    left_join(sites, by = join_by(field_key)) %>% 
    inspan(., 1999, meta$its_raw)
```

``` r
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
```

| field_type | n_otu |  stat_avg |   stat_sd |
|:-----------|------:|----------:|----------:|
| corn       |    13 | 0.8547888 | 0.1101429 |
| restored   |     1 | 0.8577461 |        NA |
| remnant    |     5 | 0.7301094 | 0.0763120 |

Indicator species stats: plant pathogens

We see the same trend as before, where more indicators are found in
cornfields, and their indicator stats are stronger. Composition at the
level of taxonomic order isn’t telling the whole story.

Plant pathogen indicators are nearly all in *Ascomycota.*

``` r
ppat_inspan %>% 
    mutate(field_type = factor(
        field_type,
        ordered = TRUE,
        levels = c("corn", "restored", "remnant")
    )) %>%
    arrange(field_type, -stat) %>% 
    kable(format = "pandoc", caption = "Indicator species of plant pathogens")
```

| otu_num  |         A |      B |      stat | p.value | field_type | primary_lifestyle | phylum        | class             | order          | family               | genus            | species                     |
|:---------|----------:|-------:|----------:|--------:|:-----------|:------------------|:--------------|:------------------|:---------------|:---------------------|:-----------------|:----------------------------|
| otu_172  | 0.9781659 | 1.0000 | 0.9890227 |  0.0010 | corn       | plant_pathogen    | Ascomycota    | Dothideomycetes   | Pleosporales   | Corynesporascaceae   | Corynespora      | Corynespora_cassiicola      |
| otu_200  | 0.9457006 | 1.0000 | 0.9724714 |  0.0015 | corn       | plant_pathogen    | Ascomycota    | Dothideomycetes   | Pleosporales   | Phaeosphaeriaceae    | Ophiosphaerella  | unidentified                |
| otu_553  | 0.9256198 | 1.0000 | 0.9620914 |  0.0060 | corn       | plant_pathogen    | Ascomycota    | Sordariomycetes   | Magnaporthales | Magnaporthaceae      | Gaeumannomyces   | NA                          |
| otu_21   | 0.9229813 | 1.0000 | 0.9607192 |  0.0010 | corn       | plant_pathogen    | Ascomycota    | Dothideomycetes   | Pleosporales   | Phaeosphaeriaceae    | Setophoma        | Setophoma_terrestris        |
| otu_1841 | 1.0000000 | 0.8000 | 0.8944272 |  0.0015 | corn       | plant_pathogen    | Ascomycota    | Dothideomycetes   | Pleosporales   | Pleosporaceae        | Curvularia       | NA                          |
| otu_432  | 0.9952517 | 0.8000 | 0.8923011 |  0.0015 | corn       | plant_pathogen    | Ascomycota    | Sordariomycetes   | Glomerellales  | Glomerellaceae       | Colletotrichum   | NA                          |
| otu_391  | 0.7493404 | 1.0000 | 0.8656445 |  0.0110 | corn       | plant_pathogen    | Ascomycota    | Dothideomycetes   | Pleosporales   | Torulaceae           | Dendryphion      | NA                          |
| otu_13   | 0.7329077 | 1.0000 | 0.8561003 |  0.0070 | corn       | plant_pathogen    | Ascomycota    | Sordariomycetes   | Glomerellales  | Plectosphaerellaceae | Plectosphaerella | Plectosphaerella_cucumerina |
| otu_796  | 0.9084249 | 0.8000 | 0.8524904 |  0.0055 | corn       | plant_pathogen    | Ascomycota    | Dothideomycetes   | Capnodiales    | Mycosphaerellaceae   | Cercospora       | NA                          |
| otu_325  | 1.0000000 | 0.6000 | 0.7745967 |  0.0060 | corn       | plant_pathogen    | Ascomycota    | Sordariomycetes   | Diaporthales   | Diaporthaceae        | Diaporthe        | NA                          |
| otu_521  | 0.9388753 | 0.6000 | 0.7505499 |  0.0220 | corn       | plant_pathogen    | Ascomycota    | Sordariomycetes   | Glomerellales  | Plectosphaerellaceae | Lectera          | NA                          |
| otu_1013 | 0.8387097 | 0.6000 | 0.7093841 |  0.0430 | corn       | plant_pathogen    | Ascomycota    | Sordariomycetes   | Xylariales     | Microdochiaceae      | Microdochium     | Microdochium_colombiense    |
| otu_2242 | 1.0000000 | 0.4000 | 0.6324555 |  0.0445 | corn       | plant_pathogen    | Ascomycota    | Sordariomycetes   | Diaporthales   | Diaporthaceae        | Phaeocytostroma  | Phaeocytostroma_ambiguum    |
| otu_332  | 0.9055118 | 0.8125 | 0.8577461 |  0.0275 | restored   | plant_pathogen    | Ascomycota    | Sordariomycetes   | Glomerellales  | Plectosphaerellaceae | Plectosphaerella | NA                          |
| otu_135  | 0.7471148 | 1.0000 | 0.8643580 |  0.0175 | remnant    | plant_pathogen    | Ascomycota    | Sordariomycetes   | Hypocreales    | Nectriaceae          | Ilyonectria      | NA                          |
| otu_1716 | 1.0000000 | 0.5000 | 0.7071068 |  0.0215 | remnant    | plant_pathogen    | Ascomycota    | Sordariomycetes   | Hypocreales    | Nectriaceae          | Volutella        | NA                          |
| otu_942  | 0.9944751 | 0.5000 | 0.7051507 |  0.0215 | remnant    | plant_pathogen    | Ascomycota    | Dothideomycetes   | Pleosporales   | Pleosporaceae        | Curvularia       | NA                          |
| otu_319  | 0.6552901 | 0.7500 | 0.7010475 |  0.0355 | remnant    | plant_pathogen    | Basidiomycota | Ustilaginomycetes | Ustilaginales  | Ustilaginaceae       | Ustilago         | Ustilago_nunavutica         |
| otu_1    | 0.4527729 | 1.0000 | 0.6728840 |  0.0490 | remnant    | plant_pathogen    | Ascomycota    | Sordariomycetes   | Hypocreales    | Nectriaceae          | Fusarium         | Fusarium_oxysporum          |

Indicator species of plant pathogens

### Wood saprotrophs

#### Trends over time

``` r
guiltime("wood_saprotroph") 
```

    ## $bm_summary
    ## 
    ## Call:
    ## lm(formula = seq_sum ~ yr_since, data = d %>% filter(region == 
    ##     "BM"))
    ## 
    ## Residuals:
    ##       1       2       3       4       5       6       7 
    ## -829.91   59.82  150.26  371.29  -18.79   66.72  200.60 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)  2891.47     277.63  10.415 0.000141 ***
    ## yr_since      -78.10      18.68  -4.182 0.008639 ** 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 423.7 on 5 degrees of freedom
    ## Multiple R-squared:  0.7777, Adjusted R-squared:  0.7332 
    ## F-statistic: 17.49 on 1 and 5 DF,  p-value: 0.008639
    ## 
    ## 
    ## $fl_summary
    ## 
    ## Call:
    ## lm(formula = seq_sum ~ yr_since, data = d %>% filter(region == 
    ##     "FL"))
    ## 
    ## Residuals:
    ##       1       2       3       4       5       6 
    ##  -165.3   414.4  -232.7  -378.8  1768.2 -1405.8 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)  
    ## (Intercept)  2498.68     956.41   2.613   0.0593 .
    ## yr_since      -21.09      35.23  -0.598   0.5818  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 1173 on 4 degrees of freedom
    ## Multiple R-squared:  0.08218,    Adjusted R-squared:  -0.1473 
    ## F-statistic: 0.3581 on 1 and 4 DF,  p-value: 0.5818
    ## 
    ## 
    ## $plot

<img src="microbial_guild_taxonomy_files/figure-gfm/wsap_guiltime-1.png" style="display: block; margin: auto;" />

Interestingly a strong negative relationship over time since restoration
($R^2_{Adj}=0.73, p<0.01$) in sharp contrast to the increasing
relationship found with soil saprotrophs. Apparently many wood
saprotrophs live in cornfield soil…let’s see:

#### Diversity

``` r
(wsap <- rerare(spe$its_raw, meta$its_raw, primary_lifestyle, "wood_saprotroph", sites))
```

    ## $seq_depth
    ## [1] 701
    ## 
    ## $zero_otu_num
    ## [1] 6
    ## 
    ## $rrfd
    ## # A tibble: 25 × 119
    ##    field_key otu_11 otu_20 otu_29 otu_39 otu_76 otu_117 otu_120 otu_130 otu_169
    ##        <dbl>  <dbl>  <dbl>  <dbl>  <dbl>  <dbl>   <dbl>   <dbl>   <dbl>   <dbl>
    ##  1         1      0     42    225    177      0      11       0     122       0
    ##  2         2     23    182      9    129    182      11       8       0       1
    ##  3         3     95    212      0     11      9       0       0       0       0
    ##  4         4      0     71     47    302      0       0       0       0       0
    ##  5         5     25     40      0    382      0       0       0      85      44
    ##  6         6    494      1    143      9      0      10       0       0       0
    ##  7         7    135    302     46     16      6       0       0       0       0
    ##  8         8     29     69     37    213      0      18       0       0       0
    ##  9         9     24     24    128    211      0       0       9      55       0
    ## 10        10      9     54     77    161      8       5     101     147       0
    ## # … with 15 more rows, and 109 more variables: otu_202 <dbl>, otu_252 <dbl>,
    ## #   otu_266 <dbl>, otu_287 <dbl>, otu_322 <dbl>, otu_329 <dbl>, otu_331 <dbl>,
    ## #   otu_333 <dbl>, otu_341 <dbl>, otu_365 <dbl>, otu_370 <dbl>, otu_397 <dbl>,
    ## #   otu_415 <dbl>, otu_437 <dbl>, otu_438 <dbl>, otu_487 <dbl>, otu_508 <dbl>,
    ## #   otu_589 <dbl>, otu_599 <dbl>, otu_606 <dbl>, otu_632 <dbl>, otu_633 <dbl>,
    ## #   otu_634 <dbl>, otu_703 <dbl>, otu_770 <dbl>, otu_786 <dbl>, otu_790 <dbl>,
    ## #   otu_793 <dbl>, otu_818 <dbl>, otu_852 <dbl>, otu_853 <dbl>, …
    ## 
    ## $rrfd_speTaxa
    ## # A tibble: 465 × 14
    ##    field_key otu_num seq_abund phylum   class order family genus species prima…¹
    ##        <dbl> <chr>       <dbl> <chr>    <chr> <chr> <chr>  <chr> <chr>   <chr>  
    ##  1         1 otu_20         42 Ascomyc… Sord… Hypo… Bione… Clon… <NA>    wood_s…
    ##  2         1 otu_29        225 Ascomyc… Sord… Hypo… Nectr… Mari… Marian… wood_s…
    ##  3         1 otu_39        177 Ascomyc… Doth… Pleo… Cucur… Pyre… uniden… wood_s…
    ##  4         1 otu_117        11 Ascomyc… Doth… Pleo… Lindg… Cloh… Clohes… wood_s…
    ##  5         1 otu_130       122 Basidio… Agar… Trec… Hydno… Subu… <NA>    wood_s…
    ##  6         1 otu_333         4 Ascomyc… Leot… Helo… Helot… Scyt… uniden… wood_s…
    ##  7         1 otu_415        30 Ascomyc… Leot… Helo… Helot… Scyt… Scytal… wood_s…
    ##  8         1 otu_508        15 Ascomyc… Doth… Pleo… Lenti… Keis… Keissl… wood_s…
    ##  9         1 otu_599         5 Ascomyc… Doth… Pleo… Didym… Para… Paraph… wood_s…
    ## 10         1 otu_633         6 Ascomyc… Doth… Pleo… Lenti… Keis… Keissl… wood_s…
    ## # … with 455 more rows, 4 more variables: field_name <chr>, region <chr>,
    ## #   field_type <ord>, yr_since <dbl>, and abbreviated variable name
    ## #   ¹​primary_lifestyle

Sequence depth is low; these aren’t abundant taxa.

``` r
(wsap_div <- calc_diversity(wsap$rrfd))
```

<div data-pagedtable="false">

<script data-pagedtable-source type="application/json">
{"columns":[{"label":["field_key"],"name":[1],"type":["dbl"],"align":["right"]},{"label":["field_name"],"name":[2],"type":["chr"],"align":["left"]},{"label":["region"],"name":[3],"type":["chr"],"align":["left"]},{"label":["field_type"],"name":[4],"type":["ord"],"align":["right"]},{"label":["yr_since"],"name":[5],"type":["dbl"],"align":["right"]},{"label":["hill_index"],"name":[6],"type":["ord"],"align":["right"]},{"label":["value"],"name":[7],"type":["dbl"],"align":["right"]}],"data":[{"1":"1","2":"BBRP1","3":"BM","4":"restored","5":"16","6":"N0","7":"16.0000000"},{"1":"1","2":"BBRP1","3":"BM","4":"restored","5":"16","6":"N1","7":"6.6404678"},{"1":"1","2":"BBRP1","3":"BM","4":"restored","5":"16","6":"N2","7":"4.8201615"},{"1":"1","2":"BBRP1","3":"BM","4":"restored","5":"16","6":"E10","7":"0.4150292"},{"1":"1","2":"BBRP1","3":"BM","4":"restored","5":"16","6":"E20","7":"0.3012601"},{"1":"2","2":"ERRP1","3":"BM","4":"restored","5":"3","6":"N0","7":"26.0000000"},{"1":"2","2":"ERRP1","3":"BM","4":"restored","5":"3","6":"N1","7":"7.8801363"},{"1":"2","2":"ERRP1","3":"BM","4":"restored","5":"3","6":"N2","7":"5.4741832"},{"1":"2","2":"ERRP1","3":"BM","4":"restored","5":"3","6":"E10","7":"0.3030822"},{"1":"2","2":"ERRP1","3":"BM","4":"restored","5":"3","6":"E20","7":"0.2105455"},{"1":"3","2":"FGC1","3":"FG","4":"corn","5":"0","6":"N0","7":"16.0000000"},{"1":"3","2":"FGC1","3":"FG","4":"corn","5":"0","6":"N1","7":"4.4628829"},{"1":"3","2":"FGC1","3":"FG","4":"corn","5":"0","6":"N2","7":"3.2036678"},{"1":"3","2":"FGC1","3":"FG","4":"corn","5":"0","6":"E10","7":"0.2789302"},{"1":"3","2":"FGC1","3":"FG","4":"corn","5":"0","6":"E20","7":"0.2002292"},{"1":"4","2":"FGREM1","3":"FG","4":"remnant","5":"NA","6":"N0","7":"17.0000000"},{"1":"4","2":"FGREM1","3":"FG","4":"remnant","5":"NA","6":"N1","7":"6.8826617"},{"1":"4","2":"FGREM1","3":"FG","4":"remnant","5":"NA","6":"N2","7":"4.2862090"},{"1":"4","2":"FGREM1","3":"FG","4":"remnant","5":"NA","6":"E10","7":"0.4048625"},{"1":"4","2":"FGREM1","3":"FG","4":"remnant","5":"NA","6":"E20","7":"0.2521299"},{"1":"5","2":"FGRP1","3":"FG","4":"restored","5":"15","6":"N0","7":"16.0000000"},{"1":"5","2":"FGRP1","3":"FG","4":"restored","5":"15","6":"N1","7":"5.2307053"},{"1":"5","2":"FGRP1","3":"FG","4":"restored","5":"15","6":"N2","7":"3.0422220"},{"1":"5","2":"FGRP1","3":"FG","4":"restored","5":"15","6":"E10","7":"0.3269191"},{"1":"5","2":"FGRP1","3":"FG","4":"restored","5":"15","6":"E20","7":"0.1901389"},{"1":"6","2":"FLC1","3":"FL","4":"corn","5":"0","6":"N0","7":"10.0000000"},{"1":"6","2":"FLC1","3":"FL","4":"corn","5":"0","6":"N1","7":"2.5774792"},{"1":"6","2":"FLC1","3":"FL","4":"corn","5":"0","6":"N2","7":"1.8518068"},{"1":"6","2":"FLC1","3":"FL","4":"corn","5":"0","6":"E10","7":"0.2577479"},{"1":"6","2":"FLC1","3":"FL","4":"corn","5":"0","6":"E20","7":"0.1851807"},{"1":"7","2":"FLC2","3":"FL","4":"corn","5":"0","6":"N0","7":"21.0000000"},{"1":"7","2":"FLC2","3":"FL","4":"corn","5":"0","6":"N1","7":"6.9456204"},{"1":"7","2":"FLC2","3":"FL","4":"corn","5":"0","6":"N2","7":"4.1509423"},{"1":"7","2":"FLC2","3":"FL","4":"corn","5":"0","6":"E10","7":"0.3307438"},{"1":"7","2":"FLC2","3":"FL","4":"corn","5":"0","6":"E20","7":"0.1976639"},{"1":"8","2":"FLREM1","3":"FL","4":"remnant","5":"NA","6":"N0","7":"21.0000000"},{"1":"8","2":"FLREM1","3":"FL","4":"remnant","5":"NA","6":"N1","7":"9.2053471"},{"1":"8","2":"FLREM1","3":"FL","4":"remnant","5":"NA","6":"N2","7":"5.9834281"},{"1":"8","2":"FLREM1","3":"FL","4":"remnant","5":"NA","6":"E10","7":"0.4383499"},{"1":"8","2":"FLREM1","3":"FL","4":"remnant","5":"NA","6":"E20","7":"0.2849251"},{"1":"9","2":"FLRP1","3":"FL","4":"restored","5":"40","6":"N0","7":"27.0000000"},{"1":"9","2":"FLRP1","3":"FL","4":"restored","5":"40","6":"N1","7":"11.5385536"},{"1":"9","2":"FLRP1","3":"FL","4":"restored","5":"40","6":"N2","7":"6.8940501"},{"1":"9","2":"FLRP1","3":"FL","4":"restored","5":"40","6":"E10","7":"0.4273538"},{"1":"9","2":"FLRP1","3":"FL","4":"restored","5":"40","6":"E20","7":"0.2553352"},{"1":"10","2":"FLRP4","3":"FL","4":"restored","5":"36","6":"N0","7":"22.0000000"},{"1":"10","2":"FLRP4","3":"FL","4":"restored","5":"36","6":"N1","7":"9.1481665"},{"1":"10","2":"FLRP4","3":"FL","4":"restored","5":"36","6":"N2","7":"6.9016025"},{"1":"10","2":"FLRP4","3":"FL","4":"restored","5":"36","6":"E10","7":"0.4158258"},{"1":"10","2":"FLRP4","3":"FL","4":"restored","5":"36","6":"E20","7":"0.3137092"},{"1":"11","2":"FLRP5","3":"FL","4":"restored","5":"35","6":"N0","7":"24.0000000"},{"1":"11","2":"FLRP5","3":"FL","4":"restored","5":"35","6":"N1","7":"6.4510404"},{"1":"11","2":"FLRP5","3":"FL","4":"restored","5":"35","6":"N2","7":"3.0193424"},{"1":"11","2":"FLRP5","3":"FL","4":"restored","5":"35","6":"E10","7":"0.2687934"},{"1":"11","2":"FLRP5","3":"FL","4":"restored","5":"35","6":"E20","7":"0.1258059"},{"1":"12","2":"FLRSP1","3":"FL","4":"restored","5":"10","6":"N0","7":"23.0000000"},{"1":"12","2":"FLRSP1","3":"FL","4":"restored","5":"10","6":"N1","7":"10.3388369"},{"1":"12","2":"FLRSP1","3":"FL","4":"restored","5":"10","6":"N2","7":"7.8043516"},{"1":"12","2":"FLRSP1","3":"FL","4":"restored","5":"10","6":"E10","7":"0.4495146"},{"1":"12","2":"FLRSP1","3":"FL","4":"restored","5":"10","6":"E20","7":"0.3393196"},{"1":"13","2":"FLRSP2","3":"FL","4":"restored","5":"10","6":"N0","7":"15.0000000"},{"1":"13","2":"FLRSP2","3":"FL","4":"restored","5":"10","6":"N1","7":"5.1239540"},{"1":"13","2":"FLRSP2","3":"FL","4":"restored","5":"10","6":"N2","7":"3.3103237"},{"1":"13","2":"FLRSP2","3":"FL","4":"restored","5":"10","6":"E10","7":"0.3415969"},{"1":"13","2":"FLRSP2","3":"FL","4":"restored","5":"10","6":"E20","7":"0.2206882"},{"1":"14","2":"FLRSP3","3":"FL","4":"restored","5":"10","6":"N0","7":"19.0000000"},{"1":"14","2":"FLRSP3","3":"FL","4":"restored","5":"10","6":"N1","7":"6.8828932"},{"1":"14","2":"FLRSP3","3":"FL","4":"restored","5":"10","6":"N2","7":"4.0193442"},{"1":"14","2":"FLRSP3","3":"FL","4":"restored","5":"10","6":"E10","7":"0.3622575"},{"1":"14","2":"FLRSP3","3":"FL","4":"restored","5":"10","6":"E20","7":"0.2115444"},{"1":"15","2":"KORP1","3":"BM","4":"restored","5":"28","6":"N0","7":"19.0000000"},{"1":"15","2":"KORP1","3":"BM","4":"restored","5":"28","6":"N1","7":"8.3328848"},{"1":"15","2":"KORP1","3":"BM","4":"restored","5":"28","6":"N2","7":"5.5465997"},{"1":"15","2":"KORP1","3":"BM","4":"restored","5":"28","6":"E10","7":"0.4385729"},{"1":"15","2":"KORP1","3":"BM","4":"restored","5":"28","6":"E20","7":"0.2919263"},{"1":"16","2":"LPC1","3":"LP","4":"corn","5":"0","6":"N0","7":"15.0000000"},{"1":"16","2":"LPC1","3":"LP","4":"corn","5":"0","6":"N1","7":"4.2439412"},{"1":"16","2":"LPC1","3":"LP","4":"corn","5":"0","6":"N2","7":"2.9199829"},{"1":"16","2":"LPC1","3":"LP","4":"corn","5":"0","6":"E10","7":"0.2829294"},{"1":"16","2":"LPC1","3":"LP","4":"corn","5":"0","6":"E20","7":"0.1946655"},{"1":"17","2":"LPREM1","3":"LP","4":"remnant","5":"NA","6":"N0","7":"23.0000000"},{"1":"17","2":"LPREM1","3":"LP","4":"remnant","5":"NA","6":"N1","7":"7.3483313"},{"1":"17","2":"LPREM1","3":"LP","4":"remnant","5":"NA","6":"N2","7":"4.3433387"},{"1":"17","2":"LPREM1","3":"LP","4":"remnant","5":"NA","6":"E10","7":"0.3194927"},{"1":"17","2":"LPREM1","3":"LP","4":"remnant","5":"NA","6":"E20","7":"0.1888408"},{"1":"18","2":"LPRP1","3":"LP","4":"restored","5":"4","6":"N0","7":"19.0000000"},{"1":"18","2":"LPRP1","3":"LP","4":"restored","5":"4","6":"N1","7":"7.6405610"},{"1":"18","2":"LPRP1","3":"LP","4":"restored","5":"4","6":"N2","7":"5.4454295"},{"1":"18","2":"LPRP1","3":"LP","4":"restored","5":"4","6":"E10","7":"0.4021348"},{"1":"18","2":"LPRP1","3":"LP","4":"restored","5":"4","6":"E20","7":"0.2866016"},{"1":"19","2":"LPRP2","3":"LP","4":"restored","5":"4","6":"N0","7":"17.0000000"},{"1":"19","2":"LPRP2","3":"LP","4":"restored","5":"4","6":"N1","7":"6.7266943"},{"1":"19","2":"LPRP2","3":"LP","4":"restored","5":"4","6":"N2","7":"4.9774728"},{"1":"19","2":"LPRP2","3":"LP","4":"restored","5":"4","6":"E10","7":"0.3956879"},{"1":"19","2":"LPRP2","3":"LP","4":"restored","5":"4","6":"E20","7":"0.2927925"},{"1":"20","2":"MBREM1","3":"BM","4":"remnant","5":"NA","6":"N0","7":"15.0000000"},{"1":"20","2":"MBREM1","3":"BM","4":"remnant","5":"NA","6":"N1","7":"7.1531922"},{"1":"20","2":"MBREM1","3":"BM","4":"remnant","5":"NA","6":"N2","7":"5.7643022"},{"1":"20","2":"MBREM1","3":"BM","4":"remnant","5":"NA","6":"E10","7":"0.4768795"},{"1":"20","2":"MBREM1","3":"BM","4":"remnant","5":"NA","6":"E20","7":"0.3842868"},{"1":"21","2":"MBRP1","3":"BM","4":"restored","5":"18","6":"N0","7":"19.0000000"},{"1":"21","2":"MBRP1","3":"BM","4":"restored","5":"18","6":"N1","7":"7.7658197"},{"1":"21","2":"MBRP1","3":"BM","4":"restored","5":"18","6":"N2","7":"4.8884434"},{"1":"21","2":"MBRP1","3":"BM","4":"restored","5":"18","6":"E10","7":"0.4087274"},{"1":"21","2":"MBRP1","3":"BM","4":"restored","5":"18","6":"E20","7":"0.2572865"},{"1":"22","2":"MHRP1","3":"BM","4":"restored","5":"7","6":"N0","7":"18.0000000"},{"1":"22","2":"MHRP1","3":"BM","4":"restored","5":"7","6":"N1","7":"9.5381111"},{"1":"22","2":"MHRP1","3":"BM","4":"restored","5":"7","6":"N2","7":"7.6076510"},{"1":"22","2":"MHRP1","3":"BM","4":"restored","5":"7","6":"E10","7":"0.5298951"},{"1":"22","2":"MHRP1","3":"BM","4":"restored","5":"7","6":"E20","7":"0.4226473"},{"1":"23","2":"MHRP2","3":"BM","4":"restored","5":"2","6":"N0","7":"16.0000000"},{"1":"23","2":"MHRP2","3":"BM","4":"restored","5":"2","6":"N1","7":"4.8992435"},{"1":"23","2":"MHRP2","3":"BM","4":"restored","5":"2","6":"N2","7":"3.1867976"},{"1":"23","2":"MHRP2","3":"BM","4":"restored","5":"2","6":"E10","7":"0.3062027"},{"1":"23","2":"MHRP2","3":"BM","4":"restored","5":"2","6":"E20","7":"0.1991748"},{"1":"24","2":"PHC1","3":"BM","4":"corn","5":"0","6":"N0","7":"16.0000000"},{"1":"24","2":"PHC1","3":"BM","4":"corn","5":"0","6":"N1","7":"5.9282632"},{"1":"24","2":"PHC1","3":"BM","4":"corn","5":"0","6":"N2","7":"4.0533601"},{"1":"24","2":"PHC1","3":"BM","4":"corn","5":"0","6":"E10","7":"0.3705164"},{"1":"24","2":"PHC1","3":"BM","4":"corn","5":"0","6":"E20","7":"0.2533350"},{"1":"25","2":"PHRP1","3":"BM","4":"restored","5":"11","6":"N0","7":"15.0000000"},{"1":"25","2":"PHRP1","3":"BM","4":"restored","5":"11","6":"N1","7":"5.3644015"},{"1":"25","2":"PHRP1","3":"BM","4":"restored","5":"11","6":"N2","7":"3.8460112"},{"1":"25","2":"PHRP1","3":"BM","4":"restored","5":"11","6":"E10","7":"0.3576268"},{"1":"25","2":"PHRP1","3":"BM","4":"restored","5":"11","6":"E20","7":"0.2564007"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>

</div>

``` r
(wasp_comp <- gudicom(wsap_div, wsap$rrfd_speTaxa, "wood_saprotroph"))
```

    ## $Hills_field_type

<img src="microbial_guild_taxonomy_files/figure-gfm/wsap_composition-1.png" style="display: block; margin: auto;" />

    ## 
    ## $Hills_yrs_since_restoration

<img src="microbial_guild_taxonomy_files/figure-gfm/wsap_composition-2.png" style="display: block; margin: auto;" />

    ## 
    ## $Composition

<img src="microbial_guild_taxonomy_files/figure-gfm/wsap_composition-3.png" style="display: block; margin: auto;" />

<div data-pagedtable="false">

<script data-pagedtable-source type="application/json">
{"columns":[{"label":["field_type"],"name":[1],"type":["ord"],"align":["right"]},{"label":["order"],"name":[2],"type":["chr"],"align":["left"]},{"label":["seq_comp"],"name":[3],"type":["dbl"],"align":["right"]}],"data":[{"1":"corn","2":"Agaricales","3":"16.113753"},{"1":"corn","2":"Hypocreales","3":"24.554290"},{"1":"corn","2":"Other (OTU<2%)","3":"7.512173"},{"1":"corn","2":"Phomatosporales","3":"16.596881"},{"1":"corn","2":"Sordariales","3":"31.215778"},{"1":"corn","2":"Trechisporales","3":"4.007124"},{"1":"restored","2":"Agaricales","3":"2.615450"},{"1":"restored","2":"Auriculariales","3":"2.408394"},{"1":"restored","2":"Chaetothyriales","3":"2.975075"},{"1":"restored","2":"Helotiales","3":"3.433536"},{"1":"restored","2":"Hymenochaetales","3":"4.586120"},{"1":"restored","2":"Hypocreales","3":"21.139855"},{"1":"restored","2":"Minutisphaerales","3":"3.423697"},{"1":"restored","2":"Other (OTU<2%)","3":"6.316131"},{"1":"restored","2":"Pezizales","3":"2.998687"},{"1":"restored","2":"Pleosporales","3":"27.876570"},{"1":"restored","2":"Sordariales","3":"7.941404"},{"1":"restored","2":"Trechisporales","3":"14.285082"},{"1":"remnant","2":"Agaricales","3":"2.103971"},{"1":"remnant","2":"Auriculariales","3":"2.156571"},{"1":"remnant","2":"Chaetothyriales","3":"12.150434"},{"1":"remnant","2":"Helotiales","3":"8.152889"},{"1":"remnant","2":"Hymenochaetales","3":"3.418953"},{"1":"remnant","2":"Hypocreales","3":"20.014026"},{"1":"remnant","2":"Minutisphaerales","3":"9.678268"},{"1":"remnant","2":"Other (OTU<2%)","3":"4.330674"},{"1":"remnant","2":"Pezizales","3":"7.784694"},{"1":"remnant","2":"Pleosporales","3":"24.669063"},{"1":"remnant","2":"Sordariales","3":"5.540458"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>

</div>

With diversity, not much jumps out.

Diversity appears high across fields and years compared with other
guilds. While *Agaric* soil saprotrophs increased strongly from corn to
remnants, they declined when characterized as wood saprotrophs.

#### Indicators

``` r
wsap_inspan <- 
    wsap$rrfd %>% 
    left_join(sites, by = join_by(field_key)) %>% 
    inspan(., 1999, meta$its_raw)
```

``` r
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
```

| field_type | n_otu | stat_avg |   stat_sd |
|:-----------|------:|---------:|----------:|
| corn       |     3 | 0.808547 | 0.0941557 |
| remnant    |     3 | 0.696997 | 0.0679090 |

Indicator species stats: wood saprotrophs

Few species show specificity or fidelity. Corn fields have a few unusual
taxa, though. Less so with remnants, and none with restored fields.

``` r
wsap_inspan %>% 
    mutate(field_type = factor(
        field_type,
        ordered = TRUE,
        levels = c("corn", "restored", "remnant")
    )) %>%
    arrange(field_type, -stat) %>% 
    kable(format = "pandoc", caption = "Indicator species of wood saprotrophs")
```

| otu_num |         A |    B |      stat | p.value | field_type | primary_lifestyle | phylum        | class           | order            | family              | genus             | species                    |
|:--------|----------:|-----:|----------:|--------:|:-----------|:------------------|:--------------|:----------------|:-----------------|:--------------------|:------------------|:---------------------------|
| otu_589 | 0.9504950 | 0.80 | 0.8720069 |  0.0020 | corn       | wood_saprotroph   | Ascomycota    | Sordariomycetes | Hypocreales      | Stachybotryaceae    | Stachybotrys      | Stachybotrys_limonispora   |
| otu_11  | 0.7280681 | 1.00 | 0.8532691 |  0.0045 | corn       | wood_saprotroph   | Ascomycota    | Sordariomycetes | Sordariales      | Chaetomiaceae       | Humicola          | Humicola_grisea            |
| otu_341 | 0.8175182 | 0.60 | 0.7003649 |  0.0375 | corn       | wood_saprotroph   | Basidiomycota | Agaricomycetes  | Agaricales       | Psathyrellaceae     | Psathyrella       | NA                         |
| otu_599 | 0.7822878 | 0.75 | 0.7659738 |  0.0480 | remnant    | wood_saprotroph   | Ascomycota    | Dothideomycetes | Pleosporales     | Didymosphaeriaceae  | Paraphaeosphaeria | Paraphaeosphaeria_michotii |
| otu_881 | 0.9655172 | 0.50 | 0.6948083 |  0.0220 | remnant    | wood_saprotroph   | Ascomycota    | Eurotiomycetes  | Chaetothyriales  | Herpotrichiellaceae | Minimelanolocus   | Minimelanolocus_asiaticus  |
| otu_970 | 0.7943262 | 0.50 | 0.6302088 |  0.0475 | remnant    | wood_saprotroph   | Ascomycota    | Dothideomycetes | Minutisphaerales | Minutisphaeraceae   | Minutisphaera     | unidentified               |

Indicator species of wood saprotrophs

### Litter saprotrophs

#### Trends over time

``` r
guiltime("litter_saprotroph") 
```

    ## $bm_summary
    ## 
    ## Call:
    ## lm(formula = seq_sum ~ yr_since, data = d %>% filter(region == 
    ##     "BM"))
    ## 
    ## Residuals:
    ##       1       2       3       4       5       6       7 
    ##   270.5  2002.2  -105.1   573.2  -498.3 -1160.7 -1081.8 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)  
    ## (Intercept)  3067.40     785.52   3.905   0.0114 *
    ## yr_since      -78.87      52.84  -1.493   0.1958  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 1199 on 5 degrees of freedom
    ## Multiple R-squared:  0.3082, Adjusted R-squared:  0.1699 
    ## F-statistic: 2.228 on 1 and 5 DF,  p-value: 0.1958
    ## 
    ## 
    ## $fl_summary
    ## 
    ## Call:
    ## lm(formula = seq_sum ~ yr_since, data = d %>% filter(region == 
    ##     "FL"))
    ## 
    ## Residuals:
    ##       1       2       3       4       5       6 
    ##  -77.48   57.80   32.87 -515.39  585.61  -83.39 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)  
    ## (Intercept)  1120.70     322.66   3.473   0.0255 *
    ## yr_since      -20.93      11.89  -1.761   0.1531  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 395.6 on 4 degrees of freedom
    ## Multiple R-squared:  0.4367, Adjusted R-squared:  0.2958 
    ## F-statistic:   3.1 on 1 and 4 DF,  p-value: 0.1531
    ## 
    ## 
    ## $plot

<img src="microbial_guild_taxonomy_files/figure-gfm/lsap_guiltime-1.png" style="display: block; margin: auto;" />

#### Diversity

``` r
(lsap <- rerare(spe$its_raw, meta$its_raw, primary_lifestyle, "litter_saprotroph", sites))
```

    ## $seq_depth
    ## [1] 297
    ## 
    ## $zero_otu_num
    ## [1] 18
    ## 
    ## $rrfd
    ## # A tibble: 25 × 124
    ##    field_key otu_18 otu_105 otu_126 otu_133 otu_147 otu_151 otu_225 otu_226
    ##        <dbl>  <dbl>   <dbl>   <dbl>   <dbl>   <dbl>   <dbl>   <dbl>   <dbl>
    ##  1         1    114       0      14       0       0      97       4       0
    ##  2         2     85      88       1       0      23       0       2      41
    ##  3         3     59      28     163       0       0       0       9       0
    ##  4         4     25       0       0       0       0       0       0       0
    ##  5         5     27       0      53       0       0       0      14       0
    ##  6         6     82       0     100       0      27       0      46       0
    ##  7         7     61       0     129       0      29       0       0       0
    ##  8         8      7       0      15     200       5       0       1       0
    ##  9         9     60       0       0       0       0       0       0       0
    ## 10        10     62      24       0       0       0      64       0       0
    ## # … with 15 more rows, and 115 more variables: otu_265 <dbl>, otu_267 <dbl>,
    ## #   otu_272 <dbl>, otu_286 <dbl>, otu_302 <dbl>, otu_326 <dbl>, otu_358 <dbl>,
    ## #   otu_393 <dbl>, otu_414 <dbl>, otu_445 <dbl>, otu_457 <dbl>, otu_484 <dbl>,
    ## #   otu_503 <dbl>, otu_542 <dbl>, otu_551 <dbl>, otu_560 <dbl>, otu_574 <dbl>,
    ## #   otu_608 <dbl>, otu_618 <dbl>, otu_623 <dbl>, otu_653 <dbl>, otu_660 <dbl>,
    ## #   otu_698 <dbl>, otu_707 <dbl>, otu_729 <dbl>, otu_732 <dbl>, otu_761 <dbl>,
    ## #   otu_766 <dbl>, otu_789 <dbl>, otu_804 <dbl>, otu_827 <dbl>, …
    ## 
    ## $rrfd_speTaxa
    ## # A tibble: 433 × 14
    ##    field_key otu_num seq_abund phylum   class order family genus species prima…¹
    ##        <dbl> <chr>       <dbl> <chr>    <chr> <chr> <chr>  <chr> <chr>   <chr>  
    ##  1         1 otu_18        114 Ascomyc… Doth… Capn… Clado… Clad… <NA>    litter…
    ##  2         1 otu_126        14 Ascomyc… Sord… Sord… Chaet… Chae… <NA>    litter…
    ##  3         1 otu_151        97 Basidio… Agar… Agar… Entol… Clit… uniden… litter…
    ##  4         1 otu_225         4 Chytrid… Rhiz… Rhiz… Rhizo… Rhiz… Rhizop… litter…
    ##  5         1 otu_265         6 Chytrid… Rhiz… Rhiz… Rhizo… Rhiz… uniden… litter…
    ##  6         1 otu_272        14 Ascomyc… Doth… Pleo… Phaeo… Neos… <NA>    litter…
    ##  7         1 otu_286         6 Ascomyc… Doth… Pleo… Phaeo… Neos… <NA>    litter…
    ##  8         1 otu_326         4 Ascomyc… Doth… Pleo… Dicty… Dict… Dictyo… litter…
    ##  9         1 otu_414         2 Ascomyc… Euro… Chae… Cyphe… Cyph… <NA>    litter…
    ## 10         1 otu_457         2 Ascomyc… Euro… Chae… Cyphe… Cyph… uniden… litter…
    ## # … with 423 more rows, 4 more variables: field_name <chr>, region <chr>,
    ## #   field_type <ord>, yr_since <dbl>, and abbreviated variable name
    ## #   ¹​primary_lifestyle

Sequencing depth of 297, perhaps too rare to justify examination.

``` r
(lsap_div <- calc_diversity(lsap$rrfd))
```

<div data-pagedtable="false">

<script data-pagedtable-source type="application/json">
{"columns":[{"label":["field_key"],"name":[1],"type":["dbl"],"align":["right"]},{"label":["field_name"],"name":[2],"type":["chr"],"align":["left"]},{"label":["region"],"name":[3],"type":["chr"],"align":["left"]},{"label":["field_type"],"name":[4],"type":["ord"],"align":["right"]},{"label":["yr_since"],"name":[5],"type":["dbl"],"align":["right"]},{"label":["hill_index"],"name":[6],"type":["ord"],"align":["right"]},{"label":["value"],"name":[7],"type":["dbl"],"align":["right"]}],"data":[{"1":"1","2":"BBRP1","3":"BM","4":"restored","5":"16","6":"N0","7":"20.0000000"},{"1":"1","2":"BBRP1","3":"BM","4":"restored","5":"16","6":"N1","7":"6.1443222"},{"1":"1","2":"BBRP1","3":"BM","4":"restored","5":"16","6":"N2","7":"3.8032596"},{"1":"1","2":"BBRP1","3":"BM","4":"restored","5":"16","6":"E10","7":"0.3072161"},{"1":"1","2":"BBRP1","3":"BM","4":"restored","5":"16","6":"E20","7":"0.1901630"},{"1":"2","2":"ERRP1","3":"BM","4":"restored","5":"3","6":"N0","7":"19.0000000"},{"1":"2","2":"ERRP1","3":"BM","4":"restored","5":"3","6":"N1","7":"7.1484238"},{"1":"2","2":"ERRP1","3":"BM","4":"restored","5":"3","6":"N2","7":"4.9827148"},{"1":"2","2":"ERRP1","3":"BM","4":"restored","5":"3","6":"E10","7":"0.3762328"},{"1":"2","2":"ERRP1","3":"BM","4":"restored","5":"3","6":"E20","7":"0.2622481"},{"1":"3","2":"FGC1","3":"FG","4":"corn","5":"0","6":"N0","7":"13.0000000"},{"1":"3","2":"FGC1","3":"FG","4":"corn","5":"0","6":"N1","7":"4.4046910"},{"1":"3","2":"FGC1","3":"FG","4":"corn","5":"0","6":"N2","7":"2.8280273"},{"1":"3","2":"FGC1","3":"FG","4":"corn","5":"0","6":"E10","7":"0.3388224"},{"1":"3","2":"FGC1","3":"FG","4":"corn","5":"0","6":"E20","7":"0.2175406"},{"1":"4","2":"FGREM1","3":"FG","4":"remnant","5":"NA","6":"N0","7":"15.0000000"},{"1":"4","2":"FGREM1","3":"FG","4":"remnant","5":"NA","6":"N1","7":"9.4393295"},{"1":"4","2":"FGREM1","3":"FG","4":"remnant","5":"NA","6":"N2","7":"7.5800464"},{"1":"4","2":"FGREM1","3":"FG","4":"remnant","5":"NA","6":"E10","7":"0.6292886"},{"1":"4","2":"FGREM1","3":"FG","4":"remnant","5":"NA","6":"E20","7":"0.5053364"},{"1":"5","2":"FGRP1","3":"FG","4":"restored","5":"15","6":"N0","7":"14.0000000"},{"1":"5","2":"FGRP1","3":"FG","4":"restored","5":"15","6":"N1","7":"7.2971107"},{"1":"5","2":"FGRP1","3":"FG","4":"restored","5":"15","6":"N2","7":"5.0581455"},{"1":"5","2":"FGRP1","3":"FG","4":"restored","5":"15","6":"E10","7":"0.5212222"},{"1":"5","2":"FGRP1","3":"FG","4":"restored","5":"15","6":"E20","7":"0.3612961"},{"1":"6","2":"FLC1","3":"FL","4":"corn","5":"0","6":"N0","7":"11.0000000"},{"1":"6","2":"FLC1","3":"FL","4":"corn","5":"0","6":"N1","7":"5.4740331"},{"1":"6","2":"FLC1","3":"FL","4":"corn","5":"0","6":"N2","7":"4.3631103"},{"1":"6","2":"FLC1","3":"FL","4":"corn","5":"0","6":"E10","7":"0.4976394"},{"1":"6","2":"FLC1","3":"FL","4":"corn","5":"0","6":"E20","7":"0.3966464"},{"1":"7","2":"FLC2","3":"FL","4":"corn","5":"0","6":"N0","7":"11.0000000"},{"1":"7","2":"FLC2","3":"FL","4":"corn","5":"0","6":"N1","7":"4.9757134"},{"1":"7","2":"FLC2","3":"FL","4":"corn","5":"0","6":"N2","7":"3.6983355"},{"1":"7","2":"FLC2","3":"FL","4":"corn","5":"0","6":"E10","7":"0.4523376"},{"1":"7","2":"FLC2","3":"FL","4":"corn","5":"0","6":"E20","7":"0.3362123"},{"1":"8","2":"FLREM1","3":"FL","4":"remnant","5":"NA","6":"N0","7":"16.0000000"},{"1":"8","2":"FLREM1","3":"FL","4":"remnant","5":"NA","6":"N1","7":"4.0846780"},{"1":"8","2":"FLREM1","3":"FL","4":"remnant","5":"NA","6":"N2","7":"2.1470925"},{"1":"8","2":"FLREM1","3":"FL","4":"remnant","5":"NA","6":"E10","7":"0.2552924"},{"1":"8","2":"FLREM1","3":"FL","4":"remnant","5":"NA","6":"E20","7":"0.1341933"},{"1":"9","2":"FLRP1","3":"FL","4":"restored","5":"40","6":"N0","7":"15.0000000"},{"1":"9","2":"FLRP1","3":"FL","4":"restored","5":"40","6":"N1","7":"11.2493145"},{"1":"9","2":"FLRP1","3":"FL","4":"restored","5":"40","6":"N2","7":"9.3135889"},{"1":"9","2":"FLRP1","3":"FL","4":"restored","5":"40","6":"E10","7":"0.7499543"},{"1":"9","2":"FLRP1","3":"FL","4":"restored","5":"40","6":"E20","7":"0.6209059"},{"1":"10","2":"FLRP4","3":"FL","4":"restored","5":"36","6":"N0","7":"13.0000000"},{"1":"10","2":"FLRP4","3":"FL","4":"restored","5":"36","6":"N1","7":"8.6776545"},{"1":"10","2":"FLRP4","3":"FL","4":"restored","5":"36","6":"N2","7":"7.0437595"},{"1":"10","2":"FLRP4","3":"FL","4":"restored","5":"36","6":"E10","7":"0.6675119"},{"1":"10","2":"FLRP4","3":"FL","4":"restored","5":"36","6":"E20","7":"0.5418277"},{"1":"11","2":"FLRP5","3":"FL","4":"restored","5":"35","6":"N0","7":"14.0000000"},{"1":"11","2":"FLRP5","3":"FL","4":"restored","5":"35","6":"N1","7":"8.0493580"},{"1":"11","2":"FLRP5","3":"FL","4":"restored","5":"35","6":"N2","7":"6.5626813"},{"1":"11","2":"FLRP5","3":"FL","4":"restored","5":"35","6":"E10","7":"0.5749541"},{"1":"11","2":"FLRP5","3":"FL","4":"restored","5":"35","6":"E20","7":"0.4687630"},{"1":"12","2":"FLRSP1","3":"FL","4":"restored","5":"10","6":"N0","7":"16.0000000"},{"1":"12","2":"FLRSP1","3":"FL","4":"restored","5":"10","6":"N1","7":"7.1964198"},{"1":"12","2":"FLRSP1","3":"FL","4":"restored","5":"10","6":"N2","7":"4.2475562"},{"1":"12","2":"FLRSP1","3":"FL","4":"restored","5":"10","6":"E10","7":"0.4497762"},{"1":"12","2":"FLRSP1","3":"FL","4":"restored","5":"10","6":"E20","7":"0.2654723"},{"1":"13","2":"FLRSP2","3":"FL","4":"restored","5":"10","6":"N0","7":"23.0000000"},{"1":"13","2":"FLRSP2","3":"FL","4":"restored","5":"10","6":"N1","7":"7.4297147"},{"1":"13","2":"FLRSP2","3":"FL","4":"restored","5":"10","6":"N2","7":"3.8012928"},{"1":"13","2":"FLRSP2","3":"FL","4":"restored","5":"10","6":"E10","7":"0.3230311"},{"1":"13","2":"FLRSP2","3":"FL","4":"restored","5":"10","6":"E20","7":"0.1652736"},{"1":"14","2":"FLRSP3","3":"FL","4":"restored","5":"10","6":"N0","7":"15.0000000"},{"1":"14","2":"FLRSP3","3":"FL","4":"restored","5":"10","6":"N1","7":"6.3748300"},{"1":"14","2":"FLRSP3","3":"FL","4":"restored","5":"10","6":"N2","7":"4.2747274"},{"1":"14","2":"FLRSP3","3":"FL","4":"restored","5":"10","6":"E10","7":"0.4249887"},{"1":"14","2":"FLRSP3","3":"FL","4":"restored","5":"10","6":"E20","7":"0.2849818"},{"1":"15","2":"KORP1","3":"BM","4":"restored","5":"28","6":"N0","7":"16.0000000"},{"1":"15","2":"KORP1","3":"BM","4":"restored","5":"28","6":"N1","7":"5.8889673"},{"1":"15","2":"KORP1","3":"BM","4":"restored","5":"28","6":"N2","7":"3.1363200"},{"1":"15","2":"KORP1","3":"BM","4":"restored","5":"28","6":"E10","7":"0.3680605"},{"1":"15","2":"KORP1","3":"BM","4":"restored","5":"28","6":"E20","7":"0.1960200"},{"1":"16","2":"LPC1","3":"LP","4":"corn","5":"0","6":"N0","7":"22.0000000"},{"1":"16","2":"LPC1","3":"LP","4":"corn","5":"0","6":"N1","7":"8.4777019"},{"1":"16","2":"LPC1","3":"LP","4":"corn","5":"0","6":"N2","7":"5.2142224"},{"1":"16","2":"LPC1","3":"LP","4":"corn","5":"0","6":"E10","7":"0.3853501"},{"1":"16","2":"LPC1","3":"LP","4":"corn","5":"0","6":"E20","7":"0.2370101"},{"1":"17","2":"LPREM1","3":"LP","4":"remnant","5":"NA","6":"N0","7":"24.0000000"},{"1":"17","2":"LPREM1","3":"LP","4":"remnant","5":"NA","6":"N1","7":"13.2489165"},{"1":"17","2":"LPREM1","3":"LP","4":"remnant","5":"NA","6":"N2","7":"9.1626675"},{"1":"17","2":"LPREM1","3":"LP","4":"remnant","5":"NA","6":"E10","7":"0.5520382"},{"1":"17","2":"LPREM1","3":"LP","4":"remnant","5":"NA","6":"E20","7":"0.3817778"},{"1":"18","2":"LPRP1","3":"LP","4":"restored","5":"4","6":"N0","7":"21.0000000"},{"1":"18","2":"LPRP1","3":"LP","4":"restored","5":"4","6":"N1","7":"8.5202849"},{"1":"18","2":"LPRP1","3":"LP","4":"restored","5":"4","6":"N2","7":"4.3258791"},{"1":"18","2":"LPRP1","3":"LP","4":"restored","5":"4","6":"E10","7":"0.4057279"},{"1":"18","2":"LPRP1","3":"LP","4":"restored","5":"4","6":"E20","7":"0.2059942"},{"1":"19","2":"LPRP2","3":"LP","4":"restored","5":"4","6":"N0","7":"15.0000000"},{"1":"19","2":"LPRP2","3":"LP","4":"restored","5":"4","6":"N1","7":"8.5366371"},{"1":"19","2":"LPRP2","3":"LP","4":"restored","5":"4","6":"N2","7":"7.1732130"},{"1":"19","2":"LPRP2","3":"LP","4":"restored","5":"4","6":"E10","7":"0.5691091"},{"1":"19","2":"LPRP2","3":"LP","4":"restored","5":"4","6":"E20","7":"0.4782142"},{"1":"20","2":"MBREM1","3":"BM","4":"remnant","5":"NA","6":"N0","7":"17.0000000"},{"1":"20","2":"MBREM1","3":"BM","4":"remnant","5":"NA","6":"N1","7":"9.8107852"},{"1":"20","2":"MBREM1","3":"BM","4":"remnant","5":"NA","6":"N2","7":"7.4343869"},{"1":"20","2":"MBREM1","3":"BM","4":"remnant","5":"NA","6":"E10","7":"0.5771050"},{"1":"20","2":"MBREM1","3":"BM","4":"remnant","5":"NA","6":"E20","7":"0.4373169"},{"1":"21","2":"MBRP1","3":"BM","4":"restored","5":"18","6":"N0","7":"29.0000000"},{"1":"21","2":"MBRP1","3":"BM","4":"restored","5":"18","6":"N1","7":"10.6234450"},{"1":"21","2":"MBRP1","3":"BM","4":"restored","5":"18","6":"N2","7":"6.9625858"},{"1":"21","2":"MBRP1","3":"BM","4":"restored","5":"18","6":"E10","7":"0.3663257"},{"1":"21","2":"MBRP1","3":"BM","4":"restored","5":"18","6":"E20","7":"0.2400892"},{"1":"22","2":"MHRP1","3":"BM","4":"restored","5":"7","6":"N0","7":"20.0000000"},{"1":"22","2":"MHRP1","3":"BM","4":"restored","5":"7","6":"N1","7":"10.5672372"},{"1":"22","2":"MHRP1","3":"BM","4":"restored","5":"7","6":"N2","7":"7.6114419"},{"1":"22","2":"MHRP1","3":"BM","4":"restored","5":"7","6":"E10","7":"0.5283619"},{"1":"22","2":"MHRP1","3":"BM","4":"restored","5":"7","6":"E20","7":"0.3805721"},{"1":"23","2":"MHRP2","3":"BM","4":"restored","5":"2","6":"N0","7":"15.0000000"},{"1":"23","2":"MHRP2","3":"BM","4":"restored","5":"2","6":"N1","7":"8.1645061"},{"1":"23","2":"MHRP2","3":"BM","4":"restored","5":"2","6":"N2","7":"6.0206812"},{"1":"23","2":"MHRP2","3":"BM","4":"restored","5":"2","6":"E10","7":"0.5443004"},{"1":"23","2":"MHRP2","3":"BM","4":"restored","5":"2","6":"E20","7":"0.4013787"},{"1":"24","2":"PHC1","3":"BM","4":"corn","5":"0","6":"N0","7":"19.0000000"},{"1":"24","2":"PHC1","3":"BM","4":"corn","5":"0","6":"N1","7":"11.0081907"},{"1":"24","2":"PHC1","3":"BM","4":"corn","5":"0","6":"N2","7":"9.2510750"},{"1":"24","2":"PHC1","3":"BM","4":"corn","5":"0","6":"E10","7":"0.5793785"},{"1":"24","2":"PHC1","3":"BM","4":"corn","5":"0","6":"E20","7":"0.4868987"},{"1":"25","2":"PHRP1","3":"BM","4":"restored","5":"11","6":"N0","7":"20.0000000"},{"1":"25","2":"PHRP1","3":"BM","4":"restored","5":"11","6":"N1","7":"7.7905937"},{"1":"25","2":"PHRP1","3":"BM","4":"restored","5":"11","6":"N2","7":"4.2094488"},{"1":"25","2":"PHRP1","3":"BM","4":"restored","5":"11","6":"E10","7":"0.3895297"},{"1":"25","2":"PHRP1","3":"BM","4":"restored","5":"11","6":"E20","7":"0.2104724"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>

</div>

``` r
(lsap_comp <- gudicom(lsap_div, lsap$rrfd_speTaxa, "litter_saprotroph"))
```

    ## $Hills_field_type

<img src="microbial_guild_taxonomy_files/figure-gfm/lsap_composition-1.png" style="display: block; margin: auto;" />

    ## 
    ## $Hills_yrs_since_restoration

<img src="microbial_guild_taxonomy_files/figure-gfm/lsap_composition-2.png" style="display: block; margin: auto;" />

    ## 
    ## $Composition

<img src="microbial_guild_taxonomy_files/figure-gfm/lsap_composition-3.png" style="display: block; margin: auto;" />

<div data-pagedtable="false">

<script data-pagedtable-source type="application/json">
{"columns":[{"label":["field_type"],"name":[1],"type":["ord"],"align":["right"]},{"label":["order"],"name":[2],"type":["chr"],"align":["left"]},{"label":["seq_comp"],"name":[3],"type":["dbl"],"align":["right"]}],"data":[{"1":"corn","2":"Agaricales","3":"6.201055"},{"1":"corn","2":"Capnodiales","3":"18.922806"},{"1":"corn","2":"Chaetosphaeriales","3":"15.982100"},{"1":"corn","2":"Chaetothyriales","3":"4.283203"},{"1":"corn","2":"Other (OTU<2%)","3":"2.796868"},{"1":"corn","2":"Pezizales","3":"3.995525"},{"1":"corn","2":"Pleosporales","3":"7.927122"},{"1":"corn","2":"Rhizophlyctidales","3":"11.315327"},{"1":"corn","2":"Sordariales","3":"28.575995"},{"1":"restored","2":"Agaricales","3":"15.743820"},{"1":"restored","2":"Capnodiales","3":"22.018195"},{"1":"restored","2":"Chaetosphaeriales","3":"10.023328"},{"1":"restored","2":"Chaetothyriales","3":"6.908759"},{"1":"restored","2":"Glomerellales","3":"2.407878"},{"1":"restored","2":"Helotiales","3":"6.459597"},{"1":"restored","2":"Hypocreales","3":"4.227017"},{"1":"restored","2":"Other (OTU<2%)","3":"2.917237"},{"1":"restored","2":"Pezizales","3":"2.273172"},{"1":"restored","2":"Pleosporales","3":"8.038609"},{"1":"restored","2":"Rhizophlyctidales","3":"7.316246"},{"1":"restored","2":"Sebacinales","3":"4.013131"},{"1":"restored","2":"Sordariales","3":"7.653012"},{"1":"remnant","2":"Agaricales","3":"9.440698"},{"1":"remnant","2":"Cantharellales","3":"9.044030"},{"1":"remnant","2":"Capnodiales","3":"10.412535"},{"1":"remnant","2":"Chaetosphaeriales","3":"5.632685"},{"1":"remnant","2":"Chaetothyriales","3":"9.698532"},{"1":"remnant","2":"Helotiales","3":"9.758033"},{"1":"remnant","2":"Hypocreales","3":"16.600555"},{"1":"remnant","2":"Other (OTU<2%)","3":"2.538675"},{"1":"remnant","2":"Pleosporales","3":"7.973027"},{"1":"remnant","2":"Rhizophlyctidales","3":"3.272511"},{"1":"remnant","2":"Sebacinales","3":"4.998017"},{"1":"remnant","2":"Sordariales","3":"3.014677"},{"1":"remnant","2":"Xylariales","3":"7.616025"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>

</div>

With no litter in cornfields, it’s perhaps not surprising to see
increasing trends across field types with this guild. Trends over time
aren’t convincing, except possibly in Fermi.

#### Indicators

``` r
lsap_inspan <- 
    lsap$rrfd %>% 
    left_join(sites, by = join_by(field_key)) %>% 
    inspan(., 1999, meta$its_raw)
```

``` r
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
```

| field_type | n_otu |  stat_avg |   stat_sd |
|:-----------|------:|----------:|----------:|
| corn       |     3 | 0.8000703 | 0.0741808 |
| remnant    |     1 | 0.6813851 |        NA |

Indicator species stats: litter saprotrophs

``` r
lsap_inspan %>% 
    mutate(field_type = factor(
        field_type,
        ordered = TRUE,
        levels = c("corn", "restored", "remnant")
    )) %>%
    arrange(field_type, -stat) %>% 
    kable(format = "pandoc", caption = "Indicator species of litter saprotrophs")
```

| otu_num  |         A |   B |      stat | p.value | field_type | primary_lifestyle | phylum          | class                 | order             | family             | genus         | species                |
|:---------|----------:|----:|----------:|--------:|:-----------|:------------------|:----------------|:----------------------|:------------------|:-------------------|:--------------|:-----------------------|
| otu_126  | 0.7845069 | 1.0 | 0.8857239 |  0.0055 | corn       | litter_saprotroph | Ascomycota      | Sordariomycetes       | Sordariales       | Chaetomiaceae      | Chaetomium    | NA                     |
| otu_358  | 0.9572650 | 0.6 | 0.7578647 |  0.0240 | corn       | litter_saprotroph | Ascomycota      | Eurotiomycetes        | Chaetothyriales   | Cyphellophoraceae  | Cyphellophora | Cyphellophora_suttonii |
| otu_1009 | 0.9541284 | 0.6 | 0.7566221 |  0.0140 | corn       | litter_saprotroph | Ascomycota      | Pezizomycetes         | Pezizales         | Pyronemataceae     | Cheilymenia   | Cheilymenia_stercorea  |
| otu_1302 | 0.9285714 | 0.5 | 0.6813851 |  0.0370 | remnant    | litter_saprotroph | Chytridiomycota | Rhizophlyctidomycetes | Rhizophlyctidales | Rhizophlyctidaceae | Rhizophlyctis | unidentified           |

Indicator species of litter saprotrophs

## AMF OTUs

Function output is verbose but retained as explained previously.

``` r
# amf_otu_summary <- amf_tax(spe_meta$amf_otu, "otu")
```

Claroideoglomeraceae differs across field types with a likelihood ratio
test result p\<0.01. Tukey’s post-hoc test with Holm correction
performed, letters on the figure show differences.

``` r
# amf_otu_summary %>% 
#     filter(family == "Claroideoglomeraceae") %>% 
#     ggplot(aes(x = field_type, y = seq_sum)) +
#     geom_boxplot(varwidth = TRUE, fill = "gray90", outlier.shape = NA) +
#     geom_beeswarm(aes(fill = region), shape = 21, size = 2, dodge.width = 0.2) +
#     annotate("text", label = c("a", "b", "ab"), x = c(1,2,3), y = rep(350, 3)) +
#     labs(x = "", y = "Sequence abundance", title = "AMF variation in field types, 97% OTU",
#          caption = "Likelihood ratio test p<0.01, Tukey's post-hoc with Holm correction at p<0.05") +
#     scale_fill_discrete_qualitative(palette = "Dark3") +
#     theme_classic()
```

Gigasporaceae increased with time since restoration by a simple linear
regression, $R^2_{adj}$ = 0.81, p \< 0.01

``` r
# amf_otu_summary %>% 
#     filter(field_type == "restored", region == "BM", family == "Gigasporaceae") %>% 
#     ggplot(aes(x = yr_since, y = seq_sum)) +
#     geom_smooth(method = "lm", linewidth = 0.4, se = FALSE) +
#     geom_point(size = 2, shape = 21, fill = "gray60") +
#     labs(x = "Years since restoration", 
#          y = "Sequence abundance", 
#          title = "Gigasporaceae abundance since restoration, 97% OTU",
#          caption = "R2Adj = 0.81, p<0.01") +
#     theme_classic()
```

# Conclusions: taxa and guilds

Little variation exists here for ITS or AMF sequences among field types,
although classes of fungi identified through ITS sequences remain to be
closely examined. It’s striking that plant pathogens decline as
restorations age while the AMF family *Gigasporaceae* increases, but
this contrast was not found in any other group of AMF and the
*Gigasporaceae* aren’t particularly abundant to begin with.
