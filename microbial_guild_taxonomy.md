Microbial data: microbial guilds and taxonomy
================
Beau Larkin

Last updated: 07 March, 2023

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
  - <a href="#amf" id="toc-amf">AMF</a>
    - <a href="#individual-families" id="toc-individual-families">Individual
      families</a>
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

    ## 
    ## [[3]]

<img src="microbial_guild_taxonomy_files/figure-gfm/its_tax_trophic_otu-2.png" style="display: block; margin: auto;" />

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
    ## remnant - corn == 0       1968.2      915.1   2.151  0.07681 . 
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
    ## (Intercept)   -58.87        NaN     NaN      NaN
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
    ## remnant - corn == 0      -1732.5      642.9  -2.695   0.0187 *
    ## remnant - restored == 0  -1218.2      538.5  -2.262   0.0598 .
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## (Adjusted p values reported -- single-step method)
    ## 
    ##     corn restored  remnant 
    ##      "b"     "ab"      "a" 
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
| corn       |    90 | 0.8255659 | 0.0998572 |
| restored   |    10 | 0.8157168 | 0.0351250 |
| remnant    |    50 | 0.7379027 | 0.0810312 |

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

| otu_num  |         A |      B |      stat | p.value | field_type | primary_lifestyle      | phylum            | class              | order                              | family                             | genus            | species                 |
|:---------|----------:|-------:|----------:|--------:|:-----------|:-----------------------|:------------------|:-------------------|:-----------------------------------|:-----------------------------------|:-----------------|:------------------------|
| otu_537  | 1.0000000 | 1.0000 | 1.0000000 |  0.0005 | corn       | soil_saprotroph        | Basidiomycota     | Agaricomycetes     | Agaricales                         | Bolbitiaceae                       | Conocybe         | Conocybe_apala          |
| otu_204  | 0.9937578 | 1.0000 | 0.9968740 |  0.0005 | corn       | NA                     | Mortierellomycota | Mortierellomycetes | Mortierellales                     | Mortierellaceae                    | NA               | NA                      |
| otu_172  | 0.9772048 | 1.0000 | 0.9885367 |  0.0015 | corn       | plant_pathogen         | Ascomycota        | Dothideomycetes    | Pleosporales                       | Corynesporascaceae                 | Corynespora      | Corynespora_cassiicola  |
| otu_188  | 0.9759492 | 1.0000 | 0.9879014 |  0.0005 | corn       | NA                     | NA                | NA                 | NA                                 | NA                                 | NA               | NA                      |
| otu_9    | 0.9753667 | 1.0000 | 0.9876066 |  0.0025 | corn       | soil_saprotroph        | Basidiomycota     | Tremellomycetes    | Cystofilobasidiales                | Mrakiaceae                         | Tausonia         | Tausonia_pullulans      |
| otu_200  | 0.9724757 | 1.0000 | 0.9861418 |  0.0005 | corn       | plant_pathogen         | Ascomycota        | Dothideomycetes    | Pleosporales                       | Phaeosphaeriaceae                  | Ophiosphaerella  | unidentified            |
| otu_59   | 0.9602305 | 1.0000 | 0.9799135 |  0.0005 | corn       | soil_saprotroph        | Mortierellomycota | Mortierellomycetes | Mortierellales                     | Mortierellaceae                    | Mortierella      | NA                      |
| otu_694  | 0.9400850 | 1.0000 | 0.9695798 |  0.0005 | corn       | NA                     | NA                | NA                 | NA                                 | NA                                 | NA               | NA                      |
| otu_553  | 0.9378783 | 1.0000 | 0.9684412 |  0.0015 | corn       | plant_pathogen         | Ascomycota        | Sordariomycetes    | Magnaporthales                     | Magnaporthaceae                    | Gaeumannomyces   | NA                      |
| otu_364  | 0.9318632 | 1.0000 | 0.9653306 |  0.0005 | corn       | NA                     | Ascomycota        | Sordariomycetes    | Sordariales                        | Lasiosphaeriaceae                  | Cladorrhinum     | NA                      |
| otu_332  | 0.9219288 | 0.8125 | 0.8654867 |  0.0380 | restored   | plant_pathogen         | Ascomycota        | Sordariomycetes    | Glomerellales                      | Plectosphaerellaceae               | Plectosphaerella | NA                      |
| otu_177  | 0.9809886 | 0.7500 | 0.8577537 |  0.0225 | restored   | NA                     | Ascomycota        | Dothideomycetes    | Pleosporales                       | NA                                 | NA               | NA                      |
| otu_817  | 1.0000000 | 0.6875 | 0.8291562 |  0.0200 | restored   | NA                     | Ascomycota        | NA                 | NA                                 | NA                                 | NA               | NA                      |
| otu_461  | 0.8351648 | 0.8125 | 0.8237545 |  0.0235 | restored   | NA                     | Ascomycota        | Dothideomycetes    | Pleosporales                       | Phaeosphaeriaceae                  | NA               | NA                      |
| otu_35   | 0.7234228 | 0.9375 | 0.8235344 |  0.0410 | restored   | animal_parasite        | Ascomycota        | Sordariomycetes    | Hypocreales                        | Clavicipitaceae                    | Metarhizium      | NA                      |
| otu_193  | 0.8307978 | 0.8125 | 0.8215979 |  0.0485 | restored   | NA                     | Basidiomycota     | Agaricomycetes     | Sebacinales                        | unidentified                       | unidentified     | unidentified            |
| otu_107  | 0.8061297 | 0.8125 | 0.8093086 |  0.0320 | restored   | NA                     | Ascomycota        | Dothideomycetes    | Pleosporales                       | NA                                 | NA               | NA                      |
| otu_114  | 0.6963432 | 0.9375 | 0.8079739 |  0.0060 | restored   | soil_saprotroph        | Mortierellomycota | Mortierellomycetes | Mortierellales                     | Mortierellaceae                    | Mortierella      | unidentified            |
| otu_33   | 0.5843320 | 1.0000 | 0.7644161 |  0.0425 | restored   | plant_pathogen         | Ascomycota        | Sordariomycetes    | Hypocreales                        | Nectriaceae                        | Fusarium         | NA                      |
| otu_10   | 0.5687968 | 1.0000 | 0.7541862 |  0.0135 | restored   | NA                     | Ascomycota        | NA                 | NA                                 | NA                                 | NA               | NA                      |
| otu_772  | 0.9272420 | 1.0000 | 0.9629340 |  0.0035 | remnant    | NA                     | Ascomycota        | Sordariomycetes    | NA                                 | NA                                 | NA               | NA                      |
| otu_629  | 0.9159892 | 1.0000 | 0.9570732 |  0.0015 | remnant    | NA                     | Ascomycota        | Leotiomycetes      | Helotiales                         | Hyaloscyphaceae                    | Microscypha      | unidentified            |
| otu_159  | 0.8185686 | 1.0000 | 0.9047478 |  0.0025 | remnant    | NA                     | Ascomycota        | Sordariomycetes    | Sordariomycetes_ord_Incertae_sedis | Sordariomycetes_fam_Incertae_sedis | Pleurophragmium  | unidentified            |
| otu_135  | 0.7768230 | 1.0000 | 0.8813757 |  0.0035 | remnant    | plant_pathogen         | Ascomycota        | Sordariomycetes    | Hypocreales                        | Nectriaceae                        | Ilyonectria      | NA                      |
| otu_854  | 1.0000000 | 0.7500 | 0.8660254 |  0.0020 | remnant    | NA                     | Ascomycota        | NA                 | NA                                 | NA                                 | NA               | NA                      |
| otu_1740 | 1.0000000 | 0.7500 | 0.8660254 |  0.0020 | remnant    | NA                     | Glomeromycota     | Glomeromycetes     | Glomerales                         | Glomeraceae                        | NA               | NA                      |
| otu_1098 | 0.9716841 | 0.7500 | 0.8536762 |  0.0055 | remnant    | NA                     | NA                | NA                 | NA                                 | NA                                 | NA               | NA                      |
| otu_1468 | 0.9332261 | 0.7500 | 0.8366119 |  0.0035 | remnant    | NA                     | Ascomycota        | Sordariomycetes    | NA                                 | NA                                 | NA               | NA                      |
| otu_140  | 0.9276552 | 0.7500 | 0.8341111 |  0.0320 | remnant    | soil_saprotroph        | Ascomycota        | Sordariomycetes    | Hypocreales                        | Stachybotryaceae                   | Striaticonidium  | Striaticonidium_cinctum |
| otu_369  | 0.6807314 | 1.0000 | 0.8250645 |  0.0150 | remnant    | unspecified_saprotroph | Ascomycota        | Sordariomycetes    | Hypocreales                        | Bionectriaceae                     | Gliomastix       | Gliomastix_roseogrisea  |

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
    ## [1] 3
    ## 
    ## $rrfd
    ## # A tibble: 25 × 248
    ##    field_key otu_2 otu_9 otu_14 otu_27 otu_37 otu_41 otu_47 otu_49 otu_55 otu_59
    ##        <dbl> <dbl> <dbl>  <dbl>  <dbl>  <dbl>  <dbl>  <dbl>  <dbl>  <dbl>  <dbl>
    ##  1         1   741     0     93      2    184     16    877      0    613      0
    ##  2         2  1452     0    453      2    146     64      0      0      0      0
    ##  3         3     7    30    409      0    250    461      0      0      0    285
    ##  4         4     0     0      0      0    279    220     22      0     19      0
    ##  5         5   116    47      0      0    122     33      0      0      0    214
    ##  6         6   145  1270    397    138     33    242      0      0      0     67
    ##  7         7   339  2265     48     35     93     49      0      0      0    237
    ##  8         8   601    38    373     71    219     62     19      0      0     25
    ##  9         9   716     0    293      3    115     43    949      0    739      7
    ## 10        10   576     0    330      0     54     71    372      0   1200      0
    ## # … with 15 more rows, and 237 more variables: otu_61 <dbl>, otu_66 <dbl>,
    ## #   otu_70 <dbl>, otu_75 <dbl>, otu_79 <dbl>, otu_88 <dbl>, otu_89 <dbl>,
    ## #   otu_100 <dbl>, otu_102 <dbl>, otu_106 <dbl>, otu_114 <dbl>, otu_132 <dbl>,
    ## #   otu_134 <dbl>, otu_140 <dbl>, otu_144 <dbl>, otu_154 <dbl>, otu_168 <dbl>,
    ## #   otu_170 <dbl>, otu_180 <dbl>, otu_182 <dbl>, otu_192 <dbl>, otu_195 <dbl>,
    ## #   otu_209 <dbl>, otu_216 <dbl>, otu_221 <dbl>, otu_223 <dbl>, otu_224 <dbl>,
    ## #   otu_234 <dbl>, otu_276 <dbl>, otu_280 <dbl>, otu_283 <dbl>, …
    ## 
    ## $rrfd_speTaxa
    ## # A tibble: 939 × 14
    ##    field_key otu_num seq_abund phylum   class order family genus species prima…¹
    ##        <dbl> <chr>       <dbl> <chr>    <chr> <chr> <chr>  <chr> <chr>   <chr>  
    ##  1         1 otu_2         741 Mortier… Mort… Mort… Morti… Mort… Mortie… soil_s…
    ##  2         1 otu_14         93 Mortier… Mort… Mort… Morti… Mort… <NA>    soil_s…
    ##  3         1 otu_27          2 Basidio… Trem… Filo… Pisku… Soli… <NA>    soil_s…
    ##  4         1 otu_37        184 Mortier… Mort… Mort… Morti… Mort… <NA>    soil_s…
    ##  5         1 otu_41         16 Mortier… Mort… Mort… Morti… Mort… Mortie… soil_s…
    ##  6         1 otu_47        877 Ascomyc… Geog… Geog… Geogl… Geog… uniden… soil_s…
    ##  7         1 otu_55        613 Basidio… Agar… Agar… Hygro… Hygr… <NA>    soil_s…
    ##  8         1 otu_75          4 Basidio… Trem… Filo… Pisku… Soli… Solico… soil_s…
    ##  9         1 otu_79         51 Mortier… Mort… Mort… Morti… Mort… <NA>    soil_s…
    ## 10         1 otu_114        66 Mortier… Mort… Mort… Morti… Mort… uniden… soil_s…
    ## # … with 929 more rows, 4 more variables: field_name <chr>, region <chr>,
    ## #   field_type <ord>, yr_since <dbl>, and abbreviated variable name
    ## #   ¹​primary_lifestyle

``` r
ssap_div <- calc_diversity(ssap$rrfd)
```

Diversity measures are stored in this data frame for further use…

``` r
(ssap_comp <- gudicom(ssap_div, ssap$rrfd_speTaxa, "soil_saprotroph"))
```

<img src="microbial_guild_taxonomy_files/figure-gfm/ssap_composition-1.png" style="display: block; margin: auto;" /><img src="microbial_guild_taxonomy_files/figure-gfm/ssap_composition-2.png" style="display: block; margin: auto;" /><img src="microbial_guild_taxonomy_files/figure-gfm/ssap_composition-3.png" style="display: block; margin: auto;" />

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
| corn       |     7 | 0.8771517 | 0.1282574 |
| restored   |     2 | 0.7819169 | 0.0063084 |
| remnant    |     2 | 0.7507356 | 0.0617005 |

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
| otu_9    | 0.9454921 | 1.0000 | 0.9723642 |  0.0035 | corn       | soil_saprotroph   | Basidiomycota     | Tremellomycetes    | Cystofilobasidiales | Mrakiaceae       | Tausonia           | Tausonia_pullulans      |
| otu_59   | 0.9348759 | 1.0000 | 0.9668898 |  0.0010 | corn       | soil_saprotroph   | Mortierellomycota | Mortierellomycetes | Mortierellales      | Mortierellaceae  | Mortierella        | NA                      |
| otu_134  | 0.8376068 | 1.0000 | 0.9152086 |  0.0045 | corn       | soil_saprotroph   | Mortierellomycota | Mortierellomycetes | Mortierellales      | Mortierellaceae  | Mortierella        | NA                      |
| otu_61   | 0.8642313 | 0.8000 | 0.8314957 |  0.0495 | corn       | soil_saprotroph   | Basidiomycota     | Agaricomycetes     | Phallales           | Phallaceae       | Phallus            | Phallus_rugulosus       |
| otu_41   | 0.6751055 | 1.0000 | 0.8216480 |  0.0200 | corn       | soil_saprotroph   | Mortierellomycota | Mortierellomycetes | Mortierellales      | Mortierellaceae  | Mortierella        | Mortierella_minutissima |
| otu_1341 | 1.0000000 | 0.4000 | 0.6324555 |  0.0420 | corn       | soil_saprotroph   | Basidiomycota     | Agaricomycetes     | Agaricales          | Entolomataceae   | Entoloma           | unidentified            |
| otu_114  | 0.6596157 | 0.9375 | 0.7863776 |  0.0180 | restored   | soil_saprotroph   | Mortierellomycota | Mortierellomycetes | Mortierellales      | Mortierellaceae  | Mortierella        | unidentified            |
| otu_2    | 0.6044382 | 1.0000 | 0.7774562 |  0.0315 | restored   | soil_saprotroph   | Mortierellomycota | Mortierellomycetes | Mortierellales      | Mortierellaceae  | Mortierella        | Mortierella_exigua      |
| otu_140  | 0.8413532 | 0.7500 | 0.7943645 |  0.0415 | remnant    | soil_saprotroph   | Ascomycota        | Sordariomycetes    | Hypocreales         | Stachybotryaceae | Striaticonidium    | Striaticonidium_cinctum |
| otu_2138 | 1.0000000 | 0.5000 | 0.7071068 |  0.0195 | remnant    | soil_saprotroph   | Ascomycota        | Leotiomycetes      | Thelebolales        | Pseudeurotiaceae | Gymnostellatospora | NA                      |

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
    ## [1] 8
    ## 
    ## $rrfd
    ## # A tibble: 25 × 153
    ##    field_key otu_1 otu_3 otu_7 otu_13 otu_16 otu_21 otu_23 otu_28 otu_33 otu_43
    ##        <dbl> <dbl> <dbl> <dbl>  <dbl>  <dbl>  <dbl>  <dbl>  <dbl>  <dbl>  <dbl>
    ##  1         1   484    35  1585      7     48      0    144      0     33     14
    ##  2         2   805   244   119    207     85      3    106      8    207     40
    ##  3         3   128   116    62    725    208    775     78     48     36     56
    ##  4         4   695     1   667      0      0      0     55      0     44      0
    ##  5         5   118   110  1180     53    107     33     37     29     45     67
    ##  6         6   489    98   496      4    271    517    145     43     17    104
    ##  7         7   444   445    35    241    295    173    102    202     85     79
    ##  8         8   718   288   544     68     45     57     30     49    189     30
    ##  9         9   783   428   372    102     22      0    293      0    105      5
    ## 10        10   507   579   681     60     12      1     37      0    149      0
    ## # … with 15 more rows, and 142 more variables: otu_53 <dbl>, otu_58 <dbl>,
    ## #   otu_65 <dbl>, otu_68 <dbl>, otu_87 <dbl>, otu_99 <dbl>, otu_135 <dbl>,
    ## #   otu_137 <dbl>, otu_153 <dbl>, otu_172 <dbl>, otu_179 <dbl>, otu_200 <dbl>,
    ## #   otu_212 <dbl>, otu_279 <dbl>, otu_285 <dbl>, otu_289 <dbl>, otu_294 <dbl>,
    ## #   otu_304 <dbl>, otu_315 <dbl>, otu_319 <dbl>, otu_325 <dbl>, otu_332 <dbl>,
    ## #   otu_383 <dbl>, otu_391 <dbl>, otu_408 <dbl>, otu_416 <dbl>, otu_425 <dbl>,
    ## #   otu_432 <dbl>, otu_504 <dbl>, otu_511 <dbl>, otu_521 <dbl>, …
    ## 
    ## $rrfd_speTaxa
    ## # A tibble: 843 × 14
    ##    field_key otu_num seq_abund phylum   class order family genus species prima…¹
    ##        <dbl> <chr>       <dbl> <chr>    <chr> <chr> <chr>  <chr> <chr>   <chr>  
    ##  1         1 otu_1         484 Ascomyc… Sord… Hypo… Nectr… Fusa… Fusari… plant_…
    ##  2         1 otu_3          35 Ascomyc… Sord… Glom… Plect… Gibe… <NA>    plant_…
    ##  3         1 otu_7        1585 Ascomyc… Doth… Pleo… Peric… Peri… <NA>    plant_…
    ##  4         1 otu_13          7 Ascomyc… Sord… Glom… Plect… Plec… Plecto… plant_…
    ##  5         1 otu_16         48 Ascomyc… Sord… Hypo… Nectr… Nect… Nectri… plant_…
    ##  6         1 otu_23        144 Ascomyc… Doth… Pleo… Pleos… Alte… <NA>    plant_…
    ##  7         1 otu_33         33 Ascomyc… Sord… Hypo… Nectr… Fusa… <NA>    plant_…
    ##  8         1 otu_43         14 Ascomyc… Sord… Hypo… Nectr… Fusa… Fusari… plant_…
    ##  9         1 otu_58         32 Ascomyc… Doth… Pleo… Phaeo… Para… <NA>    plant_…
    ## 10         1 otu_65          7 Ascomyc… Sord… Hypo… Nectr… Gibb… Gibber… plant_…
    ## # … with 833 more rows, 4 more variables: field_name <chr>, region <chr>,
    ## #   field_type <ord>, yr_since <dbl>, and abbreviated variable name
    ## #   ¹​primary_lifestyle

``` r
ppat_div <- calc_diversity(ppat$rrfd)
```

``` r
(ppat_comp <- gudicom(ppat_div, ppat$rrfd_speTaxa, "plant_pathogen", other_threshold = 1))
```

<img src="microbial_guild_taxonomy_files/figure-gfm/ppat_composition-1.png" style="display: block; margin: auto;" /><img src="microbial_guild_taxonomy_files/figure-gfm/ppat_composition-2.png" style="display: block; margin: auto;" /><img src="microbial_guild_taxonomy_files/figure-gfm/ppat_composition-3.png" style="display: block; margin: auto;" />

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
| corn       |    12 | 0.8619885 | 0.0944745 |
| restored   |     1 | 0.8370619 |        NA |
| remnant    |     5 | 0.7613326 | 0.0863075 |

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

| otu_num  |         A |    B |      stat | p.value | field_type | primary_lifestyle | phylum        | class             | order          | family               | genus            | species                     |
|:---------|----------:|-----:|----------:|--------:|:-----------|:------------------|:--------------|:------------------|:---------------|:---------------------|:-----------------|:----------------------------|
| otu_172  | 0.9768875 | 1.00 | 0.9883762 |  0.0005 | corn       | plant_pathogen    | Ascomycota    | Dothideomycetes   | Pleosporales   | Corynesporascaceae   | Corynespora      | Corynespora_cassiicola      |
| otu_200  | 0.9484198 | 1.00 | 0.9738685 |  0.0010 | corn       | plant_pathogen    | Ascomycota    | Dothideomycetes   | Pleosporales   | Phaeosphaeriaceae    | Ophiosphaerella  | unidentified                |
| otu_21   | 0.9165629 | 1.00 | 0.9573729 |  0.0005 | corn       | plant_pathogen    | Ascomycota    | Dothideomycetes   | Pleosporales   | Phaeosphaeriaceae    | Setophoma        | Setophoma_terrestris        |
| otu_553  | 0.9086682 | 1.00 | 0.9532409 |  0.0040 | corn       | plant_pathogen    | Ascomycota    | Sordariomycetes   | Magnaporthales | Magnaporthaceae      | Gaeumannomyces   | NA                          |
| otu_432  | 1.0000000 | 0.80 | 0.8944272 |  0.0005 | corn       | plant_pathogen    | Ascomycota    | Sordariomycetes   | Glomerellales  | Glomerellaceae       | Colletotrichum   | NA                          |
| otu_391  | 0.7425302 | 1.00 | 0.8617019 |  0.0105 | corn       | plant_pathogen    | Ascomycota    | Dothideomycetes   | Pleosporales   | Torulaceae           | Dendryphion      | NA                          |
| otu_13   | 0.7305373 | 1.00 | 0.8547147 |  0.0050 | corn       | plant_pathogen    | Ascomycota    | Sordariomycetes   | Glomerellales  | Plectosphaerellaceae | Plectosphaerella | Plectosphaerella_cucumerina |
| otu_796  | 0.9045521 | 0.80 | 0.8506713 |  0.0035 | corn       | plant_pathogen    | Ascomycota    | Dothideomycetes   | Capnodiales    | Mycosphaerellaceae   | Cercospora       | NA                          |
| otu_325  | 1.0000000 | 0.60 | 0.7745967 |  0.0055 | corn       | plant_pathogen    | Ascomycota    | Sordariomycetes   | Diaporthales   | Diaporthaceae        | Diaporthe        | NA                          |
| otu_1841 | 1.0000000 | 0.60 | 0.7745967 |  0.0075 | corn       | plant_pathogen    | Ascomycota    | Dothideomycetes   | Pleosporales   | Pleosporaceae        | Curvularia       | NA                          |
| otu_521  | 0.9435383 | 0.60 | 0.7524114 |  0.0245 | corn       | plant_pathogen    | Ascomycota    | Sordariomycetes   | Glomerellales  | Plectosphaerellaceae | Lectera          | NA                          |
| otu_1013 | 0.8351648 | 0.60 | 0.7078834 |  0.0495 | corn       | plant_pathogen    | Ascomycota    | Sordariomycetes   | Xylariales     | Microdochiaceae      | Microdochium     | Microdochium_colombiense    |
| otu_332  | 0.9342302 | 0.75 | 0.8370619 |  0.0225 | restored   | plant_pathogen    | Ascomycota    | Sordariomycetes   | Glomerellales  | Plectosphaerellaceae | Plectosphaerella | NA                          |
| otu_135  | 0.7621316 | 1.00 | 0.8730015 |  0.0105 | remnant    | plant_pathogen    | Ascomycota    | Sordariomycetes   | Hypocreales    | Nectriaceae          | Ilyonectria      | NA                          |
| otu_504  | 0.6911172 | 1.00 | 0.8313346 |  0.0230 | remnant    | plant_pathogen    | Ascomycota    | Dothideomycetes   | Pleosporales   | Massarinaceae        | Stagonospora     | NA                          |
| otu_319  | 0.6983745 | 0.75 | 0.7237271 |  0.0290 | remnant    | plant_pathogen    | Basidiomycota | Ustilaginomycetes | Ustilaginales  | Ustilaginaceae       | Ustilago         | Ustilago_nunavutica         |
| otu_1716 | 1.0000000 | 0.50 | 0.7071068 |  0.0180 | remnant    | plant_pathogen    | Ascomycota    | Sordariomycetes   | Hypocreales    | Nectriaceae          | Volutella        | NA                          |
| otu_1    | 0.4509032 | 1.00 | 0.6714933 |  0.0460 | remnant    | plant_pathogen    | Ascomycota    | Sordariomycetes   | Hypocreales    | Nectriaceae          | Fusarium         | Fusarium_oxysporum          |

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
    ## [1] 8
    ## 
    ## $rrfd
    ## # A tibble: 25 × 117
    ##    field_key otu_11 otu_20 otu_29 otu_39 otu_76 otu_117 otu_120 otu_130 otu_169
    ##        <dbl>  <dbl>  <dbl>  <dbl>  <dbl>  <dbl>   <dbl>   <dbl>   <dbl>   <dbl>
    ##  1         1      0     37    229    177      0      12       0     120       0
    ##  2         2     23    167     12    133    176      10       9       0       2
    ##  3         3     94    194      0      7      8       0       0       0       0
    ##  4         4      0     71     47    302      0       0       0       0       0
    ##  5         5     22     55      0    370      0       0       0      88      29
    ##  6         6    495      4    144     14      0      13       0       0       0
    ##  7         7    139    288     45     18      8       0       0       0       0
    ##  8         8     37     47     33    234      0      21       0       0       0
    ##  9         9     27     19    145    183      0       0      10      66       0
    ## 10        10     12     42     66    157      1       5     107     152       0
    ## # … with 15 more rows, and 107 more variables: otu_202 <dbl>, otu_252 <dbl>,
    ## #   otu_266 <dbl>, otu_287 <dbl>, otu_322 <dbl>, otu_329 <dbl>, otu_331 <dbl>,
    ## #   otu_333 <dbl>, otu_341 <dbl>, otu_365 <dbl>, otu_370 <dbl>, otu_397 <dbl>,
    ## #   otu_415 <dbl>, otu_437 <dbl>, otu_438 <dbl>, otu_487 <dbl>, otu_508 <dbl>,
    ## #   otu_589 <dbl>, otu_599 <dbl>, otu_606 <dbl>, otu_632 <dbl>, otu_633 <dbl>,
    ## #   otu_634 <dbl>, otu_703 <dbl>, otu_770 <dbl>, otu_786 <dbl>, otu_790 <dbl>,
    ## #   otu_793 <dbl>, otu_818 <dbl>, otu_852 <dbl>, otu_853 <dbl>, …
    ## 
    ## $rrfd_speTaxa
    ## # A tibble: 466 × 14
    ##    field_key otu_num seq_abund phylum   class order family genus species prima…¹
    ##        <dbl> <chr>       <dbl> <chr>    <chr> <chr> <chr>  <chr> <chr>   <chr>  
    ##  1         1 otu_20         37 Ascomyc… Sord… Hypo… Bione… Clon… <NA>    wood_s…
    ##  2         1 otu_29        229 Ascomyc… Sord… Hypo… Nectr… Mari… Marian… wood_s…
    ##  3         1 otu_39        177 Ascomyc… Doth… Pleo… Cucur… Pyre… uniden… wood_s…
    ##  4         1 otu_117        12 Ascomyc… Doth… Pleo… Lindg… Cloh… Clohes… wood_s…
    ##  5         1 otu_130       120 Basidio… Agar… Trec… Hydno… Subu… <NA>    wood_s…
    ##  6         1 otu_333         4 Ascomyc… Leot… Helo… Helot… Scyt… uniden… wood_s…
    ##  7         1 otu_415        27 Ascomyc… Leot… Helo… Helot… Scyt… Scytal… wood_s…
    ##  8         1 otu_508        17 Ascomyc… Doth… Pleo… Lenti… Keis… Keissl… wood_s…
    ##  9         1 otu_599         5 Ascomyc… Doth… Pleo… Didym… Para… Paraph… wood_s…
    ## 10         1 otu_633         5 Ascomyc… Doth… Pleo… Lenti… Keis… Keissl… wood_s…
    ## # … with 456 more rows, 4 more variables: field_name <chr>, region <chr>,
    ## #   field_type <ord>, yr_since <dbl>, and abbreviated variable name
    ## #   ¹​primary_lifestyle

Sequence depth is low; these aren’t abundant taxa.

``` r
wsap_div <- calc_diversity(wsap$rrfd)
```

``` r
(wasp_comp <- gudicom(wsap_div, wsap$rrfd_speTaxa, "wood_saprotroph"))
```

<img src="microbial_guild_taxonomy_files/figure-gfm/wsap_composition-1.png" style="display: block; margin: auto;" /><img src="microbial_guild_taxonomy_files/figure-gfm/wsap_composition-2.png" style="display: block; margin: auto;" /><img src="microbial_guild_taxonomy_files/figure-gfm/wsap_composition-3.png" style="display: block; margin: auto;" />

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

| field_type | n_otu |  stat_avg |   stat_sd |
|:-----------|------:|----------:|----------:|
| corn       |     3 | 0.8153075 | 0.0908325 |
| remnant    |     2 | 0.7327736 | 0.0532201 |

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

| otu_num |         A |    B |      stat | p.value | field_type | primary_lifestyle | phylum        | class           | order           | family              | genus             | species                    |
|:--------|----------:|-----:|----------:|--------:|:-----------|:------------------|:--------------|:----------------|:----------------|:--------------------|:------------------|:---------------------------|
| otu_589 | 0.9789030 | 0.80 | 0.8849420 |  0.0030 | corn       | wood_saprotroph   | Ascomycota    | Sordariomycetes | Hypocreales     | Stachybotryaceae    | Stachybotrys      | Stachybotrys_limonispora   |
| otu_11  | 0.7198085 | 1.00 | 0.8484153 |  0.0020 | corn       | wood_saprotroph   | Ascomycota    | Sordariomycetes | Sordariales     | Chaetomiaceae       | Humicola          | Humicola_grisea            |
| otu_341 | 0.8462485 | 0.60 | 0.7125651 |  0.0275 | corn       | wood_saprotroph   | Basidiomycota | Agaricomycetes  | Agaricales      | Psathyrellaceae     | Psathyrella       | NA                         |
| otu_599 | 0.7913669 | 0.75 | 0.7704059 |  0.0355 | remnant    | wood_saprotroph   | Ascomycota    | Dothideomycetes | Pleosporales    | Didymosphaeriaceae  | Paraphaeosphaeria | Paraphaeosphaeria_michotii |
| otu_881 | 0.9664430 | 0.50 | 0.6951413 |  0.0155 | remnant    | wood_saprotroph   | Ascomycota    | Eurotiomycetes  | Chaetothyriales | Herpotrichiellaceae | Minimelanolocus   | Minimelanolocus_asiaticus  |

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
    ## [1] 22
    ## 
    ## $rrfd
    ## # A tibble: 25 × 120
    ##    field_key otu_18 otu_105 otu_126 otu_133 otu_147 otu_151 otu_164 otu_225
    ##        <dbl>  <dbl>   <dbl>   <dbl>   <dbl>   <dbl>   <dbl>   <dbl>   <dbl>
    ##  1         1    126       0       9       0       0      91       0       7
    ##  2         2     87     107       2       0      14       0       0       1
    ##  3         3     44      30     176       0       0       0       0      12
    ##  4         4     32       0       0       0       0       0       0       0
    ##  5         5     26       0      51       0       0       0       0      12
    ##  6         6     88       0      95       0      22       0       0      51
    ##  7         7     62       0     124       0      37       0       0       0
    ##  8         8      9       0      18     203       7       0       0       0
    ##  9         9     60       0       0       0       0       0       0       0
    ## 10        10     56      19       0       0       0      65       0       0
    ## # … with 15 more rows, and 111 more variables: otu_226 <dbl>, otu_265 <dbl>,
    ## #   otu_267 <dbl>, otu_272 <dbl>, otu_286 <dbl>, otu_302 <dbl>, otu_326 <dbl>,
    ## #   otu_358 <dbl>, otu_393 <dbl>, otu_414 <dbl>, otu_445 <dbl>, otu_457 <dbl>,
    ## #   otu_484 <dbl>, otu_503 <dbl>, otu_542 <dbl>, otu_551 <dbl>, otu_560 <dbl>,
    ## #   otu_574 <dbl>, otu_608 <dbl>, otu_618 <dbl>, otu_623 <dbl>, otu_653 <dbl>,
    ## #   otu_660 <dbl>, otu_698 <dbl>, otu_707 <dbl>, otu_729 <dbl>, otu_732 <dbl>,
    ## #   otu_761 <dbl>, otu_766 <dbl>, otu_789 <dbl>, otu_804 <dbl>, …
    ## 
    ## $rrfd_speTaxa
    ## # A tibble: 429 × 14
    ##    field_key otu_num seq_abund phylum   class order family genus species prima…¹
    ##        <dbl> <chr>       <dbl> <chr>    <chr> <chr> <chr>  <chr> <chr>   <chr>  
    ##  1         1 otu_18        126 Ascomyc… Doth… Capn… Clado… Clad… <NA>    litter…
    ##  2         1 otu_126         9 Ascomyc… Sord… Sord… Chaet… Chae… <NA>    litter…
    ##  3         1 otu_151        91 Basidio… Agar… Agar… Entol… Clit… uniden… litter…
    ##  4         1 otu_225         7 Chytrid… Rhiz… Rhiz… Rhizo… Rhiz… Rhizop… litter…
    ##  5         1 otu_265         2 Chytrid… Rhiz… Rhiz… Rhizo… Rhiz… uniden… litter…
    ##  6         1 otu_272        23 Ascomyc… Doth… Pleo… Phaeo… Neos… <NA>    litter…
    ##  7         1 otu_286         5 Ascomyc… Doth… Pleo… Phaeo… Neos… <NA>    litter…
    ##  8         1 otu_326         2 Ascomyc… Doth… Pleo… Dicty… Dict… Dictyo… litter…
    ##  9         1 otu_414         2 Ascomyc… Euro… Chae… Cyphe… Cyph… <NA>    litter…
    ## 10         1 otu_457         3 Ascomyc… Euro… Chae… Cyphe… Cyph… uniden… litter…
    ## # … with 419 more rows, 4 more variables: field_name <chr>, region <chr>,
    ## #   field_type <ord>, yr_since <dbl>, and abbreviated variable name
    ## #   ¹​primary_lifestyle

Sequencing depth of 297, perhaps too rare to justify examination.

``` r
lsap_div <- calc_diversity(lsap$rrfd)
```

``` r
(lsap_comp <- gudicom(lsap_div, lsap$rrfd_speTaxa, "litter_saprotroph"))
```

<img src="microbial_guild_taxonomy_files/figure-gfm/lsap_composition-1.png" style="display: block; margin: auto;" /><img src="microbial_guild_taxonomy_files/figure-gfm/lsap_composition-2.png" style="display: block; margin: auto;" /><img src="microbial_guild_taxonomy_files/figure-gfm/lsap_composition-3.png" style="display: block; margin: auto;" />

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
| corn       |     2 | 0.8347634 | 0.0850886 |
| remnant    |     1 | 0.6593805 |        NA |

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

| otu_num  |         A |   B |      stat | p.value | field_type | primary_lifestyle | phylum          | class                 | order             | family             | genus         | species               |
|:---------|----------:|----:|----------:|--------:|:-----------|:------------------|:----------------|:----------------------|:------------------|:-------------------|:--------------|:----------------------|
| otu_126  | 0.8008999 | 1.0 | 0.8949301 |  0.0055 | corn       | litter_saprotroph | Ascomycota      | Sordariomycetes       | Sordariales       | Chaetomiaceae      | Chaetomium    | NA                    |
| otu_1009 | 1.0000000 | 0.6 | 0.7745967 |  0.0055 | corn       | litter_saprotroph | Ascomycota      | Pezizomycetes         | Pezizales         | Pyronemataceae     | Cheilymenia   | Cheilymenia_stercorea |
| otu_1302 | 0.8695652 | 0.5 | 0.6593805 |  0.0365 | remnant    | litter_saprotroph | Chytridiomycota | Rhizophlyctidomycetes | Rhizophlyctidales | Rhizophlyctidaceae | Rhizophlyctis | unidentified          |

Indicator species of litter saprotrophs

## AMF

Function output is verbose but retained as explained previously.

``` r
amf_otu_summary <- amf_tax(spe_meta$amf_rfy)
```

    ## ---------------------------------
    ## [1] "AMF"
    ## ---------------------------------
    ## 
    ## 
    ## family                     corn   restored   remnant
    ## ---------------------  --------  ---------  --------
    ## Glomeraceae             13953.2    13172.2   14150.8
    ## Claroideoglomeraceae      468.2     1853.2    1696.8
    ## Paraglomeraceae          1613.0     1009.8     419.8
    ## Diversisporaceae          554.4      395.4     264.2
    ## Gigasporaceae              36.3      102.7      97.5
    ## Acaulosporaceae             2.0       25.7      48.5
    ## Archaeosporaceae            0.0      137.1      12.5
    ## Ambisporaceae               0.0        0.0       1.0
    ## 
    ## ---------------------------------
    ## [1] "Compare abundances across field types with mixed model"
    ## ---------------------------------
    ## 
    ## ---------------------------------
    ## [1] "Claroideoglomeraceae"
    ## ---------------------------------
    ## Linear mixed model fit by maximum likelihood  ['lmerMod']
    ## Formula: seq_sum ~ field_type + (1 | region)
    ##    Data: amf_df %>% filter(family == test_families[i])
    ##       AIC       BIC    logLik  deviance  df.resid 
    ##  419.6787  425.7731 -204.8393  409.6787        20 
    ## Random effects:
    ##  Groups   Name        Std.Dev.
    ##  region   (Intercept)   0.0   
    ##  Residual             875.4   
    ## Number of obs: 25, groups:  region, 4
    ## Fixed Effects:
    ##  (Intercept)  field_type.L  field_type.Q  
    ##       1339.4         868.7        -629.3  
    ## optimizer (nloptwrap) convergence code: 0 (OK) ; 0 optimizer warnings; 1 lme4 warnings 
    ## ----------------------------------------------------
    ## 
    ## Linear mixed model fit by maximum likelihood  ['lmerMod']
    ## Formula: seq_sum ~ 1 + (1 | region)
    ##    Data: amf_df %>% filter(family == test_families[i])
    ##       AIC       BIC    logLik  deviance  df.resid 
    ##  423.8524  427.5090 -208.9262  417.8524        22 
    ## Random effects:
    ##  Groups   Name        Std.Dev.
    ##  region   (Intercept)    0    
    ##  Residual             1031    
    ## Number of obs: 25, groups:  region, 4
    ## Fixed Effects:
    ## (Intercept)  
    ##        1551  
    ## optimizer (nloptwrap) convergence code: 0 (OK) ; 0 optimizer warnings; 1 lme4 warnings 
    ## ----------------------------------------------------
    ## 
    ## Data: amf_df %>% filter(family == test_families[i])
    ## Models:
    ## mmod_null: seq_sum ~ 1 + (1 | region)
    ## mmod: seq_sum ~ field_type + (1 | region)
    ##           npar    AIC    BIC  logLik deviance  Chisq Df Pr(>Chisq)  
    ## mmod_null    3 423.85 427.51 -208.93   417.85                       
    ## mmod         5 419.68 425.77 -204.84   409.68 8.1737  2    0.01679 *
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
    ## Fit: lmer(formula = seq_sum ~ field_type + (1 | region), data = amf_df %>% 
    ##     filter(family == test_families[i]), REML = FALSE)
    ## 
    ## Linear Hypotheses:
    ##                         Estimate Std. Error z value Pr(>|z|)   
    ## restored - corn == 0      1385.0      448.5   3.088  0.00575 **
    ## remnant - corn == 0       1228.5      587.2   2.092  0.08943 . 
    ## remnant - restored == 0   -156.4      489.3  -0.320  0.94430   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## (Adjusted p values reported -- single-step method)
    ## 
    ##     corn restored  remnant 
    ##      "a"      "b"     "ab" 
    ## 
    ## 
    ## ---------------------------------
    ## [1] "Diversisporaceae"
    ## ---------------------------------
    ## Linear mixed model fit by maximum likelihood  ['lmerMod']
    ## Formula: seq_sum ~ field_type + (1 | region)
    ##    Data: amf_df %>% filter(family == test_families[i])
    ##       AIC       BIC    logLik  deviance  df.resid 
    ##  382.8467  388.9411 -186.4234  372.8467        20 
    ## Random effects:
    ##  Groups   Name        Std.Dev.
    ##  region   (Intercept) 117.5   
    ##  Residual             405.5   
    ## Number of obs: 25, groups:  region, 4
    ## Fixed Effects:
    ##  (Intercept)  field_type.L  field_type.Q  
    ##       393.79       -190.03         11.83  
    ## ----------------------------------------------------
    ## 
    ## Linear mixed model fit by maximum likelihood  ['lmerMod']
    ## Formula: seq_sum ~ 1 + (1 | region)
    ##    Data: amf_df %>% filter(family == test_families[i])
    ##       AIC       BIC    logLik  deviance  df.resid 
    ##  379.8286  383.4852 -186.9143  373.8286        22 
    ## Random effects:
    ##  Groups   Name        Std.Dev.
    ##  region   (Intercept) 124.2   
    ##  Residual             412.7   
    ## Number of obs: 25, groups:  region, 4
    ## Fixed Effects:
    ## (Intercept)  
    ##       393.9  
    ## ----------------------------------------------------
    ## 
    ## Data: amf_df %>% filter(family == test_families[i])
    ## Models:
    ## mmod_null: seq_sum ~ 1 + (1 | region)
    ## mmod: seq_sum ~ field_type + (1 | region)
    ##           npar    AIC    BIC  logLik deviance  Chisq Df Pr(>Chisq)
    ## mmod_null    3 379.83 383.49 -186.91   373.83                     
    ## mmod         5 382.85 388.94 -186.42   372.85 0.9819  2     0.6121
    ## ----------------------------------------------------
    ## 
    ## 
    ##   Simultaneous Tests for General Linear Hypotheses
    ## 
    ## Multiple Comparisons of Means: Tukey Contrasts
    ## 
    ## 
    ## Fit: lmer(formula = seq_sum ~ field_type + (1 | region), data = amf_df %>% 
    ##     filter(family == test_families[i]), REML = FALSE)
    ## 
    ## Linear Hypotheses:
    ##                         Estimate Std. Error z value Pr(>|z|)
    ## restored - corn == 0      -148.9      209.6  -0.710    0.754
    ## remnant - corn == 0       -268.7      272.5  -0.986    0.581
    ## remnant - restored == 0   -119.9      228.9  -0.524    0.858
    ## (Adjusted p values reported -- single-step method)
    ## 
    ##     corn restored  remnant 
    ##      "a"      "a"      "a" 
    ## 
    ## 
    ## ---------------------------------
    ## [1] "Glomeraceae"
    ## ---------------------------------
    ## Linear mixed model fit by maximum likelihood  ['lmerMod']
    ## Formula: seq_sum ~ field_type + (1 | region)
    ##    Data: amf_df %>% filter(family == test_families[i])
    ##       AIC       BIC    logLik  deviance  df.resid 
    ##  451.7504  457.8448 -220.8752  441.7504        20 
    ## Random effects:
    ##  Groups   Name        Std.Dev.
    ##  region   (Intercept)    0    
    ##  Residual             1662    
    ## Number of obs: 25, groups:  region, 4
    ## Fixed Effects:
    ##  (Intercept)  field_type.L  field_type.Q  
    ##      13758.7         139.7         718.3  
    ## optimizer (nloptwrap) convergence code: 0 (OK) ; 0 optimizer warnings; 1 lme4 warnings 
    ## ----------------------------------------------------
    ## 
    ## Linear mixed model fit by maximum likelihood  ['lmerMod']
    ## Formula: seq_sum ~ 1 + (1 | region)
    ##    Data: amf_df %>% filter(family == test_families[i])
    ##       AIC       BIC    logLik  deviance  df.resid 
    ##  449.3055  452.9621 -221.6528  443.3055        22 
    ## Random effects:
    ##  Groups   Name        Std.Dev.
    ##  region   (Intercept)    0    
    ##  Residual             1715    
    ## Number of obs: 25, groups:  region, 4
    ## Fixed Effects:
    ## (Intercept)  
    ##       13485  
    ## optimizer (nloptwrap) convergence code: 0 (OK) ; 0 optimizer warnings; 1 lme4 warnings 
    ## ----------------------------------------------------
    ## 
    ## Data: amf_df %>% filter(family == test_families[i])
    ## Models:
    ## mmod_null: seq_sum ~ 1 + (1 | region)
    ## mmod: seq_sum ~ field_type + (1 | region)
    ##           npar    AIC    BIC  logLik deviance  Chisq Df Pr(>Chisq)
    ## mmod_null    3 449.31 452.96 -221.65   443.31                     
    ## mmod         5 451.75 457.84 -220.88   441.75 1.5551  2     0.4595
    ## ----------------------------------------------------
    ## 
    ## 
    ##   Simultaneous Tests for General Linear Hypotheses
    ## 
    ## Multiple Comparisons of Means: Tukey Contrasts
    ## 
    ## 
    ## Fit: lmer(formula = seq_sum ~ field_type + (1 | region), data = amf_df %>% 
    ##     filter(family == test_families[i]), REML = FALSE)
    ## 
    ## Linear Hypotheses:
    ##                         Estimate Std. Error z value Pr(>|z|)
    ## restored - corn == 0      -781.0      851.8  -0.917    0.625
    ## remnant - corn == 0        197.5     1115.2   0.177    0.983
    ## remnant - restored == 0    978.6      929.4   1.053    0.538
    ## (Adjusted p values reported -- single-step method)
    ## 
    ##     corn restored  remnant 
    ##      "a"      "a"      "a" 
    ## 
    ## 
    ## ---------------------------------
    ## [1] "Paraglomeraceae"
    ## ---------------------------------
    ## Linear mixed model fit by maximum likelihood  ['lmerMod']
    ## Formula: seq_sum ~ field_type + (1 | region)
    ##    Data: amf_df %>% filter(family == test_families[i])
    ##       AIC       BIC    logLik  deviance  df.resid 
    ##  432.6036  438.6980 -211.3018  422.6036        20 
    ## Random effects:
    ##  Groups   Name        Std.Dev.
    ##  region   (Intercept)    0    
    ##  Residual             1134    
    ## Number of obs: 25, groups:  region, 4
    ## Fixed Effects:
    ##  (Intercept)  field_type.L  field_type.Q  
    ##     1014.188      -843.755         5.358  
    ## optimizer (nloptwrap) convergence code: 0 (OK) ; 0 optimizer warnings; 1 lme4 warnings 
    ## ----------------------------------------------------
    ## 
    ## Linear mixed model fit by maximum likelihood  ['lmerMod']
    ## Formula: seq_sum ~ 1 + (1 | region)
    ##    Data: amf_df %>% filter(family == test_families[i])
    ##       AIC       BIC    logLik  deviance  df.resid 
    ##  430.9737  434.6303 -212.4869  424.9737        22 
    ## Random effects:
    ##  Groups   Name        Std.Dev.
    ##  region   (Intercept)    0    
    ##  Residual             1189    
    ## Number of obs: 25, groups:  region, 4
    ## Fixed Effects:
    ## (Intercept)  
    ##        1036  
    ## optimizer (nloptwrap) convergence code: 0 (OK) ; 0 optimizer warnings; 1 lme4 warnings 
    ## ----------------------------------------------------
    ## 
    ## Data: amf_df %>% filter(family == test_families[i])
    ## Models:
    ## mmod_null: seq_sum ~ 1 + (1 | region)
    ## mmod: seq_sum ~ field_type + (1 | region)
    ##           npar    AIC    BIC  logLik deviance  Chisq Df Pr(>Chisq)
    ## mmod_null    3 430.97 434.63 -212.49   424.97                     
    ## mmod         5 432.60 438.70 -211.30   422.60 2.3701  2     0.3057
    ## ----------------------------------------------------
    ## 
    ## 
    ##   Simultaneous Tests for General Linear Hypotheses
    ## 
    ## Multiple Comparisons of Means: Tukey Contrasts
    ## 
    ## 
    ## Fit: lmer(formula = seq_sum ~ field_type + (1 | region), data = amf_df %>% 
    ##     filter(family == test_families[i]), REML = FALSE)
    ## 
    ## Linear Hypotheses:
    ##                         Estimate Std. Error z value Pr(>|z|)
    ## restored - corn == 0      -603.2      580.8  -1.039    0.548
    ## remnant - corn == 0      -1193.3      760.4  -1.569    0.255
    ## remnant - restored == 0   -590.1      633.7  -0.931    0.616
    ## (Adjusted p values reported -- single-step method)
    ## 
    ##     corn restored  remnant 
    ##      "a"      "a"      "a" 
    ## 
    ## 
    ## ---------------------------------
    ## [1] "Gigasporaceae"
    ## ---------------------------------
    ## Linear mixed model fit by maximum likelihood  ['lmerMod']
    ## Formula: seq_sum ~ field_type + (1 | region)
    ##    Data: amf_df %>% filter(family == test_families[i])
    ##       AIC       BIC    logLik  deviance  df.resid 
    ##  239.5879  244.3101 -114.7940  229.5879        14 
    ## Random effects:
    ##  Groups   Name        Std.Dev.
    ##  region   (Intercept)   0.0   
    ##  Residual             101.8   
    ## Number of obs: 19, groups:  region, 3
    ## Fixed Effects:
    ##  (Intercept)  field_type.L  field_type.Q  
    ##        78.85         43.25        -29.23  
    ## optimizer (nloptwrap) convergence code: 0 (OK) ; 0 optimizer warnings; 1 lme4 warnings 
    ## ----------------------------------------------------
    ## 
    ## Linear mixed model fit by maximum likelihood  ['lmerMod']
    ## Formula: seq_sum ~ 1 + (1 | region)
    ##    Data: amf_df %>% filter(family == test_families[i])
    ##       AIC       BIC    logLik  deviance  df.resid 
    ##  236.6176  239.4509 -115.3088  230.6176        16 
    ## Random effects:
    ##  Groups   Name        Std.Dev.
    ##  region   (Intercept)   0.0   
    ##  Residual             104.6   
    ## Number of obs: 19, groups:  region, 3
    ## Fixed Effects:
    ## (Intercept)  
    ##       91.68  
    ## optimizer (nloptwrap) convergence code: 0 (OK) ; 0 optimizer warnings; 1 lme4 warnings 
    ## ----------------------------------------------------
    ## 
    ## Data: amf_df %>% filter(family == test_families[i])
    ## Models:
    ## mmod_null: seq_sum ~ 1 + (1 | region)
    ## mmod: seq_sum ~ field_type + (1 | region)
    ##           npar    AIC    BIC  logLik deviance  Chisq Df Pr(>Chisq)
    ## mmod_null    3 236.62 239.45 -115.31   230.62                     
    ## mmod         5 239.59 244.31 -114.79   229.59 1.0297  2     0.5976
    ## ----------------------------------------------------
    ## 
    ## 
    ##   Simultaneous Tests for General Linear Hypotheses
    ## 
    ## Multiple Comparisons of Means: Tukey Contrasts
    ## 
    ## 
    ## Fit: lmer(formula = seq_sum ~ field_type + (1 | region), data = amf_df %>% 
    ##     filter(family == test_families[i]), REML = FALSE)
    ## 
    ## Linear Hypotheses:
    ##                         Estimate Std. Error z value Pr(>|z|)
    ## restored - corn == 0      66.381     64.756   1.025    0.553
    ## remnant - corn == 0       61.167     92.915   0.658    0.783
    ## remnant - restored == 0   -5.214     76.941  -0.068    0.997
    ## (Adjusted p values reported -- single-step method)
    ## 
    ##     corn restored  remnant 
    ##      "a"      "a"      "a" 
    ## 
    ## 
    ## ---------------------------------
    ## [1] "Test abundances with years since restoration"
    ## ---------------------------------
    ## [1] "Claroideoglomeraceae"
    ## 
    ## Call:
    ## lm(formula = seq_sum ~ yr_since, data = mod_data %>% filter(family == 
    ##     all7[i]))
    ## 
    ## Residuals:
    ##       1       2       3       4       5       6       7 
    ## -764.55  -34.43 -170.06  115.20 -629.93 -345.80 1829.57 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)  
    ## (Intercept)  2439.55     620.97   3.929   0.0111 *
    ## yr_since      -58.37      41.77  -1.398   0.2211  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 947.8 on 5 degrees of freedom
    ## Multiple R-squared:  0.2809, Adjusted R-squared:  0.1371 
    ## F-statistic: 1.953 on 1 and 5 DF,  p-value: 0.2211
    ## 
    ## [1] "Diversisporaceae"
    ## 
    ## Call:
    ## lm(formula = seq_sum ~ yr_since, data = mod_data %>% filter(family == 
    ##     all7[i]))
    ## 
    ## Residuals:
    ##       1       2       3       4       5       6       7 
    ##  -67.85 -154.22  -52.81  113.83   51.13   65.45   44.47 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)  
    ## (Intercept)  187.226     67.458   2.775   0.0391 *
    ## yr_since       5.664      4.538   1.248   0.2672  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 103 on 5 degrees of freedom
    ## Multiple R-squared:  0.2376, Adjusted R-squared:  0.08507 
    ## F-statistic: 1.558 on 1 and 5 DF,  p-value: 0.2672
    ## 
    ## [1] "Gigasporaceae"
    ## 
    ## Call:
    ## lm(formula = seq_sum ~ yr_since, data = mod_data %>% filter(family == 
    ##     all7[i]))
    ## 
    ## Residuals:
    ##       1       2       3       4       5       6       7 
    ##  97.336  10.425  -7.284 -38.600 -33.448  15.893 -44.322 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)  
    ## (Intercept)  -22.830     35.236  -0.648    0.546  
    ## yr_since       8.468      2.370   3.573    0.016 *
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 53.78 on 5 degrees of freedom
    ## Multiple R-squared:  0.7186, Adjusted R-squared:  0.6623 
    ## F-statistic: 12.77 on 1 and 5 DF,  p-value: 0.016
    ## 
    ## [1] "Glomeraceae"
    ## 
    ## Call:
    ## lm(formula = seq_sum ~ yr_since, data = mod_data %>% filter(family == 
    ##     all7[i]))
    ## 
    ## Residuals:
    ##        1        2        3        4        5        6        7 
    ##   283.75 -3281.74   -93.03   -10.71  1808.33  2086.99  -793.59 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept) 11923.55    1281.00   9.308 0.000241 ***
    ## yr_since      110.73      86.17   1.285 0.255083    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 1955 on 5 degrees of freedom
    ## Multiple R-squared:  0.2483, Adjusted R-squared:  0.09793 
    ## F-statistic: 1.651 on 1 and 5 DF,  p-value: 0.2551
    ## 
    ## [1] "Paraglomeraceae"
    ## 
    ## Call:
    ## lm(formula = seq_sum ~ yr_since, data = mod_data %>% filter(family == 
    ##     all7[i]))
    ## 
    ## Residuals:
    ##       1       2       3       4       5       6       7 
    ##   549.0  2700.7   168.0  -137.5 -1009.3 -1436.6  -834.3 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)
    ## (Intercept)  1644.09     990.18   1.660    0.158
    ## yr_since      -45.25      66.61  -0.679    0.527
    ## 
    ## Residual standard error: 1511 on 5 degrees of freedom
    ## Multiple R-squared:  0.08452,    Adjusted R-squared:  -0.09857 
    ## F-statistic: 0.4616 on 1 and 5 DF,  p-value: 0.527

From the mean sequence abundances in field types, the following families
look interesting:

- Claroideoglomeraceae: low in corn
- Paraglomeraceae: highest in corn, declines through restoration and
  remnant
- Diversisporaceae: highest in corn, declines through restoration and
  remnant
- Gigasporaceae: low in corn

These don’t all show up in linear modeling (shown next), but since this
modeling isn’t likely valid, we can explore any of these that appear
worthwhile. For now, let’s look at a boxplot for each.

``` r
spe_meta$amf_rfy %>% 
    filter(family %in% c("Claroideoglomeraceae", "Paraglomeraceae", "Diversisporaceae", "Gigasporaceae")) %>% 
    group_by(region, field_type, field_key, family) %>% 
    summarize(seq_sum = sum(seq_abund), .groups = "drop") %>% 
    ggplot(aes(x = field_type, y = seq_sum)) +
    facet_wrap(vars(family), scales = "free_y") +
    geom_boxplot(varwidth = TRUE, fill = "gray90", outlier.shape = NA) +
    geom_beeswarm(aes(fill = region), shape = 21, size = 2, dodge.width = 0.2) +
    labs(x = "", y = "Sum of sequence abundance", title = "AMF abundance in families and field types") +
    scale_fill_discrete_qualitative(palette = "Dark3") +
    theme_bw()
```

<img src="microbial_guild_taxonomy_files/figure-gfm/amf_boxplot-1.png" style="display: block; margin: auto;" />

make some comments…

### Individual families

Claroideoglomeraceae differs across field types with a likelihood ratio
test result p\<0.05. Tukey’s post-hoc test with Holm correction
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
regression, $R^2_{adj}$ = 0.66, p \< 0.05

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
