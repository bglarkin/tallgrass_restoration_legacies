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
#' # Functions 
#' Functions streamline data processing, model fitting, and results output. Functions for this 
#' script are found in a [supplemental script](microbial_guild_taxonomy_functions.R) and are loaded 
#' here for convenience. 
#+ mgt_functions
source("supporting_files/microbial_guild_taxonomy_functions.R")
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
#' 
#' ## Species metadata
#' Load taxonomy for all and guilds (called *primary lifestyle* in Fungal Traits)
#' for ITS OTUs. Replace NA values with "unidentified" to show complete numbers of 
#' unidentified groups.
meta <- list(
    its = read_csv(
        paste0(getwd(), "/clean_data/spe_ITS_metadata.csv"),
        show_col_types = FALSE
    ),
    amf = read_csv(
        paste0(getwd(), "/clean_data/spe_18S_metadata.csv"),
        show_col_types = FALSE
    )
) %>% 
    map(. %>% mutate(across(everything(), ~ replace_na(., "unidentified"))))
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
        join_spe_meta(spe$its_raw, meta$its) %>%
        write_csv(paste0(getwd(), "/clean_data/speTaxa_ITS_raw.csv")),
    its_rfy = 
        join_spe_meta(spe$its_rfy, meta$its) %>%
        write_csv(paste0(getwd(), "/clean_data/speTaxa_ITS_rfy.csv")),
    amf_raw = 
        join_spe_meta(spe$amf_raw, meta$amf) %>%
        write_csv(paste0(getwd(), "/clean_data/speTaxa_18S_raw.csv")),
    amf_rfy = 
        join_spe_meta(spe$amf_rfy, meta$amf) %>%
        write_csv(paste0( getwd(), "/clean_data/speTaxa_18S_rfy.csv" ))
)
#' 
#' # Analysis and Results
#' ## ITS sequences
#' Recall the number of OTUs recovered in each dataset. The effect of rarefying did not change
#' richness or diversity very much. 
# Number of OTUs in raw and rarefied datasets
Map(function(x) ncol(x)-1, spe[1:2])
#' 
#' ### Unassigned taxa
#' Only 21.8 percent of the ITS sequences were assigned to species. In terms of the analysis done here, 
#' its possibly more alarming that only 36 percent were assigned to primary lifestyles or guilds. 
#' This suggests that when we see guilds concentrating in certain habitats, it's possible that
#' the difference doesn't exist. This is particularly possible because we have one habitat, cornfields, 
#' which has probably been studied more than the others. 
#+ unassigned_table
meta$its %>% 
    select(-otu_num, -otu_ID) %>% 
    map(\(x) round(length(which(x == "unidentified")) / length(x) * 100, 1)) %>% 
    bind_rows() %>% 
    kable(format = "pandoc", caption = "Percent unidentified OTUs in each taxonomic group or guild")
#' ### Composition in field types
#' Function outputs are verbose, but details may be necessary later so they are displayed here.
#+ its_tax_trophic_otu,message=FALSE,fig.height=7,fig.align='center'
its_taxaGuild(spe_meta$its_rfy)
#' 
#' The top guilds are:
#' 
#' 1. Unidentified (not shown on column charts)
#' 1. plant pathogens
#' 1. soil saprotrophs
#' 1. wood saprotrophs
#' 1. dung saprotrophs
#' 1. litter saprotrophs
#' 
#' Compared with the sequence abundance in the NA group, plant pathogens and soil saprotrophs are
#' abundant enough to feel somewhat confident about in terms of coverage. 
#' 
#+ its_guilds_otu,message=FALSE
its_rfy_guilds <- its_test_taxaGuild(spe_meta$its_rfy)
#' 
#' Model tests on `field_type` are shaky due to unbalance, but are included here
#' to point out trends that we may be able to present in some better way. Trends 
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
#' a weakness that we lack replication within blocks for `field_type` in some regions. 
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
    inspan(., 1999, meta$its)
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
#' More species have higher specificity and fidelity to cornfields than to the other field types. The top ten
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
    write_csv(., paste0(getwd(), "/microbial_guild_taxonomy_files/its_inspan.csv")) %>% 
    kable(format = "pandoc", caption = "Indicator species of ITS OTUs (top 10 per field type)")
#' 
#' ### Soil saprotrophs
#' #### Trends over time
#+ ssap_guiltime,fig.width=8,fig.height=4.5,fig.align='center'
guiltime("soil_saprotroph")
#' Sequence abundance of soil saprotrophs increases over time in the Blue Mounds 
#' area ($R^2_{Adj}=0.56, p<0.05$), but
#' this appears to be leveraged by Karla Ott's property, though. With all
#' that big bluestem...maybe there is more litter and soil carbon? It will be good 
#' to look at trends in soil chemistry. 
#' 
#' #### Diversity
#+ ssap_filgu
ssap <- filgu(spe$its_rfy, meta$its, primary_lifestyle, "soil_saprotroph", sites)
#' Out of 2889 OTUs, 260 are in this group.
#' Most OTUs contain few sequences, but several range from hundreds to 25,000 sequences.
#' The 25 samples are all retained, and vary from 4000 to 16000 sequences. None are so small that 
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
#' *Agarics* generally decrease over time and *Geoglossales* increase.
#' 
#' Soil saprotrophs remain an interesting guild. 
#' 
#' #### Indicators
#+ ssap_inspan
ssap_inspan <- 
    ssap$filspe %>% 
    left_join(sites, by = join_by(field_key)) %>% 
    inspan(., 1999, meta$its)
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
    write_csv(., paste0(getwd(), "/microbial_guild_taxonomy_files/ssap_inspan.csv")) %>% 
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
ppat <- filgu(spe$its_rfy, meta$its, primary_lifestyle, "plant_pathogen", sites)
#' Out of 2889 OTUs, 159 are in this group. 
#' All samples are retained and contain 3000-16000 sequences, so none are so limited as to bias 
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
#' In the Blue Mounds area, trends in pathogen composition over time aren't obvious. Possibly
#' *Glomerales* pathogens decrease over time and *Pleosporales* increase. 
#' 
#' #### Indicators
#+ ppat_inspan
ppat_inspan <- 
    ppat$filspe %>% 
    left_join(sites, by = join_by(field_key)) %>% 
    inspan(., 1999, meta$its)
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
    write_csv(., paste0(getwd(), "/microbial_guild_taxonomy_files/ppat_inspan.csv")) %>% 
    kable(format = "pandoc", caption = "Indicator species of plant pathogens")
#' 
#' ### Wood saprotrophs
#' #### Trends over time
#+ wsap_guiltime,fig.width=8,fig.height=4.5,fig.align='center'
guiltime("wood_saprotroph") 
#' Interestingly a strong negative relationship over time since restoration ($R^2_{Adj}=0.72, p<0.01$)
#' in sharp contrast to the increasing relationship found with soil saprotrophs. Apparently many wood
#' saprotrophs live in cornfield soil...let's see:
#' 
#' #### Diversity
#+ wsap_filgu
wsap <- filgu(spe$its_rfy, meta$its, primary_lifestyle, "wood_saprotroph", sites)
#' Out of 2889 OTUs, 120 are in this group. 
#' Samples contain 800-4400 sequences. Sequence depth is low; these aren't abundant or numerous taxa.
#' Only 123 OTUs comprise this group. 
#+ wsap_div
wsap_div <- calc_diversity(wsap$filspe)
#+ wsap_composition,message=FALSE,fig.width=7,fig.height=7,fig.align='center'
wasp_comp <- gudicom(wsap_div, wsap$filspeTaxa, "wood_saprotroph", other_threshold = 3)
#' With diversity, not much jumps out. 
#' 
#' Diversity appears high across fields and years compared with other guilds.
#' While *Agaric* soil saprotrophs increased strongly from corn to remnants, 
#' they declined when characterized as wood saprotrophs.
#' 
#' Notable changes in composition are evident over time. *Tubeufiales* declines with time
#' since restoration; *Hypocreales* increases. *Pleosporales* also appear to increase, but 
#' the colors are difficult to discern. Remember to look at tabular data.  
#' 
#' #### Indicators
#+ wsap_inspan
wsap_inspan <- 
    wsap$filspe %>% 
    left_join(sites, by = join_by(field_key)) %>% 
    inspan(., 1999, meta$its)
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
lsap <- filgu(spe$its_rfy, meta$its, primary_lifestyle, "litter_saprotroph", sites)
#' Out of 2889 OTUs, 139 are in this group. 
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
    inspan(., 1999, meta$its)
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
#' ## AMF
#' Recall the number of OTUs recovered in each dataset. The effect of rarefying did not change
#' richness or diversity very much. 
Map(function(x) ncol(x)-1, spe[3:4])
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
#' 
#' From the mean sequence abundances in field types and trends over time, the following families look interesting:
#' 
#' - *Claroideoglomeraceae:* low in corn; significantly by likelihood ratio test. Declines with time in BM,
#' but this was not significant.
#' - *Paraglomeraceae:* highest in corn, declines through restoration and remnant, declines in BM and FL but 
#' likely not a significant trend
#' - *Diversisporaceae:* highest in corn, declines through restoration and remnant
#' - *Gigasporaceae:* low in corn, and also the only one with a significant change with years
#' since restoration, and this only in Blue Mounds. Gigasporaceae increase over time 
#' ($R^2_{Adj}=0.66, p<0.05$). These are rare taxa though, and I'm not sure we can really 
#' say much about them. 
#' 
#' In the next section, we will examine these families more closely.
#' 
#' ### Claroideoglomeraceae
#+ claroid_filgu
claroid <- filgu(spe$amf_rfy, meta$amf, family, "Claroideoglomeraceae", sites)
#' Out of 144 AMF OTUs, 17 map to this family. Most are low abundance across sites, but all
#' samples are retained and contain sufficient sequences to draw meaningful conclusions. 
#+ claroid_div
claroid_div <- calc_diversity(claroid$filspe)
#+ claroid_divplot,message=FALSE,results=FALSE,fig.width=7,fig.height=7,fig.align='center'
gudicom(claroid_div, claroid$filspeTaxa, "Claroideoglomeraceae", gene = "amf")
#' Little change over time, but alpha diversity in cornfields is low compared with restored and 
#' remnant fields.
#' 
#' ### Paraglomeraceae
#+ para_filgu
para <- filgu(spe$amf_rfy, meta$amf, family, "Paraglomeraceae", sites)
#' Out of 144 AMF OTUs, only 6 map to this family. Most are low abundance across sites, but all
#' samples are retained. Any interpretation here is likely to be dominated by a couple high-abundance
#' OTUs, and a couple of samples have close to zero detections. Is this real?
#+ para_div
para_div <- calc_diversity(para$filspe)
#+ para_divplot,message=FALSE,results=FALSE,fig.width=7,fig.height=7,fig.align='center'
gudicom(para_div, para$filspeTaxa, "Paraglomeraceae", gene = "amf")
#' Richness declines with time since restoration in the Blue Mounds, but with few sequences 
#' and likely non-significant correlations, I don't see doing much more with this group. 
#' 
#' ### Diversisporaceae
#+ diver_filgu
diver <- filgu(spe$amf_rfy, meta$amf, family, "Diversisporaceae", sites)
#' Out of 144 AMF OTUs, only 8 map to this family. Most are low abundance across sites, but all
#' samples are retained. Any interpretation here is likely to be dominated by a couple high-abundance
#' OTUs, and a couple of samples have close to zero detections. Is this real?
#+ diver_div
diver_div <- calc_diversity(diver$filspe)
#+ diver_divplot,message=FALSE,results=FALSE,fig.width=7,fig.height=7,fig.align='center'
gudicom(diver_div, diver$filspeTaxa, "Diversisporaceae", gene = "amf")
#' Richness declines with time since restoration in the Blue Mounds. Few sequences 
#' and likely non-significant correlations, but these taxa definitely don't like cornfields. 
#' 
#' ### Gigasporaceae
#+ giga_filgu
giga <- filgu(spe$amf_rfy, meta$amf, family, "Gigasporaceae", sites)
#' Out of 144 AMF OTUs, only 4 map to this family. Most are low abundance across sites, and 
#' only 19 samples contain these taxa. Any interpretation here is likely to be dominated by a couple high-abundance
#' OTUs, and a couple of samples have close to zero detections. Is this real? 
#' 
#' These OTUs all dropped with previous attempts to re-rarefy in the guild. Perhaps that also supports
#' that these are just too low abundance to work more with. 
#+ giga_div
giga_div <- calc_diversity(giga$filspe)
#+ giga_seq_abund_years,message=FALSE,fig.width=7,fig.height=5,fig.align='center'
giga$filspeTaxa %>% 
    filter(field_type == "restored", region %in% c("BM", "FL")) %>% 
    group_by(region, field_type, field_key, yr_since) %>% 
    summarize(seq_sum = sum(seq_abund), .groups = "drop") %>% 
    ggplot(aes(x = yr_since, y = seq_sum)) +
    facet_wrap(vars(region), scales = "free") +
    geom_smooth(aes(linetype = region), method = "lm", se = FALSE) +
    geom_point() +
    scale_linetype_manual(values = c("solid", NA), guide = "none") +
    labs(x = "Years since restoration", y = "Sum of sequences", title = "Gigasporaceae") +
    theme_bw()
#' 
#' Pity that there are so few of these AMF. It's a nice relationship. Maybe there is a natural 
#' history angle here, like an interaction between Gigasporaceae and plant pathogens, but 
#' it will be hard to argue that it matters much given the low abundance observed. 
#' 
#' # Conclusions: taxa and guilds
#' 
#' 1. Much work remains researching the natural history of taxa identified through 
#' indicator species analysis and changes in composition across field types. 
#' 1. Cornfields are weird. They harbor more and stronger indicator species, differ obviously in 
#' composition, diversity, and richness. This is not a surprise, but it is obvious. 
#' 1. Soil saprotrophs remain interesing.
#'    1. Soil saprotrophs increase with years since restoration in the Blue Mounds. 
#'    1. Richness increases from corn to restored to remnant, and sequence
#'    abundance increases with restoration age in the Blue Mounds. 
#'    *Agarics* increase strongly from corn to remnant; *Cystofilobasidiales*
#'    and *Filobasidiales* aren't found outside of cornfields. Generally, cornfield composition looks 
#'    different than the other two, but remnants do appear somewhat intermediate. *Mortierellales* appear less 
#'.   in remnants than corn or former corn fields.
#'    1. *Agarics* generally decrease over time and *Geminibasidiales* increase in composition
#'    at Blue Mounds 
#' 1. A strong decline in pathogens is seen in Blue Mounds’ restored fields. Declines are noticeable
#' in abundance, richness, Shannon's, and Simpson's diversity (although the latter three need
#' tests for significance of correlations).
#'    1. Changes in composition across field types are subtle but present, and would benefit from 
#'    additional natural history work. Composition over time in Blue Mounds doesn't change in 
#'    obvious linear progressions. 
#' 1. Wood saprotrophs decline across years in Blue Mounds in a trend that is nearly reciprocal 
#' to that seen with soil saprotrophs. Possibly a natural history angle there. 
#'    1. Diversity appears high across fields and years compared with other guilds.
#'    While *Agaric* soil saprotrophs increased strongly from corn to remnants, 
#'    they declined when characterized as wood saprotrophs.
#'    1. Notable changes in composition are evident over time. *Tubeufiales* declines with time
#'    since restoration; *Hypocreales* increases. 
#' 1. AMF generate less interest, mostly because fewer species and smaller sequence abundances 
#' prevent development of nice relationships. 
#'    1. *Claroideoglomeraceae* is probably the strongest family, with substantial differences across
#'    field types in sequence abundance. 
#'    1. *Gigasporaceae* significantly increases across the Blue Mounds series, but 
#'    sequence abundances and richness are so small in this family that further discussion 
#'    may be inappropriate or irrelevant. 
#'    
#' # Appendix: Rarefy in guilds?
#' One practice that's becoming popular is to take the raw sequence abundance data, filter it to 
#' a guild or taxonomic group, and then rarefy in that group. I did that with these guilds
#' and families, and the action reduced sequence abundances greatly without affecting 
#' what we would interpret from the data. This appendix shows some diagnostics on that process.
#' 
#' ## Diversity with ITS sequences
#' Let's graphically compare diversity metrics between datasets. *Pre-rarefied* data are rarefied
#' as a full set (all samples and OTUs), then filtered to subsets of guilds. *Post-rarefied*
#' data are filtered to subsets of guilds from raw ITS sequence abundances and then rarefied
#' within the subset. 
#' 
#' First, post-rarefied datasets are produced for the guilds of interest in this report. Then,
#' the Hill's series of diversity metrics are calculated. 
#' 
ssap_rrfd <- rerare(spe$its_raw, meta$its, primary_lifestyle, "soil_saprotroph", sites)
ppat_rrfd <- rerare(spe$its_raw, meta$its, primary_lifestyle, "plant_pathogen", sites)
wsap_rrfd <- rerare(spe$its_raw, meta$its, primary_lifestyle, "wood_saprotroph", sites)
lsap_rrfd <- rerare(spe$its_raw, meta$its, primary_lifestyle, "litter_saprotroph", sites)
#' 
#' The diversity metrics are bound to the pre-rarefied sets produced earlier in this report.
#' These data are wrangled to facilitate plotting. 
#' 
rrfd_compare <- 
    bind_rows(
        list(
            ssap_postrare = calc_diversity(ssap_rrfd$rrfd),
            ppat_postrare = calc_diversity(ppat_rrfd$rrfd),
            wsap_postrare = calc_diversity(wsap_rrfd$rrfd),
            lsap_postrare = calc_diversity(lsap_rrfd$rrfd),
            ssap_prerare = ssap_div,
            ppat_prerare = ppat_div,
            wsap_prerare = wsap_div,
            lsap_prerare = lsap_div
        ),
        .id = "guild_rrfd_key"
    ) %>% separate_wider_delim(guild_rrfd_key, "_", names = c("guild", "rrfd_step")) %>% 
    pivot_wider(names_from = rrfd_step, values_from = value)
#' 
#' Results are plottted. 
#+ rrfd_compare_fig,fig.align='center',fig.width=9,fig.height=8
ggplot(rrfd_compare %>% filter(hill_index %in% c("N0", "N1", "N2")), aes(x = prerare, y = postrare)) +
    facet_wrap(vars(hill_index, guild), scales = "free", ncol = 4) +
    geom_abline(intercept = 0, slope = 1) +
    geom_point(aes(fill = field_type), size = 2, shape = 21) +
    labs(x = "Index value from pre-rarefied abundances", y = "Index value from post-rarefied abundances",
         caption = "N0 = richness, N1 = Shannon's Diversity, N2 = Simpson's Diversity\nlsap = litter saprotrophs, ppat = plant pathogens, ssap = soil saprotrophs, wsap = wood saprotrophs") +
    scale_fill_discrete_qualitative(palette = "Harmonic") +
    theme_bw()
#' Repeating this with AMF is unnecessary...
#' 
#' # Appendix: Composition in regions
#' ## Blue Mounds
#+ taxaguild_bm,fig.align='center'
its_taxaGuild(spe_meta$its_rfy %>% filter(region == "BM"), other_threshold = 5)
#' Soil saprotrophs increase from corn to remnant. Pathogens spike in restored fields. Dung saprotrophs disappear 
#' from corn to remnant; likely because these species are adapted to labile/high nutrient environments. Litter and wood
#' saprotrophs decrease from corn to remnant...this seems odd. 
#' 
#' ## Faville Grove
#+ taxaguild_fg,fig.align='center'
its_taxaGuild(spe_meta$its_rfy %>% filter(region == "FG"), other_threshold = 5)
#' Weak trends are apparent. From corn to remnant, soil saprotrophs increase and pathogens 
#' maybe spike again in restored. 
#' 
#' ## Fermi
#+ taxaguild_fL,fig.align='center'
its_taxaGuild(spe_meta$its_rfy %>% filter(region == "FL"), other_threshold = 5)
#' Soil saprotrophs decrease in this case, which is interesting. The remnant fields aroud Fermi weren't a carbon-rich as
#' the restored fields, maybe. 
#' 
#' ## Lake Petite
#+ taxaguild_lp,fig.align='center'
its_taxaGuild(spe_meta$its_rfy %>% filter(region == "LP"), other_threshold = 5)
#' Soil saprotrophs highest in remnants, but lowest in restored. Pathogens flat. Again, wood and litter saprotrophs
#' lowest in remnant, although also missing from corn in this case. 



