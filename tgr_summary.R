#' ---
#' title: "Tallgrass Prairie Restoration Legacies, Summary"
#' author: "Beau Larkin\n"
#' date: "Last updated: `r format(Sys.time(), '%d %B, %Y')`"
#' output:
#'  pdf_document:
#'      df_print: tibble
#'      fig_caption: yes
#'      highlight: tango
#'      toc: true
#'      toc_depth: 2
#' ---
#' # Packages and Libraries
packages_needed = c("tidyverse",
                    "png",
                    "knitr",
                    "conflicted")
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
#' # Description
#' The goal is to present results and discuss whether a path forward exists. If so, we will determine
#' the strategy that presents the best story, organization, and interpretation of these results. 
#' 
#' # Data
#' For this summary, I'll pull as many objects as possible from existing files to reduce the number of 
#' interspersed code chunks. A few quick new analyses will be necessary, though. Data are loaded here. 
#+ data_sites
sites <-
    read_csv(paste0(getwd(), "/clean_data/sites.csv"), show_col_types = FALSE) %>%
    mutate(
        field_type = factor(
            field_type,
            ordered = TRUE,
            levels = c("corn", "restored", "remnant")),
        yr_since = replace(yr_since, which(field_type %in% c("remnant", "corn")), NA)) %>%
    select(-lat, -long, -yr_restore, -yr_rank) %>% 
    arrange(field_key)
#' 
#' # Methods
#' ## Sites and design
#' The survey followed an unbalanced complete block design. Corn, restored, and remnant fields are compared, with 
#' at least one of each field type in each block. I have called blocks "regions" so far.
#' We collected samples and data from four regions, shown on the map below.  
#+ fig:site_map,echo=FALSE,fig.align='center',out.width="90%",fig.cap="\\label{fig:site_map}Labels show centroids of regions used for this work. BM = Blue Mounds, FG = Faville Grove, FL = Fermilab, LP = Lake Petite."
include_graphics("site_locations_files/figure-gfm/site_map-1.png")
#' 
#' The design is unbalanced because there are more restored fields than corn or remnant. In all but one case,
#' only single corn and remnant fields were available in each region. This means that we only have
#' replication to separate field types when using the entire block design. 
#+ fields_regions_types_table,echo=FALSE
kable(table(sites$region, sites$field_type),
      format = "pandoc",
      caption = "Count of fields by type and region:\nBM = Blue Mounds, FG = Faville Grove,\nFL = Fermilab, LP = Lake Petite")
#' 
#' The rules for establishing a chronosequence are strict. We cannot call fields from all regions a chronosequence. With 
#' seven restored fields in a small geographic area, the Blue Mounds fields are our best bet for this, but we will likely 
#' have to call this a "pseudochronosequence" and avoid some inferences. Mantel tests (not shown) failed to find correlations 
#' between soil variables and pairwise distance, giving us a little confidence that we're avoiding systemic bias. 
#+ fig:site_map_bm,echo=FALSE,fig.align='center',out.width="90%",fig.cap="\\label{fig:site_map_bm}Labels show individual fields in the Blue Mounds region."
include_graphics("site_locations_files/figure-gfm/site_map_bm-1.png")
#' 
#' ## Fungal communities
#' We collected soil cores from 10 haphazardly-selected locations in each field. Genomic DNA was extracted, and the 
#' lab and bioinformatics pipeline delivered community data from ITS or 18S regions clustered as OTUs or SVs. In preliminary work,
#' inferences made with SVs were weaker but not qualitatively different. I continued with OTUs only.  
#' 
#' A few samples had failed to amplify, leaving some fields represented with 9 samples instead of 10. 
#' This unbalance, which is normally not a problem, became unacceptable for a couple of reasons. First, I had planned to use permutation tests in 
#' ordinations, and application function `how()` from package `permute` requires balance at the plot (i.e., field) level. Second, 
#' I wanted to summarize the subsamples at the field level because fields are our replicates. Any comparisons with 
#' field metadata, where we have one data point per field, would be pseudoreplicated when regressed against subsample level data. 
#' 
#' Finally, after rarefaction it was clear that several samples contained very few sequences compared with the 
#' others. I used an iterative process to remove subsamples with few sequences, choosing a rarefication cutoff near the plateau of OTU 
#' recovery at depth (see figures). This ended up being 8 subsamples for ITS and 7 for 18S datasets. These subsample sequence 
#' values were then rarefied and used as-is, or, when reporting at the field level, subsample sequence values were summed
#' and rarefied. 
#' 
#' This process didn't remove many OTUs and doesn't change any major interpretation. 
#' 
#+ fig:its_rarefaction_pre,echo=FALSE,fig.align='center',out.width="100%",fig.cap="\\label{fig:its_rarefaction_pre}Pre"
include_graphics("microbial_diagnostics_pre_files/figure-gfm/its_rarefaction_curve_fig-1.png")
#+ fig:its_rarefaction_post,echo=FALSE,fig.align='center',out.width="100%",fig.cap="\\label{fig:its_rarefaction_post}Post"
include_graphics("microbial_diagnostics_post_files/figure-gfm/its_rarefaction_curve_fig-1.png")
#' 
#' Even after this process, it's clear that fungal communities were undersampled. 
#+ fig:its_accum,echo=FALSE,fig.align='center',out.width="100%",fig.cap="\\label{fig:its_accum}ITS"
include_graphics("microbial_diagnostics_post_files/figure-gfm/its_species_accumulation_fig-1.png")
#+ fig:amf_accum,echo=FALSE,fig.align='center',out.width="100%",fig.cap="\\label{fig:amf_accum}18S"
include_graphics("microbial_diagnostics_post_files/figure-gfm/amf_species_accumulation_fig-1.png")
#' 
#' ## Plant data
#' Plant community data resulted from two independent surveys. In the Wisconsin regions, haphazard transects were established
#' and 10 meter frames placed, with percent cover estimated for all species, resulting in a dataset with plant composition. 
#' In Fermi, relevé methods were used, resulting in presence/absence data only. Plant metadata were assembled from multiple 
#' sources and include plant traits and natural history details. 
#' 
#' ## Soil data
#' In the field, soil was pooled from 10 haphazardly-selected locations, mixed, and sampled once for soil chemical analysis.
#' Soil data includes abiotic macro and micronutrients, organic carbon, and properties like pH. 
#' Average precipitation was determined for each field using PRISM climate data and is included with the soil data. 
#' 
#' ## Design/site data
#' 
#' - Field type: corn, restored, remnant
#' - Field age: years since restoration (NA for corn and remnant)
#' - Region: blocks
#' 
#' ## Response data
#' Also taken from the pooled soil in each field, one sample was taken for analysis of Water Stable Aggregates (I don't know who did these),
#' and one for quantification of microbial biomass. 

#' ## Data sources, summary
#' 
#' - Fungal genomic data, ITS and 18S, OTU clusters
#' - Plant community data, composition in Wisconsin and presence/absence in Fermi
#' - Plant traits and natural history
#' - Soil properties
#' - Fungal biomass
#' - Water stable aggregates
#' - Site metadata and design
#' 
#' # Results
#' 
#' ## Fungal communities
#+ microbial_communities.R,echo=FALSE,message=FALSE,warning=FALSE,include=FALSE
source("microbial_communities.R")
#' 
#' ### Soil fungi, field types, and regions
#' Cornfields stand apart while remnants are mixed with restored fields. As the age of restored fields
#' increases, their distance from corn and remnants grows. The high dimensionality here keeps the explanatory power
#' of axes 1 & 2 low. Clusters based on field type were significant based on a permanova (with regions as blocks and 
#' subsamples held constant within fields, the replicate). However, a pairwise post-hoc permutation failed to separate 
#' remnant fields from any other type. This was also true when field centroids were used exclusively and the 
#' permutation test adjusted to this scheme. 
#+ fig:its_samps_unified_fig-1,echo=FALSE,fig.align='center',out.width="100%",fig.cap="ITS-based fungal communities, most abundant guilds inset."
include_graphics("microbial_communities_files/figure-gfm/its_samps_unified_fig-1.png")
#' 
#+ pcoa_its_samps_permanova,echo=FALSE
print(pcoa_its_samps$permanova)
#+ pcoa_its_samps_pairwise,echo=FALSE
print(pcoa_its_samps$pairwise_contrasts)
#' 
#' ### Soil fungal guilds, field types, and regions
#' By sequence abundance across the entire dataset, the top guilds are:
#' 
#' 1. Unidentified (not shown on column charts)
#' 1. plant pathogens
#' 1. soil saprotrophs
#' 1. wood saprotrophs
#' 1. dung saprotrophs
#' 1. litter saprotrophs
#' 
#' Compared with the sequence abundance in the NA group, plant pathogens and soil saprotrophs are
#' abundant enough to feel somewhat confident about in terms of coverage. These are also the ones
#' that varied most with age in restored fields as we'll see later, and they were most 
#' conspicuous in an indicator species analysis. Trends across field types within regions
#' aren't very consistent...
#+ fig:its_guilds_regions,echo=FALSE,fig.align='center',out.width="100%",fig.cap="Abundant or dynamic guilds across field types in regions, ITS dataset."
include_graphics("microbial_communities_files/figure-gfm/its_guilds_regions_fig-1.png")
#' ## heading
#' more text
#' 
#' ### heading again
#+ fig:amf_pcoa_unified,echo=FALSE,fig.align='center',out.width="100%",fig.cap="18S-based fungal communities, most important families inset."
include_graphics("microbial_communities_files/figure-gfm/amf_samps_unified_fig-1.png")



