# Introduction
Agricultural practices in cropland, like soil tillage and fertilizer applications, profoundly alter the soil's physical and chemical properties over time. Soil carbon is greatly reduced and belowground metabolic pathways are greatly simplified, leading to attrition of microbial diversity. Restoration of abandoned cropland is increasingly popular, and usually involves seeding a great density of native plants into post-agricultural soil. Although these plants can establish, bacterial communities often fail to recover to a pre-agricultural state, instead forming novel communities [(Badger and Docherty 2022)](https://link.springer.com/10.1007/s00248-022-02150-1). We sampled soil from fields in current row-crop agriculture, active restoration of varying ages, and prairie remnants to investigate differences in fungal communities. Our goal was to learn whether fungi, like bacteria, transition to novel communities during prairie restoration, or whether post-restoration communities begin to resemble those in remnants.   

## Source Data
Raw source data were output from QIIME2. ETL of these files is accomplished in [process_data.md](process_data.md), which produces clean .csv files which are included in this here and available for use. 

## Site Data
The file [site_locations.md](site_locations.md) shows locations of regions and sites and displays associated metadata. Climate data (precipitation) are also downloaded, processed, and displayed by this script. 

Sites are tested for autocorrelation between geographic distance and fungal species or soil chemical properties in [spatial_correlation.md](spatial_correlation.md). Mild or insignificant autocorrelation was found, particularly in the Blue Mounds restored fields, suggesting that we can at least present these fields as a pseudochronosequence. 

### Additional Site Data
Site data include experimental metadata as described above. It also includes some measured data that is aggregated and the field 
level rather than containing subsamples in fields:

- [Soil abiotic properties](soil_properties.md): This script provides a quick overview of the soil abiotic property data and produces products (e.g., ordination axis values) for use in downstream analysis. A basic correlation with microbial biomass 
is also displayed. 
- [Percent water stable aggregates](soil_wsa.md): This script provides a quick overview of WSA in fields and
regions. WSA is higher in restored fields based on a mixed linear model, but isn't correlated with years 
since restoration. 
- [Microbial biomass](microbial_biomass.md) data include site-species tables derived from high-throughput sequencing and PLFA/NLFA extractions 
and analysis data which Ylva did. This report presents basic visualizations of microbial biomass inferred with PLFA/NLFA 
quantification.

## Diagnostics
In iterative workflow is used to discover the optimal number of samples to keep from each field to ensure equal 
sampling effort and an adequate representation of diversity. For details, see the associated scripts listed here.

### Workflow
1. The script `process_data.R` is run first.
1. Next, `microbial_diagnostics_pre.R` is run to investigate sequencing depth in samples 
and species accumulation in fields.
1. Then, `process_data.R` is run again, this time with the number of samples retained
per field set to the levels recommended in `microbial_diagnostics_pre.R`. 
1. Finally, `microbial_diagnostics_post.R` is run. It is very similar to the "_pre" script,
but a different file is used so that the two may be compared. 

## Microbial Species Diversity
Microbial data analyzed here include site-species tables derived from high-throughput sequencing of ITS and 18S genes 
and clustering into 97% similar OTUs.
The report [microbial_diversity.md](microbial_diversity.md) presents basic statistics and visualizations of species richness, Shannon's diversity/evenness, and Simpson's diversity/evenness in the microbial species data across field types. 

## Microbial Taxonomy and Guilds
OTUs are annotated with taxonomic information, and ITS OTUs additionally are annotated with [fungal traits](https://link.springer.com/article/10.1007/s13225-020-00466-2) "primary lifestyles" (aka guilds). 
In [this report](microbial_guild_taxonomy.md), sequence abundances
in taxonomic groups or fungal guilds are compared across field types and with time
since restoration. Indicator species analysis is also performed to identify particular species matches 
with interesting stories. 

## Microbial Communities
In [this report](microbial_communities.md), multivariate analysis is performed on three datasets: ITS (rarefied, Bray-Curtis distance), 18S (rarefied,
Bray-Curtis distance), and 18S (rarefied, UNIFRAC distance). Unconstrained ordinations are produced. Cornfields clustered away from all other field types, but separation of remnants and restored fields differed among datasets. With ITS, years since restoration looked like a strong signal among restored fields, but remnants appeared intermediate among restored fields where age is concerned. This could be evidence of a "novel microbial assemblage" as a successional endpoint in restored fields (although other scenarios are equally plausible). With both 18S datasets, remnant and restored fields clustered closer together, and away from corn. Years since restoration was still a significant signal, but less so than with ITS. 

Years since will also be tested **post-hoc** using `envfit()`. 

## Plant Composition
Plant data comprises two data sets. Quadrat surveys from 10 plots in each field were done in Wisconsin sites 
in 2016. At Fermi, we only have relevé data with presence/absence from summer 2017. Plant metadata includes 
taxonomy and life history traits, and should coverboth plant data sets. With the abundance-data sites, trait 
data are reported in percent cover. In the presence-data sites, traits are in counts of species with that trait per field. 

This [report](plant.md) produces basic visualization and diagnostic views of the plant data. Two traits matrices
are output, one for sites with abundance data (16 sites), and one for sites with presence data only (20 sites). 

## Constrained Analyses
This [report](tgr_constrained.md) attempts to find strong explanations for microbial community differences. Because the site differences 
due to region might confound such an analysis, site variables like soil abiotic properties and precipitation were first condensed into two ordination axes 
using PCA and used as covariables in the constrained analysis. The constrained analysis was done with dbRDA and tested the explanatory 
power of years since restoration, plant community axes or plant traits, and soil properties that we don't believe are due to 
regional differences and rather come from an agricultural legacy: SOM, NO<sub>3</sub>, P, and K.  

## Summary
