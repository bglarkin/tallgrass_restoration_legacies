# Introduction
Agricultural practices in cropland, like soil tillage and fertilizer applications, profoundly alter the soil's physical and chemical properties over time. Soil carbon is greatly reduced and belowground metabolic pathways are greatly simplified, leading to attrition of microbial diversity. Restoration of abandoned cropland is increasingly popular, and usually involves seeding a great density of native plants into post-agricultural soil. Although these plants can establish, bacterial communities often fail to recover to a pre-agricultural state, instead forming novel communities [(Badger and Docherty 2022)](https://link.springer.com/10.1007/s00248-022-02150-1). We sampled soil from fields in current row-crop agriculture, active restoration of varying ages, and prairie remnants to investigate differences in fungal communities. Our goal was to learn whether fungi, like bacteria, transition to novel communities during prairie restoration, or whether post-restoration communities begin to resemble those in remnants.   

## Source Data
Raw source data were output from QIIME2. ETL of these files is accomplished in [process_data.md](process_data.md), which produces clean .csv files which are included in this here and available for use. 

## Site Data
The file [site_locations.md](site_locations.md) shows locations of regions and sites and displays associated metadata. Climate data (precipitation) are also downloaded, processed, and displayed by this script. 

### Additional Site Data
The site data also include several tables which aren't wrangled in scripts, including:

- [Soil chemical properties](clean_data/soil.csv)
- Water stable aggregates
- PLFA/NLFA
- Other?

Sites are tested for autocorrelation between geographic distance and fungal species or soil chemical properties in [spatial_correlation.md](spatial_correlation.md). Mild or insignificant autocorrelation was found, particularly in the Blue Mounds restored fields, suggesting that we can at least present these fields as a pseudochronosequence. 

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
OTUs are annotated with taxonomic information, and ITS OTUs additionally are annotated with [Fungal traits](https://link.springer.com/article/10.1007/s13225-020-00466-2) "primary lifestyles" (aka guilds). 
In [this report](microbial_guild_taxonomy.md), sequence abundances
in taxonomic groups or fungal guilds are compared across field types and with time
since restoration. Indicator species analysis is also performed to identify particular species matches 
with interesting stories. 

## Microbial Communities
In [this report](microbial_communities.md), multivariate analysis is performed on three datasets: ITS (rarefied, Bray-Curtis distance), 18S (rarefied,
Bray-Curtis distance), and 18S (rarefied, UNIFRAC distance). Unconstrained ordinations are produced. Cornfields clustered away from all other field types, but separation of remnants and restored fields differed among datasets. With ITS, years since restoration looked like a strong signal among restored fields, but remnants appeared intermediate among restored fields where age is concerned. This could be evidence of a "novel microbial assemblage" as a successional endpoint in restored fields (although other scenarios are equally plausible). With both 18S datasets, remnant and restored fields clustered closer together, and away from corn. Years since restoration was still a significant signal, but less so than with ITS. 

Years since will also be tested **post-hoc** using `envfit()`. 
