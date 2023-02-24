# Introduction
Agricultural practices in cropland, like soil tillage and fertilizer applications, profoundly alter the soil's physical and chemical properties over time. Soil carbon is greatly reduced and belowground metabolic pathways are greatly simplified, leading to attrition of microbial diversity. Restoration of abandoned cropland is increasingly popular, and usually involves seeding a great density of native plants into post-agricultural soil. Although these plants can establish, bacterial communities often fail to recover to a pre-agricultural state, instead forming novel communities [(Badger and Docherty 2022)](https://link.springer.com/10.1007/s00248-022-02150-1). We sampled soil from fields in current row-crop agriculture, active restoration of varying ages, and prairie remnants to investigate changes in fungal communities. Our goal was to learn whether fungi, like bacteria, transition to novel communities during prairie restoration, or whether post-restoration communities begin to resemble those in remnants.   

## Source Data
Raw source data were output from QIIME2. ETL of these files is accomplished in [process_data.md](process_data.md), which produces clean .csv files which are included in this here and available for use. 

## Site Data
The file [site_locations.md](site_locations.md) shows locations of regions and sites and displays associated metadata. Climate data (precipitation) are also downloaded, processed, and displayed by this script. 

### Additional Site Data
The site data also include several tables which aren't wrangled in scripts, including:

- [Soil chemical properties](clean_data/soil.csv)

Sites are tested for autocorrelation between geographic distance and fungal species or soil chemical properties in [spatial_correlation.md](spatial_correlation.md). Mild or insignificant autocorrelation was found, particularly in the Blue Mounds restored fields, suggesting that we can at least present these fields as a pseudochronosequence. 

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
