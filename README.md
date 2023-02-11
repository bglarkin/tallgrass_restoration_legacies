# Introduction
Agricultural practices in cropland, like soil tillage and fertilizer applications, profoundly alter the soil's physical and chemical properties over time. Soil carbon is greatly reduced and belowground metabolic pathways are greatly simplified, leading to attrition of microbial diversity. Restoration of abandoned cropland is increasingly popular, and usually involves seeding a great density of native plants into post-agricultural soil. Although these plants can establish, bacterial communities often fail to recover to a pre-agricultural state, instead forming novel communities [(Badger and Docherty 2022)](https://link.springer.com/10.1007/s00248-022-02150-1). We sampled soil from fields in current row-crop agriculture, active restoration of varying ages, and prairie remnants to investigate changes in fungal communities. Our goal was to learn whether fungi, like bacteria, transition to novel communities during prairie restoration, or whether post-restoration communities begin to resemble those in remnants.   

## Source Data
Raw source data were output from QIIME. These are relatively large files, some of which are binaries (e.g., .xlsx) and not compatible with GitHub. To save space and avoid compatibility complications, they are not included in this repository. ETL of these files is accomplished in [process_data.md](process_data.md), which produces clean .csv files which are included in this here and available for use. **Resampled averages of sequence abundances are also produced here, and are used throughout all analyses.** Some samples failed to amplify or were cut after rarefaction, leading to difference in sampling effort across fields. Samples in fields were randomly resampled to the minimum and then averaged within fields to make comparisons legitimate. 

## Site Data
The file [site_locations.md](site_locations.md) shows locations of regions and sites and displays associated metadata. Climate data (precipitation) are also downloaded, processed, and displayed by this script. 

## Microbial Species Diversity
Microbial data analyzed here include site-species tables derived from high-throughput sequencing of ITS and 18S genes 
and clustering into 97% similar OTUs and 100% similar SVs.
The report [microbial_diversity.md](microbial_diversity.md) presents basic statistics and visualizations of species richness, Shannon's diversity/evenness, and Simpson's diversity/evenness in the microbial species data across field types. 

## Microbial Taxonomy and Guilds
Sequence clusters identified in QIIME are annotated with taxonomic information and
metadata from [FunGuild](https://github.com/UMNFuN/FUNGuild). In [this report](microbial_guild_taxonomy.md), sequence abundances
in taxonomic groups or fungal guilds are compared across field types and with time
since restoration.
