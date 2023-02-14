#' ---
#' title: "Microbial data: community differences"
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
#' Microbial data include site-species tables derived from high-throughput sequencing and 
#' clustering in QIIME by Lorinda Bullington and PLFA/NLFA data which Ylva Lekberg did. 
#' 
#' This presents basic visualizations of community differences among sites/regions
#' based on ITS data. 
#' 
#' One goal here is to see whether choosing OTU or SV clusters presents qualitatively different
#' outcomes in ordinations. We will choose one (OTUs) if they start to look similar, as they 
#' have so far.
#' 
#' Multivariate analyses/visualizations proposed:
#' - PCoA of species abundances in fields
#' - Explore/refine environment data (soilchem, site metadata)
#' - Possible constrained ordinations (LDA, db-RDA, etc.)
#' - RLQ/Fourth Corner of species, traits, and env
#' - Variation Partitioning of species vs. plants, and species vs. soilchem
#' 
#' Species distance matrices are resampled to the minimum number which successfully amplified per
#' field. This was done to equalize sampling effort. This procedure can easily be undone in 
#' the [process_data script](process_data.md)
#' 
#' # Packages and libraries
packages_needed = c("tidyverse", "vegan", "colorspace", "ape")
packages_installed = packages_needed %in% rownames(installed.packages())
#+ packages,message=FALSE
if (any(!packages_installed)) {
    install.packages(packages_needed[!packages_installed])
}
#+ libraries,message=FALSE
for (i in 1:length(packages_needed)) {
    library(packages_needed[i], character.only = T)
}
#' 
#' # Data
#' ## Sites-species tables
#' CSV files were produced in `process_data.R`
spe <- list(
    its_otu = read_csv(
        paste0(getwd(), "/clean_data/spe_ITS_otu_siteSpeMatrix_avg.csv"),
        show_col_types = FALSE
    ),
    its_sv  = read_csv(
        paste0(getwd(), "/clean_data/spe_ITS_sv_siteSpeMatrix_avg.csv"),
        show_col_types = FALSE
    ),
    amf_otu = read_csv(
        paste0(getwd(), "/clean_data/spe_18S_otu_siteSpeMatrix_avg.csv"),
        show_col_types = FALSE
    ),
    amf_sv  = read_csv(
        paste0(getwd(), "/clean_data/spe_18S_sv_siteSpeMatrix_avg.csv"),
        show_col_types = FALSE
    )
)
#' ## Sites-species tables with guild/taxonomy metadata
#' CSV files were produced in the [microbial diversity script](microbial_diversity.md)
spe_meta <- list(
    its_otu =
        read_csv(
            paste0(getwd(), "/clean_data/speGuild_ITS_otu.csv"),
            show_col_types = FALSE
        ),
    its_sv  =
        read_csv(
            paste0(getwd(), "/clean_data/speGuild_ITS_sv.csv"),
            show_col_types = FALSE
        ),
    amf_otu = 
        read_csv(
            paste0(getwd(), "/clean_data/speTaxa_18S_otu.csv"),
            show_col_types = FALSE
        ),
    amf_sv = 
        read_csv(
            paste0(getwd(), "/clean_data/speTaxa_18S_sv.csv"),
            show_col_types = FALSE
        )
)







# Ordinations ----------------------
# Set method for calculating distance in ordinations. Bray-Curtis or Ruzicka are both
# appropriate, but Bray-Curtis has produced axes with better explanatory power.
# Ruzicka is used with method="jaccard"
d_method <- "bray"
# ——————————————————————————————————

pcoa_fun <- function(data, ax, corr="none") {
    # ax = number of axes desired in the output
    p <- pcoa(data, correction = corr)
    p_vals <- data.frame(p$values) %>% 
        rownames_to_column(var = "Dim") %>% 
        mutate(Dim = as.integer(Dim))
    p_vec <- data.frame(p$vectors)
    
    if(corr == "none" | ncol(p_vals) == 6) {
        p_bstick <- ggplot(p_vals, aes(x = factor(Dim), y = Relative_eig)) + 
            geom_col(fill = "gray70", color = "gray30") + 
            geom_line(aes(x = Dim, y = Broken_stick), color = "red") +
            geom_point(aes(x = Dim, y = Broken_stick), color = "red") +
            labs(x = "Dimension", title = "PCoA Eigenvalues and Broken Stick Model") +
            theme_bw()
        p_ncomp <- with(p_vals, which(Relative_eig < Broken_stick)[1]-1)
        print(p_bstick)
    } else {
        p_bstick <- ggplot(p_vals, aes(x = factor(Dim), y = Rel_corr_eig)) + 
            geom_col(fill = "gray70", color = "gray30") + 
            geom_line(aes(x = Dim, y = Broken_stick), color = "red") +
            geom_point(aes(x = Dim, y = Broken_stick), color = "red") +
            labs(x = "Dimension", title = "PCoA Eigenvalues and Broken Stick Model") +
            theme_bw()
        p_ncomp <- with(p_vals, which(Rel_corr_eig < Broken_stick)[1]-1)
        print(p_bstick)
    }
    
    output <- list(values = p_vals[1:(ax+1), ], site_vectors = p_vec[, 1:ax],
                   comp_exceed_bs = p_ncomp)
    return(output)
}

# __Site-averaged abundances of 97% OTUs ----------------------
# Bray-Curtis used. Ruzicka may be more appropriate, but changed nothing except reducing
# the variance explained on axes. 
its_otu_savg_bc <- vegdist(data.frame(its_otu_savg, row.names = 1), method = d_method)


#PCoA
pcoa_otu_savg <- pcoa_fun(its_otu_savg_bc, ax = 2)
pcoa_otu_savg
pcoa_otu_savg_site <- 
    pcoa_otu_savg$site_vectors %>% 
    rownames_to_column(var = "site_key") %>% 
    mutate(site_key = as.integer(site_key)) %>% 
    left_join(sites, by = "site_key")
ggplot(pcoa_otu_savg_site, aes(x = Axis.1, y = Axis.2)) +
    geom_point(aes(fill = site_type, shape = region), size = 10) + 
    geom_text(aes(label = yr_since)) +
    scale_fill_discrete_qualitative(palette = "harmonic") +
    scale_shape_manual(values = c(21, 22, 23, 24)) +
    labs(x = "Axis 1 (18.9%)", y = "Axis 2 (11.6%)", title = "PCoA of site-averaged 97% OTUs", caption = "Text indicates years since restoration") +
    theme_bw() +
    guides(fill = guide_legend(override.aes = list(shape = 21)))


# __Subsample abundances of 97% OTUs ----------------------
its_otu_all_bc <- vegdist(data.frame(its_otu_all[, -c(1:2)], row.names = 1), method = d_method)


# PCoA
pcoa_otu_all <- pcoa_fun(its_otu_all_bc, ax = 2, corr = "lingoes")
pcoa_otu_all
pcoa_otu_all_site <- 
    pcoa_otu_all$site_vectors %>% 
    rownames_to_column(var = "sample_key") %>% 
    separate(sample_key, into = c("site_key", "sample"), sep = "_", remove = FALSE) %>% 
    mutate(site_key = as.numeric(site_key)) %>% 
    left_join(sites %>% select(-lat, -long, -yr_restore), by = "site_key")
ggplot(pcoa_otu_all_site, aes(x = Axis.1, y = Axis.2)) +
    geom_point(aes(fill = site_type, shape = region), size = 5) + 
    # geom_text(aes(label = yr_since)) +
    scale_fill_discrete_qualitative(palette = "harmonic") +
    scale_shape_manual(values = c(21, 22, 23, 24)) +
    labs(x = "Axis 1 (7.7%)", y = "Axis 2 (4.7%)", title = "PCoA of 97% OTUs", caption = "Text indicates years since restoration") +
    theme_bw() +
    guides(fill = guide_legend(override.aes = list(shape = 21)))
# Spiders for fun
p_spid <- cmdscale(its_otu_all_bc, k = nrow(its_otu_all)-1, eig = FALSE, add = TRUE)
plot(scores(p_spid, choices = c(1,2)), type = "n")
ordispider(p_spid, groups = factor(pcoa_otu_all_site$site_key), spiders = "centroid", label = TRUE)
points(scores(p_spid, choices = c(1,2)), col = as.numeric(factor(pcoa_otu_all_site$site_type)), pch = 16, cex = 1.3)



# Possibly show guild or taxononmy information on ordinations