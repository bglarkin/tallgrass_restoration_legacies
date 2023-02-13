# Microbial data: community differences
# Tallgrass prairie
# 2022-01-06, BL 


# Description ----------------------
# ——————————————————————————————————
# Microbial data include site-species tables derived from high-throughput sequencing 
# and PLFA/NLFA data which Ylva did. 
# 
# This presents basic visualizations of community differences among sites/regions
# based on ITS data. 
# 
# Multivariate analyses/visualizations proposed:
# - UMAP of most important/numerous guilds (in 2d or 3d)
# - Explore/refine environment data (soilchem, site metadata)
# - Possible constrained ordinations (LDA, db-RDA, etc.)
# - RLQ/Fourth Corner of species, traits, and env
# - Variation Partitioning of species vs. plants, and species vs. soilchem
# 
# - Is there an opportunity to invent something new, like a machine learning
#   version of constrained ordination (UMAP, then clustering, then test set
#   of data relating success of assignment to clusters/env variables?)
# 
# HDBSCAN clustering will be attempted on UMAP and PCoA ordination data
# Using UMAP for clustering, recommendations taken from [Using UMAP for Clustering](https://umap-learn.readthedocs.io/en/latest/clustering.html#traditional-clustering)
# The idea is to make very tightly packed clusters and increase density gradients
# - Increase the number of neighbors (but so far I found the fewest noise points with nn = 5)
# - Number of components can increase because visualization isn't the primary goal. But >10 produced more noise and fewer clusters.
# - keep min_dist very low


# Resources ------------------------
# ——————————————————————————————————
packages_needed = c("tidyverse", "uwot", "vegan", "rgl", "colorspace", "ape", "dbscan")
packages_installed = packages_needed %in% rownames(installed.packages())

if (any(! packages_installed))
    install.packages(packages_needed[! packages_installed])
for (i in 1:length(packages_needed)) {
    library(packages_needed[i], character.only = T)
} 


# Data -----------------------------
# ——————————————————————————————————
# Sites-species tables
its_otu_savg <- read_csv(paste0(getwd(), "/clean_data/spe_ITS_otu_siteSpeMatrix_avg.csv")) # "savg" is site averaged subsamples
its_otu_all  <- read_csv(paste0(getwd(), "/clean_data/spe_ITS_otu_siteSpeMatrix_allReps.csv")) %>% # "all" is OTU counts in site subsamples
    mutate(sample_key = paste(site_key, sample, sep = "_")) %>% 
    select(site_key, sample, sample_key, everything())
its_sv  <- read_csv(paste0(getwd(), "/clean_data/spe_ITS_sv_siteSpeMatrix_avg.csv"))
# Site metadata and design
# Age of remnants must be assigned to facilitate plotting/modelling
rem_age <- 50
sites   <- read_csv(paste0(getwd(), "/clean_data/site.csv")) %>% 
    filter(site_type != "oldfield") %>% 
    mutate(site_type = factor(site_type, ordered = TRUE, levels = c("corn", "restored", "remnant")),
           yr_since = replace(yr_since, which(site_type == "remnant"), rem_age)) %>% 
    glimpse()
# Unified spe-meta-guild-taxonomy tables needed to facilitate filtering of species prior to ordination
# These tables were produced in microbial_diversity.R and written to csv
spe_meta_its_otu <- read_csv(paste0(getwd(), "/clean_data/speMeta_ITS_otu.csv")) %>% glimpse()
spe_meta_its_sv  <- read_csv(paste0(getwd(), "/clean_data/speMeta_ITS_sv.csv")) %>% glimpse()


# Ordinations ----------------------
# Set method for calculating distance in ordinations. Bray-Curtis or Ruzicka are both
# appropriate, but Bray-Curtis has produced axes with better explanatory power.
# Ruzicka is used with method="jaccard"
d_method <- "bray"
# ——————————————————————————————————
# Functions for UMAP and PCoA will simplify code and reduce errors
umap_fun <- function(data, nn, nc, md, sp=1, fast=TRUE, init = "spectral"){
    u <- umap(data, n_neighbors = nn, n_components = nc, min_dist = md, spread = sp,
              init = init, verbose = FALSE, ret_nn = TRUE, fast_sgd = fast)
    u_locs_sites <- data.frame(u$embedding) 
    output <- list(locations = u_locs_sites,
                   nn = u$nn)
    return(output)
}
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
# UMAP
u_otu_savg <- umap_fun(its_otu_savg_bc, nn = 12, nc = 2, md = 0.01)
u_plot_savg_data <- 
    u_otu_savg$locations %>% 
    rownames_to_column(var = "site_key") %>% 
    mutate(site_key = as.numeric(site_key)) %>% 
    left_join(sites %>% select(-lat, -long, -yr_restore), by = "site_key")
ggplot(u_plot_savg_data, aes(x = X1, y = X2)) +
    geom_point(aes(fill = site_type, shape = region), size = 10) + 
    geom_text(aes(label = yr_since)) +
    scale_fill_discrete_qualitative(palette = "harmonic") +
    scale_shape_manual(values = c(21, 22, 23, 24)) +
    labs(x = "Axis 1", y = "Axis 2", title = "UMAP of site-averaged 97% OTUs", caption = "Text indicates years since restoration") +
    theme_bw() +
    guides(fill = guide_legend(override.aes = list(shape = 21)))
# ____UMAP and HDBSCAN unsupervised clustering --------------------------
u_otu_savg_cl <- umap_fun(its_otu_savg_bc, nn = 4, nc = 10, md = 0, fast = FALSE)
(cl_otu_savg_cl <- hdbscan(u_otu_savg_cl$locations, minPts = 2))
ggplot(u_plot_savg_data, aes(x = X1, y = X2)) +
    geom_polygon(
        data = u_plot_savg_data %>% 
            mutate(cluster = cl_otu_savg_cl$cluster) %>% 
            filter(cluster != 0) %>% 
            group_by(cluster) %>% 
            slice(chull(X1, X2)),
        aes(group = cluster),
        alpha = 0.1,
        color = "gray20",
        size = 0.2
    ) +
    geom_point(
        data = u_plot_savg_data %>% 
            mutate(cl_prob = cl_otu_savg_cl$membership_prob / max(cl_otu_savg_cl$membership_prob)),
        aes(fill = site_type, shape = region, alpha = cl_prob), size = 10) + 
    geom_text(aes(label = yr_since)) +
    scale_fill_discrete_qualitative(palette = "harmonic") +
    scale_shape_manual(values = c(21, 22, 23, 24)) +
    labs(x = "Axis 1", y = "Axis 2", 
         title = "UMAP of site-averaged 97% OTUs with HDBSCAN clustering",
         caption = "Polygons surround clusters, shading indicates cluster membership probability") +
    theme_bw() +
    guides(fill = guide_legend(override.aes = list(shape = 21)))

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
# UMAP
u_otu_all <- umap_fun(its_otu_all_bc, nn = 5, nc = 2, md = 0.4, sp = 10, init = "spca")
u_otu_plot_data <- 
    u_otu_all$locations %>% 
    rownames_to_column(var = "sample_key") %>% 
    separate(sample_key, into = c("site_key", "sample"), sep = "_", remove = FALSE) %>% 
    mutate(site_key = as.numeric(site_key)) %>% 
    left_join(sites %>% select(-lat, -long, -yr_restore), by = "site_key")
ggplot(u_otu_plot_data, aes(x = X1, y = X2)) +
    geom_point(aes(fill = site_type, shape = region), size = 5) + 
    # geom_text(aes(label = yr_since)) +
    scale_fill_discrete_qualitative(palette = "harmonic") +
    scale_shape_manual(values = c(21, 22, 23, 24)) +
    labs(x = "Axis 1", y = "Axis 2", title = "UMAP of site-averaged 97% OTUs", caption = "Text indicates years since restoration") +
    theme_bw() +
    guides(fill = guide_legend(override.aes = list(shape = 21)))
# Attempt to place convex hulls around sites
ggplot(u_otu_plot_data, aes(x = X1, y = X2)) +
    geom_polygon(
        data = u_otu_plot_data %>% group_by(site_name) %>% slice(chull(X1, X2)),
        aes(group = site_name),
        alpha = 0.1,
        color = "black",
        size = 0.3
    ) +
    geom_point(aes(fill = site_type, shape = region), size = 5) +
    scale_fill_discrete_qualitative(palette = "harmonic") +
    scale_shape_manual(values = c(21, 22, 23, 24)) +
    labs(
        x = "Axis 1",
        y = "Axis 2",
        title = "UMAP of site-averaged 97% OTUs",
        caption = "Polygons surround sites"
    ) +
    theme_bw() +
    guides(fill = guide_legend(override.aes = list(shape = 21)))
# ____UMAP and HDBSCAN unsupervised clustering --------------------------
# u_otu_all_10 <- umap_fun(its_otu_all_bc, nn = 5, nc = 6, md = 0, fast = FALSE) # Produces nearly as many clusters as sites, but usually has low noise
u_otu_all_10 <- umap_fun(its_otu_all_bc, nn = 8, nc = 12, md = 0, sp = 1, fast = FALSE, init = "spectral") # 8 nn because that's the minimum subsamples, 12 dimensions because that's half the number of sites
(cl_otu_all_10 <- hdbscan(u_otu_all_10$locations, minPts = 8))
ggplot(u_otu_plot_data, aes(x = X1, y = X2)) +
    geom_polygon(
        data = u_otu_plot_data %>% 
            mutate(cluster = cl_otu_all_10$cluster) %>% 
            filter(cluster != 0) %>% 
            group_by(cluster) %>% 
            slice(chull(X1, X2)),
        aes(group = cluster),
        alpha = 0.1,
        color = "gray20",
        size = 0.2
    ) +
    geom_point(
        data = u_otu_plot_data %>% 
            mutate(cl_prob = cl_otu_all_10$membership_prob / max(cl_otu_all_10$membership_prob)),
        aes(fill = site_type, shape = region, alpha = cl_prob), size = 5) +
    scale_fill_discrete_qualitative(palette = "harmonic") +
    scale_shape_manual(values = c(21, 22, 23, 24)) +
    # geom_text(aes(label = yr_since)) +
    labs(
        x = "Axis 1",
        y = "Axis 2",
        title = "UMAP of site-averaged 97% OTUs with HDBSCAN clustering",
        caption = "Polygons surround clusters, shading indicates cluster membership probability"
    ) +
    theme_bw() +
    guides(fill = guide_legend(override.aes = list(shape = 21)))

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
# ____PCoA and HDBSCAN unsupervised clustering ------------------
# The first 9 components exceed the broken stick model with Ruzicka distance. 
# The first 14 components exceed the broken stick with Bray-Curtis distance.
pcoa_otu_all_cl <- pcoa_fun(its_otu_all_bc, ax = 14, corr = "lingoes") 
(cl_otu_all_pcoa <- hdbscan(pcoa_otu_all_cl$site_vectors, minPts = 3))
ggplot(pcoa_otu_all_site, aes(x = Axis.1, y = Axis.2)) +
    geom_polygon(
        data = pcoa_otu_all_site %>% 
            mutate(cluster = cl_otu_all_pcoa$cluster) %>% 
            filter(cluster != 0) %>% 
            group_by(cluster) %>% 
            slice(chull(Axis.1, Axis.2)),
        aes(group = cluster),
        alpha = 0.1,
        color = "gray20",
        size = 0.2
    ) +
    geom_point(
        data = pcoa_otu_all_site %>% 
            mutate(cl_prob = cl_otu_all_pcoa$membership_prob / max(cl_otu_all_pcoa$membership_prob)),
        aes(fill = site_type, shape = region, alpha = cl_prob), size = 5) + 
    geom_text(aes(label = yr_since)) +
    scale_fill_discrete_qualitative(palette = "harmonic") +
    scale_shape_manual(values = c(21, 22, 23, 24)) +
    labs(x = "Axis 1 (7.7%)", y = "Axis 2 (4.7%)", 
         title = "PCoA of 97% OTUs with HDBSCAN clustering", 
         caption = "Polygons surround clusters, shading indicates cluster membership probability") +
    theme_bw() +
    guides(fill = guide_legend(override.aes = list(shape = 21)))

