#' ---
#' title: "Supplemental: microbial communities"
#' author: "Beau Larkin\n"
#' date: "Last updated: `r format(Sys.time(), '%d %B, %Y')`"
#' output:
#'   github_document:
#'     toc: true
#'     toc_depth: 3
#'     df_print: paged
#'     fig_width: 8
#'     fig_height: 7
#' ---
#'
#' # Description
#' Functions for the related script are stored and executed here. 
#' 
#' ## PCoA function
#+ pcoa_function
pcoa_fun <- function(s, d, env=sites, corr="none", adonis_index="bray", df_name, nperm=1999) {
    set.seed <- 397
    # Multivariate analysis
    p <- pcoa(d, correction = corr)
    p_vals <- data.frame(p$values) %>% 
        rownames_to_column(var = "Dim") %>% 
        mutate(Dim = as.integer(Dim))
    p_vec <- data.frame(p$vectors)
    # Wrangle site data
    env <- data.frame(env)
    # Global permutation test (PERMANOVA)
    gl_permtest <-
        with(env,
             adonis2(
                 d ~ field_type,
                 data = env,
                 permutations = nperm,
                 method = adonis_index,
                 add = if (corr == "none") FALSE else "lingoes",
                 strata = region
             ))
    # Pairwise post-hoc test
    group_var <- as.character(env$field_type)
    groups <- as.data.frame(t(combn(unique(group_var), m = 2)))
    contrasts <- data.frame(
        group1 = groups$V1,
        group2 = groups$V2,
        R2 = NA,
        F_value = NA,
        df1 = NA,
        df2 = NA,
        p_value = NA
    )
    for (i in seq(nrow(contrasts))) {
        group_subset <- group_var == contrasts$group1[i] | group_var == contrasts$group2[i]
        contrast_matrix <- data.frame(s[group_subset, ], row.names = 1)
        fit <- with(env[group_subset, ],
                    adonis2(
                        contrast_matrix ~ group_var[group_subset],
                        permutations = nperm,
                        method = adonis_index,
                        add = if (corr == "none") FALSE else "lingoes",
                        strata = region
                    ))
        
        contrasts$R2[i] <- round(fit$R2[1], digits = 3)
        contrasts$F_value[i] <- round(fit$F[1], digits = 3)
        contrasts$df1[i] <- fit$Df[1]
        contrasts$df2[i] <- fit$Df[2]
        contrasts$p_value[i] <- fit$`Pr(>F)`[1]
    }
    contrasts$p_value_adj <- p.adjust(contrasts$p_value, method = "fdr")
    # Diagnostic plots
    if(corr == "none" | ncol(p_vals) == 6) {
        p_bstick <- ggplot(p_vals, aes(x = factor(Dim), y = Relative_eig)) + 
            geom_col(fill = "gray70", color = "gray30") + 
            geom_line(aes(x = Dim, y = Broken_stick), color = "red") +
            geom_point(aes(x = Dim, y = Broken_stick), color = "red") +
            labs(x = "Dimension", 
                 title = paste0("PCoA Eigenvalues and Broken Stick Model (", df_name, ")")) +
            theme_bw()
        p_ncomp <- with(p_vals, which(Relative_eig < Broken_stick)[1]-1)
        eig <- round(p_vals$Relative_eig[1:2] * 100, 1)
    } else {
        p_bstick <- ggplot(p_vals, aes(x = factor(Dim), y = Rel_corr_eig)) + 
            geom_col(fill = "gray70", color = "gray30") + 
            geom_line(aes(x = Dim, y = Broken_stick), color = "red") +
            geom_point(aes(x = Dim, y = Broken_stick), color = "red") +
            labs(x = "Dimension", 
                 title = paste0("PCoA Eigenvalues and Broken Stick Model (", df_name, ")")) +
            theme_bw()
        p_ncomp <- with(p_vals, which(Rel_corr_eig < Broken_stick)[1]-1)
        eig <- round(p_vals$Rel_corr_eig[1:2] * 100, 1)
    }
    ncomp <- if(p_ncomp <= 2) {2} else {p_ncomp}
    # Ordination plot
    scores <-
        p_vec[, 1:ncomp] %>%
        rownames_to_column(var = "field_key") %>%
        mutate(field_key = as.integer(field_key)) %>%
        left_join(sites, by = "field_key") %>% 
        select(-field_name)
    # Output data
    output <- list(dataset                        = df_name,
                   components_exceed_broken_stick = p_ncomp,
                   correction_note                = p$note,
                   values                         = p_vals[1:(ncomp+1), ], 
                   eigenvalues                    = eig,
                   site_vectors                   = scores,
                   broken_stick_plot              = p_bstick,
                   permanova                      = gl_permtest,
                   pairwise_contrasts             = kable(contrasts, format = "pandoc"))
    return(output)
}
#' 
#' ## PCoA function with subsamples
#' Subsamples from fields are included here. 
#+ pcoa_samps_function
pcoa_samps_fun <- function(s, d, env=sites, corr="none", adonis_index="bray", df_name, nperm=1999) {
    set.seed <- 438
    # Multivariate analysis
    p <- pcoa(d, correction = corr)
    p_vals <- data.frame(p$values) %>% 
        rownames_to_column(var = "Dim") %>% 
        mutate(Dim = as.integer(Dim))
    p_vec <- data.frame(p$vectors)
    # Wrangle site data
    env_w <- env %>% 
        left_join(s %>% select(field_key, sample), by = join_by(field_key), multiple = "all") %>% 
        mutate(field_sample = paste(field_key, sample, sep = "_")) %>% 
        column_to_rownames(var = "field_sample") %>% 
        as.data.frame()
    # Permutation tests (PERMANOVA)
    # Fields as replicate strata with subsamples
    # Regions as blocks
    # Global test
    gl_perm_design <- with(env_w, 
                           how(within = Within(type="none"), 
                               plots  = Plots(strata=field_key, type="free"),
                               blocks = region,
                               nperm  = nperm))
    gl_permtest <- adonis2(
        d ~ field_type,
        data = env_w,
        permutations = gl_perm_design,
        method = adonis_index)
    # Pairwise post-hoc test
    group_var <- as.character(env_w$field_type)
    groups <- as.data.frame(t(combn(unique(group_var), m = 2)))
    contrasts <- data.frame(
        group1 = groups$V1,
        group2 = groups$V2,
        R2 = NA,
        F_value = NA,
        df1 = NA,
        df2 = NA,
        p_value = NA
    )
    for (i in seq(nrow(contrasts))) {
        group_subset <-
            group_var == contrasts$group1[i] |
            group_var == contrasts$group2[i]
        contrast_matrix <- s[group_subset, ]
        pw_perm_design <- with(env_w[group_subset,],
                               how(
                                   within = Within(type = "none"),
                                   plots  = Plots(strata = field_key, type = "free"),
                                   blocks = region,
                                   nperm  = nperm
                               ))
        fit <- adonis2(
            contrast_matrix ~ group_var[group_subset],
            add = if (corr == "none") FALSE else "lingoes",
            permutations = pw_perm_design,
            method = adonis_index
        )
        
        contrasts$R2[i] <- round(fit$R2[1], digits = 3)
        contrasts$F_value[i] <- round(fit$F[1], digits = 3)
        contrasts$df1[i] <- fit$Df[1]
        contrasts$df2[i] <- fit$Df[2]
        contrasts$p_value[i] <- fit$`Pr(>F)`[1]
    }
    contrasts$p_value_adj <- p.adjust(contrasts$p_value, method = "fdr")
    # Diagnostic plots
    if(corr == "none" | ncol(p_vals) == 6) {
        p_bstick <- ggplot(p_vals, aes(x = factor(Dim), y = Relative_eig)) + 
            geom_col(fill = "gray70", color = "gray30") + 
            geom_line(aes(x = Dim, y = Broken_stick), color = "red") +
            geom_point(aes(x = Dim, y = Broken_stick), color = "red") +
            labs(x = "Dimension", 
                 title = paste0("PCoA Eigenvalues and Broken Stick Model (", df_name, ")")) +
            theme_bw()
        p_ncomp <- with(p_vals, which(Relative_eig < Broken_stick)[1]-1)
        eig <- round(p_vals$Relative_eig[1:2] * 100, 1)
    } else {
        p_bstick <- ggplot(p_vals, aes(x = factor(Dim), y = Rel_corr_eig)) + 
            geom_col(fill = "gray70", color = "gray30") + 
            geom_line(aes(x = Dim, y = Broken_stick), color = "red") +
            geom_point(aes(x = Dim, y = Broken_stick), color = "red") +
            labs(x = "Dimension", 
                 title = paste0("PCoA Eigenvalues and Broken Stick Model (", df_name, ")")) +
            theme_bw()
        p_ncomp <- with(p_vals, which(Rel_corr_eig < Broken_stick)[1]-1)
        eig <- round(p_vals$Rel_corr_eig[1:2] * 100, 1)
    }
    ncomp <- if(p_ncomp <= 2) {2} else {p_ncomp}
    # Ordination plot
    scores <-
        p_vec[, 1:ncomp] %>%
        rownames_to_column(var = "field_sample") %>%
        separate_wider_delim(field_sample, delim = "_", names = c("field_key", "sample_key"), cols_remove = TRUE) %>% 
        mutate(field_key = as.integer(field_key)) %>%
        left_join(env, by = join_by(field_key))
    # Output data
    output <- list(dataset                        = df_name,
                   components_exceed_broken_stick = p_ncomp,
                   correction_note                = p$note,
                   values                         = p_vals[1:(ncomp+1), ], 
                   eigenvalues                    = eig,
                   site_vectors                   = scores,
                   broken_stick_plot              = p_bstick,
                   permanova                      = gl_permtest,
                   pairwise_contrasts             = contrasts, format = "pandoc")
    return(output)
}
#' 
#' ## PCoA function with subsamples, Blue Mounds only
#+ pcoa_samps_bm_function
pcoa_samps_bm_fun <- function(s, d, env=sites_resto_bm, corr="none", adonis_index="bray", df_name, nperm=1999) {
    set.seed <- 845
    # Multivariate analysis
    p <- pcoa(d, correction = corr)
    p_vals <- data.frame(p$values) %>% 
        rownames_to_column(var = "Dim") %>% 
        mutate(Dim = as.integer(Dim))
    p_vec <- data.frame(p$vectors)
    # Wrangle site data
    env_w <- env %>% 
        left_join(s %>% select(field_key, sample), by = join_by(field_key), multiple = "all") %>% 
        mutate(field_sample = paste(field_key, sample, sep = "_")) %>% 
        column_to_rownames(var = "field_sample")
    # Permutation test (PERMANOVA)
    p_permtest <- adonis2(d ~ field_key, data = env_w, permutations = nperm, method = adonis_index)
    # Diagnostic plots
    if(corr == "none" | ncol(p_vals) == 6) {
        p_bstick <- ggplot(p_vals, aes(x = factor(Dim), y = Relative_eig)) + 
            geom_col(fill = "gray70", color = "gray30") + 
            geom_line(aes(x = Dim, y = Broken_stick), color = "red") +
            geom_point(aes(x = Dim, y = Broken_stick), color = "red") +
            labs(x = "Dimension", 
                 title = paste0("PCoA Eigenvalues and Broken Stick Model (", df_name, ")")) +
            theme_bw()
        p_ncomp <- with(p_vals, which(Relative_eig < Broken_stick)[1]-1)
        eig <- round(p_vals$Relative_eig[1:2] * 100, 1)
    } else {
        p_bstick <- ggplot(p_vals, aes(x = factor(Dim), y = Rel_corr_eig)) + 
            geom_col(fill = "gray70", color = "gray30") + 
            geom_line(aes(x = Dim, y = Broken_stick), color = "red") +
            geom_point(aes(x = Dim, y = Broken_stick), color = "red") +
            labs(x = "Dimension", 
                 title = paste0("PCoA Eigenvalues and Broken Stick Model (", df_name, ")")) +
            theme_bw()
        p_ncomp <- with(p_vals, which(Rel_corr_eig < Broken_stick)[1]-1)
        eig <- round(p_vals$Rel_corr_eig[1:2] * 100, 1)
    }
    ncomp <- if(p_ncomp <= 2) {2} else {p_ncomp}
    # Permutation test (ENVFIT)
    # Fields as strata
    h = with(env_w, 
             how(within = Within(type="none"), 
                 plots = Plots(strata=field_key, type="free"), 
                 nperm = nperm))
    fit <- envfit(p_vec ~ yr_since, env_w, permutations = h, choices = c(1:ncomp))
    fit_sc <- scores(fit, c("vectors"))
    # Ordination plotting data
    scores <-
        p_vec[, 1:ncomp] %>%
        rownames_to_column(var = "field_sample") %>%
        separate_wider_delim(field_sample, delim = "_", names = c("field_key", "sample_key"), cols_remove = TRUE) %>% 
        mutate(field_key = as.integer(field_key)) %>%
        left_join(env, by = join_by(field_key)) %>% 
        select(-field_type)
    # Output data
    output <- list(dataset                        = df_name,
                   components_exceed_broken_stick = p_ncomp,
                   correction_note                = p$note,
                   values                         = p_vals[1:(ncomp+1), ], 
                   eigenvalues                    = eig,
                   site_vectors                   = scores,
                   broken_stick_plot              = p_bstick,
                   permanova                      = p_permtest,
                   vector_fit                     = fit,
                   vector_fit_scores              = fit_sc)
    return(output)
}