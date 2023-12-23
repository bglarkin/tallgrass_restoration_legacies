#' ---
#' title: "Supplemental: constrained analyses"
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
#' ## PCoA axes
#' Produce community axes for use in constrained analyses later. 
#+ pcoa_function
pcoa_fun <- function(s, ft=c("restored"), rg=c("BM"), method="bray", binary=FALSE, corr="none") {
    d <- vegdist(
        data.frame(
            s %>% 
                filter(field_type %in% ft, region %in% rg) %>% 
                select(-field_type, -region) %>% 
                select(field_name, where(~ is.numeric(.) && sum(.) > 0)),
            row.names = 1),
        method = method,
        binary = binary)
    p <- pcoa(d, correction = corr)
    p_vals <- data.frame(p$values) %>% 
        rownames_to_column(var = "Dim") %>% 
        mutate(Dim = as.integer(Dim))
    p_vec <- data.frame(p$vectors)
    if(corr == "none" | ncol(p_vals) == 6) {
        p_ncomp <- with(p_vals, which(Relative_eig < Broken_stick)[1]-1)
    } else {
        p_ncomp <- with(p_vals, which(Rel_corr_eig < Broken_stick)[1]-1)
    }
    ncomp <- if(p_ncomp <= 2) {2} else {p_ncomp}
    # Ordination plot
    scores <- p_vec[, 1:2]
    # Output data
    output <- list(correction_note = p$note,
                   values          = p_vals[1:(ncomp+1), ],
                   site_vectors    = scores)
    return(output)
}
#' 
#' ## Constrained analyses
#' This function performs a db-RDA on microbial data with soil/site covariables and 
#' forward selects on plant communities, agricultural nutrients, and soil carbon.
#+ dbrda_function
dbrda_fun <- function(s, pspe_pcoa="none", ft=c("restored"), rg=c("BM")) {
    fspe_bray <- vegdist(
        data.frame(
            s %>% 
                filter(field_type %in% ft, region %in% rg) %>% 
                select(-field_type, -region) %>% 
                select(field_name, where(~ is.numeric(.) && sum(.) > 0)),
            row.names = 1),
        method = "bray")
    if(is.data.frame(pspe_pcoa) == TRUE) {
        pspe_ax <- pspe_pcoa %>% 
            rownames_to_column("field_name") %>% 
            rename(plant1 = Axis.1, plant2 = Axis.2)
    }
    ptr_norm <- decostand(
        data.frame(
            ptr %>% 
                filter(field_type %in% ft, region %in% rg) %>% 
                select(-field_type, -region),
            row.names = 1),
        "normalize") %>% 
        rownames_to_column("field_name")
    soil_cov_z <- decostand(
        data.frame(
            soil %>% 
                filter(field_type %in% ft, region %in% rg) %>% 
                select(-field_type, -region, -OM, -NO3, -P, -K),
            row.names = 1),
        "standardize")
    soil_cov_sc <- data.frame(scores(rda(soil_cov_z), display = "sites")) %>% 
        rename(soil1 = PC1, soil2 = PC2) %>% 
        rownames_to_column("field_name")
    soil_expl_z <- decostand(
        data.frame(
            soil %>% 
                filter(field_type %in% ft, region %in% rg) %>% 
                select(field_name, OM, NO3, P, K, -field_type, -region),
            row.names = 1),
        "standardize") %>% 
        rownames_to_column("field_name")
    yr_z <- decostand(
        data.frame(
            sites %>% 
                filter(field_type %in% ft, region %in% rg) %>% 
                select(field_name, yr_since),
            row.names = 1),
        "standardize") %>% 
        rownames_to_column("field_name")
    # Create explanatory data frame and covariables matrix
    env <- if(is.data.frame(pspe_pcoa) == FALSE) {
        soil_cov_sc %>% 
            left_join(ptr_norm, by = join_by(field_name)) %>% 
            left_join(soil_expl_z, by = join_by(field_name)) %>% 
            left_join(yr_z, by = join_by(field_name)) %>% 
            column_to_rownames("field_name") %>% 
            drop_na()
    } else {
        soil_cov_sc %>% 
            left_join(pspe_ax, by = join_by(field_name)) %>% 
            left_join(soil_expl_z, by = join_by(field_name)) %>% 
            left_join(yr_z, by = join_by(field_name)) %>% 
            column_to_rownames("field_name") %>% 
            drop_na()
    }
    # Forward select on explanatory data with covariables
    mod_null <- dbrda(fspe_bray ~ 1, data = env, sqrt.dist = TRUE)
    mod_full <- dbrda(fspe_bray ~ ., data = env, sqrt.dist = TRUE)
    set.seed <- 395
    mod_step <- ordistep(mod_null, 
                         scope = formula(mod_full), 
                         direction = "forward", 
                         permutations = how(nperm = 1999), 
                         trace = FALSE)
    mod_r2   <- RsquareAdj(mod_step, permutations = 1999)
    mod_glax <- anova(mod_step, permutations = how(nperm = 1999))
    mod_inax <- anova(mod_step, by = "axis", permutations = how(nperm = 1999))
    # Produce plot data including borderline vars if possible
    mod_scor <- if(is.data.frame(pspe_pcoa) == FALSE) {
        scores(
            dbrda(fspe_bray ~ yr_since + forb + C4_grass, data = env, sqrt.dist = TRUE),
            choices = c(1,2),
            display = c("bp", "sites"), tidy = FALSE
        )
    } else {
        scores(
            dbrda(fspe_bray ~ yr_since + plant1 + OM, data = env, sqrt.dist = TRUE),
            choices = c(1,2),
            display = c("bp", "sites"), tidy = FALSE
        )
    }
    return(list(
        plot_data = mod_scor,
        select_mod = mod_step,
        R2 = mod_r2 %>% map(\(x) round(x, 2)),
        global_axis_test = mod_glax,
        individual_axis_test = mod_inax
    ))
} 
#' 
#' ## Plot results
#+ plot_dbrda_function
plot_dbrda <- function(site_sc, site_bp) {
    site_df <- site_sc %>% 
        data.frame() %>% 
        rownames_to_column(var = "field_name") %>% 
        left_join(sites, by = join_by(field_name))
    bp_df <- site_bp %>% 
        data.frame() %>% 
        rownames_to_column(var = "envvar") %>% 
        mutate(
            origin = 0,
            m = dbRDA2 / dbRDA1, 
            d = sqrt(dbRDA1^2 + dbRDA2^2), 
            dadd = sqrt((max(dbRDA1)-min(dbRDA2))^2 + (max(dbRDA2)-min(dbRDA2))^2)*0.1,
            labx = ((d+dadd)*cos(atan(m)))*(dbRDA1/abs(dbRDA1)), 
            laby = ((d+dadd)*sin(atan(m)))*(dbRDA1/abs(dbRDA1)))
    ggplot(site_df, aes(x = dbRDA1, y = dbRDA2)) +
        geom_text(aes(label = field_name)) +
        geom_segment(data = bp_df, 
                     aes(x = origin, xend = dbRDA1, y = origin, yend = dbRDA2), 
                     arrow = arrow(length = unit(3, "mm")),
                     color = "blue") +
        geom_text(data = bp_df, 
                  aes(x = labx, y = laby, label = envvar),
                  color = "blue") +
        theme_bw()
}