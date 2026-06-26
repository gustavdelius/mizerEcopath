#' Prepare a TMB Objective Function for Optimising Model Parameters
#'
#' This function returns a list with the data to be passed to the TMB objective
#' function. The data includes the observed catch data, the model parameters,
#' and some precomputed values that are used in the likelihood calculation.
#' The main preprocessing makes sure that we have a comprehensive set of bins
#' that cover the entire size range, even though there will not be observations
#' at all sizes. Missing observations should be interpreted as a 0 count.
#'
#' @param params A MizerParams object
#' @param species The species for which to update parameters. By default
#'   the first species in the model.
#' @param pars A named list of parameter values as returned by the TMB
#'   optimiser, containing entries for gear selectivity and mortality parameters.
#' @param data The data list as returned by [prepare_data()].
#'
#' @return A MizerParams object with updated parameters.
#' @export
update_params <- function(params, species = 1, pars, data) {
    params <- validParams(params)
    sp <- species_params(params)

    # Need to add mu_mat column if it does not exist
    # TODO: change once we have introduced a standard species param for this
    mat_idx <- colSums(outer(params@w, sp$w_mat, "<"))
    mu_mat <- ext_mort(params)[cbind(seq_len(nrow(sp)), mat_idx)]
    params <- set_species_param_default(params, "mu_mat", mu_mat)

    species <- valid_species_arg(params, species, error_on_empty = TRUE)
    if (length(species) > 1) {
        stop("Only one species can be updated at a time.")
    }

    sp_select <- sp$species == species
    sps <- sp[sp_select, ]

    # Ensure the species_params table has a D_ext column (the coefficient of the
    # external-diffusion power law d(w) = D_ext * w^(n+1)). Derive it from the
    # ext_diffusion slot / legacy d_over_g for any species that lacks it, so that
    # the per-species row assignment below has a column to write into.
    if (!"D_ext" %in% names(sp)) {
        D_ext_all <- vapply(seq_len(nrow(sp)), function(i) {
            diffusion_coefficient(params, seq_len(nrow(sp)) == i, sp$n[i])
        }, numeric(1))
        params@species_params$D_ext <- D_ext_all
        sp <- species_params(params)
        sps <- sp[sp_select, ]
    }

    gp <- params@gear_params
    gp_select <- gp$species == species
    # Only the matched gears (those present in the catch data) are updated; any
    # other model gears (e.g. a survey gear) are left untouched. Align to the
    # gear order used when the data was prepared.
    matched_idx <- which(gp_select)[match(attr(data, "gears"), gp$gear[gp_select])]
    gps <- gp[matched_idx, ]

    # Update the gear parameters
    gplist <- list()
    gpnames <- c( 'logit_l50', 'log_ratio_left', 'log_l50_right_offset', 'log_ratio_right',
                  'log_catchability')

    for (i in gpnames) gplist[[i]] <- as.numeric(pars[[i]])

    l50 <- min(data$l) + (max(data$l) - min(data$l)) * plogis(gplist$logit_l50)
    l25 <- l50 * (1 - plogis(gplist$log_ratio_left))
    l50_right <- l50 + exp(gplist$log_l50_right_offset)
    l25_right <- l50_right * (1 + exp(gplist$log_ratio_right))
    catchability <- exp(gplist$log_catchability)

    gp_res <- data.frame( l50 = l50, l25 = l25, l50_right = l50_right, l25_right = l25_right, catchability = catchability)

    gps[,'l50'] <- gp_res$l50
    gps[,'l25'] <- gp_res$l25
    gps[,'l50_right'] <- ifelse(gps$sel_func=='double_sigmoid_length', gp_res$l50_right, NA)
    gps[,'l25_right'] <- ifelse(gps$sel_func=='double_sigmoid_length', gp_res$l25_right, NA)
    gps[,'catchability'] <- gp_res$catchability

    gear_params(params)[matched_idx, ] <- gps

    # Getting the "m" value
    sps$m <- pars[["m"]]

    # recalculate the power-law mortality rate
    sps$mu_mat <- pars[["mu_mat"]]
    # Note that `mu_mat` is the mortality at the w just below w_mat
    mat_idx <- sum(params@w < sps$w_mat)
    w_mat <- params@w[mat_idx]
    ext_mort(params)[sp_select, ] <-
        pars[["mu_mat"]] * (params@w / w_mat)^sps$d

    # Update the external-diffusion power law from the optimised D_ext so that
    # the steady state below (which is diffusion-aware) uses the same diffusion
    # as the TMB objective.
    sps$D_ext <- exp(pars[["log_D_ext"]])
    params@ext_diffusion[sp_select, ] <- sps$D_ext * params@w^(sps$n + 1)

    params@species_params[sp_select, ] <- sps
    params <- setReproduction(params)

    # Calculate the new steady state

    params <- mizer::steadySingleSpecies(params, species = species)
    # Rescale it to get the observed biomass
    params <- matchBiomasses(params, species = species)
    # Set the reproduction level to zero
    rl <- numeric(length(species))
    names(rl) <- species
    params <- setBevertonHolt(params, reproduction_level = rl)

    return(params)
}
