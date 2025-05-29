#' Function that updates params object with new parameter values
#'
#' This function updates the gear selectivity parameters `l50`, `l25`,
#' the descending distance `d50` from `l50` to `l50_right`, the catchability,
#' the steepness `U` of the maturity ogive and the coefficient `M` of the
#' power-law mortality rate. It then recalculates the steady state of the model
#' and rescales it to match the observed biomass.
#'
#' @param params A MizerParams object to update.
#' @param species The species to update (index, name, or logical).
#' @param pars A named numeric vector of parameter values from `matchCatch`.
#' @param data The list of data passed to the objective function.
#' @param w_select A logical vector of weight bins used in the likelihood.
#' @return The updated MizerParams object.
#' @export
update_params <- function(params, species = 1, pars, data, w_select) {
    params <- validParams(params)

    # Ensure new selectivity columns exist
    if (!"l50_right" %in% names(params@gear_params))
        params@gear_params$l50_right <- NA_real_
    if (!"l25_right" %in% names(params@gear_params))
        params@gear_params$l25_right <- NA_real_

    sp <- species_params(params)
    # Ensure mu_mat column exists
    mat_idx <- colSums(outer(params@w, sp$w_mat, "<"))
    mu_mat_default <- ext_mort(params)[cbind(seq_len(nrow(sp)), mat_idx)]
    params <- set_species_param_default(params, "mu_mat", mu_mat_default)

    # Validate single species
    species <- valid_species_arg(params, species, error_on_empty = TRUE)
    if (length(species) > 1) stop("Only one species can be updated at a time.")

    # Subset for this species
    sp_select <- sp$species == species
    sps <- sp[sp_select, ]
    gp <- params@gear_params
    gp_select <- gp$species == species
    gps <- gp[gp_select, ]
    if (nrow(gps) > 1) stop("Only one gear per species is currently supported.")

    # Update parameters from optimizer
    if (data$use_counts) {
        # ascending selectivity
        gps$l50 <- pars["l50"]
        gps$l25 <- pars["ratio"] * pars["l50"]
        # descending distance from l50
        d50 <- pars["d50"]
        gps$l50_right <- pars["l50"] + abs(d50)
        gps$l25_right <- pars["r_right"] * gps$l50_right
    }
    # catchability
    if (data$yield_lambda > 0) {
        gps$catchability <- pars["catchability"]
    }

    # Write back gear_params
    gp[gp_select, ] <- gps
    params@gear_params <- gp

    # Update external mortality at maturity
    sps$mu_mat <- pars["mu_mat"]
    mat_idx <- sum(params@w < sps$w_mat)
    w_mat <- params@w[mat_idx]
    ext_mort(params)[sp_select, ] <- pars["mu_mat"] * (params@w / w_mat)^sps$d
    params@species_params[sp_select, ] <- sps

    # Recompute reproduction, steady state, and rescale biomass
    params <- setReproduction(params)
    params <- steadySingleSpecies(params)
    total <- sum(params@initial_n[sp_select, w_select] * params@w[w_select] * params@dw[w_select])
    factor <- data$biomass / total
    params@initial_n[sp_select, ] <- params@initial_n[sp_select, ] * factor
    params <- setBevertonHolt(params, reproduction_level = 0)

    return(params)
}
