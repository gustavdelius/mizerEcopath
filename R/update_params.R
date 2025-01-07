#' Function that updates params object with new parameter values
#'
#' This function currently updates the gear selectivity parameters `l50` and
#' `l25, the catchability, the steepness `U` of the maturity ogive and the
#' coefficient `M` of the power-law mortality rate. It then recalculates the
#' steady state of the model and rescales it to match the observed biomass.
#'
#' @param params The MizerParams object to update
#' @param species The species to update. By default the first species in the
#'   model.
#' @param pars A named numeric vector of parameter values
#' @param data The list with data as passed to the objective function
#' @param w_select A logical vector indicating which weight bins were used in
#'   the likelihood calculation
#' @return The updated MizerParams object
#' @export
update_params <- function(params, species = 1, pars, data, w_select) {
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

    gp <- params@gear_params
    gp_select <- gp$species == species
    gps <- gp[gp_select, ]
    if (nrow(gps) > 1) {
        stop("The code currently assumes that there is only a single gear for each species.")
    }

    # Update the gear parameters
    if (data$use_counts) {
        gps$l50 <- pars["l50"]
        gps$l25 <- pars["ratio"] * pars["l50"]
    }
    if (data$yield_lambda > 0) {
        gps$catchability <- pars["catchability"]
    }
    gear_params(params)[gp_select, ] <- gps

    # recalculate the power-law mortality rate
    sps$mu_mat <- pars["mu_mat"]
    # Note that `mu_mat` is the mortality at the w just below w_mat
    mat_idx <- sum(params@w < sps$w_mat)
    w_mat <- params@w[mat_idx]
    ext_mort(params)[sp_select, ] <-
        pars["mu_mat"] * (params@w / w_mat)^sps$d

    params@species_params[sp_select, ] <- sps
    params <- setReproduction(params)

    # Calculate the new steady state ----
    params <- steadySingleSpecies(params)
    # Rescale it to get the observed biomass
    total <- sum(params@initial_n[sp_select, w_select] *
                     params@w[w_select] * params@dw[w_select])
    factor <- data$biomass / total
    params@initial_n[sp_select, ] <- params@initial_n[sp_select, ] * factor
    params <- setBevertonHolt(params, reproduction_level = 0)

    return(params)
}
