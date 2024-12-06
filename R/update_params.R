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
#' @param biomass The observed biomass of the species in the range of the model
#'   specified by `w_select`.
#' @param w_select A logical vector indicating which weight bins were used in
#'   the likelihood calculation
#' @return The updated MizerParams object
#' @export
update_params <- function(params, species = 1, pars, biomass, w_select) {
    params <- validParams(params)
    species <- valid_species_arg(params, species, error_on_empty = TRUE)
    if (length(species) > 1) {
        stop("Only one species can be updated at a time.")
    }
    sp <- species_params(params)
    sp_select <- sp$species == species
    sps <- sp[sp_select, ]

    gp <- params@gear_params
    gp_select <- gp$species == species
    gps <- gp[gp_select, ]
    if (nrow(gps) > 1) {
        stop("The code currently assumes that there is only a single gear for each species.")
    }

    # Update the gear parameters
    gps$l50 <- pars["l50"]
    gps$l25 <- pars["ratio"] * pars["l50"]
    gps$catchability <- pars["catchability"]
    gear_params(params)[gp_select, ] <- gps

    # recalculate the power-law mortality rate
    # sps$M <- pars["M"]
    ext_mort(params)[sp_select, ] <- pars["M"] * params@w^sps$d

    # Update the steepness of the maturity ogive
    sps$w_mat25 <- sps$w_mat / 3^(1 / pars["U"])
    params@species_params[sp_select, ] <- sps
    params <- setReproduction(params)

    # Calculate the new steady state ----
    params <- steadySingleSpecies(params)
    # Rescale it to get the observed biomass
    total <- sum(params@initial_n[sp_select, w_select] *
                     params@w[w_select] * params@dw[w_select])
    factor <- biomass / total
    params@initial_n[sp_select, ] <- params@initial_n[sp_select, ] * factor

    return(params)
}
