#' Match the consumption of the model to the Ecopath consumption
#'
#' This function sets the metabolic respiration rate so that the consumption
#' matches the species parameter `ecopath_consumption`, while keeping the same
#' energy available for growth and reproduction. Thus the function also adjusts
#' the external encounter rate to compensate for the changed respiration rate.
#' To do this the function assumes that both the encounter rate and the metabolic
#' respiration rate are given by power laws with the same exponent `n`, so it
#' sets the species parameter `p` to the same value as `n`.
#'
#' Any of the selected species for which `ecopath_consumption` is NA will be
#' quietly ignored.
#' If for any species the production is higher than the `ecopath_consumption`,
#' then this will lead to a negative metabolic respiration rate. In this case
#' the function will issue a warning for all such species.
#'
#' @param params A MizerParams object
#' @param species A vector of species names or indices. If NULL, all species for
#'   which the species parameter `ecopath_consumption` is provided are affected.
#'
#' @return A MizerParams object with adjusted encounter and metabolic respiration
#'   rates.
#' @family match functions
#' @examples
#' params <- matchConsumption(celtic_params)
#' # The consumption now matches the observation
#' all.equal(getConsumption(params),
#'           species_params(params)$ecopath_consumption,
#'           check.attributes = FALSE)
#' # The energy available for growth and reproduction is not changed
#' all.equal(getEReproAndGrowth(params),
#'           getEReproAndGrowth(celtic_params))
#' @export
matchConsumption <- function(params, species = NULL) {
    if (!is(params, "MizerParams")) {
        stop("params must be a MizerParams object.")
    }
    if (!hasName(params@species_params, "ecopath_consumption")) {
        stop("You must provide the ecopath_consumption species parameter.")
    }
    species <- valid_species_arg(params, species)

    # Identify selected species
    sp <- params@species_params
    sp_select <- sp$species %in% species & !is.na(sp$ecopath_consumption)

    if (sum(sp_select) == 0) {
        return(params)
    }

    # Set p = n for selected species
    params@species_params$p[sp_select] <- sp$n[sp_select]

    # Calculate R = alpha * Q - P for each selected species
    total_production <- getTotalProduction(params)
    R <- sp$alpha[sp_select] * sp$ecopath_consumption[sp_select] - total_production[sp_select]

    # Warn if any species require negative metabolic respiration
    negative_species <- sp$species[sp_select][R < 0]
    if (length(negative_species) > 0) {
        warning("Negative metabolic respiration required for species: ",
                paste(negative_species, collapse = ", "), ".")
    }

    # Store old metabolic rates
    metab_old <- params@metab[sp_select, , drop = FALSE]

    # Reset metabolic rate to w^n for each selected species
    # We need a matrix of w raised to the power n for each selected species
    w <- params@w
    n_vals <- sp$n[sp_select]
    # Create a matrix of w^(n) for each species (rows = species, cols = size classes)
    w_matrix <- matrix(w, nrow = sum(sp_select), ncol = length(w), byrow = TRUE)
    n_matrix <- matrix(n_vals, nrow = sum(sp_select), ncol = length(w))
    metab_unit <- w_matrix ^ n_matrix
    params@metab[sp_select, ] <- metab_unit

    # Calculate scaling factor ks so that total metabolic respiration matches R
    # getMetabolicRespiration(params) should return a vector (one value per species)
    tot_metab_resp <- getMetabolicRespiration(params)[sp_select]
    ks <- R / tot_metab_resp

    # Update ks in species_params
    params@species_params$ks[sp_select] <- ks

    # Scale metab by ks
    # Multiply each species' w^n vector by its ks
    metab_new <- sweep(metab_unit, 1, ks, "*")
    params@metab[sp_select, ] <- metab_new

    # Increase the encounter rate to compensate: ext_encounter + (metab_new - metab_old) / alpha
    ext_en <- ext_encounter(params)
    met_diff <- metab_new - metab_old
    alpha_vec <- sp$alpha[sp_select]
    adj_mat <- sweep(met_diff, 1, alpha_vec, "/")
    ext_en[sp_select, ] <- ext_en[sp_select, ] + adj_mat
    ext_encounter(params) <- ext_en

    return(params)
}
