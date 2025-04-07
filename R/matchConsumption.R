#' Match the model’s total consumption to observed estimates
#'
#' This function sets the metabolic respiration rate so that the total consumption
#' matches the parameter `consumption_observed` for each species. The function preserves
#' the energy available for growth and reproduction by increasing external encounter
#' rates to compensate for changes in metabolic loss.
#'
#' Thus the function also adjusts the external encounter rate to compensate for the
#' changed respiration rate. To do this the function assumes that both the encounter rate
#' and the metabolic respiration rate are given by power laws with the same exponent
#' `n`, so it sets the species parameter `p` to the same value as `n`. A warning
#' is issued if the exponent `p` had to be changed for any species.
#'
#' Any of the selected species for which `consumption_observed` is NA will be
#' quietly ignored.
#'
#' If the resulting metabolic loss rate is less than 10% of the production rate,
#' the function will set the metabolic loss rate to 10% of the production rate
#' and issue a warning.
#'
#' Unless the function issues a warning that it has changed `p`, the energy
#' available for growth and reproduction is not changed and hence the steady
#' state spectra are also unchanged.
#'
#' @param params A MizerParams object
#' @param species A vector of species names or indices. If NULL,
#'   applies to all species with a provided `consumption_observed`.
#'#' @return A `MizerParams` object with adjusted encounter and metabolic
#'   respiration rates.
#' @family match functions
#' @examples
#' params <- matchConsumption(celtic_params)
#' # The consumption now matches the observation
#' all.equal(getConsumption(params),
#'           species_params(params)$consumption_observed,
#'           check.attributes = FALSE)
#' # The energy available for growth and reproduction is not changed
#' all.equal(getEReproAndGrowth(params),
#'           getEReproAndGrowth(celtic_params))
#' @export
matchConsumption <- function(params, species = NULL) {
    if (!is(params, "MizerParams")) {
        stop("params must be a MizerParams object.")
    }
    if (!hasName(params@species_params, "consumption_observed")) {
        stop("You must provide the consumption_observed species parameter.")
    }
    species <- valid_species_arg(params, species)

    # Identify selected species
    sp <- params@species_params
    sp_select <- sp$species %in% species & !is.na(sp$consumption_observed)

    if (sum(sp_select) == 0) {
        return(params)
    }

    # Set p = n for selected species
    wrong <- !is.null(sp$p[sp_select]) & !is.na(sp$p[sp_select]) &
        sp$p[sp_select] != sp$n[sp_select]
    if (any(wrong)) {
        wrong_select <- sp_select[wrong]
        warning("Exponent `p` changed for ",
                paste(sp$species[wrong_select], collapse = ", "), ".")
    }
    params@species_params$p[sp_select] <- sp$n[sp_select]

    # Calculate R = alpha * Q - P for each selected species
    total_production <- getTotalProduction(params)
    R <- sp$alpha[sp_select] * sp$consumption_observed[sp_select] - total_production[sp_select]

    # If the needed ks value is below the minimal acceptable value, set it to the
    # minimal acceptable value and issue a warning.
    problem_species <- sp$species[sp_select][R < 0.1 * total_production[sp_select]]
    if (length(problem_species) > 0) {
        warning("Perfect match to Ecopath consumption not possible for: ",
                paste(problem_species, collapse = ", "),
                " because it would lead to a low metabolic rate of less than 10% of the production rate.")
        R <- pmax(R, 0.1 * total_production[sp_select])
    }

    # Store old metabolic rates
    metab_old <- params@metab[sp_select, , drop = FALSE]

    # Reset metabolic rate to w^n for each selected species
    w <- params@w
    n_vals <- sp$n[sp_select]
    # Create a matrix of w^(n) for each species (rows = species, cols = size classes)
    w_matrix <- matrix(w, nrow = sum(sp_select), ncol = length(w), byrow = TRUE)
    n_matrix <- matrix(n_vals, nrow = sum(sp_select), ncol = length(w))
    metab_unit <- w_matrix ^ n_matrix
    params@metab[sp_select, ] <- metab_unit

    # Calculate scaling factor ks so that total metabolic respiration matches R
    # getMetabolicRespiration(params) returns a vector (one value per species)
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
