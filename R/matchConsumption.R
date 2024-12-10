#' Match the consumption of the model to the Ecopath consumption
#'
#' This function sets the metabolic respiration rate so that the consumption
#' matches the species parameter `ecopath_consumption`, while keeping the same
#' energy available for growth and reproduction. Thus the function also adjusts
#' the external encounter rate to compensate for the changed respiration rate.
#' To be able to do this the function needs to assume that both the encounter
#' rate and the metabolic respiration rate are given by power laws with the same
#' exponent `n`.
#'
#' If the production for a species is higher than the consumption in
#' `ecopath_consumption`, then this will lead to a negative metabolic
#' respiration rate. In this case the function will issue a warning.
#'
#' @param params A MizerParams object
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
    species <- valid_species_arg(params, species, error_on_empty = TRUE)
    if (length(species) > 1) {
        for (i in seq_along(species)) {
            params <- matchConsumption(params, i)
        }
        return(params)
    }

    # Select the desired species
    sp <- params@species_params
    sp_select <- sp$species == species
    sps <- sp[sp_select, ]
    if (is.na(sps$ecopath_consumption)) {
        stop("The ecopath_consumption parameter for species ", species, " is NA.")
    }
    if (sps$p != sps$n) {
        stop("The encounter rate and metabolic respiration rate must have the same exponent.")
    }

    # According to Ecopath, R = alpha * Q - P
    R = sps$alpha * sps$ecopath_consumption - getTotalProduction(params)[sp_select]
    if (R < 0) {
        warning("Negative metabolic respiration required for species ", species, ".")
    }
    metab_old <- params@metab[sp_select, ]
    # Set metabolic rate with unit coefficient
    params@metab[sp_select, ] <- params@w ^ sps$n
    # Rescale the metabolic rate
    ks <- R / getMetabolicRespiration(params)[sp_select]
    params@species_params[sp_select, "ks"] <- ks
    params@metab[sp_select, ] <- ks * params@w ^ sps$n
    # Increase the encounter rate to compensate
    ext_encounter(params)[sp_select, ] <- ext_encounter(params)[sp_select, ] +
        (metab(params)[sp_select, ] - metab_old) / sps$alpha

    return(params)
}
