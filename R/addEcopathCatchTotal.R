#' Sets up gear parameters from Ecopath catch data
#'
#' Currently this sets up only a single gear for the total catch. The
#' `yield_observed` gear parameter is set to the total catch for each species.
#' The `catchability` is set to the ratio of the observed yield to the observed
#' biomass. The selectivity of the gear is assumed to be sigmoidal with 50% of
#' fish selected at maturity size and 25% selected at 90% of maturity size
#' (measured in weight). Fishing is turned on with an initial effort of 1. The
#' result will therefore not be at steady state yet.
#'
#' @param params A MizerParams object
#' @param ecopath_catch The Ecopath Catch data frame
#' @export
addEcopathCatchTotal <- function(params, ecopath_catch) {
    sp <- params@species_params
    if (!hasName(sp, "ecopath_groups") ||
        !hasName(sp, "biomass_observed")) {
        stop("You must use `addEcopathParams()` first.")
    }
    if (!hasName(ecopath_catch, "TotalCatch..t.km..year.") ||
        !hasName(ecopath_catch, "Group.name")) {
        stop("The ecopath_catch argument is invalid. It must be a data frame with columns 'TotalCatch..t.km..year.' and 'Group.name'.")
    }
    gp <- data.frame(
        species = sp$species,
        gear = "total",
        sel_func = "sigmoid_length",
        l50 = w2l(sp$w_mat, params),
        l25 = w2l(sp$w_mat, params) * 0.9,
        # catchability and yield_observed will be extracted from Ecopath later
        catchability = 0,
        yield_observed = 0
    )
    # Extract total yield for each species from Ecopath
    for (i in seq_len(nrow(sp))) {
        for (group in sp$ecopath_groups[[i]]) {
            yield <- ecopath_catch$TotalCatch..t.km..year.[ecopath_catch$Group.name == group]
            if (length(yield) == 0) {
                warning("No catch data found for group ", group)
            }
            gp$yield_observed[i] <- sum(gp$yield_observed[i], yield, na.rm = TRUE)
        }
    }
    # Set catchability
    gp$catchability <- gp$yield_observed / sp$biomass_observed
    # Set gear parameters
    gear_params(params) <- validGearParams(gp, sp)
    # Turn on fishing
    initial_effort(params) <- 1

    return(params)
}
