#' Add Ecopath catch data to gear parameters
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
        l25 = w2l(sp$w_mat, params) * 0.8,
        catchability = 1,
        yield_observed = 0
    )
    for (i in seq_len(nrow(sp))) {
        for (group in sp$ecopath_groups[[i]]) {
            yield <- ecopath_catch$TotalCatch..t.km..year.[ecopath_catch$Group.name == group]
            if (length(yield) == 0) {
                warning("No catch data found for group ", group)
            }
            gp$yield_observed[i] <- sum(gp$yield_observed[i], yield, na.rm = TRUE)
        }
    }
    gp$catchability <- gp$yield_observed / sp$biomass_observed
    gear_params(params) <- gp
    initial_effort(params) <- 1

    return(params)
}
