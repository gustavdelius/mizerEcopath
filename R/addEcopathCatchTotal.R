#' Add Ecopath catch data to gear parameters
#'
#' @export
addEcopathCatchTotal <- function(params, catch) {
    sp <- p@species_params
    if (!hasName(sp, "ecopath_groups") ||
        !hasName(sp, "biomass_observed")) {
        stop("You must use `addEcopathParams()` first.")
    }
    if (!hasName(catch, "TotalCatch..t.km..year.") ||
        !hasName(catch, "Group.name")) {
        stop("The catch argument is invalid. It must be a data frame with columns 'TotalCatch..t.km..year.' and 'Group.name'.")
    }
    gp <- data.frame(
        species = sp$species,
        gear = "total",
        sel_func = "sigmoid_length",
        l50 = w2l(sp$w_mat, p),
        l25 = w2l(sp$w_mat, p) * 0.8,
        catchability = 1,
        yield_observed = 0
    )
    for (i in seq_len(nrow(sp))) {
        for (group in sp$ecopath_groups[[i]]) {
            yield <- catch$TotalCatch..t.km..year.[catch$Group.name == group]
            if (length(yield) == 0) {
                warning("No catch data found for group ", group)
            }
            if (is.na(yield)) {
                warning("Missing catch data for group ", group)
            }
            gp$yield_observed[i] <- gp$yield_observed[i] + yield
        }
    }
    gp$catchability <- gp$yield_observed / sp$biomass_observed
    gear_params(params) <- gp
    initial_effort(params) <- 1

    return(params)
}
