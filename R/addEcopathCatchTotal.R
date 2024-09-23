#' Add Ecopath catch data to gear parameters
#'
#' @export
addEcopathCatchTotal <- function(params, species_params,
                                 groups_to_species, catch) {
    sp <- validSpeciesParams(species_params)
    gp <- data.frame(
        species = sp$species,
        gear = "total",
        sel_func = "sigmoid_length",
        l50 = w2l(sp$w_mat, p),
        l25 = w2l(sp$w_mat, p) * 0.8,
        catchability = 1,
        yield_observed = 0
    )
    for (species in gp$species) {
        for (group in groups_to_species[[species]]) {
            yield <- catch$TotalCatch..t.km..year.[catch$Group.name == group]
            gp$yield_observed[gp$species == species] <-
                gp$yield_observed[gp$species == species] + yield
        }
    }
    gp$catchability <- gp$yield_observed / sp$biomass_observed
    gear_params(params) <- gp
    initial_effort(params) <- 1

    return(params)
}
