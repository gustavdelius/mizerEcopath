#' Set initial abundances to single-species steady state abundances
#'
#' `r lifecycle::badge("experimental")`
#'
#' This is a thin wrapper around [mizer::steadySingleSpecies()] so that
#' all `mizerEcopath` code uses the implementation provided by the
#' installed version of `mizer`.
#'
#' @inheritParams mizer::steadySingleSpecies
#' @export
steadySingleSpecies <- function(params,
                                species = NULL,
                                keep = c("egg", "biomass", "number"),
                                ...) {
    mizer::steadySingleSpecies(
        params = params,
        species = species,
        keep = keep,
        ...
    )
}
