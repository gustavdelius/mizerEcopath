#' Calculate age at maturity from von Bertalanffy growth parameters
#'
#' This is not a good way to determine the age at maturity because the von
#' Bertalanffy growth curve is not reliable for larvae and juveniles. However
#' this was used in previous versions of mizer and is supplied for
#' backwards compatibility.
#'
#' Uses the age at maturity that is implied by the von Bertalanffy growth curve
#' specified by the `w_inf`, `k_vb`, `t0`, `a` and `b` parameters in the
#' species_params data frame.
#'
#' If any of `k_vb` is missing for a species, the function returns NA for that
#' species. Default values of `b = 3` and `t0 = 0` are used if these are
#' missing. If `w_inf` is missing, `w_max` is used instead.
#'
#' @param object A MizerParams object or a species_params data frame
#' @return A named vector. The names are the species names and the values are
#'   the ages at maturity.
#' @concept helper
age_mat_vB <- function(object) {
    if (is(object, "MizerParams")) {
        sp <- object@species_params
    } else {
        if (!is.data.frame(object)) {
            stop("The first argument must be either a MizerParams object or a species_params data frame.")
        }
        sp <- validSpeciesParams(object)
    }
    sp <- set_species_param_default(sp, "t0", 0)
    sp <- set_species_param_default(sp, "b", 3)
    sp <- set_species_param_default(sp, "k_vb", NA)
    sp <- set_species_param_default(sp, "w_inf", sp$w_max)

    a_mat <- -log(1 - (sp$w_mat / sp$w_inf) ^ (1/sp$b)) / sp$k_vb + sp$t0
    names(a_mat) <- sp$species
    a_mat
}
