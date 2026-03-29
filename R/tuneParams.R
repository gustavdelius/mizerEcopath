#' Launch shiny gadget for tuning parameters
#'
#' The function opens a shiny gadget, an interactive web page. This page has a
#' side panel with controls for various model parameters and a main panel with
#' tabs for various diagnostic plots.
#'
#' @details This gadget is meant for tuning a model to steady state. It is not
#' meant for tuning the dynamics of the model. That should be done in a second
#' step using functions like `setRmax()` or `changeResource()`.
#'
#' There is an "Instructions" button near the top left of the gadget that gives
#' you a quick overview of the user interface.
#'
#' After you click the "Return" button in the side panel, the function will
#' return the parameter object in the state at that time, with `Rmax` set to
#' `Inf` and `erepro` set to the value it had after the last run to steady
#' state.
#'
#' At any time the gadget allows the user to download the current params object
#' as an .rds file via the "Download" button.
#'
#' @inheritParams tuningGadget
#' @param match A character vector. Determines which quantities should be
#'   matched to observations each time the "steady" button is pressed. Possible
#'   entries are `"growth"` (using [matchGrowth()]), `"biomass"` (using
#'   [matchBiomasses()]) and `"yield"` (using [matchYields()]).
#'
#' @return The tuned MizerParams object
#' @md
#' @export
tuneParams <- function(params,
                       controls = c("predation",
                                   "fishing",
                                   "reproduction",
                                   "other",
                                   "interaction",
                                   "resource"),
                       tabs = c("Spectra",
                                "Abundance",
                                "Growth",
                                "Repro",
                                "Catch",
                                "Diet",
                                "Death",
                                "Resource",
                                "Rates",
                                "Sim"),
                       match = c("none"),
                       preserve = c("erepro", "reproduction_level", "R_max"),
                       return_app = FALSE,
                       ...) {

    tuningGadget(
        params = params,
        controls = controls,
        tabs = tabs,
        match_choices = c("growth", "biomass", "yield"),
        match_selected = match,
        action_label = HTML("<u>s</u>teady"),
        action_tooltip = "Find steady state. Keyboard shortcut: s",
        action_key = 83,
        action_handler = tuneParams_steady_handler,
        prepare_params_hook = tuneParams_prepare_hook,
        finalise_params_hook = tuneParams_finalise_hook,
        preserve = preserve,
        return_app = return_app,
        ...
    )
}

#' @noRd
tuneParams_prepare_hook <- function(p) {
    # Make sure every species has a gear by adding a "no gear" gear
    # This will be removed again at the end
    sp <- species_params(p)
    gp <- gear_params(p)
    missing_species <- setdiff(sp$species, gp$species)
    if (length(missing_species) > 0) {
        gp_missing <- data.frame(species = missing_species,
                                 gear = "no gear",
                                 catchability = 0,
                                 sel_func = "sigmoid_length",
                                 l25 = 10,
                                 l50 = 15)
        gp_missing <- validGearParams(gp_missing, sp)
        gear_params(p) <- dplyr::bind_rows(gp, gp_missing)
    }
    return(p)
}

#' @noRd
tuneParams_finalise_hook <- function(p) {
    # Remove the "no gear" gear that was added at the beginning
    gp <- gear_params(p)
    gp <- gp[gp$gear != "no gear", ]
    gear_params(p) <- gp
    p
}

#' @noRd
tuneParams_steady_handler <- function(p, input, session, params, params_old,
                                      logs, flags, trigger_update, ...) {
    tuneParams_run_steady(p, params = params,
                          params_old = params_old,
                          logs = logs, session = session, input = input)
}


#' Tune Growth
#'
#' A simplified instance of `tuneParams()` that is useful for tuning growth
#' rates.
#' @inheritParams tuneParams
#' @return The tuned MizerParams object
#' @export
tuneGrowth <- function(params, match = "biomass") {
    match <- match.arg(match)
    # Want to include tab appropriate for matched quantity
    tabname <- match
    substr(tabname, 1, 1) <- toupper(substr(match, 1, 1))
    if (tabname == "Yield") {
        tabname <- "Catch"
    }
    tuneParams(params, controls = c("growth"),
               tabs = c("Growth", tabname, "Spectra"),
               match = match)
}
