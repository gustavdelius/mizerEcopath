#' Launch shiny gadget for matching mizer model to Ecopath data
#'
#' @description `r lifecycle::badge("experimental")`
#'
#' The function opens a shiny gadget, an interactive web page. This page has a
#' side panel with controls for various model parameters and a main panel with
#' tabs for various diagnostic plots.
#'
#' @details This gadget is meant for tuning a model to match Ecopath data.
#'
#' There is an "Instructions" button near the top left of the gadget that gives
#' you a quick overview of the user interface.
#'
#' After you click the "Return" button in the side panel, the function will
#' return the parameter object in the state at that time, with `Rmax` set to
#' `Inf` and `erepro` set to the value it had after the last match to Ecopath.
#'
#' @param params MizerParams object to tune. If missing, the gadget tries to
#'   recover information from log files left over from aborted previous runs.
#' @param diet An optional diet matrix to display in the Diet tab. If NULL,
#'   the diet tab will not show observed diet data.
#' @inheritParams matchCatch
#' @inheritParams matchDiet
#' @param controls A character vector of names of input parameter control
#'   sections that should be displayed in the sidebar. See
#'   [tuningGadget()] for details.
#' @param tabs A character vector of names of the tabs that should be displayed
#'   in the main section. See [tuningGadget()] for details.
#' @param match A character vector. Determines which quantities should be
#'   matched to observations each time the "match" button is pressed.
#' @param preserve Specifies whether the `reproduction_level` should be
#'   preserved or the maximum reproduction rate `R_max` or the reproductive
#'   efficiency `erepro` (Default).
#' @param return_app Boolean. For testing purposes only.
#' @param ... Other params needed by individual tabs.
#'
#' @return The tuned MizerParams object
#' @md
#' @export
tuneEcopath <- function(params, catch = NULL, diet = NULL,
                        controls = c("ecopathGrowth",
                                     "ecopathDiffusion",
                                     "ecopathFishing",
                                     "ecopathReproduction",
                                     "ecopathOther",
                                     "ecopathMatch"),
                        tabs = c("ecopathSurvey",
                                 "ecopathSpectra",
                                 "ecopathCatch",
                                 "ecopathGrowth",
                                 "ecopathRepro",
                                 "ecopathDiet",
                                 "ecopathDeath"),
                        match = c("growth", "yield", "catch", "consumption"),
                        preserve = c("erepro", "reproduction_level", "R_max"),
                        return_app = FALSE,
                        ...) {

    tuningGadget(
        params = params,
        controls = controls,
        tabs = tabs,
        match_choices = c("growth", "yield", "catch", "consumption"),
        match_selected = match,
        action_label = HTML("<u>m</u>atch"),
        action_tooltip = "Match ecopath parameters. Keyboard shortcut: m",
        action_key = 77,
        action_handler = function(p, input, session, params, params_old,
                                  logs, flags, trigger_update, ...) {
            ecopath_match_handler(p, catch, params = params,
                                 params_old = params_old,
                                 logs = logs, session = session, input = input)
            rm(list = ls(flags), pos = flags)
            trigger_update(runif(1))
        },
        prepare_params_hook = ecopath_prepare_hook,
        update_species_hook = ecopath_update_species_hook,
        tab_extra_args = list(catch = catch, diet = diet),
        preserve = preserve,
        return_app = return_app,
        ...
    )
}

# Ecopath-specific prepare_params hook ----
ecopath_prepare_hook <- function(p) {
    no_sp <- nrow(p@species_params)
    p <- set_species_param_default(p, "d", p@species_params$n - 1)
    p <- set_species_param_default(p, "yield_lambda", 0)
    p <- set_species_param_default(p, "production_lambda", 0)
    # Determine gonad proportion
    current <- getGonadicProduction(p) / getProduction(p)
    p <- set_species_param_default(p, "gonad_proportion", current)
    # Determine mu_mat
    mat_idx <- colSums(outer(p@w, p@species_params$w_mat, "<"))
    mu_mat <- ext_mort(p)[cbind(seq_len(no_sp), mat_idx)]
    p <- set_species_param_default(p, "mu_mat", mu_mat)
    p <- set_species_param_default(p, "d_over_g", 0)
    p <- set_species_param_default(p, "spawning_mu", 0.5)
    p <- set_species_param_default(p, "spawning_kappa", 5)
    p <- set_species_param_default(p, "annuli_min_age", 0)
    p <- set_species_param_default(p, "annuli_date", 0)

    p <- steadySingleSpecies(p)
    p <- matchBiomasses(p)
    return(p)
}

# Ecopath-specific update_species hook ----
ecopath_update_species_hook <- function(sp, p) {
    p <- matchBiomasses(p, species = sp)
    p
}

# Ecopath-specific match handler ----
ecopath_match_handler <- function(p, catch, params, params_old, logs,
                                  session, input) {

    tryCatch({
        # Debugging code to check biomass preservation
        pb <- matchBiomasses(p)
        if (!isTRUE(all.equal(getBiomass(p, use_cutoff = TRUE),
                              getBiomass(pb, use_cutoff = TRUE)))) {
            stop("Biomass has changed before matchGrowth")
        }
        if ("growth" %in% input$match) {
            p <- matchGrowth(p, species = input$sp, keep = "biomass")
            pb <- matchBiomasses(p)
            if (!isTRUE(all.equal(getBiomass(p, use_cutoff = TRUE),
                                  getBiomass(pb, use_cutoff = TRUE)))) {
                stop("Biomass has changed after matchGrowth")
            }
        }
        if ("catch" %in% input$match) {
            p <- matchCatch(p, species = input$sp, catch = catch,
                            production_lambda = 10^input$production_lambda,
                            yield_lambda = 10^input$yield_lambda)
            pb <- matchBiomasses(p)
            if (!isTRUE(all.equal(getBiomass(p, use_cutoff = TRUE),
                                  getBiomass(pb, use_cutoff = TRUE)))) {
                stop("Biomass has changed after matchCatch")
            }
        }
        if ("yield" %in% input$match) {
            p <- matchYield(p, species = input$sp, keep = "biomass")
            pb <- matchBiomasses(p)
            if (!isTRUE(all.equal(getBiomass(p, use_cutoff = TRUE),
                                  getBiomass(pb, use_cutoff = TRUE)))) {
                stop("Biomass has changed after matchYield")
            }
        }
        if ("consumption" %in% input$match) {
            p <- matchConsumption(p, species = input$sp)
            pb <- matchBiomasses(p)
            if (!isTRUE(all.equal(getBiomass(p, use_cutoff = TRUE),
                                  getBiomass(pb, use_cutoff = TRUE)))) {
                stop("Biomass has changed after matchConsumption")
            }
        }

        # Update the reactive params objects
        params_old(p)
        tuneParams_add_to_logs(logs, p, params)
    },
    error = error_fun)
}
