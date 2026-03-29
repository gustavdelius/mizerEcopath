# Environment to store hooks set by tuningGadget()
tuning_env <- new.env(parent = emptyenv())

#' Prepare the params object for tuning
#'
#' This function is called once when a tuning gadget is started. It prepares
#' the params object for tuning by setting defaults for the species parameters
#' that the tuning gadget uses but which may not be set in the params object.
#'
#' @param p The params object
#' @param hook An optional function `(p) -> p` called after the shared
#'   preparation to apply package-specific defaults.
#' @return The prepared params object
#' @keywords internal
#' @export
prepare_params <- function(p, hook = NULL) {
    p@species_params$species <- as.character(p@species_params$species)
    rownames(p@species_params) <- p@species_params$species
    p <- set_species_param_default(p, "a", 0.006)
    p <- set_species_param_default(p, "b", 3)
    p <- set_species_param_default(p, "k_vb", NA)
    p <- set_species_param_default(p, "t0", 0)
    p <- set_species_param_default(p, "w_mat25",
                                   p@species_params$w_mat/(3^(1/10)))

    # Apply package-specific hook
    if (!is.null(hook)) {
        p <- hook(p)
    }

    return(p)
}


#' Finalise the params object after tuning
#'
#' This function is called when a params object is downloaded or when the done
#' button is pressed. It clears the undo attribute and restores the
#' Beverton-Holt reproduction settings.
#'
#' @param p The params object
#' @return The finalised params object
#' @keywords internal
#' @export
finalise_params <- function(p) {
    # Clear attribute that was only needed for the undo functionality
    attr(p, "changes") <- NULL

    # Set reproduction
    if ("tuneParams_old_repro_level" %in% names(p@species_params)) {
        p <- setBevertonHolt(p, reproduction_level =
                                 p@species_params$tuneParams_old_repro_level)
        p@species_params$tuneParams_old_repro_level <- NULL
    }
    if ("tuneParams_old_R_max" %in% names(p@species_params)) {
        p <- setBevertonHolt(p, R_max =  p@species_params$tuneParams_old_R_max)
        p@species_params$tuneParams_old_R_max <- NULL
    }
    if ("tuneParams_old_erepro" %in% names(p@species_params)) {
        p <- setBevertonHolt(p, erepro =  p@species_params$tuneParams_old_erepro)
        p@species_params$tuneParams_old_erepro <- NULL
    }
    p
}


#' Update the species parameters
#'
#' This function is called when a species parameter is changed. It calculates the
#' steady state for the changed species and updates the params object.
#'
#' If an `update_species_hook` was set by [tuningGadget()], it is called after
#' `steadySingleSpecies()` to allow package-specific post-processing (e.g.,
#' `matchBiomasses()` in mizerEcopath).
#'
#' @param sp The species to update
#' @param p The params object
#' @param params The reactive params object
#' @param params_old The reactive params object before the change
#' @keywords internal
#' @export
tuneParams_update_species <- function(sp, p, params, params_old) {
    # wrap the code in trycatch so that when there is a problem we can
    # simply stay with the old parameters
    tryCatch({
        # The spectrum for the changed species is calculated with new
        # parameters but in the context of the original community
        p_old <- params_old()
        p@initial_n <- p_old@initial_n

        p <- steadySingleSpecies(p)

        # Apply update_species_hook if set by tuningGadget()
        hook <- tuning_env$update_species_hook
        if (!is.null(hook)) {
            p <- hook(sp, p)
        }

        # Update the reactive params object
        tuneParams_update_params(p, params)
    },
    error = function(e) {
        error_fun(e)
        tuneParams_update_params(p, params)
    })
}


#' Run to steady state
#'
#' This function is called when the user clicks the "Steady" button in
#' [tuneParams()]. It runs the model to steady state and updates the params
#' object.
#'
#' @param p The params object
#' @param params The reactive params object
#' @param params_old The reactive params object before the change
#' @param logs The logs object
#' @param session The Shiny session object
#' @param input The Shiny input object
#' @param return_sim Whether to return the simulation object
#' @keywords internal
#' @export
tuneParams_run_steady <- function(p, params, params_old, logs, session, input,
                                  return_sim = FALSE) {

    tryCatch({
        # Create a Progress object
        progress <- shiny::Progress$new(session)
        on.exit(progress$close())

        if ("biomass" %in% input$match) {
            p <- matchBiomasses(p)
        }
        if ("growth" %in% input$match) {
            p <- matchGrowth(p, keep = "biomass")
            sp <- p@species_params[p@species_params$species == input$sp, ]
            updateSliderInput(session, "gamma", value = sp$gamma)
            updateSliderInput(session, "h", value = sp$h)
            updateSliderInput(session, "ks", value = sp$ks)
            updateSliderInput(session, "k", value = sp$k)
        }
        if ("yield" %in% input$match) {
            p <- matchYield(p, keep = "biomass")
            gp_idx <- which(p@gear_params$species == input$sp &
                                p@gear_params$gear == input$gear)
            catchability <- p@gear_params[gp_idx, "catchability"]
            updateSliderInput(session, "catchability",
                              value = catchability)
        }

        # Run to steady state
        if (return_sim) {
            # This is for the "Steady" tab where we want to show the
            # evolution of biomass over time during the run to steady
            # to diagnose eventual problems.
            return(mizer::steady(p, t_max = 100, tol = 1e-2,
                          return_sim = TRUE,
                          progress_bar = progress))
        }
        p <- mizer::steady(p, t_max = 100, tol = 1e-2,
                           progress_bar = progress)

        # Update the reactive params objects
        params_old(p)
        tuneParams_add_to_logs(logs, p, params)
    },
    error = error_fun)
}


#' Update the abundance of a species
#'
#' This function is called when the abundance of a species is changed directly
#' i.e., not when it is changed as a consequence of a parameter change.
#' It will make the change permanent by also saving it in params_old.
#'
#' @param p The params object
#' @param sp The species to update
#' @param params The reactive params object
#' @param params_old The reactive params object before the change
#' @concept helper
#' @export
tuneParams_update_abundance <- function(p, sp, params, params_old) {

    # We need to update `params_old()` because otherwise the change
    # will not persist past the next parameter change.
    p_old <- isolate(params_old())
    p_old@initial_n[sp, ] <- p@initial_n[sp, ]
    params_old(p_old)

    # Update the reactive params object
    tuneParams_update_params(p, params)
}

#' Update the params object
#'
#' This function is called whenever the params object needs to be updated, unless
#' you also need to write it to the logs, in which case call
#' `tuneParams_add_to_logs()` instead.
#' It indicates that the params have changed and updates the reactive params object.
#'
#' @param p The params object
#' @param params The reactive params object
#' @keywords internal
#' @export
tuneParams_update_params <- function(p, params) {

    # indicate that the params have changed. This will be used in the Undo
    # functionality.
    if (is.null(attr(p, "changes"))) {
        attr(p, "changes") <- 1
        # Now that a change has taken place, there is certainly something to undo
        shinyjs::enable("undo")
        shinyjs::enable("undo_all")
    } else {
        attr(p, "changes") <- attr(p, "changes") + 1
    }

    params(p)
}

#' Add the params object to the logs
#'
#' This function is called when the params object needs to be written to the logs.
#' It updates the params object and writes it to the logs.
#'
#' @param logs The logs object
#' @param p The params object
#' @param params The reactive params object
#' @keywords internal
#' @export
tuneParams_add_to_logs <- function(logs, p, params) {

    # Clear attribute used in undo functionality
    attr(p, "changes") <- NULL
    # update params object
    params(p)

    # Save params object to disk
    time = format(Sys.time(), "_%Y_%m_%d_at_%H_%M_%S")
    file = paste0(tempdir(), "/mizer_params", time, ".rds")
    saveRDS(p, file = file)
    # Update logs
    if (logs$idx < length(logs$files)) {
        file.remove(logs$files[(logs$idx + 1):length(logs$files)])
    }
    logs$files <- append(logs$files[min(1, logs$idx):logs$idx], file)
    logs$idx <- logs$idx + 1
    shinyjs::disable("redo")
    if (logs$idx > 1) {
        shinyjs::enable("undo")
        shinyjs::enable("undo_all")
    } else {
        shinyjs::disable("undo")
        shinyjs::disable("undo_all")
    }
}

#' Error function
#'
#' This function is called when an error occurs. It shows a modal dialog with the error message.
#'
#' @param e The error object
#' @keywords internal
#' @export
error_fun <- function(e) {
    showModal(modalDialog(
        title = "Invalid parameters",
        HTML(paste0("These parameter do not lead to an acceptable steady state. ",
                    "Please choose other values.<br>",
                    "The error message was:<br>", e)),
        easyClose = TRUE
    ))}

#' Convert the tab name given by the user to lower case
#'
#' This function is used to convert the tab name given by the user to lower case,
#' because the names of the tab functions will always start with lower case.
#'
#' @param tab The tab name
#' @return The tab name in lower case
#' @keywords internal
#' @export
tab_name <- function(tab) {
    tabname <- tab
    substr(tabname, 1, 1) <- tolower(substr(tab, 1, 1))
    tabname
}

#' Return the title for the tab
#'
#' This function is used to return the title for the tab. This is either defined
#' by the tab author or otherwise is the tab name supplied by the user.
#'
#' @param tab The tab name
#' @return The tab title
#' @keywords internal
#' @export
tab_title <- function(tab) {
    tabname <- tab_name(tab)
    title_var <- paste0(tabname, "TabTitle")
    if (!is.null(title <- get0(title_var))) {
        if (!is.string(title)) {
            stop(title_var, "should contain a string with the title for the tab")
        }
        return(title)
    }
    tab
}
