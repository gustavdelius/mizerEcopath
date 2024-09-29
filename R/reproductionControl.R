#' Controlling the reproduction parameters in the tuning gadget
#'
#' @param input Reactive holding the inputs
#' @param output Reactive holding the outputs
#' @param session Shiny session
#' @param params Reactive value holding updated MizerParams object
#' @param params_old Reactive value holding non-updated MizerParams object
#' @param flags Environment holding flags to skip certain observers
#' @param ... Unused
reproductionControl <- function(input, output, session, params, params_old,
                                flags, ...) {
    observeEvent(
        list(input$w_mat, input$wfrac, input$gonad_proportion, input$m),
        {
            p <- params()
            sp <- input$sp
            if (!identical(sp, flags$sp_old_repro)) {
                flags$sp_old_repro <- sp
                return()
            }
            # Update w_max in respone to change in gonad_proportion
            ratio <- input$gonad_proportion /
                p@species_params[sp, "gonad_proportion"]
            if (ratio != 1) {
                w_maxratio <- ratio ^ (1 / (p@species_params[sp, "n"] -
                                                p@species_params[sp, "m"]))
                p@species_params[sp, "w_max"] <- p@species_params[sp, "w_max"] *
                    w_maxratio
            }

            # Update slider min/max so that they are a fixed proportion of the
            # parameter value
            updateSliderInput(session, "w_mat",
                              min = signif(input$w_mat / 2, 2),
                              max = signif(input$w_mat * 1.5, 2))
            updateSliderInput(session, "gonad_proportion",
                              min = signif(input$gonad_proportion / 2, 2),
                              max = signif(input$gonad_proportion * 1.5, 2))

            p@species_params[sp, "w_mat25"]   <- input$w_mat * input$wfrac
            p@species_params[sp, "w_mat"]   <- input$w_mat
            p@species_params[sp, "gonad_proportion"]   <- input$gonad_proportion
            p@species_params[sp, "m"]     <- input$m

            p <- setReproduction(p)
            tuneParams_update_species(sp, p, params, params_old)
        },
        ignoreInit = TRUE)
}

#' @rdname reproductionControl
#' @param The MizerParams object currently being tuned.
#' @param input Reactive holding the inputs
#' @return A tagList with sliders for the exponents
reproductionControlUI <- function(p, input) {
    sp <- p@species_params[input$sp, ]
    tagList(
        tags$h3(tags$a(id = "reproduction"), "Reproduction"),
        sliderInput("w_mat", "w_mat", value = sp$w_mat,
                    min = signif(sp$w_mat / 2, 2),
                    max = signif(sp$w_mat * 1.5, 2)),
        sliderInput("wfrac", "w_mat25/w_mat", value = sp$w_mat25/sp$w_mat,
                    min = 0.01,
                    max = 1,
                    step = 0.01),
        sliderInput("gonad_proportion", "gonad_proportion",
                    value = sp$gonad_proportion,
                    min = signif(sp$gonad_proportion / 2, 2),
                    max = signif(sp$gonad_proportion * 1.5, 2)),
        # m must always be larger than n
        sliderInput("m", "m", value = sp$m,
                    min = signif(sp$n, 2),
                    max = signif(1.5, 2),
                    step = 0.01)
    )
}

