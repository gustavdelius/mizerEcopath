#' A single slider that adjusts both `h` and `gamma`
#'
#' @param input Reactive holding the inputs
#' @param output Reactive holding the outputs
#' @param session Shiny session
#' @param params Reactive value holding updated MizerParams object
#' @param params_old Reactive value holding non-updated MizerParams object
#' @param flags Environment holding flags to skip certain observers
#' @param ... Unused
growthControl <- function(input, output, session, params, params_old, flags,
                                    ...) {
    observeEvent(input$gamma, {
        p <- params()
        sp <- input$sp
        sp_sel <- p@species_params$species == sp
        if (!identical(sp, flags$sp_old_pred)) {
            flags$sp_old_pred <- sp
            return()
        }
        factor <- input$gamma / p@species_params[sp, "gamma"]
        # adjust gamma
        updateSliderInput(session, "gamma",
                          min = signif(input$gamma / 2, 3),
                          max = signif(input$gamma * 1.5, 3))
        p@species_params[sp, "gamma"] <- input$gamma
        p <- setSearchVolume(p)

        # adjust h
        p@species_params[sp, "h"] <- p@species_params[sp, "h"] * factor
        p <- setMaxIntakeRate(p)

        # adjust external encounter
        ext_encounter(p)[sp_sel, ] <- ext_encounter(p)[sp_sel, ] * factor

        tuneParams_update_species(sp, p, params, params_old)
    },
    ignoreInit = TRUE,
    ignoreNULL = TRUE
    )
}

#' @rdname growthControl
#' @param params The MizerParams object currently being tuned.
#' @param input Reactive holding the inputs
growthControlUI <- function(p, input) {
    sp <- p@species_params[input$sp, ]
    tagList(
        tags$h3(tags$a(id = "predation"), "Predation"),
        popify(sliderInput("gamma", "Search volume coefficient 'gamma'",
                    value = sp$gamma,
                    min = signif(sp$gamma / 2, 3),
                    max = signif(sp$gamma * 1.5, 3),
                    step = sp$gamma / 50, ticks = FALSE),
               title = "Adjusting predation parameters",
               content= "This slider adjusts the coefficient `h` of the maximum intake rate and the external encounter rate at the same time as `gamma` to keep the feeding level at the smallest size unchanged.")
        )
}
