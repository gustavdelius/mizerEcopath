#' Adjusting predation parameters for growth
#'
#' This control provides a single slider for the search volume coefficient
#' `gamma`. Each time `gamma` is changed, the maximum intake rate coefficient
#' `h` is also adjusted by the same factor. This maintains the same growth
#' rate at small sizes while allowing for tuning of the overall growth
#' potential.
#'
#' @inheritParams abundanceControl
#' @family gadget controls
#' @export
growthControl <- function(input, output, session, params, params_old, flags,
                                    ...) {
    observeEvent(input$gamma, {
        p <- params()
        sp <- input$sp
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

        tuneParams_update_species(sp, p, params, params_old)
    },
    ignoreInit = TRUE,
    ignoreNULL = TRUE
    )
}

#' @rdname growthControl
#' @inheritParams abundanceControlUI
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
               content= "This slider adjusts the coefficient `h` of the maximum intake rate at the same time as `gamma` to keep the feeding level at the smallest size unchanged.")
        )
}
