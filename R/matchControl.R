#' Controlling the matching
#'
#' @param input Reactive holding the inputs
#' @param output Reactive holding the outputs
#' @param session Shiny session
#' @param params Reactive value holding updated MizerParams object
#' @param params_old Reactive value holding non-updated MizerParams object
#' @param flags Environment holding flags to skip certain observers
#' @param ... Unused
matchControl <- function(input, output, session, params, params_old,
                         flags, ...) {
    observeEvent(
        list(input$production_lambda, input$yield_lambda),
        {
            p <- params()
            sp <- input$sp
            if (!identical(sp, flags$sp_old_match)) {
                flags$sp_old_match <- sp
                return()
            }

            # change species parameters
            p@species_params[[sp, "production_lambda"]] <- input$production_lambda
            p@species_params[[sp, "yield_lambda"]] <- input$yield_lambda
            updateSliderInput(session, "production_lambda",
                              min = signif(input$production_lambda - 2, 3),
                              max = signif(input$production_lambda + 2, 3))
            updateSliderInput(session, "yield_lambda",
                              min = signif(input$yield_lambda - 2, 3),
                              max = signif(input$yield_lambda + 2, 3))

            tuneParams_update_species(sp, p, params, params_old)
        },
        ignoreInit = TRUE)
}

#' @rdname matchControl
#'
#' @param params The MizerParams object currently being tuned.
#' @param input Reactive holding the inputs
#' @return A tagList with sliders for the exponents
matchControlUI <- function(params, input) {
    sp <- params@species_params[input$sp, ]
    tagList(
        tags$h3(tags$a(id = "match"), "Matching strength"),
        sliderInput("production_lambda", "Production",
                    value = sp[["production_lambda"]],
                    min = sp[["production_lambda"]] - 2,
                    max = sp[["production_lambda"]] + 2,
                    step = 0.1),
        sliderInput("yield_lambda", "Yield",
                    value = sp[["yield_lambda"]],
                    min = sp[["yield_lambda"]] - 2,
                    max = sp[["yield_lambda"]] + 2,
                    step = 0.1)
    )
}
