#' Ecopath-specific other parameter control
#'
#' This control adjusts the external mortality at 1g and the metabolic
#' rate. If the external mortality at 1g is changed, then the
#' external mortality at all other sizes is scaled by the same factor.
#'
#' @param input Reactive holding the inputs
#' @param output Reactive holding the outputs
#' @param session Shiny session
#' @param params Reactive value holding updated MizerParams object
#' @param params_old Reactive value holding non-updated MizerParams object
#' @param flags Environment holding flags to skip certain observers
#' @param ... Unused
#' @family gadget controls
#' @export
ecopathOtherControl <- function(input, output, session, params, params_old,
                         flags, ...) {
    observe({
        req(input$z_ext)
        p <- isolate(params())
        sp <- isolate(input$sp)
        if (!identical(sp, flags$sp_old_other)) {
            flags$sp_old_other <- sp
            return()
        }
        z_ext <- ext_mort(p)[input$sp, 1] / p@w[1]^p@species_params[sp, "d"]

        if (z_ext != input$z_ext) {
            updateSliderInput(session, "z_ext",
                              min = signif(input$z_ext / 2, 2),
                              max = signif((input$z_ext + 0.1) * 1.5, 2))
            if (z_ext > 0) {
                ext_mort(p)[sp, ] <- ext_mort(p)[sp, ] * (input$z_ext / z_ext)
            } else {
                ext_mort(p)[sp, ] <- input$z_ext *
                    p@w ^ p@species_params[sp, "d"]
            }
        }

        p@species_params[sp, "z_ext"] <- input$z_ext
        p <- setMetabolicRate(p, reset = TRUE)
        tuneParams_update_species(sp, p, params, params_old)
    })
}

#' @rdname ecopathOtherControl
#' @param params The MizerParams object currently being tuned.
#' @param input Reactive holding the inputs
#' @return A tagList with sliders for the exponents
ecopathOtherControlUI <- function(params, input) {
    sp <- params@species_params[input$sp, ]
    tagList(
        tags$h3(tags$a(id = "ext_mort"), "Mort"),
        sliderInput("z_ext", "External mortality at 1g",
                    value = sp$z_ext,
                    min = signif(sp$z_ext / 2, 2),
                    max = signif((sp$z_ext + 0.1) * 1.5, 2),
                    step = 0.05)
    )
}

ecopathOtherControlTitle <- "Other"
ecopathOtherControlDescription <-
    "Adjust the external mortality at 1g (scaling the external mortality at all sizes by the same factor) and the metabolic rate."
