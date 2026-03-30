#' Ecopath-specific reproduction control
#'
#' This control adjusts reproduction parameters using length-based inputs:
#' *   **l_mat**: Length at 50% maturity.
#' *   **l_mat25/l_mat**: Ratio between length at 25% and 50% maturity.
#' *   **l_repro_max**: Maximum length for reproduction.
#' *   **m**: Allometric exponent for reproductive effort.
#'
#' @inheritParams ecopathOtherControl
#' @family gadget controls
#' @export
ecopathReproductionControl <- function(input, output, session, params, params_old,
                                flags, ...) {
    observeEvent(
        list(input$l_mat, input$lfrac,
             input$m, input$l_repro_max),
        {
            p <- params()
            sp <- input$sp
            sps <- p@species_params[sp, ]
            if (!identical(sp, flags$sp_old_repro)) {
                flags$sp_old_repro <- sp
                return()
            }

            updateSliderInput(session, "l_mat",
                              min = signif(input$l_mat / 2, 2),
                              max = signif(input$l_mat * 1.5, 2))
            updateSliderInput(session, "l_repro_max",
                              min = signif(input$l_repro_max / 2, 2),
                              max = signif(input$l_repro_max * 1.5, 2))

            l_mat_25 <- input$l_mat * input$lfrac
            w_mat_25 <- sps$a * (l_mat_25 ^ sps$b)
            w_mat <- sps$a * (input$l_mat ^ sps$b)
            w_repro_max <- sps$a * (input$l_repro_max ^ sps$b)
            p@species_params[sp, "w_mat25"]   <- w_mat_25
            p@species_params[sp, "w_mat"]   <- w_mat
            p@species_params[sp, "w_repro_max"] <- w_repro_max
            p@species_params[sp, "m"]     <- input$m

            p <- setReproduction(p)
            tuneParams_update_species(sp, p, params, params_old)
        },
        ignoreInit = TRUE)
}

#' @rdname ecopathReproductionControl
#' @param params The MizerParams object currently being tuned.
#' @param input Reactive holding the inputs
#' @return A tagList with sliders for the exponents
ecopathReproductionControlUI <- function(params, input) {
    sps <- params@species_params[input$sp, ]
    l_mat <- (sps$w_mat / sps$a) ^ (1 / sps$b)
    l_mat25 <- (sps$w_mat25 / sps$a) ^ (1 / sps$b)
    l_repro_max <- (sps$w_repro_max / sps$a) ^ (1 / sps$b)
    tagList(
        tags$h3(tags$a(id = "reproduction"), "Reproduction"),
        sliderInput("l_mat", "l_mat", value = l_mat,
                    min = signif(l_mat / 2, 2),
                    max = signif(l_mat * 1.5, 2)),
        sliderInput("lfrac", "l_mat25/l_mat", value = l_mat25/l_mat,
                    min = 0.01,
                    max = 1,
                    step = 0.01),
        sliderInput("l_repro_max", "l_repro_max",
                    value = l_repro_max,
                    min = signif(l_repro_max / 2, 1),
                    max = signif(l_repro_max * 1.5, 1)),
        sliderInput("m", "m", value = sps$m,
                    min = signif(sps$n, 2),
                    max = signif(1.5, 2),
                    step = 0.01)
    )
}

ecopathReproductionControlTitle <- "Reproduction"
