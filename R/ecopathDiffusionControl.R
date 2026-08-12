#' Ecopath-specific diffusion control
#'
#' This control adjusts parameters related to age-based diffusion and cohort
#' dynamics in an Ecopath-aligned mizer model. It includes settings for:
#' *   **Mean time/concentration of spawning**: Controls the seasonal timing and
#'     duration of reproduction.
#' *   **Annulus date/age**: Specifies when and at what age annual growth
#'     markers (annuli) are formed.
#' *   **Diffusion strength**: Adjusts `D_ext`, the coefficient of the
#'     external-diffusion power law `d(w) = D_ext * w^(n+1)`, which determines
#'     how much individuals of the same age can vary in size.
#'
#' @inheritParams ecopathOtherControl
#' @family gadget controls
#' @export
ecopathDiffusionControl <- function(input, output, session, params, params_old,
                             flags, ...) {
    observe({
        req(input$D_ext, input$spawning_mu, input$spawning_kappa,
            input$annuli_min_age, input$annuli_date)
        p <- isolate(params())
        sp <- isolate(input$sp)
        if (!identical(sp, flags$sp_old_diffusion)) {
            flags$sp_old_diffusion <- sp
            return()
        }

        updateSliderInput(session, "spawning_mu",
                          value = input$spawning_mu)
        updateSliderInput(session, "spawning_kappa",
                          value = input$spawning_kappa,
                          min = signif(input$spawning_kappa / 2, 2),
                          max = signif((input$spawning_kappa + 0.1) * 1.5, 2))
        updateSliderInput(session, "annuli_date",
                          value = input$annuli_date)
        updateSliderInput(session, "annuli_min_age",
                          value = input$annuli_min_age)
        updateSliderInput(session, "D_ext",
                          value = input$D_ext,
                          min = signif(input$D_ext / 2, 2),
                          max = signif((input$D_ext + 0.1) * 1.5, 2))

        # `recalculate = FALSE` records the new values, so that a later
        # recalculation does not undo them, without recalculating the external
        # diffusion rate, which is set explicitly below.
        sp_new <- p@species_params
        sp_new[[sp, "spawning_mu"]] <- input$spawning_mu
        sp_new[[sp, "spawning_kappa"]] <- input$spawning_kappa
        sp_new[[sp, "annuli_min_age"]] <- input$annuli_min_age
        sp_new[[sp, "annuli_date"]] <- input$annuli_date
        sp_new[[sp, "D_ext"]] <- input$D_ext
        species_params(p, recalculate = FALSE) <- sp_new
        # Update the external-diffusion power law d(w) = D_ext * w^(n+1). We let
        # mizer build it from the new `D_ext` rather than writing the power law
        # by hand, so that it is bin-averaged exactly as `setExtDiffusion()`
        # does when `second_order_w(p)[["bin_average"]]` is TRUE. Only the row
        # of the selected species is taken across, leaving any other species
        # whose diffusion was set by hand untouched.
        ed <- ext_diffusion(p)
        ed[sp, ] <- ext_diffusion(setExtDiffusion(p, reset = TRUE))[sp, ]
        ext_diffusion(p) <- ed
        tuneParams_update_species(sp, p, params, params_old)
    })
}

#' @rdname ecopathDiffusionControl
#' @param params The MizerParams object currently being tuned.
#' @param input Reactive holding the inputs
#' @return A tagList with sliders for the exponents
ecopathDiffusionControlUI <- function(params, input) {
    sp <- params@species_params[input$sp, ]
    tagList(
        sliderInput("spawning_mu", "Mean time of spawning",
                    value = sp$spawning_mu,
                    min = 0,
                    max = 0.99,
                    step = 0.01),
        sliderInput("spawning_kappa", "Concentration of spawning",
                    value = sp$spawning_kappa,
                    min = signif(sp$spawning_kappa / 2, 2),
                    max = signif((sp$spawning_kappa + 0.1) * 1.5, 2),
                    step = 0.01),
        sliderInput("annuli_date", "Annulus date",
                    value = sp$annuli_date,
                    min = 0,
                    max = 0.99,
                    step = 0.01),
        sliderInput("annuli_min_age", "Min age for annulus",
                    value = sp$annuli_min_age,
                    min = -1,
                    max = 3,
                    step = 0.01),
        sliderInput("D_ext", "Diffusion strength",
                    value = sp$D_ext,
                    min = signif(sp$D_ext / 2, 2),
                    max = signif((sp$D_ext + 0.1) * 1.5, 2),
                    step = 0.0001)
    )
}

ecopathDiffusionControlTitle <- "Cohorts"
ecopathDiffusionControlDescription <-
    "Control the timing and concentration of spawning, the annulus formation and the growth diffusion strength used when fitting age-at-size data."
