#' Controlling the diffusion rate
#'
#'
#' @param input Reactive holding the inputs
#' @param output Reactive holding the outputs
#' @param session Shiny session
#' @param params Reactive value holding updated MizerParams object
#' @param params_old Reactive value holding non-updated MizerParams object
#' @param flags Environment holding flags to skip certain observers
#' @param ... Unused
diffusionControl <- function(input, output, session, params, params_old,
                             flags, ...) {
    observe({
        req(input$d_over_g, input$spawning_mu, input$spawning_kappa,
            input$annuli_min_age, input$annuli_date)
        p <- isolate(params())
        sp <- isolate(input$sp)
        if (!identical(sp, flags$sp_old_diffusion)) {
            flags$sp_old_diffusion <- sp
            return()
        }

        # Update slider min/max so that they are a fixed proportion of the
        # parameter value
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
        updateSliderInput(session, "d_over_g",
                          value = input$d_over_g,
                          min = signif(input$d_over_g / 2, 2),
                          max = signif((input$d_over_g + 0.1) * 1.5, 2))

        p@species_params[[sp, "spawning_mu"]] <- input$spawning_mu
        p@species_params[[sp, "spawning_kappa"]] <- input$spawning_kappa
        p@species_params[[sp, "annuli_min_age"]] <- input$annuli_min_age
        p@species_params[[sp, "annuli_date"]] <- input$annuli_date
        p@species_params[[sp, "d_over_g"]] <- input$d_over_g
        tuneParams_update_species(sp, p, params, params_old)
    })
}

#' @rdname diffusionControl
#' @param params The MizerParams object currently being tuned.
#' @param input Reactive holding the inputs
#' @return A tagList with sliders for the exponents
diffusionControlUI <- function(params, input) {
    sp <- params@species_params[input$sp, ]
    tagList(
        tags$h3(tags$a(id = "cohorts"), "Cohorts"),
        sliderInput("spawning_mu", "Mean time of spawning",
                    value = sp$spawning_mu ,
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
        sliderInput("d_over_g", "Diffusion strength",
                    value = sp$d_over_g,
                    min = signif(sp$d_over_g / 2, 2),
                    max = signif((sp$d_over_g + 0.1) * 1.5, 2),
                    step = 0.01)
    )
}
