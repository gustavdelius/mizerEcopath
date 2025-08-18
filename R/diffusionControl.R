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
        req(input$d_over_g, input$spawning_t0, input$spawning_kappa)
        p <- isolate(params())
        sp <- isolate(input$sp)
        if (!identical(sp, flags$sp_old_diffusion)) {
            flags$sp_old_diffusion <- sp
            return()
        }

        # Update slider min/max so that they are a fixed proportion of the
        # parameter value
        updateSliderInput(session, "spawning_t0",
                          value = input$spawning_t0, # this will trigger the other observer
                          min = signif(input$spawning_t0 / 2, 2),
                          max = signif((input$spawning_t0 + 0.1) * 1.5, 2))
        updateSliderInput(session, "spawning_kappa",
                          value = input$spawning_kappa, # this will trigger the other observer
                          min = signif(input$spawning_kappa / 2, 2),
                          max = signif((input$spawning_kappa + 0.1) * 1.5, 2))
        updateSliderInput(session, "d_over_g",
                          value = input$d_over_g, # this will trigger the other observer
                          min = signif(input$d_over_g / 2, 2),
                          max = signif((input$d_over_g + 0.1) * 1.5, 2))

        p@species_params[[sp, "spawning_t0"]] <- input$spawning_t0
        p@species_params[[sp, "spawning_kappa"]] <- input$spawning_kappa
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
        sliderInput("spawning_t0", "Mean time of spawning",
                    value = sp$spawning_t0 ,
                    min = signif(sp$spawning_t0 / 2, 2),
                    max = signif((sp$spawning_t0 + 0.1) * 1.5, 2),
                    step = 0.01),
        sliderInput("spawning_kappa", "Concentration of spawning",
                    value = sp$spawning_kappa,
                    min = signif(sp$spawning_kappa / 2, 2),
                    max = signif((sp$spawning_kappa + 0.1) * 1.5, 2),
                    step = 0.01),
        sliderInput("d_over_g", "Diffusion strength",
                    value = sp$d_over_g,
                    min = signif(sp$d_over_g / 2, 2),
                    max = signif((sp$d_over_g + 0.1) * 1.5, 2),
                    step = 0.01)
    )
}
