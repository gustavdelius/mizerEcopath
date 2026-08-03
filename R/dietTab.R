#' Diet tab for tuning gadget
#'
#' This tab displays the diet composition for the selected species and
#' the prey available to it. It includes:
#' *   **Diet plot**: A plotly visualization showing what the species
#'     is eating by prey category.
#' *   **Prey availability plot**: Shows the feeding kernel, number
#'     density, and biomass density of the available prey for a
#'     given predator size, see [plotPreyAvailability()].
#' *   **Predator size slider**: Allows the user to select the predator
#'     size for which prey availability is shown.
#'
#' @inheritParams biomassTab
#' @family gadget tabs
#' @export
dietTab <- function(input, output, session, params, logs, ...) {
    
    # Plot diet ----
    output$plot_diet <- renderPlotly({
        req(input$sp)
        plotDietX(params(), input$sp, xtrans = input$xtrans)
    })
    
    # Plot prey ----
    output$plot_prey <- renderPlotly({
        plotlyPreyAvailability(params(), species = req(input$sp),
                               pred_size = 10^req(input$pred_size))
    })
    
    # Prey size slider ----
    output$pred_size_slider <- renderUI({
        p <- isolate(params())
        sp <- which.max(p@species_params$species == input$sp)
        sliderInput("pred_size", "log_10 predator size",
                    value = signif(log10(p@species_params$w_mat[sp]), 2),
                    min = signif(log10(p@species_params$w_min[sp]), 2),
                    max = signif(log10(p@species_params$w_max[sp]), 2),
                    step = 0.2,
                    width = "80%",
                    animate = animationOptions(loop = TRUE))
    })
}

#' @rdname dietTab
dietTabUI <- function(...) {
    tagList(
        plotlyOutput("plot_diet"),
        radioButtons("xtrans", "x-axis scale:",
                     choices = c("log10", "identity"),
                     selected = "log10", inline = TRUE),
        plotlyOutput("plot_prey"),
        uiOutput("pred_size_slider")
    )
}