#' Rates tab for tuning gadget
#'
#' This tab provides a detailed look at the energy budget and mortality
#' rates for the selected species. It includes:
#' *   **Energy budget plot**: Visualizes how much energy is allocated
#'     to growth vs reproduction, compared with metabolic costs.
#' *   **Mortality rates plot**: Shows the death rate, split into the
#'     contributions from each predator species, from fishing and from
#'     external mortality, see [mizerExperimental::plotDeath()].
#'
#' @inheritParams biomassTab
#' @family gadget tabs
#' @export
ratesTab <- function(input, output, session, params, logs, ...) {
    # Plot growth rates ----
    output$plotGrowth <- renderPlotly({
        req(input$sp)
        plotEnergyBudget(params(),
                         species = input$sp,
                         logarithmic = (input$axis == "Logarithmic"))
    })

    # Plot death rates ----
    output$plotDeath <- renderPlotly({
        req(input$sp)
        plotDeath(params(), species = input$sp, proportion = FALSE,
                  xtrans = ifelse(input$axis == "Logarithmic",
                                  "log10", "identity"))
    })
}

#' @rdname ratesTab
ratesTabUI <- function(...) {
    tagList(
        radioButtons("axis", "x-axis scale:",
                     choices = c("Logarithmic", "Normal"),
                     selected = "Logarithmic", inline = TRUE),
        plotlyOutput("plotGrowth"),
        plotlyOutput("plotDeath")
    )
}
