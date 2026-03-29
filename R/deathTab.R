#' Death tab for tuning gadget
#'
#' This tab displays the breakdown of mortality for the selected species.
#' It shows:
#' *   **Mortality plot**: A plotly visualization showing different
#'     sources of mortality (predation, fishing, background) across
#'     the size spectrum.
#' *   **Controls**: Options to display mortality as either a rate or
#'     a proportion, and to switch the x-axis between log and identity scales.
#'
#' @inheritParams biomassTab
#' @family gadget tabs
#' @export
deathTab <- function(input, output, session, params, logs, ...) {
    
    # Plot predators ----
    output$plot_pred <- renderPlotly({
        req(input$sp)
        plotDeath(params(), species = input$sp, 
                  proportion = input$death_prop == "Proportion",
                  xtrans = input$death_xtrans)
    })
}

#' @rdname deathTab
deathTabUI <- function(...) {
    tagList(
        plotlyOutput("plot_pred"),
        radioButtons("death_prop", "Show",
                     choices = c("Proportion", "Rate"),
                     selected = "Proportion",
                     inline = TRUE),
        radioButtons("death_xtrans", "x-axis scale:",
                     choices = c("log10", "identity"),
                     selected = "log10", inline = TRUE)
    )
}