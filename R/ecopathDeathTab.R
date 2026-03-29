#' Ecopath death tab
#'
#' This tab provides diagnostic plots for mortality and production. It
#' allows the user to explore:
#' *   **Mortality density**: Breakdown of various mortality sources by size
#'     and species.
#' *   **Production vs Species**: Comparison of biomass produced across
#'     all species in the model.
#'
#' @param input Reactive holding the inputs
#' @param output Reactive holding the outputs
#' @param session Shiny session
#' @param params Reactive value holding updated MizerParams object
#' @param logs Environment holding the log of steady states.
#' @param diet A diet matrix to match the diet of the model to.
#' @param ... Unused
#' @family gadget tabs
#' @export
ecopathDeathTab <- function(input, output, session, params, logs,
                     diet = NULL, ...) {
    # Plot mortality
    output$plot_mort <- renderPlotly({
        req(input$sp)
        p <- params()
        if (!is.null(diet)) {
            p <- setFeedingLevels(params = p, f = 0.6, f_c = 0.2)
            p <- matchDiet(p, diet)
        }
        plotlyDeathX(p, species = input$sp,
                     proportion = input$death_prop == "Proportion",
                     xtrans = input$death_xtrans,
                     xvar = input$death_xvar)
    })
    # Plot Production
    output$plot_prod <- renderPlotly({
        plotProductionVsSpecies(params())
    })
}

#' @rdname ecopathDeathTab
ecopathDeathTabUI <- function(...) {
    tagList(
        plotlyOutput("plot_mort"),
        fluidRow(
            column(
                width = 4,
                radioButtons("death_prop", "Show:",
                             choices = c("Rate", "Proportion"),
                             selected = "Rate",
                             inline = TRUE)
            ),
            column(
                width = 4,
                radioButtons("death_xvar", "Show Size as:",
                             choices = c("Length", "Weight"),
                             selected = "Length", inline = TRUE)
            ),
            column(
                width = 4,
                radioButtons("death_xtrans", "x-axis scale:",
                             choices = c("log10", "identity"),
                             selected = "identity", inline = TRUE)
            )
        ),
        plotlyOutput("plot_prod")
    )
}
