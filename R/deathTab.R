#' Serve tab with death plots
#'
#' @param input Reactive holding the inputs
#' @param output Reactive holding the outputs
#' @param session Shiny session
#' @param params Reactive value holding updated MizerParams object
#' @param logs Environment holding the log of steady states.
#' @param diet A diet matrix to match the diet of the model to.
#' @param ... Unused
deathTab <- function(input, output, session, params, logs,
                     diet = NULL, ...) {
    # Plot mortality
    output$plot_mort <- renderPlotly({
        req(input$sp)
        p <- params()
        if (!is.null(diet)) {
            p <- matchDiet(p, diet)
        }
        p <- catchSelectivity(p, catch)
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

#' @rdname deathTab
deathTabUI <- function(...) {
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
