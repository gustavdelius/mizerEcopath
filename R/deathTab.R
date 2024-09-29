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

    # Plot predators ----
    output$plot_pred <- renderPlotly({
        req(input$sp)
        p <- params()
        if (!is.null(diet)) {
            p <- matchDiet(p, diet)
        }
        plotDeath(p, species = input$sp,
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
