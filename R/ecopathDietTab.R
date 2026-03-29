#' Ecopath diet tab
#'
#' @inheritParams ecopathDeathTab
ecopathDietTab <- function(input, output, session, params, logs,
                    diet = NULL, ...) {

    # Plot diet ----
    output$plot_diet <- renderPlotly({
        req(input$sp)
        p <- params()
        if (!is.null(diet)) {
            p<-setFeedingLevels(params=p, f=0.6, f_c=0.2)
            p <- matchDiet(p, diet_matrix = diet)
            p<-steadySingleSpecies(p)
        }
        plotDiet(p, input$sp, xtrans = input$xtrans)
    })

}

#' @rdname ecopathDietTab
ecopathDietTabUI <- function(...) {
    tagList(
        plotlyOutput("plot_diet"),
        radioButtons("xtrans", "x-axis scale:",
                     choices = c("log10", "identity"),
                     selected = "log10", inline = TRUE)
    )
}
