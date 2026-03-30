#' Ecopath catch tab
#'
#' This tab provides diagnostic plots for comparison with catch data:
#' *   **Total yield plot**: A bar-like plot showing circles (model)
#'     and squares (observed) for each species caught by the gear.
#' *   **Catch size distribution**: Detailed size breakdown of
#'     yield for the selected species and gear.
#'
#' @inheritParams ecopathDeathTab
#' @param catch Data frame holding binned observed catch data.
#' @family gadget tabs
#' @export
ecopathCatchTab <- function(input, output, session, params, logs,
                     catch = NULL, ...) {

    if (!is.null(catch)) {
        assert_that(
            is.data.frame(catch),
            "catch" %in% names(catch),
            "species" %in% names(catch),
            "gear" %in% names(catch),
            all(c("length", "dl") %in% names(catch)) |
                all(c("weight", "dw") %in% names(catch))
        )
    }

    # Catch density for selected species ----
    output$plotCatchDist <- renderPlotly({
        plotlyYieldVsSize(params(), species = req(input$sp),
                          gear = input$gear,
                          catch = catch, x_var = input$catch_x)
    })

    # Total yield by species ----
    output$plotTotalYield <- renderPlot({
        plotYieldVsSpecies(params(), gear = input$gear) +
            theme(text = element_text(size = 18))
    })
}

#' @rdname ecopathCatchTab
ecopathCatchTabUI <- function(...) {
    tagList(
        plotlyOutput("plotCatchDist"),
        radioButtons("catch_x", "Show size in:",
                     choices = c("Weight", "Length"),
                     selected = "Length", inline = TRUE),
        plotOutput("plotTotalYield")
    )
}

ecopathCatchTabTitle <- "Catch"
