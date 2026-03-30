#' Ecopath growth tab
#'
#' This tab displays growth-related information for the model:
#' *   **Growth curves**: Plot size against age for either all species
#'     or the selected species, optionally with `size_at_age` data.
#' *   **Consumption rate**: Plot comparing model consumption with
#'     Ecopath-estimated consumption values.
#'
#' @inheritParams ecopathDeathTab
#' @param size_at_age A data frame with columns 'species', 'age' and 'size'
#'   giving the size of individuals of each species at each age.
#' @family gadget tabs
#' @export
ecopathGrowthTab <- function(input, output, session, params, logs,
                      size_at_age = NULL, ...) {
    # Help button ----
    help_steps <- data.frame(
        element = c(NA),
        intro = c("This still needs to be written.")
    )
    observeEvent(
        input$growth_help,
        introjs(session, options = list(
            steps = help_steps)
        )
    )

    # Click ----
    observeEvent(input$growth_click, {
        if (!is.null(input$growth_click$panelvar1) &&
            input$growth_click$panelvar1 != input$sp) {
            updateSelectInput(session, "sp",
                              selected = input$growth_click$panelvar1)
        }
    })
    # Double Click ----
    observeEvent(input$growth_dblclick, {
        if (!is.null(input$growth_dblclick$panelvar1) &&
            input$growth_dblclick$panelvar1 != input$sp) {
            updateSelectInput(session, "sp",
                              selected = input$growth_dblclick$panelvar1)
        }
        if (input$all_growth == "All") {
            updateRadioButtons(session, "all_growth",
                               selected = "Selected species")
        } else {
            updateRadioButtons(session, "all_growth",
                               selected = "All")
        }
    })


    # Plot growth curves ----
    output$plotGrowthCurve <- renderPlot({
        p <- params()
        if (input$all_growth == "All") {
            plotGrowthCurves(p, species_panel = TRUE) +
                theme(text = element_text(size = 16))
        } else {
            plotGrowthCurves(p, species = input$sp, size_at_age = size_at_age) +
                theme(text = element_text(size = 16))
        }
    })

    # Plot feeding level ----
    output$plot_consumption <- renderPlotly({
        plotConsumptionVsSpecies(params())
    })
}

#' @rdname ecopathGrowthTab
ecopathGrowthTabUI <- function(...) {
    tagList(
        radioButtons("all_growth", "Show:",
                     choices = c("All", "Selected species"),
                     selected = "All", inline = TRUE),
        plotOutput("plotGrowthCurve",
                   click = "growth_click",
                   dblclick = "growth_dblclick"),
        plotlyOutput("plot_consumption")
    )
}

ecopathGrowthTabTitle <- "Growth"
