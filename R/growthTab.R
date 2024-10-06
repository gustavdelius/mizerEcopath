#' Growth tab for tuning gadget
#'
#' @inheritParams deathTab
#' @param size_at_age A data frame with columns 'species', 'age' and 'size'
#'   giving the size of individuals of each species at each age.
growthTab <- function(input, output, session, params, logs,
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
    # See https://shiny.rstudio.com/articles/plot-interaction-advanced.html
    observeEvent(input$growth_click, {
        if (!is.null(input$growth_click$panelvar1) &&
            input$growth_click$panelvar1 != input$sp) {
            updateSelectInput(session, "sp",
                              selected = input$growth_click$panelvar1)
        }
    })
    # Double Click ----
    # See https://shiny.rstudio.com/articles/plot-interaction-advanced.html
    observeEvent(input$growth_dblclick, {
        if (!is.null(input$growth_dblclick$panelvar1) &&
            input$growth_dblclick$panelvar1 != input$sp) {
            updateSelectInput(session, "sp",
                              selected = input$growth_dblclick$panelvar1)
        }
        # Toggle between "All" and "Selected species"
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

#' @rdname growthTab
#' @inheritParams biomassTabUI
growthTabUI <- function(...) {
    tagList(
        # popify was messing up the layout and wasn't working anyway.
        # popify(div(radioButtons("all_growth", "Show:",
        #                     choices = c("All", "Selected species"),
        #                     selected = "All", inline = TRUE),
        #            style = "width: 200px; margin: auto;"),
        #        title = "Switch views",
        #        content = "Select whether to view growth curves for all species or just for the selected species. You can also toggle this by double-clicking on the plot. Single-clicking on the plot changes the selected species without changing the view."),
        radioButtons("all_growth", "Show:",
                     choices = c("All", "Selected species"),
                     selected = "All", inline = TRUE),
        plotOutput("plotGrowthCurve",
                   click = "growth_click",
                   dblclick = "growth_dblclick"),
        plotlyOutput("plot_consumption")
    )
}
