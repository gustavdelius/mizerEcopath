#' Ecopath survey tab
#'
#' This tab provides diagnostic plots for comparing the model with
#' survey-style observations:
#' *   **Catch size distribution**: Comparison of the model catch size
#'     spectrum with the provided `catch` observations.
#' *   **Age-at-length**: Visualization of age-at-length likelihood if
#'     `age_at_length` data is provided.
#' *   **Spawning density**: Plot of the von Mises spawning distribution
#'     based on the `spawning_mu` and `spawning_kappa` parameters.
#'
#' @inheritParams ecopathDeathTab
#' @param age_at_length A data frame with age-at-length observations.
#' @param catch Data frame holding binned observed catch data.
#' @keywords internal
#' @family gadget tabs
#' @export
ecopathSurveyTab <- function(input, output, session, params, logs,
                      age_at_length = NULL, catch = NULL, ...) {

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

    # Help button ----
    help_steps <- data.frame(
        element = c(NA),
        intro = c("This still needs to be written.")
    )
    observeEvent(
        input$survey_help,
        introjs(session, options = list(
            steps = help_steps)
        )
    )

    # Plots ----
    output$plotSurveyCatchDist <- renderPlotly({
        plotlyYieldVsSize(params(), species = req(input$sp),
                          gear = input$gear,
                          catch = catch, x_var = "Length")
    })
    output$plotSurveyAge <- renderPlot({
        p <- params()
        plotAgeLikelihood(p, species = input$sp, age_at_length = age_at_length) +
            theme(text = element_text(size = 16))
    })
    output$plotSpawningDensity <- renderPlot({
        mu <- input$spawning_mu
        kappa <- input$spawning_kappa

        dates <- seq(0, 1, length.out = 200)
        density <- spawning_density(dates, mu, kappa)

        plot(dates * 12, density, type = "l", lwd = 2, col = "lightblue4",
             xlab = "Month", ylab = "Spawning probability density",
             main = "Von Mises Spawning Distribution")
    })
}

#' @rdname ecopathSurveyTab
#' @keywords internal
ecopathSurveyTabUI <- function(...) {
    tagList(
        plotlyOutput("plotSurveyCatchDist"),
        plotOutput("plotSurveyAge"),
        plotOutput("plotSpawningDensity")
    )
}
