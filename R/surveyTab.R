#' Survey tab for tuning gadget
#'
#' Server logic for the survey tab in the tuning gadget. Renders a plot of the
#' survey number at length and age at length data together with the
#' corresponding model predictions.
#'
#' @inheritParams deathTab
#' @param age_at_length A data frame with age-at-length observations for the
#'   selected species. See `preprocess_length_at_age()` for expected columns.
#' @param catch Data frame holding binned observed catch data. The data can be
#'   binned either into length bins or weight bins. In the former case the data
#'   frame should have columns \code{length} and \code{dl} holding the start of
#'   the size bins in cm and the width of the size bins in cm respectively. In
#'   the latter case the data frame should have columns \code{weight} and
#'   \code{dw} holding the start of the size bins in grams and the width of the
#'   size bins in grams. The data frame also needs to have the columns
#'   \code{species} (the name of the species), \code{gear} (the name of the
#'   gear) and \code{catch} (the number of individuals of a particular species
#'   caught by a particular gear in a size bin).
#' @return Invisibly returns `NULL`. Used for its side-effect of rendering plots.
#' @keywords internal
surveyTab <- function(input, output, session, params, logs,
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
}

#' @rdname surveyTab
#' @inheritParams biomassTabUI
#' @return UI elements for the survey tab.
#' @keywords internal
surveyTabUI <- function(...) {
    tagList(
        plotlyOutput("plotSurveyCatchDist"),
        plotOutput("plotSurveyAge")
    )
}
