#' Age tab for tuning gadget
#'
#' Server logic for the age-at-length tab in the tuning gadget. Renders
#' a residual heatmap comparing observed and simulated annuli counts
#' across length classes, via `plotAge()`.
#'
#' @inheritParams deathTab
#' @param age_at_length A data frame with age-at-length observations for the
#'   selected species. See `preprocess_length_at_age()` for expected columns.
#' @return Invisibly returns `NULL`. Used for its side-effect of rendering plots.
#' @keywords internal
ageTab <- function(input, output, session, params, logs,
                      age_at_length = NULL, ...) {
    # Help button ----
    help_steps <- data.frame(
        element = c(NA),
        intro = c("This still needs to be written.")
    )
    observeEvent(
        input$age_help,
        introjs(session, options = list(
            steps = help_steps)
        )
    )

    # Plot mean age ----
    output$plotAge <- renderPlot({
        p <- params()
        plotAgeLikelihood(p, species = input$sp, age_at_length = age_at_length) +
            theme(text = element_text(size = 16))
    })
}

#' @rdname ageTab
#' @inheritParams biomassTabUI
#' @return UI elements for the age tab.
#' @keywords internal
ageTabUI <- function(...) {
    tagList(
        plotOutput("plotAge")
    )
}
