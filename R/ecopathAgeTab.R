#' Ecopath age tab
#'
#' @inheritParams ecopathDeathTab
#' @param age_at_length A data frame with age-at-length observations.
#' @keywords internal
ecopathAgeTab <- function(input, output, session, params, logs,
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

#' @rdname ecopathAgeTab
#' @keywords internal
ecopathAgeTabUI <- function(...) {
    tagList(
        plotOutput("plotAge")
    )
}
