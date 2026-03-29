#' Ecopath age tab
#'
#' This tab displays the age distribution of the selected species. It
#' allows the user to compare the model's age-at-length distribution
#' with provided observations (`age_at_length`) to verify growth
#' and cohort dynamics.
#'
#' @inheritParams ecopathDeathTab
#' @param age_at_length A data frame with age-at-length observations.
#' @keywords internal
#' @family gadget tabs
#' @export
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
