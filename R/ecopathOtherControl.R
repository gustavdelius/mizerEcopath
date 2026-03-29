#' Ecopath-specific other parameter control
#'
#' This control adjusts the external mortality at maturity and the metabolic
#' rate. If the external mortality at maturity size is changed, then the
#' external mortality at all other sizes is scaled by the same factor.
#'
#' @param input Reactive holding the inputs
#' @param output Reactive holding the outputs
#' @param session Shiny session
#' @param params Reactive value holding updated MizerParams object
#' @param params_old Reactive value holding non-updated MizerParams object
#' @param flags Environment holding flags to skip certain observers
#' @param ... Unused
#' @family gadget controls
#' @export
ecopathOtherControl <- function(input, output, session, params, params_old,
                         flags, ...) {
    observe({
        req(input$mu_mat)
        p <- isolate(params())
        sp <- isolate(input$sp)
        if (!identical(sp, flags$sp_old_other)) {
            flags$sp_old_other <- sp
            return()
        }
        mat_idx <- sum(p@w < p@species_params[sp, "w_mat"])
        mu_mat <- ext_mort(p)[input$sp, mat_idx]

        if (mu_mat != input$mu_mat) {
            updateSliderInput(session, "mu_mat",
                              min = signif(input$mu_mat / 2, 2),
                              max = signif((input$mu_mat + 0.1) * 1.5, 2))
            if (mu_mat > 0) {
                ext_mort(p)[sp, ] <- ext_mort(p)[sp, ] * (input$mu_mat / mu_mat)
            } else {
                ext_mort(p)[sp, ] <- input$mu_mat *
                    (p@w / p@w[mat_idx]) ^ p@species_params[sp, "d"]
            }
        }

        p@species_params[sp, "mu_mat"] <- input$mu_mat
        p <- setMetabolicRate(p, reset = TRUE)
        tuneParams_update_species(sp, p, params, params_old)
    })
}

#' @rdname ecopathOtherControl
#' @param params The MizerParams object currently being tuned.
#' @param input Reactive holding the inputs
#' @return A tagList with sliders for the exponents
ecopathOtherControlUI <- function(params, input) {
    sp <- params@species_params[input$sp, ]
    tagList(
        tags$h3(tags$a(id = "other"), "Other"),
        tags$h3(tags$a(id = "ext_mort"), "Mort"),
        sliderInput("mu_mat", "External mortality at maturity size",
                    value = sp$mu_mat,
                    min = signif(sp$mu_mat / 2, 2),
                    max = signif((sp$mu_mat + 0.1) * 1.5, 2),
                    step = 0.05)
    )
}
