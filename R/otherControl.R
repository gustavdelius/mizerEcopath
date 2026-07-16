#' Other parameters control for tuning gadget
#'
#' This control provides sliders for species parameters that do not fit into
#' the other categories. It includes:
#' *   **ks and p**: Metabolic rate coefficient and exponent. Changing `p`
#'     automatically adjusts `ks` to keep the metabolic rate at maturity
#'     size constant.
#' *   **k**: Activity rate coefficient.
#' *   **z_ext**: External mortality at 1g. Changing this scales
#'     the external mortality at all other sizes by the same factor.
#' *   **alpha**: Assimilation efficiency.
#'
#' @inheritParams abundanceControl
#' @family gadget controls
#' @export
otherControl <- function(input, output, session, params, params_old,
                         flags, ...) {
    observe({
        req(input$alpha, input$ks, input$k, input$z_ext)
        p <- isolate(params())
        sp <- isolate(input$sp)
        if (!identical(sp, flags$sp_old_other)) {
            flags$sp_old_other <- sp
            return()
        }
        # determine external mortality at 1g
        z_ext <- ext_mort(p)[input$sp, 1] / p@w[1]^p@species_params[sp, "d"]
        
        # Update slider min/max so that they are a fixed proportion of the
        # parameter value
        updateSliderInput(session, "ks",
                          min = signif(input$ks / 2, 2),
                          max = signif((input$ks + 0.1) * 1.5, 2))
        updateSliderInput(session, "k",
                          min = signif(input$k / 2, 2),
                          max = signif((input$k + 0.1) * 1.5, 2))
        
        if (z_ext != input$z_ext) {
            updateSliderInput(session, "z_ext",
                              min = signif(input$z_ext / 2, 2),
                              max = signif((input$z_ext + 0.1) * 1.5, 2))
            # re-calculate ext_mort.
            ext_mort(p)[sp, ] <- ext_mort(p)[sp, ] * (input$z_ext / z_ext)
        }

        p@species_params[sp, "alpha"] <- input$alpha
        p@species_params[sp, "ks"]    <- input$ks
        p@species_params[sp, "k"]     <- input$k
        p <- setMetabolicRate(p)
        tuneParams_update_species(sp, p, params, params_old)
    })

    observeEvent(
        input$p,
        {
            p <- params()
            sp <- input$sp

            # change ks so that metabolic rate at maturity stays the same
            p@species_params[[sp, "ks"]] <- p@species_params[[sp, "ks"]] *
                p@species_params[[sp, "w_mat"]] ^
                (p@species_params[[sp, "p"]] - input$p)
            p@species_params[[sp, "p"]] <- input$p
            ks <- p@species_params[[sp, "ks"]]
            updateSliderInput(session, "ks",
                              value = ks, # this will trigger the other observer
                              min = signif(ks / 2, 2),
                              max = signif((ks + 0.1) * 1.5, 2))
            updateSliderInput(session, "p",
                              min = signif(input$p - 0.1, 2),
                              max = signif(input$p + 0.1, 2))
        },
        ignoreInit = TRUE)
}

#' @rdname otherControl
#' @inheritParams abundanceControlUI
otherControlUI <- function(p, input) {
    sp <- p@species_params[input$sp, ]
    # determine external mortality at 1g
    z_ext <- ext_mort(p)[input$sp, 1] / p@w[1]^sp$d
    tagList(
        sliderInput("ks", "Coefficient of standard metabolism 'ks'",
                    value = sp$ks,
                    min = signif(sp$ks / 2, 2),
                    max = signif((sp$ks + 0.1) * 1.5, 2),
                    step = 0.05),
        sliderInput("p", "Exponent of metabolism 'p'",
                     value = sp[["p"]],
                     min = sp[["p"]] - 0.1, max = sp[["p"]] + 0.1, 
                     step = 0.005),
        sliderInput("k", "Coefficient of activity 'k'",
                    value = sp$k,
                    min = signif(sp$k / 2, 2),
                    max = signif((sp$k + 0.1) * 1.5, 2),
                    step = 0.01),
        tags$h3(tags$a(id = "ext_mort"), "Mort"),
        sliderInput("z_ext", "External mortality at 1g",
                    value = z_ext,
                    min = signif(z_ext / 2, 2),
                    max = signif((z_ext + 0.1) * 1.5, 2),
                    step = 0.05),
        sliderInput("alpha", "Assimilation efficiency 'alpha'",
                    value = sp$alpha,
                    min = 0,
                    max = 1)
    )
}

otherControlTitle <- "Other"
