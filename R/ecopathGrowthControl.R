#' Ecopath-specific encounter control
#'
#' This control allows for manual adjustment of the external encounter rate
#' coefficient `Eiw` for a species. Increasing this value increases the
#' available food from outside the model, thereby boosting growth without
#' changing the internal predation dynamics.
#'
#' @inheritParams ecopathOtherControl
#' @family gadget controls
#' @export
ecopathGrowthControl <- function(input, output, session, params, params_old, flags, ...) {

    observeEvent(input$Eiw, {
        p <- params()
        sp <- input$sp
        sp_sel <- p@species_params$species == sp

        # Skip if species selection changed
        if (!identical(sp, flags$sp_old_Eiw)) {
            flags$sp_old_Eiw <- sp
            return()
        }

        # Update ext_encounter matrix
        ext_enc <- getExtEncounter(p)
        sps <- species_params(p)
        w_bins <- w(p)

        # Vectorized update across size bins
        ext_enc[sp_sel, ] <- input$Eiw * w_bins^sps$n[sp_sel]

        # Save back to params
        p <- setExtEncounter(p, ext_encounter = ext_enc)

        #Save new Eiw in species params
        sps$Eiw[sp_sel]<-input$Eiw

        species_params(p)<-sps

        # Update reactive params object
        tuneParams_update_species(sp, p, params, params_old)
    },
    ignoreInit = TRUE,
    ignoreNULL = TRUE)

}

#' @rdname ecopathGrowthControl
#' @param p The MizerParams object currently being tuned.
#' @param input Reactive holding the inputs
ecopathGrowthControlUI <- function(p, input) {
    sp <- p@species_params[input$sp, ]
    tagList(
        tags$h3("Encounter"),

        popify(
            sliderInput("Eiw", "External encounter coefficient 'Eiw'",
                        value = sp$Eiw,
                        min = 0,
                        max = 200,
                        step = 0.01,
                        ticks = FALSE),
            title = "External encounter",
            content = "This slider scales the external encounter rate across all size bins while keeping the allometric exponent n constant."
        )
    )
}
