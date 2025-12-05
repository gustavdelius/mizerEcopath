#' A single slider that adjusts both `h` and `gamma`
#'
#' @param input Reactive holding the inputs
#' @param output Reactive holding the outputs
#' @param session Shiny session
#' @param params Reactive value holding updated MizerParams object
#' @param params_old Reactive value holding non-updated MizerParams object
#' @param flags Environment holding flags to skip certain observers
#' @param ... Unused
growthControl <- function(input, output, session, params, params_old, flags, ...) {

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

        # Update reactive params object
        tuneParams_update_species(sp, p, params, params_old)
    },
    ignoreInit = TRUE,
    ignoreNULL = TRUE)

}

#' @rdname growthControl
#' @param params The MizerParams object currently being tuned.
#' @param input Reactive holding the inputs
growthControlUI <- function(p, input) {
    sp <- p@species_params[input$sp, ]
    tagList(
        tags$h3("Encounter"),

        popify(
            sliderInput("Eiw", "External encounter coefficient 'Eiw'",
                        value = 1,  # starting value, could use initial_Eiw
                        min = 0,
                        max = 200,
                        step = 0.01,
                        ticks = FALSE),
            title = "External encounter",
            content = "This slider scales the external encounter rate across all size bins while keeping the allometric exponent n constant."
        )
    )
}
