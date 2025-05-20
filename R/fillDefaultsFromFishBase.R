#' Fill species parameters from FishBase
#'
#' Uses FishBase data to populate traits (a, b, Lmax, w_max, w_mat, age_mat, l_mat).
#'
#' @param species_params A data frame with at least a column for species name and scientific name.
#' @param scientific_name_col The column name in `species_params` that holds the scientific names.
#' @param overwrite Logical. If TRUE, replaces existing values. Default is FALSE.
#' @param verbose Logical. If TRUE, prints messages. Default is TRUE.
#' @return A `species_params` data frame with missing parameters filled in where possible.
#' @export
fillDefaultsFromFishBase <- function(species_params,
                                     scientific_name_col = "Scientific_name",
                                     overwrite = FALSE,
                                     verbose = TRUE) {

    stopifnot(is.data.frame(species_params))
    stopifnot(scientific_name_col %in% names(species_params))

    if (!requireNamespace("rfishbase", quietly = TRUE)) {
        stop("The 'rfishbase' package is required but not installed. Please install it with install.packages('rfishbase') to use this function.")
    }

    sci_names <- species_params[[scientific_name_col]]

    # Get SpecCode mapping
    spec_codes <- rfishbase::species(sci_names, fields = c("Species", "SpecCode"))
    names(spec_codes)[1] <- scientific_name_col  # Rename to match input frame

    # Merge SpecCode into species_params
    species_params <- dplyr::left_join(species_params, spec_codes, by = scientific_name_col)

    # Pull traits with SpecCode
    length_weight <- rfishbase::estimate(sci_names)
    maturity <- rfishbase::maturity(sci_names)
    species_info <- rfishbase::species(sci_names, fields = c("SpecCode", "Length"))

    # Summarise maturity by SpecCode (not Species!)
    maturity_summary <- maturity |>
        dplyr::filter(!is.na(tm), !is.na(Lm)) |>
        dplyr::group_by(SpecCode) |>
        dplyr::summarise(
            age_mat = median(tm, na.rm = TRUE),
            l_mat = median(Lm, na.rm = TRUE),
            .groups = "drop"
        )

    # Merge all FishBase traits by SpecCode
    fish_traits <- length_weight |>
        dplyr::inner_join(species_info, by = "SpecCode") |>
        dplyr::left_join(maturity_summary, by = "SpecCode") |>
        dplyr::mutate(
            w_max = a * Length^b,
            w_mat = a * l_mat^b
        )

    # Join FishBase traits into main species_params
    species_params <- dplyr::left_join(
        species_params,
        fish_traits,
        by = "SpecCode",
        suffix = c("", ".fb")
    )

    # Fill columns conditionally
    cols <- c("a", "b", "Length", "l_mat", "age_mat", "w_max", "w_mat")
    for (col in cols) {
        filled <- paste0(col, ".fb")
        if (filled %in% names(species_params)) {
            species_params[[col]] <- ifelse(
                overwrite | is.na(species_params[[col]]),
                species_params[[filled]],
                species_params[[col]]
            )
            species_params[[filled]] <- NULL
        }
    }

    if (verbose) {
        n_filled <- sum(!is.na(species_params$w_max))
        message("FishBase values filled for ", n_filled, " species.")
    }

    return(species_params)
}
