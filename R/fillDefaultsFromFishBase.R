#' Fill Species Traits Using FishBase
#'
#' Populates key life-history traits in a data frame of species using FishBase, via the `rfishbase` package.
#'
#' This function is typically used early in model setup, to populate a table of species names with biological parameters needed for mizer models. It supports cases where only species names and scientific names are initially known.
#'
#' @param species_params A data frame containing at least a column with scientific names for each species. Other
#' columns (e.g. common names or user-defined identifiers) may be included and will be preserved. This is
#' typically a minimal table created at the start of model construction (e.g., from a species list), not yet a
#' full `species_params` object.
#' @param scientific_name_col A string giving the column name in `species_params` that holds the scientific names.
#'   Defaults to `"Scientific_name"`.
#' @param overwrite Logical (default = `FALSE`). If `TRUE`, existing values in the data frame are replaced by FishBase values.
#'   If `FALSE`, only missing values are filled.
#' @param verbose Logical (default = `TRUE`). If `TRUE`, prints a summary of filled species.
#'
#' @return The same data frame, augmented with the following columns where data are available:
#' \itemize{
#'   \item \code{a}, \code{b}: Length–weight allometric coefficients
#'   \item \code{Length}: Maximum length (cm)
#'   \item \code{w_max}: Maximum weight (calculated from Length, a, and b)
#'   \item \code{l_mat}: Length at maturity (median)
#'   \item \code{age_mat}: Age at maturity (median)
#'   \item \code{w_mat}: Weight at maturity (from l_mat, a, and b)
#'   item \code{t0}: Hypothetical age at length zero (von Bertalanffy growth, median)
#' }
#'
#' @details
#' This is a convenience function that queries FishBase through the `rfishbase` package.
#' Traits are pulled from the `estimate()`, `maturity()`, and `species()` tables. Maturity traits
#' are summarised by median within FishBase `SpecCode` groups. Derived weights are calculated from
#' length–weight relationships.
#'
#' Existing columns in the input will only be overwritten if `overwrite = TRUE`.
#'
#' @note Requires the `rfishbase` package. Install it with `install.packages("rfishbase")`.
#'
#' @seealso [rfishbase::estimate()], [rfishbase::maturity()], [rfishbase::species()]
#'
#' @examples
#' \dontrun{
#' # Minimal example with only scientific names
#' species_df <- data.frame(
#'   Scientific_name = c("Merluccius merluccius", "Scomber scombrus")
#' )
#' enriched <- fillDefaultsFromFishBase(species_df)
#'
#' # Optional: include user-defined species labels
#' species_df <- data.frame(
#'   species = c("Hake", "Mackerel"),
#'   Scientific_name = c("Merluccius merluccius", "Scomber scombrus")
#' )
#' enriched <- fillDefaultsFromFishBase(species_df)
#' }
#'
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

    # Summarise t0 (von Bert growth) by species
    vbg_med <- rfishbase::popgrowth(sci_names, fields = c("SpecCode", "to")) |>
        dplyr::filter(!is.na(to)) |>
        dplyr::group_by(SpecCode) |>
        dplyr::summarise(t0 = median(to, na.rm = TRUE), .groups = "drop")

    # Merge all FishBase traits by SpecCode
    fish_traits <- length_weight |>
        dplyr::inner_join(species_info, by = "SpecCode") |>
        dplyr::left_join(maturity_summary, by = "SpecCode") |>
        dplyr::mutate(
            w_max = a * Length^b,
            w_mat = a * l_mat^b
        ) |>
    dplyr::left_join(vbg_med, by = "SpecCode")

    # Join FishBase traits into main species_params
    species_params <- dplyr::left_join(
        species_params,
        fish_traits,
        by = "SpecCode",
        suffix = c("", ".fb")
    )

    # Fill columns conditionally
    cols <- c("a", "b", "Length", "l_mat", "age_mat", "w_max", "w_mat", "t0")
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
