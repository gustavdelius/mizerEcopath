#' Sets up gear parameters from scaled catch data
#'
#' Uses scaled catch data from landings and survey to set up gear parameters
#' in a MizerParams object. Each gear has a sigmoidal selectivity with 50% of
#' individuals selected at maturity weight, and 25% selected at 90% of maturity
#' weight (measured in length).
#'
#' `yield_observed` is set from the data and `catchability` is initialised to 1.
#' Use `matchCatch()` to refine catchability after running this function.
#'
#' Survey effort is turned on with an initial effort of 1.
#'
#' @param params A MizerParams object.
#' @param landings A data frame with columns `species`, `gear`, and (if
#'   `fishing_dead_biomass` is not provided) `biomass`, which is summed per
#'   species and gear to compute `yield_observed`.
#' @param fishing_dead_biomass Optional data frame with columns `species`,
#'   `gear`, `total_weight_dead_gear_per_area`. If provided, this is used as
#'   `yield_observed` instead of summing `biomass` from `landings`.
#' @param survey A data frame with columns `species` and `gear`, used only to
#'   determine which gear-species combinations to add. `yield_observed` is set
#'   to a negligible value (1e-6) for all survey gears.
#'
#' @return A MizerParams object with updated gear parameters.
#' @export
addCatch <- function(params, landings, survey, fishing_dead_biomass = NULL) {
    sp <- params@species_params

    # Determine l50 base: use l_mat if available, otherwise convert w_mat
    l_mat_vals <- if (hasName(sp, "l_mat")) {
        setNames(sp$l_mat, sp$species)
    } else {
        setNames(w2l(sp$w_mat, params), sp$species)
    }

    if (!is.null(fishing_dead_biomass)) {
        landing_yield <- fishing_dead_biomass |>
            rename(yield_observed = total_weight_dead_gear_per_area) |>
            select(species, gear, yield_observed)
    } else {
        landing_yield <- landings |>
            group_by(species, gear) |>
            summarise(yield_observed = sum(biomass, na.rm = TRUE),
                      .groups = "drop")
    }

    survey_yield <- survey |>
        distinct(species, gear) |>
        mutate(yield_observed = 1e-6)

    # Create initial gear parameter templates
    gp_landings <- create_gear_df(landings, l_mat_vals, landing_yield)
    gp_survey <- create_gear_df(survey, l_mat_vals, survey_yield)

    # Combine landings and survey
    gp <- bind_rows(gp_landings, gp_survey) |>
        mutate(catchability = 1)

    # Set gear parameters
    gear_params(params) <- validGearParams(gp, sp)

    # Survey effort is always on
    initial_effort(params)["survey"] <- 1

    return(params)
}

# Builds a gear_params data frame from a source data frame and pre-computed
# yield_observed values. `yield_df` must have columns `species`, `gear`,
# `yield_observed`. Missing combinations are filled with 0.
create_gear_df <- function(data, l_mat_vals, yield_df) {
    unique(data[, c("gear", "species")]) |>
        data.frame(stringsAsFactors = FALSE) |>
        mutate(
            l50 = l_mat_vals[species],
            l25 = l50 * 0.9,
            sel_func = "sigmoid_length",
            catchability = 0
        ) |>
        left_join(yield_df, by = c("species", "gear")) |>
        mutate(yield_observed = coalesce(yield_observed, 0)) |>
        select(species, gear, sel_func, l50, l25, catchability, yield_observed)
}
