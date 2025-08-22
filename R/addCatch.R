#' Sets up gear parameters from scaled catch data
#'
#' Uses scaled catch data from landings and survey to set up gear parameters
#' in a MizerParams object. Each gear has a sigmoidal selectivity with 50% of
#' individuals selected at maturity weight, and 25% selected at 90% of maturity
#' weight (measured in length).
#'
#' `yield_observed` is set from the data, and `catchability` is estimated as
#' yield divided by observed biomass. This initial estimate is usually too low,
#' since not all biomass is fully selected by the gear. Use `matchCatch()` to
#' refine catchability after running this function.
#'
#' Fishing is turned on with an initial effort of 1.
#'
#' @param params A MizerParams object.
#' @param landings A data frame with columns `species`, `gear`, `biomass`.
#' @param survey A data frame with columns `species`, `gear`, `biomass`.
#' @param step sets effort depending on step, landing effort set to 0 if step =
#' 1 and 1 if step = 2 or 3.
#'
#' @return A MizerParams object with updated gear parameters and effort turned on.
#' @export
addCatch <- function(params, landings, survey,step) {
    sp <- params@species_params

    if (!hasName(sp, "ecopath_groups") || !hasName(sp, "biomass_observed")) {
        stop("You must use `addEcopathParams()` before calling `addCatch()`.")
    }

    # Helper to create gear param structure
    create_gear_df <- function(gear_names) {
        bind_rows(lapply(gear_names, function(gear) {
            data.frame(
                species = sp$species,
                gear = gear,
                sel_func = "sigmoid_length",
                l50 = w2l(sp$w_mat, params),
                l25 = w2l(sp$w_mat, params) * 0.9,
                catchability = 0,
                yield_observed = 0,
                stringsAsFactors = FALSE
            )
        }))
    }

    # Create initial gear parameter templates
    gp_landings <- create_gear_df(unique(landings$gear))
    gp_survey <- create_gear_df("survey")

    # Summarize biomass for landings
    landings_summary <- landings %>%
        group_by(species, gear) %>%
        summarise(yield = sum(biomass, na.rm = TRUE), .groups = "drop")

    gp_landings <- gp_landings %>%
        left_join(landings_summary, by = c("species", "gear")) %>%
        mutate(yield_observed = coalesce(yield, 0)) %>%
        select(-yield)

    # Summarize biomass for survey
    survey_summary <- survey %>%
        group_by(species, gear) %>%
        summarise(yield = sum(biomass, na.rm = TRUE), .groups = "drop")

    gp_survey <- gp_survey %>%
        left_join(survey_summary, by = c("species", "gear")) %>%
        mutate(yield_observed = coalesce(yield, 0)) %>%
        select(-yield)

    # Combine landings and survey
    gp <- bind_rows(gp_landings, gp_survey)

    # Safely match biomass observed
    gp <- gp %>%
        left_join(sp %>% select(species, biomass_observed), by = "species") %>%
        mutate(
            catchability = 1
        ) %>%
        select(-biomass_observed)

    # Set gear parameters
    gear_params(params) <- validGearParams(gp, sp)

    # Set fishing effort depending on step
    if (step == 1) {
        initial_effort(params)[unique(landings$gear)] <- 1e-13
    } else if (step %in% c(2, 3)) {
        initial_effort(params)[unique(landings$gear)] <- 1
    }

    # Survey effort is always on
    initial_effort(params)["survey"] <- 1

    return(params)
}
