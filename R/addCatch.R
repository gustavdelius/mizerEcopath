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
addCatch <- function(params, landings, fishing_dead_biomass, survey, step) {
    sp <- params@species_params

    create_landing_gear_df <- function(landings,fishing_dead_biomass) {
        # Keep only unique gear-species combinations
        data_unique <- unique(landings[, c("gear", "species")])
        # Build the dataframe
        df <- data.frame(
            species = data_unique$species,
            gear = data_unique$gear)%>%
            left_join(sp, by="species")%>%
            left_join(fishing_dead_biomass, by=c("species", "gear"))%>%
            mutate(sel_func = "sigmoid_length",
                   l50 = l_mat,
                   l25 = l_mat * 0.9,
                   catchability = 0,
                   yield_observed = total_weight_dead_gear_per_area)%>%
            select(species, gear, sel_func,l50,l25,catchability, yield_observed)
        return(df)
    }

    create_survey_gear_df <- function(survey) {
        # Keep only unique gear-species combinations
        data_unique <- unique(survey[, c("gear", "species")])
        # Build the dataframe
        df <- data.frame(
            species = data_unique$species,
            gear = data_unique$gear)%>%
            left_join(sp, by=c("species"="species"))%>%
            mutate(sel_func = "sigmoid_length",
                   l50 = l_mat,
                   l25 = l_mat * 0.9,
                   catchability = 0,
                   yield_observed = 1e-6)%>%
            select(species, gear, sel_func,l50,l25,catchability, yield_observed)
        return(df)
    }

    # Create initial gear parameter templates
    gp_landings <- create_landing_gear_df(landings, fishing_dead_biomass)
    gp_survey <- create_survey_gear_df(survey)

    # Combine landings and survey
    gp <- bind_rows(gp_landings, gp_survey)

    # Safely match biomass observed
    gp <- gp %>%
        mutate(
            catchability = 1)

    # Set gear parameters
    gear_params(params) <- validGearParams(gp, sp)

    # Set fishing effort depending on step
    if (step == 1) {
        initial_effort(params)[unique(landings$gear)] <- 0
        gear_params(params)$catchability<-1e-7
    } else if (step %in% c(2, 3)) {
        initial_effort(params)[unique(landings$gear)] <- 1
        gear_params(params)$catchability<-1
    }

    # Survey effort is always on
    initial_effort(params)["survey"] <- 1

    return(params)
}
