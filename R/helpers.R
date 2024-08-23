

#' Reduce Ecopath diet matrix to mizer species
#'
#' Given an Ecopath diet matrix and a dictionary that maps Ecopath groups to
#' mizer species, this function returns a diet matrix for the mizer species.
#'
#' @param ecopath_diet The Ecopath diet matrix, as exported by the Ecopath
#'   software.
#' @param dict A named list where the names are mizer species and the values
#'  are vectors of Ecopath groups.
#'
#' @return A matrix with dimnames `predator` and `prey`, where the row names are
#'   the species names of the predators and the column names are the species
#'   names of the prey and an extra column named "other". The matrix entries
#'   give for each predator the proportion of its diet that comes from each prey
#'   species or from other ecosystem components.
#' @export
reduceEcopathDiet <- function(ecopath_diet, dict) {
    ecopath_diet[is.na(ecopath_diet)] <- 0
    species <- names(dict)
    no_sp <- length(species)
    # The second column of ecopath diet matrix has the group names.
    ecopath_groups <- ecopath_diet[, 2]
    # Keep only the columns that correspond to predators
    ecopath_diet_reduced <- ecopath_diet[, 3:ncol(ecopath_diet)]
    # Get the indices of the predators from the column names
    # Just need to strip off the leading X introduced by read.csv
    preds <- sapply(names(ecopath_diet_reduced),
                    function(x) substr(x, 2, nchar(x)))
    preds <- as.integer(preds)
    # Note that the rows are the prey and the columns are the predators
    rownames(ecopath_diet_reduced) <- ecopath_groups
    colnames(ecopath_diet_reduced) <- ecopath_groups[preds]
    # Keep only predators that correspond to our species
    selected_groups <- unlist(dict)
    ignored_groups <- setdiff(ecopath_groups, selected_groups)
    ecopath_diet_reduced <- ecopath_diet_reduced[, selected_groups]
    # Create mizer diet matrix
    dm <- array(0, dim = c(no_sp, no_sp + 1),
                dimnames = list(predator = species,
                                prey = c(species, "other")))
    for (predator in species) {
        for (pred_group in dict[[predator]]) {
            for (prey in species) {
                for (prey_group in dict[[prey]]) {
                    dm[predator, prey] <- dm[predator, prey] +
                        ecopath_diet_reduced[prey_group, pred_group]
                }
            }
            # Add consumption from all other groups
            dm[predator, "other"] <- dm[predator, "other"] +
                sum(ecopath_diet_reduced[ignored_groups, pred_group])
        }
    }
    # Convert to proportions
    dm <- dm / rowSums(dm)
    return(dm)
}

#' Add Ecopath parameters to species parameters
#'
#' Uses a dictionary for translating between Ecopath groups and mizer species
#' to calculate the biomass, consumption and production rates for each species
#' based on the Ecopath parameters and adds these to the species parameters.
#'
#' The `dict` list must have the same names as the species in `species_params`
#' and the values must be vectors of Ecopath groups that should be combined to
#' give the mizer species. This is usually used to combine Ecopath stanzas.
#' These vectors of groups are added to the `species_params` data frame as the
#' column `ecopath_groups`.
#'
#' @param species_params A data frame with mizer species parameters
#' @param dict A named list where the names are mizer species and the values
#'   are vectors of Ecopath groups.
#' @param ecopath_params A data frame with Ecopath parameters for each group
#'   as exported by the Ecopath software.
#'
#' @return The mizer species parameter data frame with the added columns
#'  `biomass_observed`, `ecopath_groups`, `ecopath_consumption` and
#'  `ecopath_production`.
#' @export
addEcopathParams <- function(species_params, dict, ecopath_params) {
    sp <- validSpeciesParams(species_params)
    ecopath_params <- ecopath_params |>
        filter(!is.na(...1))
    sp$biomass_observed <- 0
    sp$ecopath_production <- 0
    sp$ecopath_consumption <- 0
    for (species in sp$species) {
        sp[species, "ecopath_groups"] <- dict[[species]]
        for (group in dict[[species]]) {
            estimates <- ecopath_params[ecopath_params$`Group name` == group, ]
            biomass <- estimates$`Biomass (t/km²)`
            sp[species, "biomass_observed"] <-
                sp[species, "biomass_observed"] + biomass

            consumption <- estimates$`Consumption / biomass (/year)` * biomass
            sp[species, "ecopath_consumption"] <-
                sp[species, "ecopath_consumption"] + consumption

            production <- estimates$`Production / consumption (/year)` * consumption
            sp[species, "ecopath_production"] <-
                sp[species, "ecopath_production"] + production
        }
    }
    return(sp)
}

#' Make model non-interacting
#'
#' Adds the predation mortality to the external mortality and the predation
#' encounter to the external encounter and then sets the interaction matrix to
#' zero and switches off the interaction with the resource, so that the model
#' becomes non-interacting.
#'
#' @param params A MizerParams object
#' @return The modified MizerParams object
#' @export
makeNoninteracting <- function(params) {

    # Put predation mortality into external mortality
    ext_mort(params) <- ext_mort(params) + getPredMort(params)

    # Put predation encounter into external encounter
    # We make the assumption that all of the encounter rate is from
    # predation and ext_encounter. We need to do this because we do not have
    # a way to calculate the encounter specifically from predation. There is
    # no getPredEncounter() function.
    ext_encounter(params) <- getEncounter(params)

    # Set the interaction matrix to zero
    interaction_matrix(params)[] <- 0
    species_params(params)$interaction_resource <- 0

    return(params)
}
