

#' Reduce Ecopath diet matrix to mizer species
#'
#' Given an Ecopath diet matrix and a dictionary that maps Ecopath groups to
#' mizer species, this function returns a diet matrix for the mizer species.
#'
#' @param ecopath_diet The Ecopath diet matrix, as exported by the Ecopath
#'   software.
#' @param groups_to_species A named list where the names are mizer species and the values
#'  are vectors of Ecopath groups.
#'
#' @return A matrix with dimnames `predator` and `prey`, where the row names are
#'   the species names of the predators and the column names are the species
#'   names of the prey and an extra column named "other". The matrix entries
#'   give for each predator the proportion of its diet that comes from each prey
#'   species or from other ecosystem components.
#' @export
reduceEcopathDiet <- function(ecopath_diet, groups_to_species) {
    ecopath_diet[is.na(ecopath_diet)] <- 0
    species <- names(groups_to_species)
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
    selected_groups <- unlist(groups_to_species)
    ignored_groups <- setdiff(ecopath_groups, selected_groups)
    ecopath_diet_reduced <- ecopath_diet_reduced[, selected_groups]
    # Create mizer diet matrix
    dm <- array(0, dim = c(no_sp, no_sp + 1),
                dimnames = list(predator = species,
                                prey = c(species, "other")))
    for (predator in species) {
        for (pred_group in groups_to_species[[predator]]) {
            for (prey in species) {
                for (prey_group in groups_to_species[[prey]]) {
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
#' Determines the biomass, consumption and production rates for each species
#' based on the Ecopath parameters and adds these to the species parameters.
#'
#' Ecopath works with "groups" where mizer works with "species". The
#' `groups_to_species` parameter is a named list that maps Ecopath groups to
#' mizer species. The names must be the same as the species in `species_params`
#' and the values must be the names of the corresponding Ecopath groups.
#'
#' In case the Ecopath model has groups that are split into stanzas, then the
#' biomass, production and consumption of these stanzas need to be added
#' together to give the values for the corresponding mizer species. In that case
#' the value in `groups_to_species` should be a vector with the names of the
#' stanzas that need to be combined.
#'
#' @param species_params A data frame with mizer species parameters
#' @param groups_to_species A named list where the names are mizer species and
#'   the values are vectors of Ecopath groups.
#' @param ecopath_params A data frame with Ecopath parameters for each group
#'   as exported by the Ecopath software.
#'
#' @return The mizer species parameter data frame with the added columns
#'  `biomass_observed`, `ecopath_consumption` and `ecopath_production`.
#' @export
addEcopathParams <- function(species_params, ecopath_params,
                             groups_to_species) {
    sp <- validSpeciesParams(species_params)
    ecopath_params <- validEcopathParams(ecopath_params, groups_to_species)
    # Remove rows that are just header rows to the stanza groups
    ecopath_params <- ecopath_params[!is.na(ecopath_params$X)]

    sp$biomass_observed <- 0
    sp$ecopath_production <- 0
    sp$ecopath_consumption <- 0
    for (species in sp$species) {
        for (group in groups_to_species[[species]]) {
            estimates <- ecopath_params[ecopath_params$Group.name == group, ]
            biomass <- estimates$Biomass..t.km..
            sp[species, "biomass_observed"] <-
                sp[species, "biomass_observed"] + biomass

            consumption <- estimates$Consumption...biomass...year. * biomass
            sp[species, "ecopath_consumption"] <-
                sp[species, "ecopath_consumption"] + consumption

            production <- estimates$Production...consumption...year. * consumption
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

#' Validate Ecopath parameter data frame
#'
#' Checks that the Ecopath parameter data frame has the required columns and
#' that the group names are unique and that all groups in the
#' `groups_to_species` list are included in the data frame.
#'
#' @param ecopath_params A data frame with Ecopath parameters for each group
#'   as exported by the Ecopath software.
#' @param groups_to_species A named list where the names are mizer species and
#'   the values are vectors of the Ecopath groups/stanzas making up the species.
#'
#' @return The validated Ecopath parameter data frame
#' @export
validEcopathParams <- function(ecopath_params, groups_to_species) {
    if (!is.data.frame(ecopath_params)) {
        stop("ecopath_params must be a data frame.")
    }
    required_cols <- c("X", "Group.name", "Biomass..t.km..",
                       "Consumption...biomass...year.",
                       "Production...consumption...year.")
    if (!all(required_cols %in% names(ecopath_params))) {
        stop("ecopath_params must have columns ",
             paste(required_cols, collapse = ", "))
    }

    # Check that the group names are unique
    if (length(unique(ecopath_params$Group.name)) != nrow(ecopath_params)) {
        stop("Group names in ecopath_params must be unique.")
    }

    # Check that all groups are included
    required_groups <- groups_to_species |> unlist()
    if (!all(required_groups %in% ecopath_params$Group.name)) {
        stop("Not all groups in groups_to_species are included in ecopath_params.")
    }

    return(ecopath_params)
}
