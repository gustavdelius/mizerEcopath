

#' Reduce Ecopath diet matrix to mizer species
#'
#' Given an Ecopath diet matrix and a dictionary that maps Ecopath groups to
#' mizer species, this function returns a diet matrix for the mizer species.
#'
#' @param species_params A data frame with mizer species parameters
#' @param ecopath_diet The Ecopath diet matrix, as exported by the Ecopath
#'   software.
#'
#' @return A matrix with dimnames `predator` and `prey`, where the row names are
#'   the species names of the predators and the column names are the species
#'   names of the prey and an extra column named "other". The matrix entries
#'   give for each predator the proportion of its diet that comes from each prey
#'   species or from other ecosystem components.
#' @export
reduceEcopathDiet <- function(species_params, ecopath_diet) {
    ecopath_diet[is.na(ecopath_diet)] <- 0
    sp <- validSpeciesParams(species_params)
    species <- sp$species
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
    selected_groups <- unlist(sp$ecopath_groups)
    ignored_groups <- setdiff(ecopath_groups, selected_groups)
    ecopath_diet_reduced <- ecopath_diet_reduced[, selected_groups]
    # Create mizer diet matrix
    dm <- array(0, dim = c(no_sp, no_sp + 1),
                dimnames = list(predator = species,
                                prey = c(species, "other")))
    for (i in seq_along(species)) {
        for (pred_group in sp$ecopath_groups[[i]]) {
            for (j in seq_along(species)) {
                for (prey_group in sp$ecopath_groups[[j]]) {
                    dm[i, j] <- dm[i, j] +
                        ecopath_diet_reduced[prey_group, pred_group]
                }
            }
            # Add consumption from all other groups
            dm[i, "other"] <- dm[i, "other"] +
                sum(ecopath_diet_reduced[ignored_groups, pred_group])
        }
    }
    # Convert to proportions
    dm <- dm / rowSums(dm)
    return(dm)
}

#' Add Ecopath parameters to species parameters
#'
#' Determines the biomass, consumption and production rates for each species in
#' the ecopath model and adds these to the species parameters.
#'
#' Ecopath works with "groups" where mizer works with "species". The
#' `species_to_groups` parameter is a named list that maps Ecopath groups to
#' mizer species. The names must be the same as the species in `species_params`
#' and the values must be the names of the corresponding Ecopath groups.
#'
#' In case the Ecopath model has groups that are split into stanzas, then the
#' biomass, production and consumption of these stanzas need to be added
#' together to give the values for the corresponding mizer species. In that case
#' the value in `species_to_groups` should be a vector with the names of the
#' stanzas that need to be combined.
#'
#' The biomasses are added to the species_params data frame in the
#' `biomass_observed` column. The consumption rates are put into a
#' `ecopath_consumption` column and the production rates are put into a
#' `ecopath_production` column. The names of the Ecopath groups associated to
#' each mizer species are put ino the `ecopath_groups` column. This column is a
#' list column so that it can store a vector of groups in the case where a
#' species is made up of several groups.
#'
#' @param species_params A data frame with mizer species parameters
#' @param species_to_groups A named list where the names are mizer species and
#'   the values are vectors of Ecopath groups, see Details below.
#' @param ecopath_params A data frame with Ecopath parameters for each group
#'   as exported by the Ecopath software.
#'
#' @return The mizer species parameter data frame with the added columns
#'  `biomass_observed`, `ecopath_consumption`, `ecopath_production` and
#'  `ecopath_groups`.
#' @export
addEcopathParams <- function(species_params, ecopath_params,
                             species_to_groups) {
    sp <- validSpeciesParams(species_params)
    ecopath_params <- validEcopathParams(ecopath_params, species_to_groups)

    sp$biomass_observed <- 0
    sp$ecopath_production <- 0
    sp$ecopath_consumption <- 0
    sp$ecopath_groups <- vector("list", nrow(sp))
    for (i in seq_len(nrow(sp))) {
        species <- sp$species[i]
        sp$ecopath_groups[[i]] <- species_to_groups[[species]]
        for (group in species_to_groups[[species]]) {
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
#' `species_to_groups` list are included in the data frame.
#'
#' @param ecopath_params A data frame with Ecopath parameters for each group
#'   as exported by the Ecopath software.
#' @param species_to_groups A named list where the names are mizer species and
#'   the values are vectors of the Ecopath groups/stanzas making up the species.
#'
#' @return The validated Ecopath parameter data frame
#' @export
validEcopathParams <- function(ecopath_params, species_to_groups) {
    if (!is.data.frame(ecopath_params)) {
        stop("ecopath_params must be a data frame.")
    }
    # Sometimes some columns have different names
    wrong <- colnames(ecopath_params) == "Biomass..t.km.2."
    if (any(wrong)) {
        colnames(ecopath_params)[wrong] <- "Biomass..t.km.."
    }
    required_cols <- c("X", "Group.name", "Biomass..t.km..",
                       "Consumption...biomass...year.",
                       "Production...consumption...year.")
    if (!all(required_cols %in% names(ecopath_params))) {
        stop("ecopath_params must have columns ",
             paste(required_cols, collapse = ", "))
    }

    # Remove rows that are just header rows to the stanza groups
    ecopath_params <- ecopath_params[!is.na(ecopath_params$X), ]
    # Check that the group names are now unique
    if (length(unique(ecopath_params$Group.name)) != nrow(ecopath_params)) {
        stop("Group names in ecopath_params must be unique.")
    }

    # Check that all groups are included
    required_groups <- species_to_groups |> unlist()
    if (!all(required_groups %in% ecopath_params$Group.name)) {
        stop("Not all groups in species_to_groups are included in ecopath_params.")
    }

    return(ecopath_params)
}
