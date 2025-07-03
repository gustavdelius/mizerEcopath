#' Reduce Ecopath diet matrix to mizer species
#'
#' Aggregates a group-level Ecopath diet matrix to the species level defined in a mizer model,
#' using a mapping dictionary from Ecopath groups (e.g., juvenile/adult stanzas) to a single mizer species.
#' Produces a species-by-species diet matrix, including an \code{"other"} category to account
#' for diet contributions from Ecopath groups not included in the model. The contributions of each stanza
#' are first summed, each scaled by its consumption rate (Q = B × Q/B), and the result
#' is then normalised to proportions.
#'
#' This produces a species-by-species diet matrix, including an \code{"other"} category to account
#' for diet contributions from Ecopath groups not included in the model.
#'
#' This function assumes that the supplied Ecopath diet matrix contains proportions (not absolute rates).
#'
#' @param species_params A data frame with mizer species parameters.
#' @param ecopath_diet The Ecopath diet matrix, as exported by the Ecopath software.
#'   Rows correspond to prey groups, columns to predator groups.
#' @param ecopath_params The  basic parameters data frame, containing at least biomass and
#' consumption/biomass for each group.
#'
#' @return A matrix with dimnames `predator` and `prey`, where the row names are
#'   the species names of the predators and the column names are the species
#'   names of the prey, plus an extra column "other" for groups not included in the model.
#'   The matrix entries give, for each predator, the proportion of its total consumption
#'   that comes from each prey species or from other ecosystem components.
#'
#' @export
#'
#' @examples
#' # Generate a reduced diet matrix
#' diet_matrix <- reduceEcopathDiet(
#'     species_params = species_params_example,
#'     ecopath_diet = ecopath_diet_example,
#'     ecopath_params = ecopath_params_example
#' )
#'
#' # View first few rows of the result
#' head(diet_matrix)
#'
reduceEcopathDiet <- function(species_params, ecopath_diet, ecopath_params)
{
    ecopath_diet[is.na(ecopath_diet)] <- 0
    sp <- validSpeciesParams(species_params)
    species <- sp$species
    no_sp <- length(species)
    ecopath_groups <- ecopath_diet[, 2]
    ecopath_diet_reduced <- ecopath_diet[, 3:ncol(ecopath_diet)]
    ecopath_params <- validEcopathParams(
        ecopath_params,
        species_to_groups = setNames(sp$ecopath_groups, sp$species)
    )
    # Compute absolute consumption Q = B * (Q/B) for each stanza
    Q <- with(ecopath_params,
              `Biomass (t/km²)` * `Consumption / biomass (/year)`)
    names(Q) <- ecopath_params$`Group name`
    preds <- sapply(names(ecopath_diet_reduced), function(x) substr(x,
                                                                    2, nchar(x)))
    preds <- as.integer(preds)
    rownames(ecopath_diet_reduced) <- ecopath_groups
    colnames(ecopath_diet_reduced) <- ecopath_groups[preds]
    selected_groups <- unlist(sp$ecopath_groups)
    ignored_groups <- setdiff(ecopath_groups, selected_groups)
    ecopath_diet_reduced <- ecopath_diet_reduced[, selected_groups]
    dm <- array(0, dim = c(no_sp, no_sp + 1), dimnames = list(predator = species,
                                                              prey = c(species, "other")))
    for (i in seq_along(species)) {
        for (pred_group in sp$ecopath_groups[[i]]) {
            for (j in seq_along(species)) {
                for (prey_group in sp$ecopath_groups[[j]]) {
                    dm[i, j] <- dm[i, j] +
                        Q[pred_group] * ecopath_diet_reduced[prey_group, pred_group]
                }
            }
            dm[i, "other"] <- dm[i, "other"] +
                Q[pred_group] * sum(ecopath_diet_reduced[ignored_groups, pred_group])
        }
    }
    dm <- dm/rowSums(dm)
    return(dm)
}

#' Add Ecopath parameters to species parameters
#'
#' Determines the biomass, consumption and production rates for each species in
#' the ecopath model and adds these to the species parameters. Ecopath works
#' with "groups" where mizer works with "species". The `species_to_groups`
#' parameter is a named list that maps Ecopath groups to mizer species. The
#' names must be the same as the species in `species_params` and the values must
#' be the names of the corresponding Ecopath groups. Species in the mizer model
#' that are not mapped to groups in the Ecopath model are left unchanged.
#'
#' In case the Ecopath model has groups that are split into stanzas, then the
#' value in `species_to_groups` should be a vector with the names of the stanzas
#' that need to be combined. The biomass, production and consumption of these
#' stanzas are added together to give the values for the corresponding mizer
#' species.
#'
#' The species mapping from `species_to_groups` is stored in the `ecopath_groups`
#' column of the species_params data frame. This column is a list column so that
#' it can store a vector of groups in the case where a species is made up of
#' several groups.
#'
#' The biomasses are added to the species_params data frame in the
#' `biomass_observed` column. The consumption rates are put into a
#' `consumption_observed` column and the production rates are put into a
#' `production_observed` column. The names of the Ecopath groups associated to
#' each mizer species are put into the `ecopath_groups` column. This column is a
#' list column so that it can store a vector of groups in the case where a
#' species is made up of several groups.
#'
#' If `biomass_observed`, `production_observed` or `consumption_observed` columns
#' already have values differing from the Ecopath values, a warning is issued
#' and the values are overwritten.
#'
#' If a species in `species_params` is not found in `species_to_groups`, a
#' message is issued and the species is skipped.
#'
#' @param species_params A data frame with mizer species parameters
#' @param species_to_groups A named list where the names are mizer species and
#'   the values are vectors of Ecopath groups, see Details below.
#' @param ecopath_params A data frame with Ecopath parameters for each group
#'   as exported by the Ecopath software.
#'
#' @return The mizer species parameter data frame with the added columns
#'  `biomass_observed`, `consumption_observed`, `production_observed` and
#'  `ecopath_groups`.
#' @export
addEcopathParams <- function(species_params, ecopath_params,
                             species_to_groups = list()) {
    # Validate species parameters
    sp <- validGivenSpeciesParams(species_params)

    # Validate Ecopath parameters
    ecopath_params <- validEcopathParams(ecopath_params, species_to_groups)

    # Ensure necessary columns exist in species_params
    required_columns <- c("biomass_observed", "production_observed", "consumption_observed", "ecopath_groups")
    for (col in required_columns) {
        if (!col %in% colnames(sp)) {
            sp[[col]] <- if (col == "ecopath_groups") vector("list", nrow(sp)) else NA_real_
        }
    }

    # Indices of mapped species
    mapped_species_indices <- which(sp$species %in% names(species_to_groups))
    unmapped_species_indices <- setdiff(seq_len(nrow(sp)), mapped_species_indices)

    # Check if any Ecopath-related columns already have non-default values for mapped species
    already_populated <- sp[mapped_species_indices, ] %>%
        select(all_of(required_columns)) %>%
        select(-ecopath_groups) %>%
        summarise(across(everything(), ~ any(!is.na(.) & . != 0))) %>%
        unlist()

    if (any(already_populated)) {
        mapped_species <- sp$species[mapped_species_indices]
        warning(
            paste0(
                "Ecopath-related columns already contain non-default values for the following mapped species: ",
                paste(mapped_species, collapse = ", "), ". ",
                "These values will be overwritten based on the Ecopath data. ",
                "Unmapped species will retain their pre-existing values or be given defaults."
            )
        )
    }

    # Only overwrite values for mapped species
    for (i in mapped_species_indices) {
        species <- sp$species[i]
        sp$ecopath_groups[[i]] <- species_to_groups[[species]]

        # Initialise the values to zero for accumulation
        sp$biomass_observed[i] <- 0
        sp$consumption_observed[i] <- 0
        sp$production_observed[i] <- 0

        for (group in species_to_groups[[species]]) {
            # Extract Ecopath estimates for the group
            estimates <- ecopath_params[ecopath_params$`Group name` == group, ]

            if (nrow(estimates) > 0) {
                # Accumulate biomass, consumption, and production
                biomass <- estimates$`Biomass (t/km²)`
                sp$biomass_observed[i] <- sp$biomass_observed[i] + biomass

                consumption <- estimates$`Consumption / biomass (/year)` * biomass
                sp$consumption_observed[i] <- sp$consumption_observed[i] + consumption

                production <- estimates$`Production / consumption (/year)` * consumption
                sp$production_observed[i] <- sp$production_observed[i] + production
            }
        }
    }

    # Issue a warning if there are unmapped species
    if (length(unmapped_species_indices) > 0) {
        unmapped_species <- sp$species[unmapped_species_indices]
        warning(
            paste0(
                "The following species were not found in species_to_groups and were skipped: ",
                paste(unmapped_species, collapse = ", "), "."
            )
        )
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
#' The function also deals with the fact that the column names in the Ecopath
#' data frame can be different from the expected names if it was loaded in with
#' `read.csv` instead of `readr::read_csv`.
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
    column_mappings <- list(
        "Biomass..t.km.." = "Biomass (t/km²)",
        "Biomass..t.km." = "Biomass (t/km²)",
        "Biomass..t.km.2." = "Biomass (t/km²)",
        "Consumption...biomass...year." = "Consumption / biomass (/year)",
        "Production...consumption...year." = "Production / consumption (/year)",
        "X" = "...1",
        "Group.name" = "Group name"
    )

    # Rename columns dynamically based on the mappings
    for (old_name in names(column_mappings)) {
        new_name <- column_mappings[[old_name]]
        if (hasName(ecopath_params, old_name)) {
            ecopath_params <- ecopath_params |> rename(!!new_name := !!sym(old_name))
        }
    }

    required_cols <- unique(column_mappings)
    if (!all(required_cols %in% names(ecopath_params))) {
        stop("ecopath_params must have columns ",
             paste(required_cols, collapse = ", "))
    }

    # Remove rows that are just header rows to the stanza groups
    ecopath_params <- ecopath_params[!is.na(ecopath_params$`...1`), ]
    # Check that the group names are now unique
    if (length(unique(ecopath_params$`Group name`)) != nrow(ecopath_params)) {
        stop("Group names in ecopath_params must be unique.")
    }

    # Check that all groups are included
    required_groups <- species_to_groups |> unlist()
    missing_groups <- setdiff(required_groups, ecopath_params$`Group name`)
    if (length(missing_groups) > 0) {
        stop("The following groups in species_to_groups are not included in ecopath_params: ",
             missing_groups)
    }

    return(ecopath_params)
}
