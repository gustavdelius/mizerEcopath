#' Sets up gear parameters from Ecopath catch data
#'
#' Uses the Catch data frame exported by Ecopath to set up gear parameters for a
#' simplified "total" gear across all fleets. The `yield_observed` gear parameter
#' is set to the total catch for each species, and `catchability` is initially set
#' to the ratio of yield to biomass. This provides a starting point for fishing
#' pressure, but does not yet ensure steady-state or matched yield.
#'
#' The selectivity function is specified via the `sel_func` argument:
#'
#' * If `sel_func = "sigmoid_length"` (default), gear selectivity is a classic
#'   sigmoidal function based on body length. The 50% selection point is set at
#'   maturity size (`w_mat`), and 25% selection occurs at 90% of maturity size
#'   (measured in weight).
#'
#' * If `sel_func = "double_sigmoid_length"`, a dome-shaped selectivity curve is
#'   used. This is the product of two sigmoidal functions: one rising (as above)
#'   and one falling. The descending limb reaches 50% selectivity at the species'
#'   maximum length (`Length`) and 25% selectivity at 110% of that length. This
#'   pattern reflects fishing that avoids very large individuals (e.g. due to
#'   gear selectivity or behaviour).
#'
#' In both cases, fishing is switched on with an initial effort of 1. However, the
#' result is not guaranteed to be at steady state, and catchability will generally
#' need further adjustment. You should call `matchCatch()` to refine catchability
#' and match the observed Ecopath yields.
#'
#' @details
#' The `sel_func` argument specifies the form of the length-based selectivity function:
#'
#' * `"sigmoid_length"` (default): Classic sigmoidal selectivity based on length, with 50% selectivity at maturity length (`w_mat`) and 25% at 90% of that length (measured in weight).
#'
#' * `"double_sigmoid_length"`: A bell-shaped (dome-shaped) selectivity curve constructed by multiplying two sigmoid functions. The first sigmoid rises to 50% selection at maturity length (`w_mat`) and reaches 25% selection at 90% of that length. The second sigmoid falls, reaching 50% selection at the maximum observed length (`Length`) and 25% at 110% of that length. This mimics fisheries that avoid very large individuals (e.g. due to escape or avoidance behaviour).
#'
#' When `sel_func = "double_sigmoid_length"`, the additional columns `l50_right` and `l25_right` are automatically added to the gear parameters, and are calculated from `species_params$Length`.

#'
#' @param params A [`MizerParams`] object created with species parameters, typically including output from `addEcopathParams()`.
#' @param ecopath_catch The Ecopath Catch data frame
#' @param sel_func A string indicating the selectivity function to use. Either `"sigmoid_length"` (default) or `"double_sigmoid_length"`.

#' @return A MizerParams object with gear parameters set up and fishing effort
#'   switched on.
#' @seealso [matchCatch()], [gear_params()]
#' @export
addEcopathCatchTotal <- function(params, ecopath_catch, sel_func = "sigmoid_length") {
    sp <- params@species_params
    if (!hasName(sp, "ecopath_groups") ||
        !hasName(sp, "biomass_observed")) {
        stop("You must use `addEcopathParams()` first.")
    }
    ecopath_catch <- validEcopathCatch(ecopath_catch, sp)
    gp <- data.frame(
        species = sp$species,
        gear = "total",
        sel_func = sel_func,
        l50 = w2l(sp$w_mat, params),
        l25 = w2l(sp$w_mat, params) * 0.9,
        # catchability and yield_observed will be extracted from Ecopath later
        catchability = 0,
        yield_observed = 0
    )

    # If double sigmoid is requested, add right-hand selectivity
    if (sel_func == "double_sigmoid_length") {
        gp$l50_right <- sp$Length
        gp$l25_right <- sp$Length * 1.1
    }

    # Extract total yield for each species from Ecopath
    for (i in seq_len(nrow(sp))) {
        for (group in sp$ecopath_groups[[i]]) {
            yield <- ecopath_catch$`TotalCatch (t/km²/year)`[ecopath_catch$`Group name` == group]
            if (length(yield) == 0) {
                warning("No catch data found for group ", group)
            }
            gp$yield_observed[i] <- sum(gp$yield_observed[i], yield, na.rm = TRUE)
        }
    }
    # Set catchability
    gp$catchability <- gp$yield_observed / sp$biomass_observed
    # Set gear parameters
    gear_params(params) <- validGearParams(gp, sp)
    # Turn on fishing
    initial_effort(params) <- 1

    return(params)
}

#' Validate Ecopath Catch data frame
#'
#' Checks that the Ecopath Catch data frame has the required columns.
#'
#' The function also deals with the fact that the column names in the Ecopath
#' data frame can be different from the expected names if it was loaded in with
#' `read.csv` instead of `readr::read_csv`.
#'
#' @param ecopath_catch A data frame with Ecopath catch for each group
#'   as exported by the Ecopath software.
#' @param species_params The species parameters data frame
#'
#' @return The validated Ecopath catch data frame
#' @export
validEcopathCatch <- function(ecopath_catch, species_params) {
    if (!is.data.frame(ecopath_catch)) {
        stop("ecopath_catch must be a data frame.")
    }
    # Sometimes some columns have different names
    column_mappings <- list(
        "TotalCatch..t.km..year." = "TotalCatch (t/km²/year)",
        "Group.name" = "Group name"
    )

    # Rename columns dynamically based on the mappings
    for (old_name in names(column_mappings)) {
        new_name <- column_mappings[[old_name]]
        if (hasName(ecopath_catch, old_name)) {
            ecopath_catch <- ecopath_catch |> rename(!!new_name := !!sym(old_name))
        }
    }

    required_cols <- unique(column_mappings)
    if (!all(required_cols %in% names(ecopath_catch))) {
        stop("ecopath_catch must have columns ",
             paste(required_cols, collapse = ", "))
    }

    return(ecopath_catch)
}
