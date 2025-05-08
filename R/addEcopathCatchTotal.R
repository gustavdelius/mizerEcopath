#' Sets up gear parameters from Ecopath catch data
#'
#' Uses the Catch data frame exported by Ecopath to set up gear parameters.
#' Currently this sets up only a single gear for the total catch. The
#' `yield_observed` gear parameter is set to the total catch for each species.
#' The `catchability` is set to the ratio of the observed yield to the observed
#' biomass. The selectivity of the gear is set to be sigmoidal with 50% of
#' fish selected at maturity size and 25% selected at 90% of maturity size
#' (measured in weight). Fishing is turned on with an initial effort of 1. The
#' result will therefore not be at steady state yet.
#'
#' The catchability set by this function for each gear will be an
#' underestimation of the catchability needed to produce the observed yield
#' because setting it to the ratio of the observed yield to the observed biomass
#' neglects that some of the biomass will be at sizes that are not fully
#' selected by the gear. You will need to call `matchCatch()` to adjust the
#' catchability to match the observed yield.
#'
#' @param params A MizerParams object
#' @param ecopath_catch The Ecopath Catch data frame
#' @return A MizerParams object with gear parameters set up and fishing effort
#'   switched on.
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
