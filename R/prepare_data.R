#' Prepare a TMB Objective Function for Optimizing Model Parameters
#'
#' This function returns a list with the data to be passed to the TMB objective
#' function. The data includes the observed catch data, the model parameters,
#' and some precomputed values that are used in the likelihood calculation.
#' The main preprocessing makes sure that we have a comprehensive set of bins
#' that cover the entire size range, even though there will not be observations
#' at all sizes. Missing observations should be interpreted as a 0 count.
#'
#' @param params A MizerParams object
#' @param species The species for which the objective function is to be prepared.
#' @param catch A data frame containing the observed binned catch data. It must
#'   contain the following columns:
#'   * `length`: The start of each bin.
#'   * `dl`: The width of each bin.
#'   * `count`: The observed count for each bin.
#' @param yield_lambda A parameter that controls the strength of the penalty for
#'   deviation from the observed yield.
#' @param pars A named list of starting values for the parameters that will be
#'   optimized.
#'
#' @return The objective function
#' @export
prepare_data <- function(params, species = 1, catch, yield_lambda) {

    # Validate MizerParams object
    params <- validParams(params)
    species <- valid_species_arg(params, species, error_on_empty = TRUE)
    if (length(species) > 1) {
        stop("Only one species can be updated at a time.")
    }
    sp <- species_params(params)
    sp_select <- sp$species == species
    sps <- sp[sp_select, ]

    gp <- params@gear_params
    gp_select <- gp$species == species
    gps <- gp[gp_select, ]
    if (nrow(gps) > 1) {
        stop("The code currently assumes that there is only a single gear for each species.")
    }

    # Validate catch data frame
    catch <- valid_catch(catch, species)

    # Extract observed bin starts, ends, and counts
    observed_bins <- data.frame(
        bin_start = catch$length,
        bin_end = catch$length + catch$dl,
        count = catch$count)
    # Add empty bins at either end to ensure that the full range is covered
    # if (min(catch$length) > 2) {
    #     observed_bins <- rbind(observed_bins,
    #                            data.frame(bin_start = 1,
    #                                       bin_end = min(catch$length),
    #                                       count = 0))
    # }
    # max_idx <- which.max(catch$length)
    # max_length <- catch$length[max_idx] + catch$dl[max_idx]
    # l_max <- (sps$w_max / sps$a)^(1/sps$b)
    # if (l_max - max_length > 1) {
    #     observed_bins <- rbind(observed_bins,
    #                            data.frame(bin_start = max_length,
    #                                       bin_end = l_max,
    #                                       count = 0))
    # }

    # Create a comprehensive set of bin edges covering all observed bins
    bin_edges <- sort(unique(c(observed_bins$bin_start, observed_bins$bin_end)))

    # Define full bins covering the observed range
    full_bins <- data.frame(
        bin_start = bin_edges[-length(bin_edges)],
        bin_end = bin_edges[-1],
        count = 0  # Initialize counts to zero
    )

    # Map observed counts to the corresponding bins in full_bins
    bin_key <- paste0(full_bins$bin_start, "_", full_bins$bin_end)
    observed_bin_key <- paste0(observed_bins$bin_start, "_", observed_bins$bin_end)
    matched_indices <- match(observed_bin_key, bin_key)
    full_bins$count[matched_indices] <- observed_bins$count

    # Extract counts and compute total observations
    counts <- full_bins$count

    # When we calculate the likelihood, we will need to have values of the
    # modelled catch density at all bin boundaries. We will use interpolation
    # for that purpose
    # Collect all unique bin boundaries for interpolation
    l_bin_boundaries <- unique(c(full_bins$bin_start, full_bins$bin_end))
    # Give this in terms of weights
    w_bin_boundaries <- sps$a * l_bin_boundaries^sps$b

    # Precompute bin widths
    w_bin_widths <- diff(w_bin_boundaries)

    # Interpolate EReproAndGrowth and repro_prop to all bin boundaries
    EReproAndGrowth <-
        approx(w(params), getEReproAndGrowth(params)[sp_select, ],
               xout = w_bin_boundaries)$y
    EReproAndGrowth[EReproAndGrowth < 0] <- 0
    if (anyNA(EReproAndGrowth)) {
        stop("Interpolation of EReproAndGrowth failed.")
    }
    repro_prop <-
        approx(w(params), repro_prop(params)[sp_select, ],
               xout = w_bin_boundaries)$y
    repro_prop[repro_prop > 1] <- 1
    if (anyNA(repro_prop)) {
        stop("Interpolation of repro_prop failed.")
    }

    # Calculate biomass above 2cm
    biomass <- getBiomass(params, min_w = min(w_bin_boundaries))[sp_select]

    # Prepare data
    data <- list(
        counts = counts,
        bin_widths = w_bin_widths,
        bin_boundaries = w_bin_boundaries,
        bin_boundary_lengths = l_bin_boundaries,
        yield = gps$yield_observed,
        biomass = biomass,
        EReproAndGrowth = EReproAndGrowth,
        repro_prop = repro_prop,
        w_mat = sps$w_mat,
        d = sps$d,
        yield_lambda = yield_lambda
    )
    return(data)
}

#' @export
valid_catch <- function(catch, species) {
    # Allow "catch" as an alternative name to "count"
    if ("catch" %in% names(catch)) {
        catch$count <- catch$catch
    }
    if (!all(c('length', 'dl', 'count') %in% names(catch))) {
        stop("Data frame 'catch' must contain columns 'length', 'dl', and 'count'.")
    }
    # If this contains data for several species, extract the desired species
    if ("species" %in% names(catch)) {
        catch <- catch[catch$species == species, ]
    }
    if ("gear" %in% names(catch) && length(unique(catch$gear)) > 1) {
        stop("The code currently assumes that there is only a single gear for each species.")
    }
    if (nrow(catch) == 0) {
        stop("No catch data for species ", species)
    }
    return(catch)
}
