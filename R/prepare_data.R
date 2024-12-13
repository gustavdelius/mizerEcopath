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
#' @param species The species for which the data is to be prepared. By default
#'   the first species in the model.
#' @param catch A data frame containing the observed binned catch data. It must
#'   contain the following columns:
#'   * `length`: The start of each bin.
#'   * `dl`: The width of each bin.
#'   * `count`: The observed count for each bin.
#' @param yield_lambda A parameter that controls the strength of the penalty for
#'   deviation from the observed yield.
#'
#' @return A list with the data to be passed to the TMB objective function. If
#'   there is no catch data for the species, the function returns NULL.
#' @export
prepare_data <- function(params, species = 1, catch, yield_lambda = 1) {

    # Validate MizerParams object and extract data for the selected species ----
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

    # Validate catch data frame and extract data for the selected species ----
    catch <- valid_catch(catch, species)

    if (nrow(catch) == 0) {
        return(NULL)
    }

    max_length <- max(catch$length + catch$dl)
    max_weight <- sps$a * max_length^sps$b
    if (max_weight > sps$w_max) {
        stop("For ", species, " you have observed catches of larger weight than the `w_max` that you specified.")
    }

    # Fill in missing zero counts ----

    # Extract observed bin starts, ends, and counts
    observed_bins <- data.frame(
        bin_start = catch$length,
        bin_end = catch$length + catch$dl,
        count = catch$count)
    # Add empty bins at either end. This will have an effect only when the
    # catch data is very poor and would be matched by curves that are still
    # large at the end of the observation interval.
    if (min(catch$length) > 2) {
        observed_bins <- rbind(observed_bins,
                               data.frame(bin_start = 1,
                                          bin_end = min(catch$length),
                                          count = 0))
    }
    max_idx <- which.max(catch$length)
    max_length <- catch$length[max_idx] + catch$dl[max_idx]
    observed_bins <- rbind(observed_bins,
                           data.frame(bin_start = max_length,
                                      bin_end = max_length + catch$dl[max_idx],
                                      count = 0))

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

    # Extract counts, bin boundaries and widths ----
    counts <- full_bins$count
    l_bin_boundaries <- unique(c(full_bins$bin_start, full_bins$bin_end))
    w_bin_boundaries <- sps$a * l_bin_boundaries^sps$b
    w_bin_widths <- diff(w_bin_boundaries)

    w <- w(params)
    # Select subset of w that totally contains the observed range
    w_min <- max(w[w <= min(w_bin_boundaries)])
    w_max <- min(w[w >= max(w_bin_boundaries)])
    w_select <- w >= w_min & w <= w_max
    w <- w[w >= w_min & w <= w_max]
    dw <- dw(params)[w_select]
    l <- (w / sps$a) ^ (1 / sps$b)

    N <- initialN(params)[sp_select, w_select]
    biomass <- sum(N * w * dw)

    growth <- getEGrowth(params)[sp_select, w_select]

    # Precompute weights for interpolation
    weight_list <- precompute_weights(w_bin_boundaries, w)

    # The w_mat relevant for calculating mortality is the w just below it
    w_mat_idx <- sum(params@w < sps$w_mat)
    w_mat <- params@w[w_mat_idx]

    # Prepare data list for TMB ----
    data <- list(
        counts = counts,
        bin_index = weight_list$bin_index,
        f_index = weight_list$f_index,
        coeff_fj = weight_list$coeff_fj,
        coeff_fj1 = weight_list$coeff_fj1,
        dw = dw,
        w = w,
        l = l,
        yield = gps$yield_observed,
        biomass = biomass,
        growth = growth,
        w_mat = w_mat,
        d = sps$d,
        yield_lambda = yield_lambda
    )
    return(data)
}

#' Precompute weights for integration of density over observed bins
#'
#' We have a set of weight bins with boundaries `w_bin_boundaries` and need
#' an efficient way to integrate a probability density to determine a
#' probability for each bin. The probability density is available at the set
#' of weights given in the vector `w`. We want to use linear interpolation
#' for the density between these values. Because the values in `w` do not
#' align with the values in `w_bin_boundaries` we need to split each bin into
#' segments. In this function we want to
#' precompute the weights with which we need to add up the density values to
#' approximate the integral over each bin.
precompute_weights <- function(w_bin_boundaries, w) {
    # Precompute overlaps and weights
    num_bins <- length(w_bin_boundaries) - 1
    num_w <- length(w)

    # Initialize lists to store precomputed data
    bin_index <- c()
    f_index <- c()
    coeff_fj <- c()
    coeff_fj1 <- c()

    for (i in 1:num_bins) { # Loop over bins
        bin_start <- w_bin_boundaries[i]
        bin_end <- w_bin_boundaries[i + 1]

        for (j in 1:(num_w - 1)) { # Loop over w
            x0 <- w[j]
            x1 <- w[j + 1]

            # Check for overlap
            segment_start <- max(x0, bin_start)
            segment_end <- min(x1, bin_end)

            if (segment_start < segment_end) {
                delta_x <- segment_end - segment_start
                dx_j <- x1 - x0

                # Calculate interpolation weights
                w0_k <- (segment_start - x0) / dx_j  # Weight at segment_start
                w1_k <- (segment_end - x0) / dx_j    # Weight at segment_end

                # Coefficients for f(j) and f(j+1)
                coeff_fj_k = delta_x * (1 - (w0_k + w1_k)/2)
                coeff_fj1_k = delta_x * (w0_k + w1_k)/2

                # Store the precomputed data
                bin_index <- c(bin_index, i - 1)  # Zero-based indexing for C++
                f_index <- c(f_index, j - 1)      # Zero-based indexing for C++
                coeff_fj <- c(coeff_fj, coeff_fj_k)
                coeff_fj1 <- c(coeff_fj1, coeff_fj1_k)
            }
        }
    }

    return(list(
        bin_index = bin_index,
        f_index = f_index,
        coeff_fj = coeff_fj,
        coeff_fj1 = coeff_fj1
    ))
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

    return(catch)
}
