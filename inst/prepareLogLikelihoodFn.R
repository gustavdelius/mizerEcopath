#' Prepare a Log-Likelihood Function for Optimizing Model Parameters
#'
#' This function processes the observed count data in bins and precomputes
#' data structures needed to calculate the log-likelihood efficiently. It returns
#' a function that can compute the log-likelihood given the probability density
#' function (PDF) values at specified lengths. The returned function can be used
#' with different model-generated PDFs to optimize model parameters by maximizing
#' the log-likelihood.
#'
#' @param df A data frame containing the observed binned count data. It must contain
#'   the following columns:
#'   \itemize{
#'     \item \code{length}: The start of each bin.
#'     \item \code{dl}: The width of each bin.
#'     \item \code{count}: The observed count for each bin.
#'   }
#'
#' @return A function that computes the log-likelihood given two arguments:
#'   \describe{
#'     \item{\code{pdf_lengths}}{A numeric vector of the lengths at which the PDF
#'     is evaluated.}
#'     \item{\code{pdf_values}}{A numeric vector of PDF values at those lengths,
#'     which must correspond to \code{pdf_lengths}.}
#'   }
#'   The returned function calculates the log-likelihood for the provided PDF values
#'   over the observed binned data.
#'
#' @examples
#' # Example usage of prepare_log_likelihood:
#' observed_data <- data.frame(
#'   length = c(1.0, 3.0, 6.0, 8.0),  # Bin starts
#'   dl = c(1.0, 1.5, 0.5, 2.0),      # Bin widths (bins of different sizes)
#'   count = c(10, 20, 15, 5)         # Observed counts
#' )
#'
#' # Prepare the log-likelihood function
#' log_likelihood_fn <- prepare_log_likelihood(observed_data)
#'
#' # Use the log-likelihood function with model-generated PDFs
#' pdf_lengths <- seq(1, 10, by = 0.1)
#' pdf_values <- dnorm(pdf_lengths, mean = 5, sd = 2)  # Example PDF (Normal distribution)
#'
#' # Calculate the log-likelihood
#' log_likelihood_value <- log_likelihood_fn(pdf_lengths, pdf_values)
#' print(log_likelihood_value)
#'
#' @seealso \code{\link[stats]{dnorm}}, \code{\link[base]{optim}}
#'
#' @export
prepare_log_likelihood <- function(df) {

    # Validate data frame
    if (!all(c('length', 'dl', 'count') %in% names(df))) {
        stop("Data frame 'df' must contain columns 'length', 'dl', and 'count'.")
    }

    # Extract observed bin starts, ends, and counts
    observed_bins <- data.frame(
        bin_start = df$length,
        bin_end = df$length + df$dl,
        count = df$count
    )

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
    N <- sum(counts)

    # Precompute any data-dependent constants in the log-likelihood formula
    data_log_likelihood_constant <- lgamma(N + 1) - sum(lgamma(counts + 1))

    # Collect all unique bin boundaries for interpolation
    all_bin_boundaries <- unique(c(full_bins$bin_start, full_bins$bin_end))

    # Precompute bin widths
    bin_widths <- full_bins$bin_end - full_bins$bin_start

    # Return a function that computes the log-likelihood given a PDF
    function(pdf_lengths, pdf_values) {
        # pdf_lengths: Numeric vector of lengths at which the PDF is evaluated
        # pdf_values: Numeric vector of PDF values at those lengths

        # Validate inputs
        if (length(pdf_lengths) != length(pdf_values)) {
            stop("pdf_lengths and pdf_values must be vectors of the same length.")
        }

        # Sort the PDF data by lengths
        pdf_order <- order(pdf_lengths)
        pdf_lengths <- pdf_lengths[pdf_order]
        pdf_values <- pdf_values[pdf_order]

        # Interpolate the PDF at all bin boundaries at once
        interpolated_pdf <- approx(pdf_lengths, pdf_values,
                                   xout = all_bin_boundaries, rule = 2)$y

        # Create a named vector for easy lookup
        names(interpolated_pdf) <- all_bin_boundaries

        # Get the densities at bin starts and ends
        density_start <- interpolated_pdf[as.character(full_bins$bin_start)]
        density_end <- interpolated_pdf[as.character(full_bins$bin_end)]

        # Approximate the integral over the bin using the trapezoidal rule
        bin_probs <- ((density_start + density_end) / 2) * bin_widths

        # Ensure probabilities are positive and normalize them
        if (any(bin_probs <= 0)) {
            stop("All bin probabilities must be positive.")
        }
        bin_probs <- bin_probs / sum(bin_probs)

        # Compute the log-likelihood using the multinomial distribution formula
        LL <- data_log_likelihood_constant + sum(counts * log(bin_probs))

        return(LL)
    }
}
