test_that("prepare_log_likelihood works correctly", {
    # Sample observed data for testing
    observed_data <- data.frame(
        length = c(1.0, 3.0, 6.0, 8.0),  # Bin starts
        dl = c(1.0, 1.5, 0.5, 2.0),      # Bin widths (bins of different sizes)
        count = c(10, 20, 15, 5)         # Observed counts
    )

    # Prepare the log-likelihood function
    log_likelihood_fn <- prepare_log_likelihood(observed_data)

    # Test the function with a basic PDF (Normal distribution)
    pdf_lengths <- seq(1, 10, by = 0.5)
    pdf_values <- dnorm(pdf_lengths, mean = 5, sd = 2)  # Example PDF

    # Verify the log-likelihood output is of type "double" (i.e., numeric)
    expect_type(log_likelihood_fn(pdf_lengths, pdf_values), "double")

    # Verify the result is finite (not NaN or infinite)
    log_likelihood_value <- log_likelihood_fn(pdf_lengths, pdf_values)
    expect_true(is.finite(log_likelihood_value))
})

test_that("prepare_log_likelihood handles missing counts correctly", {
    # Observed data with some missing counts (i.e., bins with zero counts)
    observed_data_missing <- data.frame(
        length = c(1.0, 4.0, 7.0),  # Some bins are missing
        dl = c(1.0, 1.5, 0.5),      # Bin widths
        count = c(10, 0, 15)        # One bin has zero counts
    )

    # Prepare the log-likelihood function
    log_likelihood_fn <- prepare_log_likelihood(observed_data_missing)

    # Use a basic PDF for testing
    pdf_lengths <- seq(1, 10, by = 0.5)
    pdf_values <- dnorm(pdf_lengths, mean = 5, sd = 2)

    # Verify that the function still runs and gives a numeric result
    expect_type(log_likelihood_fn(pdf_lengths, pdf_values), "double")

    # Verify that the log-likelihood is finite
    log_likelihood_value <- log_likelihood_fn(pdf_lengths, pdf_values)
    expect_true(is.finite(log_likelihood_value))
})

test_that("prepare_log_likelihood handles incorrect input data", {
    # Test for missing columns in the data frame
    observed_data_incorrect <- data.frame(
        start = c(1.0, 3.0, 6.0),  # Wrong column name
        width = c(1.0, 1.5, 0.5),  # Wrong column name
        total = c(10, 20, 15)      # Wrong column name
    )

    # Expect an error due to incorrect column names
    expect_error(prepare_log_likelihood(observed_data_incorrect),
                 "Data frame 'df' must contain columns 'length', 'dl', and 'count'.")

    # Test for mismatched pdf_lengths and pdf_values
    observed_data_correct <- data.frame(
        length = c(1.0, 3.0, 6.0),
        dl = c(1.0, 1.5, 0.5),
        count = c(10, 20, 15)
    )

    log_likelihood_fn <- prepare_log_likelihood(observed_data_correct)

    pdf_lengths <- seq(1, 10, by = 0.5)
    pdf_values <- dnorm(pdf_lengths, mean = 5, sd = 2)

    # Test with incorrect length of pdf_values
    expect_error(log_likelihood_fn(pdf_lengths, pdf_values[-1]),
                 "pdf_lengths and pdf_values must be vectors of the same length.")
})

test_that("prepare_log_likelihood handles edge cases", {
    # Test with a data frame with zero counts in all bins
    observed_data_zero_counts <- data.frame(
        length = c(1.0, 3.0, 6.0),
        dl = c(1.0, 1.5, 0.5),
        count = c(0, 0, 0)  # Zero counts in all bins
    )

    log_likelihood_fn <- prepare_log_likelihood(observed_data_zero_counts)

    pdf_lengths <- seq(1, 10, by = 0.5)
    pdf_values <- dnorm(pdf_lengths, mean = 5, sd = 2)

    # The log-likelihood should handle zero counts without error
    log_likelihood_value <- log_likelihood_fn(pdf_lengths, pdf_values)

    # Check that the result is a finite number
    expect_true(is.finite(log_likelihood_value))

    # Test for large inputs
    large_data <- data.frame(
        length = seq(1, 1000, by = 10),
        dl = rep(10, 100),
        count = rpois(100, lambda = 50)  # Random large counts
    )

    log_likelihood_fn_large <- prepare_log_likelihood(large_data)

    pdf_lengths_large <- seq(1, 1000, by = 1)
    pdf_values_large <- dnorm(pdf_lengths_large, mean = 500, sd = 100)

    # Test that it runs for large inputs and returns a numeric result
    log_likelihood_value_large <- log_likelihood_fn_large(pdf_lengths_large, pdf_values_large)
    expect_type(log_likelihood_value_large, "double")
    expect_true(is.finite(log_likelihood_value_large))
})
