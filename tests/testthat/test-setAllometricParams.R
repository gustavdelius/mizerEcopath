test_that("juvenile biomass spectrum slope is near -0.2 for default parameters", {
    # Set up a simple species parameter data frame for a single species
    sp <- data.frame(
        species = "TestFish",
        w_max = 1000,        # maximum weight
        w_mat = 100,         # weight at maturity
        n = 0.7,             # default consumption exponent
        a = 0.01,            # dummy parameter for length-weight relationship
        b = 3,               # dummy parameter for length-weight relationship
        age_mat = 5,
        Length = 50
    )

    # Create the model parameters with newAllometricParams
    p <- newAllometricParams(sp, no_w = 200)

    # Get the weight grid from MizerParams
    w <- p@w

    # Extract the abundance distribution
    N <- as.numeric(initialN(p))

    # Calculate biomass density as abundance multiplied by weight
    B <- N * w

    # Define the juvenile range (weights less than w_mat)
    juvenile <- w < p@species_params$w_mat[p@species_params$species == "TestFish"]

    # Log-transform the values and fit a linear model
    log_w <- log(w[juvenile])
    log_B <- log(B[juvenile])
    fit <- lm(log_B ~ log_w)
    slope <- coef(fit)["log_w"]

    # Check that the fitted slope is near -0.2 (with a tolerance of 0.1)
    expect_equal(as.numeric(slope), -0.2, tolerance = 0.1)
})

test_that("mortality exponent d is forced to equal n - 1", {
    # Create a species parameter data frame with an intentionally wrong d value
    sp <- data.frame(
        species = "TestFish",
        w_max = 1000,
        w_mat = 100,
        n = 0.65,            # non-default n
        a = 0.01,
        b = 3,
        age_mat = 5,
        Length = 50,
        d = 42               # bogus d value that should be overridden
    )

    # Create the model parameters
    p <- newAllometricParams(sp, no_w = 200)

    # Extract the enforced mortality exponent from the model
    actual_d <- p@species_params$d[p@species_params$species == "TestFish"]
    expected_d <- p@species_params$n[p@species_params$species == "TestFish"] - 1

    # Test that the forced relationship holds
    expect_equal(actual_d, expected_d)
})
