test_that("setFeedingLevels works", {
    # Setup: Create a base model that does not meet the requirements
    sp_params <- species_params(NS_params)
    sp_params$n <- 0.6
    sp_params$p <- 0.7
    params <- newAllometricParams(sp_params, no_w = 200)
    params@species_params$h <- 2
    params@species_params$ks <- 2

    ### Test Prerequisites / Errors ------------------------------------------
    #Check with interacting model NS
    expect_error(setFeedingLevels(NS_params, f = 0.6, f_c = 0.2),
                 regexp = "This function only works for models where all encounter is external encounter")

    # Check exponents n and p values are the same
    expect_error(setFeedingLevels(params, f = 0.6, f_c = 0.2),
                 regexp = "Exponents n and p must be equal")

    params@species_params$n <- 0.7

    # Check that h is Inf
    expect_error(setFeedingLevels(params, f = 0.6, f_c = 0.2),
                 regexp = "h must be Inf before calling this function")

    params@species_params$h <- Inf

    # Check that ks is 0
    expect_error(setFeedingLevels(params, f = 0.6, f_c = 0.2),
                 regexp = "ks must be 0 before calling this function")

    params@species_params$ks <- 0

    # Check feeding level is between 0 and 1
    expect_error(setFeedingLevels(params, f = 1.1),
                 "Feeding level must be positive and strictly less than 1.")

    # Check feeding level is between 0 and 1
    expect_error(setFeedingLevels(params, f = 0.6, f_c = -1),
                 "Critical feeding level must be positive and strictly less than 1.")

    # Check feeding level is more than critical feeding level
    expect_error(setFeedingLevels(params, f = 0.4, f_c = 0.6),
                 regexp = "Critical feeding level must be less than the feeding level ")

    # Check for Allometric rates error (should trigger because n value has beeen changed above)
    expect_error(setFeedingLevels(params, f = 0.6, f_c = 0.2),
                 regexp = "This function only works for models made up of allometric rates.")

    # Make functioning proper model
    sp_params <- species_params(NS_params)
    sp_params$n <- 0.7
    sp_params$p <- 0.7
    params <- newAllometricParams(sp_params, no_w = 200)

    ### Test Outputs ------------------------------------------
    params_changed <- setFeedingLevels(params, f = 0.7, f_c = 0.3)

    ## Compare rates before and after test
    params_rates <- getRates(params)
    params_changed_rates <- getRates(params_changed)

    expect_equal(params_rates$e_growth, params_changed_rates$e_growth)
    expect_equal(params_rates$mort, params_changed_rates$mort)
    expect_equal(params_rates$rdd, params_changed_rates$rdd)

    # Check feeding levels
    f <- getFeedingLevel(params_changed)
    f_wanted <- f
    f_wanted[] <- 0.7
    expect_equal(f, f_wanted)

    f_c <- getCriticalFeedingLevel(params_changed)
    f_wanted <- f_c
    f_wanted[] <- 0.3
    expect_equal(f_c, f_wanted)
    })
