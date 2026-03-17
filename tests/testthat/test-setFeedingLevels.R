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
    expect_error(setFeedingLevels(NS_params, feeding_level = 0.6, critical_feeding_level = 0.2),
                regexp = "This function only works for models where all encounter is external encounter")

    # Check exponents n and p values are the same
    expect_error(setFeedingLevels(params, feeding_level = 0.6, critical_feeding_level = 0.2),
        regexp = "Exponents n and p must be equal")

    params@species_params$n <- 0.7

    # Check that h is Inf
    expect_error(setFeedingLevels(params, feeding_level = 0.6, critical_feeding_level = 0.2),
                 regexp = "h must be Inf before calling this function")

    params@species_params$h <- Inf

    # Check that ks is 0
    expect_error(setFeedingLevels(params, feeding_level = 0.6, critical_feeding_level = 0.2),
                 regexp = "ks must be 0 before calling this function")

    params@species_params$ks <- 0

    # Check feeding level is between 0 and 1
    expect_error(setFeedingLevels(params, feeding_level = 1.1),
                 "Feeding level must be positive and strictly less than 1.")

    # Check feeding level is between 0 and 1
    expect_error(setFeedingLevels(params, feeding_level = 0.6, critical_feeding_level = -1),
                 "Critical feeding level must be positive and strictly less than 1.")

    # Check feeding level is more than critical feeding level
    expect_error(setFeedingLevels(params, feeding_level = 0.4, critical_feeding_level = 0.6),
                 regexp = "Critical feeding level must be less than the feeding level ")

    # Check for Allometric rates error (should trigger because n value has beeen changed above)
    expect_error(setFeedingLevels(params, feeding_level = 0.6, critical_feeding_level = 0.2),
                 regexp = "This function only works for models made up of allometric rates.")

    # Make functioning proper model
    sp_params <- species_params(NS_params)
    sp_params$n <- 0.7
    sp_params$p <- 0.7
    params <- newAllometricParams(sp_params, no_w = 200)

    ### Test Outputs ------------------------------------------
    params_changed <- setFeedingLevels(params, feeding_level = 0.7, critical_feeding_level = 0.3)

    ## Check that Eriw remains the same
    Eriw_params <- getEReproAndGrowth(params)
    Eriw_params_changed <- getEReproAndGrowth(params_changed)
    diff <- Eriw_params - Eriw_params_changed

    # Test that all differences are within [-1e-10, 1e-10]
    expect_true(all(diff >= -1e-10 & diff <= 1e-10),
                info = paste0("Max difference: ", max(abs(diff))))

    ## Compare rates before and after test
    params_rates <- getRates(params)
    params_changed_rates <- getRates(params_changed)

    # Check that mortality remains the same
    params_rates$mort
    params_changed_rates$mort
    diff <- params_rates$mort - params_changed_rates$mort
    expect_true(all(diff >= -1e-10 & diff <= 1e-10),
                info = paste0("Max difference: ", max(abs(diff))))

    # Check that rdd is the same
    params_rates$rdd
    params_changed_rates$rdd
    diff <- params_rates$rdd - params_changed_rates$rdd
    expect_true(all(diff >= -1e-7 & diff <= 1e-7),
                info = paste0("Max difference: ", max(abs(diff))))

    # Check that e growth is the same
    params_rates$e_growth
    params_changed_rates$e_growth
    diff <- params_rates$e_growth - params_changed_rates$e_growth
    expect_true(all(diff >= -1e-11 & diff <= 1e-11),
                info = paste0("Max difference: ", max(abs(diff))))
    })
