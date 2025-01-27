# test-matchCatch.R
library(testthat)
library(mizer)

test_that("matchCatch runs without error with valid inputs for a single species", {
    # Run the function for Hake
    result <- matchCatch(celtic_params, species = "Hake", catch = celtic_catch)
    # Running it again should not make further changes except on time-stamp
    result2 <- matchCatch(result, species = "Hake", catch = celtic_catch)
    result2@time_modified <- result@time_modified
    expect_equal(result, result2)
})

test_that("matchCatch throws error if catch data is missing required columns", {
    # Remove 'dl' column from catch data
    bad_catch_data <- celtic_catch
    bad_catch_data$dl <- NULL

    expect_error(
        matchCatch(celtic_params, species = "Hake", catch = bad_catch_data),
        "must contain columns"
    )
})

test_that("matchCatch can handle multiple species", {
        multi_species <- c("Hake", "Herring")
        result <- matchCatch(celtic_params, species = multi_species, catch = celtic_catch)

        expect_s4_class(result, "MizerParams")
        gp <- gear_params(result)
        # Expect that changes have been applied to both species
        expect_true(all(multi_species %in% gp$species))
})

test_that("matchCatch preserves biomass for the adjusted species", {
    # The documentation suggests biomass remains unchanged for the adjusted species
    biomass_before <- getBiomass(celtic_params)
    result <- matchCatch(celtic_params, species = "Hake", catch = celtic_catch)
    biomass_after <- getBiomass(result)

    # Check that biomass difference is small for Hake
    expect_equal(biomass_before["Hake"], biomass_after["Hake"], tolerance = 1e-3)
})

test_that("matchCatch updates gear selectivity parameters for the selected species", {
    gp_before <- gear_params(celtic_params)

    result <- matchCatch(celtic_params, species = "Hake", catch = celtic_catch)
    gp_after <- gear_params(result)

    sp_idx <- which(gp_after$species == "Hake")
    # Expect that some key gear parameters have changed
    expect_false(all(gp_before[sp_idx, c("l50", "l25", "catchability")] ==
                         gp_after[sp_idx, c("l50", "l25", "catchability")]))
})

test_that("matchCatch respects the yield_lambda parameter", {
    # Run with default yield_lambda
    result_default <- matchCatch(celtic_params, species = "Hake", catch = celtic_catch)
    gp_default <- gear_params(result_default)

    # Run with altered yield_lambda
    result_mod <- matchCatch(celtic_params, species = "Hake", catch = celtic_catch, yield_lambda = 10)
    gp_mod <- gear_params(result_mod)

    sp_idx <- which(gp_mod$species == "Hake")
    # Different yield_lambda should lead to different parameter estimates
    expect_false(all(gp_default[sp_idx, c("l50", "l25", "catchability")] ==
                         gp_mod[sp_idx, c("l50", "l25", "catchability")]))
})

test_that("matchCatch handles invalid species inputs gracefully", {
    # Species not in params
    expect_warning(
        matchCatch(celtic_params, species = c("Cod", "NotASpecies"),
                   catch = celtic_catch),
        "The following species do not exist: NotASpecies."
    )

    # Empty species vector
    expect_error(
        matchCatch(celtic_params, species = character(0), catch = celtic_catch),
        "No species have been selected."
    )
})
