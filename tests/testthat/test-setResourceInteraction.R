test_that("setResourceInteraction preserves resource levels and modifies encounters", {
    # Create a simple test params object
    params <- newMultispeciesParams(
        species_params = data.frame(
            species = "species1",
            w_inf = 100,
            interaction_resource = 0.5,
            stringsAsFactors = FALSE
        )
    ) |> suppressWarnings() |> suppressMessages()

    # Store initial values
    initial_resource <- resource_level(params)
    initial_ext_encounter <- params@ext_encounter
    initial_interaction <- params@species_params$interaction_resource

    # Apply the function
    params_new <- setResourceInteraction(params)

    # Test resource level is preserved
    expect_equal(resource_level(params_new), initial_resource)

    # Test interaction_resource increased
    expect_true(params_new@species_params$interaction_resource > initial_interaction)

    # Test external encounter rate decreased
    expect_true(all(params_new@ext_encounter <= initial_ext_encounter))
})

test_that("getResourceEncounterRate calculates encounter rates correctly", {
    # Create test params with known values
    params <- newMultispeciesParams(
        species_params = data.frame(
            species = "species1",
            w_inf = 100,
            interaction_resource = 1,
            stringsAsFactors = FALSE
        )
    )

    # Get encounter rates
    encounter_rates <- getResourceEncounterRate(params)

    # Test basic properties
    expect_true(is.matrix(encounter_rates))
    expect_equal(dim(encounter_rates)[1], nrow(params@species_params))

    # Test scaling with interaction_resource
    params2 <- params
    params2@species_params$interaction_resource <- 2
    encounter_rates2 <- getResourceEncounterRate(params2)

    # Should scale linearly with interaction_resource
    expect_equal(encounter_rates2, 2 * encounter_rates)
})

test_that("setResourceInteraction maintains steady state", {
    params <- newMultispeciesParams(
        species_params = data.frame(
            species = "species1",
            w_inf = 100,
            interaction_resource = 0.5,
            stringsAsFactors = FALSE
        )
    )

    # Run to steady state
    params <- steady(params)
    initial_n <- params@initial_n
    initial_n_pp <- params@initial_n_pp

    # Apply setResourceInteraction
    params_new <- setResourceInteraction(params)

    # Run again to steady state
    params_new <- steady(params_new)

    # Test that abundances remain the same
    expect_equal(params_new@initial_n, initial_n, tolerance = 1e-10)
    expect_equal(params_new@initial_n_pp, initial_n_pp, tolerance = 1e-10)
})
