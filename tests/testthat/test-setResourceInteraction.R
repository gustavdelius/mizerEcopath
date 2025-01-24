test_that("setResourceInteraction preserves resource levels and modifies encounters", {
    # Create a simple test params object
    params <- newMultispeciesParams(
        species_params = data.frame(
            species = "species1",
            w_inf = 100,
            interaction_resource = 0.5,
            stringsAsFactors = FALSE
        ),
    ) |> suppressWarnings() |> suppressMessages() |>
        setResource(resource_level = 0.5)

    power_encounter <- params@w ^ params@species_params$n
    ext_encounter(params) <- matrix(power_encounter, nrow = 1)

    # Store initial values
    initial_resource <- resource_level(params)
    initial_ext_encounter <- params@ext_encounter
    initial_interaction <- params@species_params$interaction_resource

    # Apply the function
    params_new <- setResourceInteraction(params)

    # Test resource level is preserved
    expect_equal(resource_level(params_new), initial_resource)

    # Test total encounter rate is preserved
    expect_equal(getEncounter(params_new), getEncounter(params))
    # Test that growth and mortality rates are preserved
    expect_equal(getEGrowth(params_new), getEGrowth(params))
    expect_equal(getMort(params_new), getMort(params))

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
    ) |> suppressWarnings() |> suppressMessages()

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
