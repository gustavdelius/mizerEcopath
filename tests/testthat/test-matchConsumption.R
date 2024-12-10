test_that("matchConsumption throws error for non-MizerParams input", {
    expect_error(matchConsumption("not_params"),
                 "params must be a MizerParams object.")
})

test_that("matchConsumption behaves correctly when ecopath_consumption missing", {
    # Create a copy of celtic_params with some NAs in ecopath_consumption
    # The corresponding species should be ignored
    params_no_ecopath <- celtic_params
    params_no_ecopath@species_params$ecopath_consumption[1] <- NA
    result <- matchConsumption(params_no_ecopath, species = 1)
    # The result should be the same as the input
    expect_identical(result, params_no_ecopath)

    # Remove the ecopath_consumption column
    params_no_ecopath@species_params$ecopath_consumption <- NULL
    expect_error(matchConsumption(params_no_ecopath),
                 "You must provide the ecopath_consumption species parameter.")
})

test_that("matchConsumption throws error if exponent mismatch", {
    # Create a scenario where one species has p != n
    params_mismatch <- celtic_params
    # Assume we have at least one species; set mismatch for the first selected species
    params_mismatch@species_params$p[1] <- params_mismatch@species_params$n[1] + 0.1
    expect_error(matchConsumption(params_mismatch,
                                  species = params_mismatch@species_params$species[1]),
                 "The encounter rate and metabolic respiration rate must have the same exponent.")
})

test_that("matchConsumption works with single species", {
    # Select a single species, ensure it runs and returns a MizerParams object
    single_sp <- celtic_params@species_params$species[1]
    result <- matchConsumption(celtic_params, species = single_sp)
    # And doing it again makes no further changes
    result2 <- matchConsumption(celtic_params, species = single_sp)
    result2@time_modified <- result@time_modified
    expect_identical(result, result2)
})

test_that("matchConsumption works with multiple species and adjusts correctly", {
    # Run for multiple species (all species in celtic_params)
    all_spp <- celtic_params@species_params$species
    result <- matchConsumption(celtic_params, species = all_spp)
    expect_s4_class(result, "MizerParams")
    # Check that ks was updated
    expect_true(all(!is.na(result@species_params$ks)))
    # Check that metab and ext_encounter arrays were updated
    expect_true(all(dim(result@metab) == dim(celtic_params@metab)))
    expect_true(all(dim(ext_encounter(result)) == dim(ext_encounter(celtic_params))))
})

test_that("matchConsumption warns when negative metabolic respiration required", {
    # Create a scenario where production > ecopath_consumption for a species
    params_negative <- celtic_params
    # Assume the first species can be manipulated
    sp_idx <- 1
    # Make production greater than ecopath consumption by reducing ecopath_consumption
    params_negative@species_params$ecopath_consumption[sp_idx] <-
        getTotalProduction(params_negative)[sp_idx] / params_negative@species_params$alpha[sp_idx] * 0.5

    expect_warning(matchConsumption(params_negative,
                                    species = params_negative@species_params$species[sp_idx]),
                   "Negative metabolic respiration required for species")
})

test_that("matchConsumption updates ks correctly", {
    # Run matchConsumption
    sp_name <- celtic_params@species_params$species[1]
    orig_ks <- celtic_params@species_params$ks[1]
    result <- matchConsumption(celtic_params, species = sp_name)

    # After matchConsumption, ks should have changed
    new_ks <- result@species_params$ks[result@species_params$species == sp_name]
    expect_false(is.na(new_ks))
    expect_false(new_ks == orig_ks)
})

test_that("matchConsumption preserves energy for growth and reproduction", {
    # Before running matchConsumption, record the EReproAndGrowth
    orig_ERG <- getEReproAndGrowth(celtic_params)
    # Run matchConsumption for all species
    result <- matchConsumption(celtic_params)
    new_ERG <- getEReproAndGrowth(result)
    # They should be equal (within floating-point tolerance)
    expect_equal(new_ERG, orig_ERG, tolerance = 1e-8)
})

test_that("matchConsumption matches consumption to ecopath_consumption", {
    # Run matchConsumption for all species
    result <- matchConsumption(celtic_params)
    # Check that model consumption matches ecopath_consumption
    model_consumption <- unname(getConsumption(result))
    expected_consumption <- result@species_params$ecopath_consumption
    expect_equal(model_consumption, expected_consumption, tolerance = 1e-8)
})
