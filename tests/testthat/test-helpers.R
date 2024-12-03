# Sample data for tests
species_params <- data.frame(
    species = c("Species1", "Species2", "Species3"),
    w_max = c(100, 200, 300),  # Required by validSpeciesParams()
    biomass_observed = c(NA, NA, NA),
    ecopath_production = c(NA, NA, NA),
    ecopath_consumption = c(NA, NA, NA),
    stringsAsFactors = FALSE
)

ecopath_params <- data.frame(
    X = c(1, 2, 3, 4),
    Group.name = c("Group1", "Group2", "Stanza1", "Stanza2"),
    Biomass..t.km.. = c(10, 20, 5, 15),
    Consumption...biomass...year. = c(2, 1.5, 3, 4),
    Production...consumption...year. = c(0.5, 0.6, 0.7, 0.8),
    stringsAsFactors = FALSE
)

species_to_groups <- list(
    Species1 = "Group1",
    Species2 = c("Stanza1", "Stanza2")
)

test_that("Function adds Ecopath parameters correctly", {
    expect_warning(result <- addEcopathParams(species_params, ecopath_params, species_to_groups),
                   "The following species were not found in species_to_groups and were skipped: Species3\\.")

    # Check added columns
    expect_true(all(c("biomass_observed", "ecopath_consumption", "ecopath_production", "ecopath_groups") %in% colnames(result)))

    # Check values for Species1
    expect_equal(result$biomass_observed[1], 10)
    expect_equal(result$ecopath_consumption[1], 20)
    expect_equal(result$ecopath_production[1], 10)

    # Check values for Species2
    expect_equal(result$biomass_observed[2], 20)  # 5 + 15
    expect_equal(result$ecopath_consumption[2], 75)  # (3*5) + (4*15)
    expect_equal(result$ecopath_production[2], 58.5)  # (0.7*15) + (0.8*60)

    # Check unmapped species
    expect_true(is.na(result$biomass_observed[3]))
    expect_true(is.na(result$ecopath_consumption[3]))
    expect_true(is.na(result$ecopath_production[3]))
    expect_null(result$ecopath_groups[[3]])
})

test_that("Function warns about existing non-default values", {
    species_params_with_values <- species_params
    species_params_with_values$biomass_observed <- c(1, 2, NA)

    expect_warning(
        addEcopathParams(species_params_with_values, ecopath_params, species_to_groups),
        "Ecopath-related columns already contain non-default values"
    ) |> suppressWarnings()
})

test_that("Function handles multiple stanzas correctly", {
    result <- addEcopathParams(species_params, ecopath_params, species_to_groups) |>
        suppressWarnings()
    expect_equal(result$biomass_observed[2], 20)  # Sum of Stanza1 and Stanza2 biomass
    expect_equal(result$ecopath_consumption[2], 75)  # Correct total consumption
    expect_equal(result$ecopath_production[2], 58.5)  # Correct total production
})

test_that("Function works with missing Ecopath parameters", {
    incomplete_ecopath_params <- ecopath_params[-1, ]  # Remove first group

    expect_error(
        addEcopathParams(species_params, incomplete_ecopath_params, species_to_groups),
        "Not all groups in species_to_groups are included in ecopath_params."
    )
})

test_that("Function initializes missing columns in species_params", {
    species_params_no_cols <- species_params[, c("species", "w_max"), drop = FALSE]
    result <- addEcopathParams(species_params_no_cols, ecopath_params, species_to_groups) |>
        suppressWarnings()

    expect_true(all(c("biomass_observed", "ecopath_consumption", "ecopath_production", "ecopath_groups") %in% colnames(result)))

    # Check that unmapped species have NA or NULL values
    expect_true(is.na(result$biomass_observed[3]))
    expect_true(is.na(result$ecopath_consumption[3]))
    expect_true(is.na(result$ecopath_production[3]))
    expect_null(result$ecopath_groups[[3]])
})

test_that("Function does not overwrite values for unmapped species", {
    # Set existing values for unmapped species
    species_params_with_values <- species_params
    species_params_with_values$biomass_observed[3] <- 100
    species_params_with_values$ecopath_production[3] <- 200
    species_params_with_values$ecopath_consumption[3] <- 300
    species_params_with_values$ecopath_groups[[3]] <- "ExistingGroup"

    result <- addEcopathParams(species_params_with_values, ecopath_params, species_to_groups) |>
        suppressWarnings()

    # Check that values for unmapped species are unchanged
    expect_equal(result$biomass_observed[3], 100)
    expect_equal(result$ecopath_production[3], 200)
    expect_equal(result$ecopath_consumption[3], 300)
    expect_equal(result$ecopath_groups[[3]], "ExistingGroup")
})
