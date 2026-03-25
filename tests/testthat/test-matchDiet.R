test_that("matchDiet throws error for non-MizerParams input", {
    expect_error(matchDiet("not_params"),
                 "params must be a MizerParams object.")
})

test_that("matchDiet with getDietMatrix leaves diet unchanged", {
    # Apply matchDiet once to obtain a result in a known state
    dm <- getDietMatrix(celtic_params)
    result1 <- matchDiet(celtic_params, diet_matrix = dm)
    # Apply matchDiet again using its own diet matrix
    dm1 <- getDietMatrix(result1)
    result2 <- matchDiet(result1, diet_matrix = dm1)
    result2@time_modified <- result1@time_modified
    expect_equal(result1, result2)
})

test_that("matchDiet aggregates non-species prey columns into other", {
    dm <- getDietMatrix(celtic_params)
    sp <- celtic_params@species_params$species

    # Split the first non-species column into two artificial prey columns;
    # the result should be identical to passing the original diet matrix.
    non_sp_cols <- setdiff(colnames(dm), sp)
    if (length(non_sp_cols) > 0) {
        col <- non_sp_cols[[1]]
        dm_split <- cbind(
            dm[, sp, drop = FALSE],
            exotic_prey_A = dm[, col] * 0.6,
            exotic_prey_B = dm[, col] * 0.4,
            dm[, non_sp_cols[-1], drop = FALSE]
        )
        result_original <- matchDiet(celtic_params, diet_matrix = dm)
        result_split    <- matchDiet(celtic_params, diet_matrix = dm_split)
        result_split@time_modified <- result_original@time_modified
        expect_equal(result_original, result_split)
    } else {
        skip("No non-species columns in getDietMatrix output to test aggregation")
    }
})

test_that("matchDiet normalises rows so scale of input does not matter", {
    dm <- getDietMatrix(celtic_params)
    # Scale all rows by arbitrary species-specific constants; the result
    # should be identical because rows are normalised internally.
    scale_factors <- seq(0.5, 2, length.out = nrow(dm))
    dm_scaled <- dm * scale_factors  # recycles row-wise (column-major order)
    result_original <- matchDiet(celtic_params, diet_matrix = dm)
    result_scaled   <- matchDiet(celtic_params, diet_matrix = dm_scaled)
    result_scaled@time_modified <- result_original@time_modified
    expect_equal(result_original, result_scaled)
})
