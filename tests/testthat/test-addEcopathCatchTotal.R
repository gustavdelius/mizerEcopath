# tests/testthat/test-addEcopathCatchTotal.R
test_that("addEcopathCatchTotal sets gear params with sigmoid_length", {
    sp <- data.frame(
        species          = "TestFish",
        w_max            = 1e4,
        w_mat            = 100,
        a                = 0.01,
        b                = 3,
        age_mat          = 2,
        Length           = 50,
        biomass_observed = 500,
        stringsAsFactors = FALSE
    )

    ecop <- data.frame(
        `...1`                            = 1,
        `Group name`                      = "GroupA",
        `Biomass (t/km²)`                 = 500,
        `Consumption / biomass (/year)`   = 5,
        `Production / consumption (/year)`= 1,
        stringsAsFactors = FALSE
    )

    sp <- addEcopathParams(sp, ecop, list(TestFish = "GroupA"))
    params <- newAllometricParams(sp)

    catch_df <- data.frame(
        `Group name`              = "GroupA",
        `TotalCatch (t/km²/year)` = 100,
        stringsAsFactors          = FALSE
    )
    p <- addEcopathCatchTotal(params, catch_df, sel_func = "sigmoid_length")

    gp <- gear_params(p)
    expect_equal(unique(gp$sel_func), "sigmoid_length")
    expect_true(all(c("l50", "l25") %in% names(gp)))
    expect_false(any(c("l50_right", "l25_right") %in% names(gp)))
})

test_that("addEcopathCatchTotal sets gear params with double_sigmoid_length", {
    sp <- data.frame(
        species          = "TestFish",
        w_max            = 1e4,
        w_mat            = 100,
        a                = 0.01,
        b                = 3,
        age_mat          = 2,
        Length           = 50,
        biomass_observed = 500,
        stringsAsFactors = FALSE
    )

    ecop <- data.frame(
        `...1` = 1,
        `Group name` = "GroupA",
        `Biomass (t/km²)` = 500,
        `Consumption / biomass (/year)` = 5,
        `Production / consumption (/year)` = 1,
        stringsAsFactors = FALSE
    )

    sp <- addEcopathParams(sp, ecop, list(TestFish = "GroupA"))
    params <- newAllometricParams(sp)

    catch_df <- data.frame(
        `Group name`              = "GroupA",
        `TotalCatch (t/km²/year)` = 100,
        stringsAsFactors          = FALSE
    )

    p2 <- addEcopathCatchTotal(params, catch_df, sel_func = "double_sigmoid_length")
    gp2 <- gear_params(p2)

    expect_equal(unique(gp2$sel_func), "double_sigmoid_length")
    expect_true(all(c("l50_right", "l25_right") %in% names(gp2)))
    expect_equal(gp2$l50_right, 50)
    expect_equal(gp2$l25_right, 50 * 1.1)
})
