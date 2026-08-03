params <- celtic_params

test_that("plotPreyAvailability returns the correct data", {
    df <- plotPreyAvailability(params, species = "Cod", pred_size = 1000,
                               return_data = TRUE)

    expect_named(df, c("w", "value", "Type"))
    expect_equal(nrow(df), 4 * length(params@w_full))
    expect_setequal(df$Type,
                    c("Feeding kernel", "Number density",
                      "Number density in log w", "Biomass density in log w"))
    expect_true(all(is.finite(df$value)))
    # Each curve is a probability density in log w
    dx <- log(params@w_full[2]) - log(params@w_full[1])
    integrals <- tapply(df$value, df$Type, function(value) sum(value) * dx)
    expect_equal(as.vector(integrals), rep(1, 4))
})

test_that("plotPreyAvailability takes the interaction with the prey into account", {
    p <- params
    # Cod does not interact with anything, so there is no prey at all
    p@interaction["Cod", ] <- 0
    species_params(p)$interaction_resource[[3]] <- 0
    df <- plotPreyAvailability(p, species = "Cod", pred_size = 1000,
                               return_data = TRUE)
    kernel <- df[df$Type == "Feeding kernel", ]
    expect_true(all(is.finite(kernel$value)))
    # while switching on the resource gives prey at small sizes
    species_params(p)$interaction_resource[[3]] <- 1
    df <- plotPreyAvailability(p, species = "Cod", pred_size = 1000,
                               return_data = TRUE)
    numbers <- df[df$Type == "Number density", ]
    expect_true(all(is.finite(numbers$value)))
    expect_true(sum(numbers$value) > 0)
})

test_that("plotPreyAvailability defaults to the maturity size", {
    df <- plotPreyAvailability(params, species = "Cod", return_data = TRUE)
    expected <- plotPreyAvailability(params, species = "Cod",
                                     pred_size = params@species_params["Cod", "w_mat"],
                                     return_data = TRUE)
    expect_identical(df, expected)
})

test_that("plotPreyAvailability checks its arguments", {
    expect_error(suppressWarnings(plotPreyAvailability(params, species = "Yeti")))
    expect_error(plotPreyAvailability(params, species = c("Cod", "Hake")),
                 "one species at a time")
})

test_that("plotPreyAvailability returns a plot", {
    expect_s3_class(plotPreyAvailability(params, species = "Cod"), "ggplot")
    expect_s3_class(plotlyPreyAvailability(params, species = "Cod"), "plotly")
})
