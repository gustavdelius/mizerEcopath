test_that("setFeedingLevel works", {
    expect_error(setFeedingLevel(NS_params, 1.1),
                 "only works for models where all encounter is external.")
    params <- makeNoninteracting(NS_params)
    expect_error(setFeedingLevel(params, c(1, 2)),
                 "The length of feeding_level vector")
    expect_error(setFeedingLevel(params, -1),
                 "must be positive")
    expect_error(setFeedingLevel(params, 1),
                 "must be positive and strictly less than 1")
    expect_error(setFeedingLevel(params, NA),
                 "feeding_level is not a numeric or integer vector")

    g <- getEGrowth(params)
    comment(g) <- "set manually"
    p <- setFeedingLevel(params, 0.3)
    desired <- getFeedingLevel(params) # matrix with right dimension
    desired[] <- 0.3
    comment(desired) <- "set manually"
    expect_equal(getFeedingLevel(p), desired)
    expect_equal(getEGrowth(p), g)
    expect_equal(given_species_params(p)$f0, rep(0.3, nrow(p@species_params)))

    p <- setFeedingLevel(params, 0)
    desired[] <- 0
    expect_equal(getFeedingLevel(p), desired)
    expect_equal(getEGrowth(p), g)
})
