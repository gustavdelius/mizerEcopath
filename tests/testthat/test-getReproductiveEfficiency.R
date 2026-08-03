params <- celtic_params

test_that("getReproductiveEfficiency works for MizerParams", {
    efficiency <- getReproductiveEfficiency(params)

    expect_length(efficiency, nrow(params@species_params))
    expect_named(efficiency, params@species_params$species)
    # It is the product of erepro and the survival of density dependence
    expect_equal(efficiency,
                 unname(params@species_params$erepro *
                            (1 - getReproductionLevel(params))),
                 ignore_attr = "names")
})

test_that("getReproductiveEfficiency is the offspring biomass per investment", {
    # The efficiency is the offspring biomass production expressed as a
    # fraction of the energy the females invest into reproduction. The identity
    # holds whichever quadrature the model uses, because getGonadicProduction()
    # discretises the reproduction integral the same way getRDI() does.
    expect_false(isTRUE(params@second_order_w[["bin_average"]]))
    expect_equal(unname(getReproductiveEfficiency(params)),
                 unname(getOffspringProduction(params) /
                            (getGonadicProduction(params) / 2)))

    p <- params
    second_order_w(p) <- TRUE
    expect_equal(unname(getReproductiveEfficiency(p)),
                 unname(getOffspringProduction(p) /
                            (getGonadicProduction(p) / 2)))
})

test_that("getReproductiveEfficiency returns NA for species not reproducing", {
    p <- params
    p@species_params$erepro[[1]] <- 0
    expect_equal(getRDI(p)[[1]], 0)

    efficiency <- getReproductiveEfficiency(p)
    expect_true(is.na(efficiency[[1]]))
    expect_equal(efficiency[-1], getReproductiveEfficiency(params)[-1])
})

test_that("getReproductiveEfficiency works for MizerSim", {
    sim <- project(params, t_max = 2, t_save = 1, progress_bar = FALSE)
    efficiency <- getReproductiveEfficiency(sim)

    expect_true(is.ArrayTimeBySpecies(efficiency))
    expect_equal(dim(efficiency), c(3, nrow(params@species_params)))
    expect_equal(dimnames(efficiency)$sp, params@species_params$species)
    # The first saved time is the initial state, so it must agree with the
    # value calculated from the MizerParams object
    expect_equal(as.numeric(efficiency[1, ]),
                 unname(getReproductiveEfficiency(params)))
})
