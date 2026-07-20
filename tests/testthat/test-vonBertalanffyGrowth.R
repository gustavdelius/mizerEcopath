# test-vonBertalanffyGrowth.R
#
# Tests for vonBertalanffyGrowth(), which returns a species-by-size array of
# von Bertalanffy growth rates computed from the w_inf, k_vb, a and b species
# parameters.
#
#   - returns a species-by-size array with the right dims and dimnames
#   - agrees with the closed-form A w^n - B w
#   - agrees with the length-based dl/dt = k_vb (L_inf - l)
#   - is zero above w_inf (never negative)
#   - uses w_max when w_inf is missing and b = 3 when b is missing
#   - returns NA for a species with missing k_vb

sp_single <- data.frame(
    species = "TestFish",
    w_max = 1000,
    w_mat = 100,
    n = 0.7,
    a = 0.01,
    b = 3,
    age_mat = 5,
    Length = 50,
    k_vb = 0.3,
    w_inf = 1000
)

test_that("vonBertalanffyGrowth returns a species-by-size array", {
    p <- newAllometricParams(sp_single)
    g <- vonBertalanffyGrowth(p)
    expect_true(is.array(g))
    expect_equal(dim(g), c(nrow(p@species_params), length(p@w)))
    expect_equal(dimnames(g)$sp, p@species_params$species)
})

test_that("vonBertalanffyGrowth agrees with the closed form A w^n - B w", {
    p <- newAllometricParams(sp_single)
    g <- vonBertalanffyGrowth(p)
    sp <- p@species_params
    w <- p@w
    for (i in seq_len(nrow(sp))) {
        nn <- 1 - 1 / sp$b[i]
        A <- sp$b[i] * sp$k_vb[i] * sp$w_inf[i] ^ (1 / sp$b[i])
        B <- sp$b[i] * sp$k_vb[i]
        expected <- pmax(A * w ^ nn - B * w, 0)
        expect_equal(unname(g[i, ]), expected)
    }
})

test_that("vonBertalanffyGrowth agrees with length-based dl/dt", {
    p <- newAllometricParams(sp_single)
    g <- vonBertalanffyGrowth(p)
    sp <- p@species_params
    w <- p@w
    below <- w < sp$w_inf[1]
    l <- (w / sp$a[1]) ^ (1 / sp$b[1])
    Linf <- (sp$w_inf[1] / sp$a[1]) ^ (1 / sp$b[1])
    dwdl <- sp$a[1] * sp$b[1] * l ^ (sp$b[1] - 1)
    dwdt <- dwdl * sp$k_vb[1] * (Linf - l)
    expect_equal(g[1, below], dwdt[below], ignore_attr = TRUE)
})

test_that("vonBertalanffyGrowth is zero above w_inf and never negative", {
    p <- newAllometricParams(sp_single)
    g <- vonBertalanffyGrowth(p)
    sp <- p@species_params
    above <- p@w > sp$w_inf[1]
    expect_true(all(g[1, above] == 0))
    expect_false(any(g < 0))
})

test_that("vonBertalanffyGrowth uses w_max and b = 3 defaults", {
    sp <- sp_single
    sp$w_inf <- NULL
    sp$b <- NULL
    p <- newAllometricParams(sp)
    # remove b so the default is exercised
    p@species_params$b <- NULL
    g <- vonBertalanffyGrowth(p)
    w <- p@w
    w_inf <- p@species_params$w_max[1]
    nn <- 1 - 1 / 3
    A <- 3 * p@species_params$k_vb[1] * w_inf ^ (1 / 3)
    B <- 3 * p@species_params$k_vb[1]
    expect_equal(unname(g[1, ]), pmax(A * w ^ nn - B * w, 0))
})

test_that("vonBertalanffyGrowth returns NA when k_vb is missing", {
    sp <- sp_single
    sp$k_vb <- NULL
    p <- newAllometricParams(sp)
    g <- vonBertalanffyGrowth(p)
    expect_true(all(is.na(g[1, ])))
})
