p_allometric <- newAllometricParams(
    data.frame(
        species = "TestFish",
        w_max = 1000,
        w_mat = 100,
        n = 0.7,
        a = 0.01,
        b = 3,
        age_mat = 5,
        Length = 50
    )
)

test_that("alignResource returns a MizerParams object", {
    p <- alignResource(p_allometric)
    expect_s4_class(p, "MizerParams")
})

test_that("alignResource resource follows power law with correct lambda below cutoff", {
    p <- alignResource(p_allometric)
    lambda <- p@resource_params$lambda
    cutoff <- p@resource_params$w_pp_cutoff
    w_full <- p@w_full
    N_R <- initialNResource(p)
    below <- w_full <= cutoff & N_R > 0
    # N_R[w] / w^(-lambda) should be constant (= kappa) throughout
    coeff <- N_R[below] / w_full[below]^(-lambda)
    expect_equal(unname(coeff), rep(unname(coeff[1]), sum(below)), tolerance = 1e-10)
})

test_that("alignResource resource is tangent to fish community abundance", {
    p <- alignResource(p_allometric)
    N_R <- initialNResource(p)
    total <- colSums(initialN(p))
    fish_sel <- p@w_full >= p@w[1]
    # The minimum log-difference should be exactly zero (tangent point)
    log_diff <- log(N_R[fish_sel]) - log(total)
    expect_equal(min(log_diff), 0, tolerance = 1e-10)
})

test_that("alignResource resource is zero above w_pp_cutoff", {
    p <- alignResource(p_allometric)
    cutoff <- p@resource_params$w_pp_cutoff
    N_R <- initialNResource(p)
    above <- p@w_full > cutoff
    expect_true(all(N_R[above] == 0))
})

test_that("alignResource updates kappa in resource_params", {
    p <- alignResource(p_allometric)
    lambda <- p@resource_params$lambda
    cutoff <- p@resource_params$w_pp_cutoff
    w_full <- p@w_full
    N_R <- initialNResource(p)
    below <- w_full <= cutoff & N_R > 0
    expected_kappa <- unname(N_R[below][1] / w_full[below][1]^(-lambda))
    expect_equal(p@resource_params$kappa, expected_kappa, tolerance = 1e-10)
})

test_that("alignResource sets initialNResource equal to resource_capacity", {
    p <- alignResource(p_allometric)
    expect_equal(initialNResource(p), resource_capacity(p))
})

test_that("alignResource respects custom w_pp_cutoff argument", {
    w_cutoff <- p_allometric@w[50]
    p <- alignResource(p_allometric, w_pp_cutoff = w_cutoff)
    expect_equal(p@resource_params$w_pp_cutoff, w_cutoff)
    N_R <- initialNResource(p)
    expect_true(all(N_R[p@w_full > w_cutoff] == 0))
})
