sp_single <- data.frame(
    species = "TestFish",
    w_max = 1000,
    w_mat = 100,
    n = 0.7,
    a = 0.01,
    b = 3,
    age_mat = 5,
    Length = 50
)

# newAllometricParams --------------------------------------------------------

test_that("newAllometricParams returns a MizerParams object", {
    p <- newAllometricParams(sp_single)
    expect_s4_class(p, "MizerParams")
})

test_that("newAllometricParams sets default exponents correctly", {
    p <- newAllometricParams(sp_single)
    sp <- p@species_params
    expect_equal(sp$n, 0.7)
    expect_equal(sp$p, 0.7)
    expect_equal(sp$d, sp$n - 1)
})

test_that("newAllometricParams fills missing d with n - 1", {
    sp <- sp_single[, setdiff(names(sp_single), "d")]
    p <- newAllometricParams(sp)
    expect_equal(p@species_params$d, p@species_params$n - 1)
})

test_that("newAllometricParams sets default alpha to 0.8", {
    p <- newAllometricParams(sp_single)
    expect_equal(p@species_params$alpha, 0.8)
})

test_that("newAllometricParams respects user-supplied alpha", {
    sp <- sp_single
    sp$alpha <- 0.6
    p <- newAllometricParams(sp)
    expect_equal(p@species_params$alpha, 0.6)
})

test_that("newAllometricParams sets ks and z0 to zero", {
    p <- newAllometricParams(sp_single)
    expect_equal(p@species_params$ks, 0)
    expect_equal(p@species_params$z0, 0)
})

test_that("newAllometricParams switches off all interactions", {
    p <- newAllometricParams(sp_single)
    expect_true(all(interaction_matrix(p) == 0))
    expect_equal(p@species_params$interaction_resource, 0)
})

test_that("newAllometricParams switches off satiation (h = Inf)", {
    p <- newAllometricParams(sp_single)
    expect_equal(p@species_params$h, Inf)
    expect_true(all(is.infinite(intake_max(p))))
})

test_that("newAllometricParams uses resource_constant dynamics", {
    p <- newAllometricParams(sp_single)
    expect_equal(p@resource_dynamics, "resource_constant")
})

test_that("newAllometricParams produces allometric rates", {
    p <- newAllometricParams(sp_single)
    expect_true(isAllometric(p))
})

test_that("newAllometricParams works with multi-species input", {
    sp <- species_params(NS_params)
    p <- newAllometricParams(sp, no_w = 100)
    expect_s4_class(p, "MizerParams")
    expect_equal(nrow(p@species_params), nrow(sp))
    expect_true(isAllometric(p))
})

test_that("newAllometricParams juvenile biomass spectrum slope is near -0.2", {
    p <- newAllometricParams(sp_single)
    w <- p@w
    N <- as.numeric(initialN(p))
    B <- N * w
    juvenile <- w < p@species_params$w_mat
    fit <- lm(log(B[juvenile]) ~ log(w[juvenile]))
    slope <- coef(fit)["log(w[juvenile])"]
    expect_equal(as.numeric(slope), -0.2, tolerance = 0.1)
})

test_that("juvenile biomass spectrum slope is near -0.2 for default parameters", {
    # Set up a simple species parameter data frame for a single species
    sp <- data.frame(
        species = "TestFish",
        w_max = 1000,        # maximum weight
        w_mat = 100,         # weight at maturity
        n = 0.7,             # default consumption exponent
        a = 0.01,            # dummy parameter for length-weight relationship
        b = 3,               # dummy parameter for length-weight relationship
        age_mat = 5,
        Length = 50
    )

    # Create the model parameters with newAllometricParams
    p <- newAllometricParams(sp, no_w = 200)

    # Get the weight grid from MizerParams
    w <- p@w

    # Extract the abundance distribution
    N <- as.numeric(initialN(p))

    # Calculate biomass density as abundance multiplied by weight
    B <- N * w

    # Define the juvenile range (weights less than w_mat)
    juvenile <- w < p@species_params$w_mat[p@species_params$species == "TestFish"]

    # Log-transform the values and fit a linear model
    log_w <- log(w[juvenile])
    log_B <- log(B[juvenile])
    fit <- lm(log_B ~ log_w)
    slope <- coef(fit)["log_w"]

    # Check that the fitted slope is near -0.2 (with a tolerance of 0.1)
    expect_equal(as.numeric(slope), -0.2, tolerance = 0.1)
})


# isAllometric ---------------------------------------------------------------

test_that("isAllometric returns TRUE for newAllometricParams output", {
    p <- newAllometricParams(sp_single)
    expect_true(isAllometric(p))
})

test_that("isAllometric returns FALSE for non-allometric model", {
    expect_false(isAllometric(NS_params))
})

test_that("isAllometric returns FALSE when encounter rate is perturbed", {
    p <- newAllometricParams(sp_single)
    ee <- ext_encounter(p)
    ee[1, ] <- ee[1, ] * seq(0.5, 1.5, length.out = ncol(ee))
    ext_encounter(p) <- ee
    expect_false(isAllometric(p))
})

test_that("isAllometric errors on zero mortality", {
    p <- newAllometricParams(sp_single)
    ext_mort(p)[] <- 0
    expect_snapshot(isAllometric(p), error = TRUE)
})
