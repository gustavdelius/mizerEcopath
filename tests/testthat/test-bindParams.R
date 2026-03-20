# bindParams ----
test_that("bindParams combines species correctly", {
    p <- newTraitParams(no_sp = 4)
    sp <- p@species_params$species
    p1 <- removeSpecies(p, sp[3:4])
    p2 <- removeSpecies(p, sp[1:2])
    pc <- bindParams(p1, p2)

    expect_identical(pc@species_params$species, sp)
    expect_equal(dim(pc@psi), c(4L, length(p@w)))
    expect_equal(pc@psi[1:2, ], p1@psi, ignore_attr = TRUE)
    expect_equal(pc@psi[3:4, ], p2@psi, ignore_attr = TRUE)
    expect_equal(pc@maturity[1:2, ], p1@maturity, ignore_attr = TRUE)
    expect_equal(pc@maturity[3:4, ], p2@maturity, ignore_attr = TRUE)
    expect_equal(pc@initial_n[1:2, ], p1@initial_n, ignore_attr = TRUE)
    expect_equal(pc@initial_n[3:4, ], p2@initial_n, ignore_attr = TRUE)
    expect_equal(pc@intake_max[1:2, ], p1@intake_max, ignore_attr = TRUE)
    expect_equal(pc@search_vol[1:2, ], p1@search_vol, ignore_attr = TRUE)
    expect_equal(pc@metab[1:2, ], p1@metab, ignore_attr = TRUE)
    expect_equal(pc@mu_b[1:2, ], p1@mu_b, ignore_attr = TRUE)
    expect_equal(pc@ext_encounter[1:2, ], p1@ext_encounter, ignore_attr = TRUE)
    expect_equal(pc@diffusion[1:2, ], p1@diffusion, ignore_attr = TRUE)
    expect_equal(pc@ft_pred_kernel_e[1:2, ], p1@ft_pred_kernel_e, ignore_attr = TRUE)
    expect_equal(pc@ft_mask[1:2, ], p1@ft_mask, ignore_attr = TRUE)
    expect_identical(names(dimnames(pc@mu_b)), names(dimnames(p1@mu_b)))
    expect_identical(names(dimnames(pc@ft_pred_kernel_e)),
                     names(dimnames(p1@ft_pred_kernel_e)))
    expect_identical(pc@w_min_idx[1:2], p1@w_min_idx)
    expect_identical(pc@A[1:2], p1@A)
    validObject(pc)
})

test_that("bindParams sets interaction matrix to 0", {
    p <- newTraitParams(no_sp = 4)
    sp <- p@species_params$species
    p1 <- removeSpecies(p, sp[3:4])
    p2 <- removeSpecies(p, sp[1:2])
    pc <- bindParams(p1, p2)

    expect_true(all(pc@interaction == 0))
    expect_identical(dimnames(pc@interaction)$predator, sp)
    expect_identical(dimnames(pc@interaction)$prey, sp)
})

test_that("bindParams combines data frames correctly", {
    p <- newTraitParams(no_sp = 4)
    sp <- p@species_params$species
    p1 <- removeSpecies(p, sp[3:4])
    p2 <- removeSpecies(p, sp[1:2])
    # add a column only in p2 to test NA-filling
    p2@species_params$extra <- "x"
    pc <- bindParams(p1, p2)

    expect_identical(nrow(pc@species_params), 4L)
    expect_identical(pc@species_params$species, sp)
    expect_true("extra" %in% names(pc@species_params))
    expect_true(all(is.na(pc@species_params$extra[1:2])))
    expect_identical(pc@species_params$extra[3:4], c("x", "x"))
})

test_that("bindParams accepts a list of params", {
    p <- newTraitParams(no_sp = 4)
    sp <- p@species_params$species
    p1 <- removeSpecies(p, sp[3:4])
    p2 <- removeSpecies(p, sp[1:2])
    pc_dots <- bindParams(p1, p2)
    pc_list <- bindParams(list(p1, p2))
    expect_identical(pc_dots@species_params$species,
                     pc_list@species_params$species)
    expect_identical(pc_dots@psi, pc_list@psi)
})

test_that("bindParams combines gears correctly", {
    p <- newTraitParams(no_sp = 4)
    sp <- p@species_params$species
    p1 <- removeSpecies(p, sp[3:4])
    p2 <- removeSpecies(p, sp[1:2])
    p2 <- renameGear(p2, setNames("gear2", dimnames(p2@selectivity)[[1]]))
    pc <- bindParams(p1, p2)

    all_gears <- c(dimnames(p1@selectivity)[[1]], "gear2")
    expect_identical(dimnames(pc@selectivity)[[1]], all_gears)
    expect_identical(dimnames(pc@catchability)[[1]], all_gears)
    expect_true(all(names(pc@initial_effort) %in% all_gears))
    # selectivity for gear2 is 0 for species 1:2
    expect_true(all(pc@selectivity["gear2", 1:2, ] == 0))
    # original gear selectivity is 0 for species 3:4
    orig_gear <- dimnames(p1@selectivity)[[1]]
    expect_true(all(pc@selectivity[orig_gear, 3:4, ] == 0))
})

test_that("bindParams combines pred_kernel correctly", {
    p <- newTraitParams(no_sp = 4)
    p <- setPredKernel(p, pred_kernel = getPredKernel(p))
    sp <- p@species_params$species
    p1 <- removeSpecies(p, sp[3:4])
    p2 <- removeSpecies(p, sp[1:2])
    pc <- bindParams(p1, p2)

    expect_equal(dim(pc@pred_kernel)[1], 4L)
    expect_equal(pc@pred_kernel[1:2, , ], p1@pred_kernel, ignore_attr = TRUE)
    expect_equal(pc@pred_kernel[3:4, , ], p2@pred_kernel, ignore_attr = TRUE)
})

test_that("bindParams errors on duplicate species", {
    p1 <- newSingleSpeciesParams()
    p2 <- newSingleSpeciesParams()
    expect_error(bindParams(p1, p2), "Duplicate species found")
})

test_that("bindParams errors on mismatched w", {
    p1 <- newSingleSpeciesParams(no_w = 50)
    p2 <- newSingleSpeciesParams(species_name = "sp2", no_w = 80)
    expect_error(bindParams(p1, p2), "identical `w` slots")
})

test_that("bindParams errors on mismatched resource_params", {
    p1 <- newSingleSpeciesParams()
    p2 <- newSingleSpeciesParams(species_name = "sp2", kappa = 1e5)
    # force same w grid
    p2@w <- p1@w
    p2@w_full <- p1@w_full
    expect_error(bindParams(p1, p2), "identical `resource_params`")
})

test_that("bindParams errors with fewer than two params", {
    p1 <- newSingleSpeciesParams()
    expect_error(bindParams(p1), "at least two")
    expect_error(bindParams(list(p1)), "at least two")
})
