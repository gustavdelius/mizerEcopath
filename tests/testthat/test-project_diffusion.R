test_that("project_diffusion handles emigration correctly", {
    params <- celtic_params
    species <- "Cod"

    # Set initial condition
    n_init <- rep(0, length(w))
    n_init[2] <- 1  # Point source at small size
    initialN(params)[species, ] <- n_init
    # Solve the PDE
    dt <- 0.05
    nsteps <- 1

    n_hist <- project_diffusion(params, species, dt = dt, nsteps = nsteps)

    # Now add emigration
    params2 <- params
    emigration(params2)[] <- 0.01
    n_hist2 <- project_diffusion(params2, species, dt = dt, nsteps = nsteps)
})
