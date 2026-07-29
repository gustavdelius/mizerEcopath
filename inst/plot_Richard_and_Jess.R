# Plot the steady state in Richard and Jess model
params <- psr
plotYieldVsSpecies(params)
plotBiomassVsSpecies(params)
plotlySpectra(params, power = 2, resource = FALSE)
plotlyFMort(params)
plotGrowthCurves(params, species_panel = TRUE)
plotlyFeedingLevel(params, include_critical = TRUE)
plotlyDiet(params, species = "Megrim")
