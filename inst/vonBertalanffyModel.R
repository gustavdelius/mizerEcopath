# Script to fit von Bertalanffy growth parameters and match catch data
library(dplyr)
library(mizer)
library(mizerEcopath)

Jess_sp <- celtic_params@species_params |>
    select(species, w_mat, w_max,
           #biomass_cutoff, #biomass_observed,
           pred_kernel_type, beta, sigma,
           kernel_exp, kernel_l_l, kernel_u_l, kernel_l_r, kernel_u_r)
James_sp <- readRDS("inst/James_sp.rds")
sp <- James_sp |>
    select(species, a, b, w_inf, k_vb, t0, w_mat, biomass_observed,
           production_observed, consumption_observed) |>
    inner_join(Jess_sp, by = "species")

saveRDS(sp, "inst/sp.rds")

# Load pre-packaged species parameters for Celtic Sea model
sp <- readRDS("inst/sp.rds")

# Initialize von Bertalanffy params using mizerEcopath helper
p <- newVonBertalanffyParams(sp)

# Attach gear parameters matching the selected species
gp <- gear_params(celtic_params) |>
    filter(species %in% sp$species)

#gp$l50_right <- gp$l50 + 5
#gp$l25_right <- gp$l50 + 10
#gp$sel_func[gp$gear == "Gillnet"] <- "double_sigmoid_length"

gp$yield_observed <- gp$yield_observed / 3
gp$yield_observed[gp$species == "Herring"] <- 0
gear_params(p) <- gp
initial_effort(p) <- 1

# Bring single-species models to steady state and set Beverton-Holt stock-recruitment
p <- steadySingleSpecies(p) |>
    setBevertonHolt()
p <- matchBiomasses(p)

# Fit selectivity/catchability to match observed Celtic Sea catches
p <- matchCatch(p, catch = celtic_catch, production_lambda = 0,
                yield_lambda = 1)

p <- tuneEcopath(p, catch = celtic_catch, diet = reduced_dm)

# Plot resulting yield comparison
plotYieldVsSpecies(p)
plotBiomassVsSpecies(p)
plotlySpectra(p, power = 2, resource = FALSE)
plotlyFMort(p)
p@species_params$z_ext
gps <- gear_params(p)
View(gps[gps$gear == "Gillnet", ])
View(gps[gps$species == "Herring", ])

pb <- p
p <- pb
p <- matchDiet(p, reduced_dm)
warnings()

plotSpectra(p, power = 2)
p@species_params$D_ext
plot(getMort(p), log_y = TRUE)
p@gear_params$catchability
