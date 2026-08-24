# Build a Celtic Sea model from Jess's data: Ecopath biomasses, production and
# consumption together with her gear-resolved commercial catch size
# distributions and her survey size distributions.
#
# The accompanying narrative is in vignettes/Jess.qmd.
#
# Stages 1-6 take a couple of minutes. Stage 8 takes about half an hour and its
# result is recorded at the bottom of the file.

library(dplyr)
library(ggplot2)
library(mizer)
library(mizerExperimental)
library(mizerEcopath)

source("inst/calibration_helpers.R")   # set_rl(), cap_rl(), setFeedingLevelInteracting()
source("inst/multigear_helpers.R")     # the multi-gear versions of the rest

load("inst/data_Jess.rda")

# Gear setup ----

# Every gear, survey and commercial alike, is *allowed* a dome-shaped
# selectivity; the fit decides which ones need one. The survey trawl
# under-samples the largest fish, and in this data set the commercial fleets for
# cod, mackerel, horse mackerel and sole do too. Forcing the selectivity to be
# monotone makes the size spectrum itself carry the fall-off, which it can only
# do by inventing mortality. Where a gear's data show no dome, the fit pushes
# l50_right out to thousands of centimetres and the curve is a sigmoid again.
gp$l50_right <- NA_real_
gp$l25_right <- NA_real_
gp$sel_func <- "double_sigmoid_length"
l_max_gear <- sp$l_max[match(gp$species, sp$species)]
gp$l50_right <- pmax(0.6 * l_max_gear, gp$l50 * 1.05)
gp$l25_right <- pmax(0.9 * l_max_gear, gp$l50_right * 1.2)

# Weighting of the size-distribution likelihood. Each gear contributes a
# multinomial negative log likelihood *per fish*, so the number of fish behind a
# distribution does not set its weight - the number of gears does. The survey is
# a single gear while the commercial fishery is split over up to four, so
# without this the commercial data would outvote the survey four to one for some
# species and one to one for others. Giving the commercial gears of a species
# 1/n each puts the two data sources on an equal footing for every species.
n_com <- gp |> filter(gear != "survey") |> count(species, name = "n_com")
gp <- gp |>
    left_join(n_com, by = "species") |>
    mutate(catch_dist_weight = if_else(gear == "survey", 1, 1 / n_com)) |>
    select(-n_com)

# Stage 1: single-species scaffolding ----

p <- newVonBertalanffyParams(sp)
gear_params(p) <- gp
initial_effort(p) <- 1
p <- steadySingleSpecies(p) |> setBevertonHolt()
p <- matchBiomasses(p)

# Stage 2: size distributions, yields and production ----

# `production_observed` comes with the data and is kept: it is what ties the
# external mortality to something observed rather than to the model's own
# default. `z_ext_lim = 20` lifts the default cap of 5, which monkfish needs.
pm <- matchCatch(p, catch = catch, z_ext_lim = 20)

# Blue whiting's production can only be reached with a juvenile spectrum
# steeper than the default community slope allows, so its cap is relaxed.
pm <- matchCatch(pm, catch = catch, species = "Blue whiting",
                 z_ext_lim = 20, lambda = 2.5)

# Stage 3: satiation ----

# Consumption, growth and the whole steady state are unchanged; what changes is
# how strongly growth responds to a change in food.
pm <- setFeedingLevel(pm, feeding_level = 0.6)

# Stage 4: making the diet matrix affordable ----

# matchDiet() pays for explicit predation out of the external mortality, and the
# diet matrix demands more predation on herring, blue whiting, mackerel and
# especially horse mackerel than they carry at the sizes where the predation
# falls. `production_observed` is data we have been asked to keep, so we cannot
# buy them more mortality; instead we move the excess predation into the `other`
# column, which leaves every predator's total consumption untouched and only
# changes who it is attributed to.
fit <- fit_diet_matrix(pm, reduced_dm, kappa = 0.9)
dm <- fit$dm
print(round(fit$lambda, 4))

# Stage 5: species interactions ----

pd <- matchDiet(pm, dm)
ps <- steady(pd, tol = 1e-10)

# Stage 6: the plankton resource ----

psr <- alignResource(ps)
resource_params(psr)$w_pp_cutoff <- 1
initialNResource(psr)[w_full(psr) > 1] <- 0
comment(psr@cc_pp) <- NULL
psr <- setResourceInteraction(psr, resource_dynamics = "resource_semichemostat",
                              tol = 1e-2)
psr <- steady(psr, tol = 1e-10)
resource_level(psr) <- 0.5
psr <- steady(psr, tol = 1e-12, t_max = 200)

# Stage 7: response knobs ----

# Herring eats almost nothing but the resource, and at f = 0.6 its growth
# response to being fished down holds its yield peak 14% above its FMSY even at
# the reproduction floor, so the reproduction level cannot reach it. Raising the
# feeding level damps that response without moving the steady state: at f = 0.85
# the peak at the floor is 0.78 of the target, which leaves the bisection room.
psr <- setFeedingLevelInteracting(psr, c(Herring = 0.85))

# Stage 8: reproduction levels ----

# Targets: FMSY where the data give one, the model's own current Fbar
# otherwise. Fbar is the biomass-weighted mean fishing mortality over mature
# individuals, the multi-gear analogue of the fully-selected F.
sps <- species_params(psr)$species
fb <- vapply(sps, function(s) Fbar(psr, s), numeric(1))
Ft <- setNames(species_params(psr)$FMSY, sps)
Ft[is.na(Ft)] <- fb[is.na(Ft)]
m_t <- Ft / fb

if (FALSE) {   # about twenty minutes per sweep
    params <- setBevertonHolt(psr, reproduction_level = 0.5)
    for (sweep in 1:4) {
        res <- parallel::mclapply(sps, function(s) {
            r <- tune_rl(params, s, m_t[[s]], steps = 6, tol = 0.10)
            r$species <- s
            r
        }, mc.cores = 6)
        for (r in res) params <- set_rl(params, r$species, r$rl)
        print(do.call(rbind, lapply(res, as.data.frame)))
    }
    # Verify the peaks on a wider grid than the sweeps use: the sweeps' grid
    # stops at 2.6 times the target and truncates three species' peaks.
    do.call(rbind, parallel::mclapply(sps, function(s) {
        y <- yield_vs_mult(params, s, mult_grid(m_t[[s]], 0.08, 6, 14))
        p <- peak_mult(y, m_t[[s]])
        data.frame(species = s, F_peak = p[["ratio"]] * Ft[[s]],
                   F_target = Ft[[s]], ratio = p[["ratio"]], edge = p[["edge"]])
    }, mc.cores = 6))
}

# The recorded result of those sweeps. Eight of the twelve peaks land within
# 11% of target. Haddock is pinned at the reproduction ceiling with its peak at
# half its target; whiting is pinned at the floor with no maximum anywhere in
# the scanned range; herring and megrim do not converge across sweeps - herring
# because its stage-7 feeding level was chosen before the other species were
# tuned, megrim because its yield curve is bimodal. See the "The sweeps do not
# converge" section of vignettes/Jess.qmd.
params <- setBevertonHolt(psr, reproduction_level = c(
    Herring  = 0.005,  Cod            = 0.005,  Megrim = 0.3038,
    Monkfish = 0.1016, Haddock        = 0.9995, Whiting = 0.005,
    Hake     = 0.005,  `Blue whiting` = 0.3193, Plaice = 0.7882,
    Mackerel = 0.2152, Sole           = 0.2917, `Horse mackerel` = 0.3031))

saveParams(params, "inst/params_final_Jess.rds")
