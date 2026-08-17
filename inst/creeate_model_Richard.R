
library(dplyr)
library(ggplot2)
library(mizer)
library(mizerEcopath)

load("inst/data_Richard_survey.rda")

# Build MizerParams ----
#p <- newVonBertalanffyParams(sp)
p <- newAllometricParams(sp)
gear_params(p) <- gp
initial_effort(p) <- 1

p <- steadySingleSpecies(p) |>
    setBevertonHolt() |>
    matchBiomasses()

# Anchor the mortality. The fit is free to move the external mortality
# coefficient `z_ext`, and with only size distributions and a yield to constrain
# it, it can wander. In steady state mortality and production are two sides of
# the same coin, so asking the fit to stay near the default model's production
# keeps the mortality near the default too. This is a regularisation term, not
# an observation.
species_params(p)$production_observed <- getSomaticProduction(p)

# Set fishing
source("inst/catchability_Richard.R")
plot(ArraySpeciesBySize(getFMortGear(p)["commercial", , ], params = p))

# Match ----
pm <- tuneEcopath(p, catch = survey_Richard, diet = reduced_dm, match = c())

# Feeding levels ----

# Raising the feeding level damps how strongly growth responds to a change in
# food, which weakens the compensation that holds a species' yield peak above
# the fishing mortality we want it to peak at.
#
# These four were chosen by measuring where every species' peak lands relative
# to its target in the *finished* model, across six configurations. Ratio of
# peak to target, wide scan so that nothing is censored at the grid edge:
#
#                     all 0.6   Her .9   these four   +cannib   Her.9+cannib0
#   Herring            >3.00     1.58       1.19        1.49        >3.00
#   Cod                 1.51     1.52       1.51        1.40         1.27
#   Megrim              1.24     1.38       1.33        1.48         2.48
#   Haddock             1.02     1.14       1.17        1.59         2.53
#   Whiting             1.64     2.01       1.69        1.07         1.16
#   Hake                1.10     1.09       1.07        1.02         1.00
#   Blue whiting        1.04     1.04       1.04        1.06         1.02
#   Plaice              1.05     1.05       1.04        1.04         0.99
#   Sole                1.07     1.07       1.07        1.07         1.07
#   mean |ratio - 1|    0.41     0.32      *0.234*      0.246        0.62
#
# Two things this table settles. Raising Cod and Whiting does nothing for Cod
# and Whiting -- they are unmoved by any feeding level, because their peaks are
# held up by cannibalism rather than by food (see below). They are raised anyway
# because the *community* result is best with them raised: Herring reaches 1.19
# only when they are, and 1.58 when it is raised alone. The knob is not a
# per-species dial.
#
# And weakening the cannibalism, which does reach Cod and Whiting, buys nothing
# overall: it improves them and costs Herring, Megrim and Haddock, netting out
# no better than the feeding levels alone (0.246 against 0.234). Since it also
# reallocates a diet entry the Ecopath data does report into anonymous external
# mortality, it is not worth paying for a tie. See `rescaleInteraction()`, and
# the Cannibalism section of `create_model_from_Richard_and_Jess_new.R`, if you
# want to revisit that.
#
# NB `setFeedingLevel()` assigns positionally; the names are for our benefit.
f <- setNames(rep(0.6, nrow(sp)), sp$species)
f[["Herring"]] <- 0.9
f[["Cod"]]     <- 0.9
f[["Megrim"]]  <- 0.9
f[["Whiting"]] <- 0.9

Er_before <- getEReproAndGrowth(pm)
pm <- setFeedingLevel(pm, feeding_level = unname(f))
# Growth must be untouched: the size-distribution fit above depends on it.
stopifnot(max(abs(getEReproAndGrowth(pm) - Er_before) /
              pmax(abs(Er_before), 1e-30)) < 1e-8)

# Diet ----

# Unlike the allometric model, no production raises are needed: every species is
# already feasible. Check rather than assume, since this depends on the feeding
# levels set above.
feas <- getDietFeasibility(pm, reduced_dm)
stopifnot(all(feas$mort_feasible))
# The remaining `enc_feasible = FALSE` entries (Megrim, Whiting, Hake) are the
# milder problem: the diet matrix gives them more food from modelled species
# than the model gives them in total, and the shortfall is clipped at zero. It
# affects a very small number of large individuals, so `steady()` absorbs it.

pd <- matchDiet(pm, reduced_dm)
ps <- steady(pd)

# Resource ----

psr <- addSemichemostatResource(ps, resource_level = 0.5, w_pp_cutoff = 1)

# No cannibalism reallocation here -- see the feeding-level section above for
# why, and note that Cod is eaten only by Cod in this model and Whiting mostly
# by itself, so the mechanism is present; it is the cost/benefit that fails.

# Reproduction ----

# Recorded from `tuneReproductionLevel(psr, sweeps = 2, criterion = "slope")`,
# which takes around ten minutes. Re-derive it rather than trusting these if
# anything above changes.
params <- setBevertonHolt(psr, reproduction_level = c(
    Herring        = 0.1111, Cod     = 0.0100, Megrim = 0.9253,
    Haddock        = 0.8674, Whiting = 0.0100, Hake   = 0.5817,
    `Blue whiting` = 0.7281, Plaice  = 0.7281, Sole   = 0.7645))

# Where the peaks end up. Use the wide scan: the default `mult` stops at 1.5, and
# several species here sit beyond it, so the default silently reports 1.5 with
# `bracketed = FALSE` rather than a measurement.
# getYieldPeaks(params, mult = seq(0.25, 3, 0.25))
#
# The yield-peak criterion is only partly met, and by a wider margin than in the
# allometric model, which lands every species between 0.97 and 1.09. Here four
# species sit between 1.17 and 1.69. That is worth stating plainly rather than
# tuning away: the von Bertalanffy model pins growth to the observed growth
# curve, leaving the fit less freedom to arrange a spectrum that also responds
# to fishing the way single-species reference points imply. Note too that the
# targets themselves are extrapolations for the worst offenders -- Herring's
# FMSY is 70 times its fitted F and Blue whiting's 26 times -- so a large miss
# there is as much a comment on the criterion as on the model.

# Save ----
saveParams(params, "inst/params_final_Richard_and_Jess_new_vB.rds")
