# Create a Celtic Sea model from Richard's biomasses, yields and catch size
# distributions together with Jess's survey size distributions, using
# `newVonBertalanffyParams()` as the starting point so that every species grows
# along its own von Bertalanffy curve while in steady state.
#
# This is the von Bertalanffy counterpart of
# `inst/create_model_from_Richard_and_Jess.R`, which starts from
# `newAllometricParams()`. The data and the sequence of stages are the same; the
# growth model is not, and several of the tuning knobs in the later stages turn
# out to push in the opposite direction. The accompanying vignette
# `vignettes/Richard_and_Jess_vB.qmd` explains each stage and reports the scans
# behind the numbers chosen here.
library(dplyr)
library(ggplot2)
library(mizer)
library(mizerExperimental)
library(mizerEcopath)

load("inst/data_Richard_and_Jess.rda")

# Local helpers -------------------------------------------------------------

# Change a species' feeding level in a fully interacting model without moving
# the steady state, the diet matrix, or the mortality it inflicts on its prey.
#
# Multiply the predator's search volume *and* its external encounter rate - the
# two additive parts of the encounter rate - by c = (1 - f_old)/(1 - f_new), and
# set its maximum intake so that the feeding level becomes f_new. Then
#   encounter'   = c E                        (so the feeding level is f_new)
#   consumption' = (1 - f') c E = (1 - f) E   (unchanged)
#   consumption of each prey, and the predation mortality it inflicts,
#                ~ (1 - f') c = (1 - f)       (unchanged)
# so nothing in the steady state moves. What moves is the elasticity of
# consumption to food, dlog(C)/dlog(E) = 1 - f. Unlike `setFeedingLevel()` this
# can be applied after the interactions are in place, which means the feeding
# level can be retuned without rebuilding the model.
setFeedingLevelInteracting <- function(params, f_new) {
    f_old <- getFeedingLevel(params)
    sp <- species_params(params)$species
    if (is.null(names(f_new))) names(f_new) <- sp
    for (s in names(f_new)) {
        i <- match(s, sp)
        cc <- (1 - f_old[i, ]) / (1 - f_new[[s]])
        sv <- search_vol(params); sv[i, ] <- sv[i, ] * cc
        search_vol(params) <- sv
        ee <- ext_encounter(params); ee[i, ] <- ee[i, ] * cc
        ext_encounter(params) <- ee
        E <- getEncounter(params)
        im <- intake_max(params)
        im[i, ] <- E[i, ] * (1 - f_new[[s]]) / f_new[[s]]
        intake_max(params) <- im
    }
    params
}

# Weaken the predation of `pred` on `prey` by the factor lambda and hand back
# what that removes: the lost food to the predators' external encounter, the lost
# predation mortality to the prey's external mortality. This is
# `makeInteracting()` run in reverse, and it leaves growth, mortality, the size
# spectrum and the catch size distributions exactly where they were; what changes
# is that the reallocated part of the food web can no longer respond to fishing.
# The allometric build uses the special case pred == prey (cannibalism).
#
# The ordering of the two restorations is not cosmetic: predation mortality
# carries a factor (1 - f) and so depends on the encounter rate, so measuring
# the mortality shortfall before the encounter has been restored would attribute
# a spurious external mortality to every prey species.
set_link <- function(params, pred, prey, lambda) {
    mort_old <- getPredMort(params)
    enc_old  <- getEncounter(params)
    inter <- interaction_matrix(params)
    inter[pred, prey] <- inter[pred, prey] * lambda
    interaction_matrix(params) <- inter
    ext_encounter(params) <- ext_encounter(params) +
        (enc_old - getEncounter(params))
    ext_mort(params) <- ext_mort(params) + (mort_old - getPredMort(params))
    params
}

# Stage 1: the von Bertalanffy starting model ---------------------------------

# The survey size distributions of Megrim, Blue whiting and Plaice fall away at
# large sizes more steeply than any monotone selectivity curve can produce, so
# for these three the survey gear is given a dome-shaped (double sigmoid)
# selectivity. This is a change to the *observation* model only: the survey gear
# has catchability 1e-12 and `yield_weight = 0`, so it removes no fish and its
# shape has no effect on the dynamics. The commercial selectivity, which does
# set the fishing mortality, is left monotone; letting it be domed too buys
# almost nothing and fits 2 cm wide knife-bands (see the vignette).
domed <- c("Megrim", "Blue whiting", "Plaice")
gp$l50_right <- NA_real_
gp$l25_right <- NA_real_
dome_sel <- gp$gear == "survey" & gp$species %in% domed
gp$sel_func[dome_sel] <- "double_sigmoid_length"
l_max_dome <- sp$l_max[match(gp$species[dome_sel], sp$species)]
gp$l50_right[dome_sel] <- 0.75 * l_max_dome
gp$l25_right[dome_sel] <- 0.95 * l_max_dome

p <- newVonBertalanffyParams(sp)
gear_params(p) <- gp
initial_effort(p) <- 1

p <- steadySingleSpecies(p) |>
    setBevertonHolt()
p <- matchBiomasses(p)

# We don't want the optimiser to wander too far from the default mortalities.
# We do this by matching to the production in the default model. In steady
# state, mortality is determined by production.
species_params(p)$production_observed <- getSomaticProduction(p)

# Stage 2: match the size distributions ---------------------------------------

pm <- matchCatch(p, catch = catch)

# Stage 3: satiation ----------------------------------------------------------

# The von Bertalanffy model already carries a metabolic loss - the `B w` term of
# the growth equation - so unlike the allometric build there is no need to
# invent one via a critical feeding level. What is still missing is satiation:
# `h = Inf`, so the fish eat in strict proportion to the food available.
# `setFeedingLevel()` gives every species a size-independent feeding level by
# lowering `h` and raising the encounter rate together, leaving consumption -
# and hence growth, the size spectrum and every Ecopath ratio - unchanged.
pm <- setFeedingLevel(pm, feeding_level = 0.6)

# Stage 4: mortality enough to carry the diet matrix --------------------------

# Imposing the diet matrix converts part of the external mortality into explicit
# predation mortality, and for Cod and Hake the diet matrix demands more
# predation than their total mortality provides, which would need a negative
# external mortality.
#
# Hake gets more total mortality - by the steady-state identity of stage 2, more
# production. A factor of 1.5 is enough; it was found by scanning upwards until
# the `matchDiet()` warning disappeared.
species_params(pm)["Hake", "production_observed"] <-
    1.5 * species_params(pm)["Hake", "production_observed"]
pm <- matchCatch(pm, catch = catch, species = "Hake")

# For Cod the whole of the demand is cannibalism, because nothing else in the
# model eats Cod. Rather than inflate Cod's mortality to accommodate it, we
# weaken the cannibalism entry of the diet matrix to 20% of its Ecopath value
# and hand the difference to the `other` column. That leaves Cod's fitted
# mortality untouched *and* gives Cod a lower yield peak in stage 8, which is
# where it is needed - see the scan in the vignette. Below 0.2 the yield peak
# rises steeply again (1.76 at lambda = 0), and above about 0.25 the clipping of
# Cod's external mortality returns, so 0.2 is close to the best available.
lambda_cod <- 0.2
dm <- reduced_dm
dm["Cod", "other"] <- dm["Cod", "other"] + (1 - lambda_cod) * dm["Cod", "Cod"]
dm["Cod", "Cod"]   <- lambda_cod * dm["Cod", "Cod"]

# Stage 5: species interactions from the diet matrix --------------------------

pd <- matchDiet(pm, dm)
# Residual warnings about negative external *encounter* (Megrim, Whiting, Hake)
# are clipped at zero; they affect only a few large individuals, which is why we
# recompute the steady state afterwards.
ps <- steady(pd, tol = 1e-10)

# Stage 6: the plankton resource ---------------------------------------------

psr <- alignResource(ps)
resource_params(psr)$w_pp_cutoff <- 1
initialNResource(psr)[w_full(psr) > 1] <- 0
comment(psr@cc_pp) <- NULL
psr <- setResourceInteraction(psr,
                              resource_dynamics = "resource_semichemostat",
                              tol = 1e-2)
psr <- steady(psr, tol = 1e-10)
resource_level(psr) <- 0.5
psr <- steady(psr, tol = 1e-12, t_max = 200)

# Stage 7: the response knobs -------------------------------------------------

# Six of the nine species can be brought onto their FMSY by the reproduction
# level alone. Herring, Cod and Whiting cannot, and each needs a different knob.
# All three knobs below leave the steady state, the biomasses, the size
# distributions and the yields exactly where they are.
allsp <- species_params(psr)$species

# Herring gets 96% of its food from the resource, which makes its growth response
# to fishing enormous and puts its yield peak at 2.6 x FMSY even with no density
# dependence in reproduction. Raising its feeding level to 0.94 damps that
# response - the elasticity of consumption to food is 1 - f - and brings the peak
# onto target.
f <- setNames(rep(0.6, 9), allsp)
f[["Herring"]] <- 0.94
psr <- setFeedingLevelInteracting(psr, f)

# For Cod and Whiting the feeding level does almost nothing and the compensation
# runs through their predation on the other modelled fish. Reallocating that
# predation to the inert external terms removes it. Whiting needs a mild
# reduction (its fish diet falls from 10.1% to 6.6% of its consumption); Cod
# needs nearly all of it (17.8% to 0.18%), which is the largest concession this
# build makes to the FMSY targets - see the vignette.
psr <- set_link(psr, "Whiting", allsp, 0.65)
psr <- set_link(psr, "Cod", allsp, 0.01)

saveRDS(psr, "inst/params_vB_stage7_Richard_and_Jess.rds")

# Stage 8: reproduction levels ----------------------------------------------

# Choose each species' reproduction level so that the *global* maximum of its
# yield-vs-F curve sits at its FMSY (at the fitted F for Plaice, which has no
# FMSY). Several of these curves are bimodal, so the cheap two-point "which side
# of the peak" test used in the allometric build is not safe here and we locate
# the global maximum of a whole curve instead.
#
# Tuning one species shifts the curves of the others, so the sweep is iterated.
# It is run as a parallel Jacobi iteration - all nine species tuned from the same
# starting point, then all the reproduction levels applied at once - because a
# sequential sweep takes several hours. Three iterations take about 15 minutes on
# six cores, and iterations 2 and 3 give identical answers.
source("inst/tune_repro_level_vB.R")
library(parallel)

Ft <- F_targets(psr)
params <- psr
for (sweep in 1:3) {
    res <- mclapply(names(Ft), function(s) {
        r <- tune_peak(params, s, Ft[[s]], steps = 7, tol = 0.10)
        r$species <- s
        r
    }, mc.cores = 6)
    for (r in res) params <- set_rl(params, r$species, r$rl)
    print(do.call(rbind, lapply(res, as.data.frame)))
}

# The tuning above uses a deliberately cheap yield curve. Recompute the peaks on
# a finer grid with a much tighter convergence tolerance - this is the number to
# quote, and for Cod it differs from the cheap one by 0.24 because 60 years is
# not long enough for the longest-lived species to re-equilibrate.
MULT_FINE <- exp(seq(log(0.15), log(2.6), length.out = 15))
peaks <- t(vapply(names(Ft), function(s) {
    y <- suppressMessages(getYieldVsF(params, s, gear = "commercial",
                                      F_range = Ft[[s]] * MULT_FINE,
                                      tol = 1e-3, t_max = 200))
    peak_ratio(y[order(y$F), ], Ft[[s]])
}, numeric(2)))
print(cbind(F_target = Ft, F_peak = peaks[, "ratio"] * Ft, peaks))

# Save ----
saveRDS(params, "inst/params_vB_final_Richard_and_Jess.rds")

# The alternative model of the "Cod" section of the vignette: moving Cod's
# commercial l50 from the fitted 58.3 cm to 52 cm brings all nine species within
# 11% of their FMSY, at the cost of about 100 log-likelihood units on Cod's own
# catch-at-length. Rebuild from stage 1 with
#   gp$l50[gp$species == "Cod" & gp$gear == "commercial"] <- 52
#   gp$l25[gp$species == "Cod" & gp$gear == "commercial"] <- 45
# and hold them fixed in the stage-2 fit with
#   map = list(m = factor(NA), l50 = factor(c(NA, NA)), l25 = factor(c(NA, NA)))
# for Cod, then repeat stages 4 to 8 unchanged. Saved as
# inst/params_vB_l50_variant_Richard_and_Jess.rds.
