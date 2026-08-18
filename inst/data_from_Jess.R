# Combine Ecopath biomasses, and yields with
# with Jess's survey and catch size distributions
library(dplyr)
library(ggplot2)
library(mizer)
library(mizerEcopath)

survey_Jess <- survey_length_distribution

catch_Jess <- catch_distribution

catch <- rbind(catch_Jess, survey_Jess)
l_max <- catch |>
    group_by(species) |>
    summarise(l_max = max(length + dl) * 1.1)

# Species params ----
ICES_FMSY <- readr::read_csv("inst/ICES_FMSY.csv") |>
    group_by(SpeciesName) |>
    summarise(FMSY = mean(FMSY, na.rm = TRUE)) |>
    rename(Scientific_name = SpeciesName)
Jess_sp <- celtic_params@species_params |>
    left_join(ICES_FMSY) |>
    select(species, Scientific_name, a, b, w_mat, age_mat,
           pred_kernel_type, beta, sigma,
           kernel_exp, kernel_l_l, kernel_u_l,
           kernel_l_r, kernel_u_r, FMSY)
James_sp <- readRDS("inst/James_sp.rds")
sp <- James_sp |>
    select(species, w_inf, k_vb, t0,
           biomass_cutoff, biomass_observed,
           production_observed, consumption_observed) |>
    inner_join(Jess_sp, by = "species")

rownames(sp) <- sp$species
sp["Cod", "kernel_l_l"] <- 3.8

sp$l_inf <- NA
sp["Hake", c("k_vb", "l_inf", "t0")] <- c(0.178, 110, 0)
sp["Hake", "w_inf"] <- NA
sp["Red gurnard", c("k_vb", "l_inf", "l_mat", "age_mat", "t0")]  <-
    c(0.21, 42.4, 27.59, 3.7, -1.21)
sp["Red gurnard", c("w_inf", "w_mat")] <- NA

sp <- sp |>
    inner_join(l_max)

sp$D_ext <- 1

# Gear params ----

gp <- gear_params(celtic_params) |>
    filter(species %in% sp$species)
gp$yield_weight <- 1
gp$yield_weight[gp$gear == "survey"] <- 0

# Save data
save(sp, gp, catch, file = "inst/data_Jess.rda")
