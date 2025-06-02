#devtools::document()
#TMB::compile("src/objective_function.cpp")
#dyn.load(TMB::dynlib("src/objective_function"))
#devtools::load_all()


# --------------------------- 1. Load packages ---------------------------
library(rfishbase)         # for FishBase look-ups
library(dplyr)
library(mizer)
library(mizerExperimental)
library(ggplot2)
library(RColorBrewer)
devtools::load_all()

# --------------------------- 2. Pick species ---------------------------
# Select species for testing
species_test <- c("Hake", "Mackerel", "Blue whiting")

# Assign  colours based on number of species
n_species <- length(species_test)
cols <- if (n_species < 3) {
    # Use base colours if 1 or 2 species
    c("firebrick", "steelblue", "darkgreen")[1:n_species]
} else {
    RColorBrewer::brewer.pal(n = n_species, name = "Dark2")
}

# Full lookup table with common and scientific names
species_map_all <- data.frame(
    species         = c("Hake", "Mackerel", "Blue whiting"),
    Scientific_name = c("Merluccius merluccius",
                        "Scomber scombrus",
                        "Micromesistius poutassou"),
    stringsAsFactors = FALSE
)

# Filter just the species you want to test
species_map <- species_map_all %>%
    filter(species %in% species_test)

# ------- 3. Pull FishBase defaults (w_max, w_mat, a, b, age_mat, …) -----
fb <- fillDefaultsFromFishBase(species_map, overwrite = FALSE, verbose = FALSE)

sp <- fb %>%
    select(species, w_max, w_mat, a, b, age_mat, Length) %>%
    mutate(
        n = 0.7,
        p = 0.7,
        d = -0.3,
        alpha = 0.8
    )

# -------- 4. Minimal Ecopath “Basic estimates” table --------------------
ecopath_basic_all <- data.frame(
    `...1`                             = 1:3,
    `Group name`                       = c("Blue whiting", "Hake", "Mackerel"),
    `Biomass (t/km²)`                  = c(0.444, 0.260, 15.653),
    `Consumption / biomass (/year)`    = c(6.666, 3.529, 1.730),
    `Production / consumption (/year)` = c(0.165, 0.312, 0.376),
    check.names = FALSE
)

ecopath_basic <- ecopath_basic_all %>%
    filter(`Group name` %in% species_test)

# ------------- 5. Map model species to Ecopath groups -------------------
species_to_groups <- as.list(fb$species)
names(species_to_groups) <- fb$species

# ------------- 6. Prepare species_params & add Ecopath info -------------
sp <- sp %>%
    mutate(
        w_max = a * Length^b,
        w_mat = w_mat,
        w_repro_max = w_max,
        n     = 0.7,
        p     = 0.7,
        alpha = 0.8
    )

sp <- addEcopathParams(sp, ecopath_basic, species_to_groups)

# -------- 7. Build a non-interacting allometric model -------------------
params <- newAllometricParams(sp, no_w = 200)   # steady state

# -------- 8. Add hand-made total catch (t km-2 year-1) ------------------
catch_df <- data.frame(
    `Group name`                       = sp$species,
    `TotalCatch (t/km²/year)`          = c(0.3, 0.1, 0.5),
    check.names = FALSE
)

# -------- 9. Add single-sigmoid length selectivity gear ----------------
params_ss <- addEcopathCatchTotal(
    params,
    catch_df,
    sel_func = "sigmoid_length"
)

# -------- 10. Plot original single-sigmoid selectivity -----------------
gp <- gear_params(params_ss)
plot(NA, xlim = c(0, 1.6 * max(gp$l50, na.rm = TRUE)), ylim = c(0, 1),
     xlab = "Length (cm)", ylab = "Selectivity",
     main = "Single-sigmoid selectivity (gear = 'total')")

i <- 0
for (sp_row in split(gp, rownames(gp))) {
    i   <- i + 1
    l50 <- sp_row$l50
    l25 <- sp_row$l25
    k   <- -log(3) / (l25 - l50)
    L   <- seq(0, 1.6 * l50, length.out = 200)
    S   <- 1 / (1 + exp(-k * (L - l50)))
    lines(L, S, col = cols[i], lwd = 2)
}
legend("bottomright", legend = gp$species, lwd = 2, col = cols, bty = "n")

# -------- 11. Add double-sigmoid selectivity gear ----------------------
params_dome <- addEcopathCatchTotal(
    params,
    catch_df,
    sel_func = "double_sigmoid_length"
)

gp2 <- gear_params(params_dome)
plot(NA, xlim = c(0, 1.6 * max(gp2$l50_right, na.rm = TRUE)), ylim = c(0, 1),
     xlab = "Length (cm)", ylab = "Selectivity",
     main = "Double-sigmoid (dome-shaped) selectivity")

i <- 0
for (sp_row in split(gp2, rownames(gp2))) {
    i   <- i + 1
    if (anyNA(sp_row$l50_right)) next  # <-- skip if optimisation not yet run
    kL <- -log(3) / (sp_row$l25 - sp_row$l50)
    kR <-  log(3) / (sp_row$l25_right - sp_row$l50_right)
    L  <- seq(0, 1.6 * sp_row$l50_right, length.out = 300)
    S_left  <- 1 / (1 + exp(-kL * (L - sp_row$l50)))
    S_right <- 1 / (1 + exp( kR * (L - sp_row$l50_right)))
    lines(L, S_left * S_right, col = cols[i], lwd = 2)
}
legend("bottomright", legend = gp2$species, lwd = 2, col = cols, bty = "n")

# --------------------------- 12. Load catch data ------------------------
path_to_fmp_data <- "C:/Documents/York academic/Projects/R projects/Mizer_FMP_EwE/Data"
landings_total <- readRDS(file.path(path_to_fmp_data, "catch_with_observer.rds"))

landings_total <- landings_total %>%
    mutate(dl = 1,
           gear = ifelse(gear == "commercial", "total", gear),
           species = ifelse(species == "Horse Mackerel", "Horse mackerel", species)) %>%
    filter(gear == "total") %>%
    group_by(species, length) %>%
    summarise(dl = first(dl), count = sum(catch), .groups = "drop") %>%
    filter(species %in% species_test) %>%
    arrange(species, length)

ggplot(landings_total, aes(x = length, y = count, colour = species)) +
    geom_line(linewidth = 1) +
    facet_wrap(~species, scales = "free_y") +
    labs(x = "Length (cm)", y = "Catch count", title = "Observed catch size distributions") +
    theme_minimal()
# ------------------- 13. Run matchCatch with real data -----------------
# --- Single-sigmoid: prep before matching catch ---
params_ss_mc <- params_ss |>
    matchGrowth() |>
    steadySingleSpecies() |>
    matchBiomasses() |>
    matchCatch(catch = landings_total)

# --- Double-sigmoid: prep before matching catch and change penalties ---
params_ds_mc <- params_dome |>
    matchGrowth() |>
    steadySingleSpecies() |>
    matchBiomasses() |>
    matchCatch(
        catch = landings_total,
        yield_lambda      = 0.25,   # was 1
        production_lambda = 0.25    # was 1
    )

# ------------------- 14. Compare before and after params ---------------
# Single-sigmoid
gp_ss_before <- gear_params(params_ss)[, c("l25", "l50")]
gp_ss_after  <- gear_params(params_ss_mc)[, c("l25", "l50")]

cat("\n--- Single-sigmoid gear parameters ---\n")
cat("Before matchCatch:\n")
print(gp_ss_before)
cat("\nAfter matchCatch:\n")
print(gp_ss_after)

# Double-sigmoid
gp_ds_before <- gear_params(params_dome)[, c("l25", "l50", "l50_right", "l25_right")]
gp_ds_after  <- gear_params(params_ds_mc)[, c("l25", "l50", "l50_right", "l25_right")]

cat("\n--- Double-sigmoid gear parameters ---\n")
cat("Before matchCatch:\n")
print(gp_ds_before)
cat("\nAfter matchCatch:\n")
print(gp_ds_after)

# ------------------- 15. Plot: Single-sigmoid before vs after ----------
plot(NA, xlim = c(0, 200), ylim = c(0, 1),
     xlab = "Length (cm)", ylab = "Selectivity",
     main = "Single-sigmoid: before (dotted) vs after (solid)")
i <- 0
for (sp in rownames(gp_ss_before)) {
    i <- i + 1
    col <- cols[i]

    # Before
    l25_b <- gp_ss_before[sp, "l25"]
    l50_b <- gp_ss_before[sp, "l50"]
    k_b   <- -log(3) / (l25_b - l50_b)
    L     <- seq(0, 1.6 * l50_b, length.out = 300)
    S_b   <- 1 / (1 + exp(-k_b * (L - l50_b)))
    lines(L, S_b, col = col, lwd = 2, lty = 2)

    # After
    l25_a <- gp_ss_after[sp, "l25"]
    l50_a <- gp_ss_after[sp, "l50"]
    k_a   <- -log(3) / (l25_a - l50_a)
    L     <- seq(0, 1.6 * l50_a, length.out = 300)
    S_a   <- 1 / (1 + exp(-k_a * (L - l50_a)))
    lines(L, S_a, col = col, lwd = 2, lty = 1)
}
# Remove ", total" to align with species_test
legend_labels <- gsub(", total", "", rownames(gp_ss_before))
legend("bottomright", legend = legend_labels, col = cols, lwd = 2, bty = "n")

# ------------------- 16. Plot: Double-sigmoid before vs after ----------
plot(NA, xlim = c(0, 1.6 * max(gp_ds_before$l50_right, na.rm = TRUE)), ylim = c(0, 1),
     xlab = "Length (cm)", ylab = "Selectivity",
     main = "Double-sigmoid: before (dotted) vs after (solid)")
i <- 0
for (sp in rownames(gp_ds_before)) {
    i <- i + 1
    col <- cols[i]
    if (anyNA(gp_ds_after[sp, c("l50_right", "l25_right")])) next  # <-- skip if not defined

    # Before
    L     <- seq(0, 1.6 * gp_ds_before[sp, "l50_right"], length.out = 300)
    kL_b  <- -log(3) / (gp_ds_before[sp, "l25"] - gp_ds_before[sp, "l50"])
    kR_b  <-  log(3) / (gp_ds_before[sp, "l25_right"] - gp_ds_before[sp, "l50_right"])
    S_left_b  <- 1 / (1 + exp(-kL_b * (L - gp_ds_before[sp, "l50"])))
    S_right_b <- 1 / (1 + exp( kR_b * (L - gp_ds_before[sp, "l50_right"])))
    lines(L, S_left_b * S_right_b, col = col, lwd = 2, lty = 2)

    # After
    L     <- seq(0, 1.6 * gp_ds_after[sp, "l50_right"], length.out = 300)
    kL_a  <- -log(3) / (gp_ds_after[sp, "l25"] - gp_ds_after[sp, "l50"])
    kR_a  <-  log(3) / (gp_ds_after[sp, "l25_right"] - gp_ds_after[sp, "l50_right"])
    S_left_a  <- 1 / (1 + exp(-kL_a * (L - gp_ds_after[sp, "l50"])))
    S_right_a <- 1 / (1 + exp( kR_a * (L - gp_ds_after[sp, "l50_right"])))
    lines(L, S_left_a * S_right_a, col = col, lwd = 2, lty = 1)
}

# Clean species names for legend
legend_labels <- gsub(", total", "", rownames(gp_ds_before))
legend("bottomright", legend = legend_labels, col = cols, lwd = 2, bty = "n")

# -------------- 17. Sensitivity to starting values (revised) --------------

run_matchCatch_with_start <- function(params, species, catch, jitter_factors) {
    sp_row <- species_params(params) %>% filter(species == !!species)
    gp_row <- gear_params(params) %>% filter(species == !!species)

    par0 <- list(
        l50          = gp_row$l50 * jitter_factors["l50"],
        ratio        = (gp_row$l25 / gp_row$l50) * jitter_factors["ratio"],
        d50          = (gp_row$l50_right - gp_row$l50) * jitter_factors["d50"],
        mu_mat       = ifelse(is.na(sp_row$mu_mat),
                              ext_mort(params)[sp_row$species == species, which.min(abs(w(params) - sp_row$w_mat))],
                              sp_row$mu_mat * jitter_factors["mu_mat"]),
        catchability = pmax(gp_row$catchability * jitter_factors["catchability"], 1e-8),
        r_right      = (gp_row$l25_right / gp_row$l50_right) * jitter_factors["r_right"]
    )

    # Replace any non-finite/NA starting value with 1 just to avoid an immediate crash
    par0 <- lapply(par0, function(x) ifelse(is.finite(x), x, 1))

    # Prepare the TMB data as in matchCatch()
    data_obj <- prepare_data(params, species = species, catch = catch,
                             yield_lambda = 0.25, production_lambda = 0.25)
    obj <- TMB::MakeADFun(data = data_obj,
                          parameters = par0,
                          DLL = "mizerEcopath",
                          silent = TRUE)

    # Run optimizer inside try() so we catch R-level errors
    opt <- try(
        nlminb(start = obj$par, objective = obj$fn, gradient = obj$gr,
               lower = rep(-Inf, length(obj$par)),
               upper = rep( Inf, length(obj$par))),
        silent = TRUE
    )
    if (inherits(opt, "try-error")) {
        return(NULL)  # a true R error, so we can’t evaluate at all
    }
    # If opt$objective is NaN or Inf, that means we started outside the valid domain
    if (!is.finite(opt$objective)) {
        return(NULL)
    }
    list(opt = opt, par = opt$par, objective = opt$objective)
}

set.seed(42)
n_runs <- 30

# (Example “broad” jitter — we’ll fix these ranges next)
jitter_grid <- replicate(n_runs, {
    c(l50          = runif(1, 0.8, 1.2),
      ratio        = runif(1, 0.8, 1.2),
      d50          = runif(1, 0.8, 1.2),
      mu_mat       = runif(1, 0.2, 2),
      catchability = runif(1, 0.2, 2),
      r_right      = runif(1, 0.8, 1.2))
}, simplify = FALSE)

# Run all 30 attempts
results <- lapply(jitter_grid, function(jit) {
    run_matchCatch_with_start(params_ds_mc, species = "Hake",
                              catch = landings_total,
                              jitter_factors = jit)
})

# Filter only the truly finite–objective runs
finite_runs <- Filter(Negate(is.null), results)

cat("Total attempts:         ", length(results), "\n")
cat("Finite-objective runs:  ", length(finite_runs), "\n\n")

if (length(finite_runs) > 0) {
    objective_values <- sapply(finite_runs, function(x) x$objective)
    print(summary(objective_values))
    hist(objective_values,
         main = "Objective values for truly valid starts",
         xlab = "Negative log‐likelihood",
         col = "skyblue", breaks = 10)
    abline(v = min(objective_values), col = "red", lwd = 2)
} else {
    cat("No finite-objective runs to summarize.\n")

    # --------------------- Inspect outlier starting points ---------------------

    # Print default parameter values from unperturbed model
    sp_row <- species_params(params_ds_mc) %>% filter(species == "Hake")
    gp_row <- gear_params(params_ds_mc) %>% filter(species == "Hake")

    default_vals <- list(
        l50          = gp_row$l50,
        ratio        = gp_row$l25 / gp_row$l50,
        d50          = gp_row$l50_right - gp_row$l50,
        mu_mat       = ifelse(is.na(sp_row$mu_mat),
                              ext_mort(params_ds_mc)[sp_row$species == "Hake", which.min(abs(w(params_ds_mc) - sp_row$w_mat))],
                              sp_row$mu_mat),
        catchability = gp_row$catchability,
        r_right      = gp_row$l25_right / gp_row$l50_right
    )

    cat("\n--- Default parameter values (unperturbed) ---\n")
    print(default_vals)

    # Determine IQR and threshold for 'outlier' objective values
    q1 <- quantile(objective_values, 0.25)
    q3 <- quantile(objective_values, 0.75)
    iqr <- q3 - q1
    threshold <- q3 + 2 * iqr

    cat("\n--- Outlier threshold (objective > ", round(threshold), ") ---\n", sep = "")
    outlier_indices <- which(objective_values > threshold)

    if (length(outlier_indices) == 0) {
        cat("No outlier runs found.\n")
    } else {
        for (i in outlier_indices) {
            cat("\n⚠️  Outlier run", i, " - Objective =", round(objective_values[i]), "\n")
            cat("Jittered starting values:\n")
            print(jitter_grid[[i]])
        }
    }

}
