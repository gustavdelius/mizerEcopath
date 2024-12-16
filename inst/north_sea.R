library(mizerEcopath)
library(rfishbase)
library(dplyr)
library(readr)

## Prepare species parameters ----

species_to_latin <- c(
    "Sprat"   = "Sprattus sprattus",
    "Sandeel" = "Ammodytes tobianus",
    "N.pout"  = "Trisopterus esmarkii",
    "Herring" = "Clupea harengus",
    "Dab"     = "Limanda limanda",
    "Whiting" = "Merlangius merlangus",
    "Sole"    = "Solea solea",
    "Gurnard" = "Eutrigla gurnardus",
    "Plaice"  = "Pleuronectes platessa",
    "Haddock" = "Melanogrammus aeglefinus",
    "Cod"     = "Gadus morhua",
    "Saithe"  = "Pollachius virens"
)
sp <- data.frame(species = names(species_to_latin),
                 latin_name = species_to_latin)


# Set maximum sizes from observed sizes

catch <- readRDS("inst/extdata/ns_catch.rds")
# Select "total" gear only
catch <- catch[catch$gear == "total", ]
max_size <- catch |>
    group_by(species) |>
    summarise(l_max = max(length) * 1.2)
missing <- !(sp$species %in% max_size$species)
max_size_fishbase <- rfishbase::species(sp$latin_name[missing]) |>
    select(latin_name = Species, l_max = Length) |>
    left_join(select(sp, species, latin_name),
              by = "latin_name")
max_size <- bind_rows(max_size, max_size_fishbase) |>
    select(species, l_max)

# Add length-weight parameters
length_weight <- estimate(sp$latin_name, fields = c("Species", "a", "b"))

sp <- sp |>
    left_join(length_weight, by = c("latin_name" = "Species")) |>
    left_join(max_size) |>
    mutate(w_max = a * l_max ^ b)

# Add maturity parameters
maturity_tbl <- rfishbase::maturity(species_to_latin)
median_maturity <- maturity_tbl |>
    group_by(Species) |>
    filter(!is.na(tm), !is.na(Lm)) |>
    summarise(age_mat = median(tm),
              l_mat = median(Lm))
sp <- sp |>
    left_join(median_maturity, by = c("latin_name" = "Species")) |>
    mutate(w_mat = a * l_mat ^ b)

comment(sp$l_mat) <- "Median of `Lm` over all observations on the 'maturity' table on FishBase that had both `Lm` and `tm`."
comment(sp$age_mat) <- "Median of `tm` over all observations on the 'maturity' table on FishBase that had both `Lm` and `tm`."
comment(sp$w_mat) <- "Calculated from `l_mat` using weight-length parameters `a` and `b`."

# Add gonadic proportion
sp$gonad_proportion <- 0.2

## Add Ecopath species parameters ----

basic_estimates <-
    read_csv("inst/extdata/North Sea-Basic estimates.csv")

# Dictionary between species and ecopath groups
species_to_groups <- list(
    "Sprat"   = "Sprat",
    "Sandeel" = "Sandeels",
    "N.pout"  = "Norway pout",
    "Herring" = c("Herring (juvenile 0, 1)", "Herring (adult)"),
    "Dab"     = "Dab",
    "Whiting" = c("Juvenile Whiting (0-1, 0-20cm)", "Whiting (adult)"),
    "Sole"    = "Sole",
    "Gurnard" = "Gurnards",
    "Plaice"  = "Plaice",
    "Haddock" = c("Juvenile Haddock (0-1, 0-20cm)", "Haddock (adult)"),
    "Cod"     = c("Juvenile Cod(0-2, 0-40cm)", "Cod (adult)"),
    "Saithe"  = c("Juvenile Saithe (0-3, 0-40cm)", "Saithe (adult)")
)

sp <- addEcopathParams(sp, basic_estimates, species_to_groups)

## Create allometric model ----
p <- newAllometricParams(sp)

sp <- p@species_params
no_sp <- nrow(sp) # Number of species

## Set up gear params ----
ecopath_catch <- read_csv("inst/extdata/North Sea-Catch.csv")

p <- addEcopathCatchTotal(p, ecopath_catch)

# Get new steady state with desired biomass and growth
p <- p |> steadySingleSpecies() |> calibrateBiomass() |> matchGrowth() |>
    matchBiomasses() |> steadySingleSpecies()

## Match catch ----
p <- matchCatch(p, catch = catch)
## Match consumption ----
p <- matchConsumption(p)
# Fix the species with unrealistic reproductive efficiency
p <- tuneEcopath(p, catch = catch)


## Aggregate ecopath diet matrix ----
ecopath_diet <- read.csv("inst/extdata/North Sea-Diet composition.csv")
dm <- reduceEcopathDiet(sp, ecopath_diet)
