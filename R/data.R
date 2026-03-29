#' Example species parameters from a few Celtic Sea species.
#'
#' A curated subset of species parameters derived from the Celtic Sea model of Lauria et al. (2016),
#' as well as arbitrary power-law exponents and stanza mapping dictionary.
#'
#' @format A data frame with 12 rows and 18 columns. Key columns include:
#' \describe{
#'   \item{species}{Common name of the species}
#'   \item{w_max}{Maximum body weight (g)}
#'   \item{w_mat}{Maturity weight (g)}
#'   \item{age_mat}{Age at maturity (years)}
#'   \item{a, b}{Length-weight relationship parameters}
#'   \item{n, p, d}{Allometric exponents for consumption, metabolism, and mortality}
#'   \item{alpha}{Assimilation efficiency}
#'   \item{ecopath_groups}{List column mapping Ecopath stanza groups}
#'   \item{biomass_observed}{Observed biomass (t/km^2) from Ecopath}
#'   \item{consumption_observed}{Observed consumption (t/km^2/year) from Ecopath}
#'   \item{production_observed}{Observed production (t/km^2/year) from Ecopath}
#' }
#' @name species_params_example
#' @docType data
#' @keywords datasets
"species_params_example"

#' Celtic Sea MizerParams object
#'
#' A 12-species MizerParams object for the Celtic Sea, built from the Ecopath
#' model of Lauria et al. (2016). Used as the main example object in tests and
#' vignettes.
#'
#' @format A MizerParams object with 12 species.
#' @name celtic_params
#' @docType data
#' @keywords datasets
"celtic_params"

#' Celtic Sea catch size distributions
#'
#' Observed catch size distributions from DATRAS survey data for Celtic Sea
#' species, used to calibrate gear selectivity via [matchCatch()].
#'
#' @format A data frame with columns:
#' \describe{
#'   \item{species_code}{ICES species code}
#'   \item{length}{Length class (cm)}
#'   \item{gear}{Survey gear identifier}
#'   \item{catch}{Observed catch count}
#'   \item{Scientific_name}{Scientific name}
#'   \item{English_name}{Common English name}
#'   \item{species}{Species name matching MizerParams}
#'   \item{dl}{Width of the length bin (cm)}
#' }
#' @name celtic_catch
#' @docType data
#' @keywords datasets
"celtic_catch"

#' North Sea 3-species model: initial state
#'
#' A 3-species MizerParams object for the North Sea before calibration.
#'
#' @format A MizerParams object with 3 species (Cod, Haddock, Sprat).
#' @name ns_3_spp_model_initial
#' @docType data
#' @keywords datasets
"ns_3_spp_model_initial"

#' North Sea 3-species model: trial state
#'
#' A 3-species MizerParams object for the North Sea at an intermediate
#' calibration stage.
#'
#' @format A MizerParams object with 3 species (Cod, Haddock, Sprat).
#' @name ns_3_spp_model_trial
#' @docType data
#' @keywords datasets
"ns_3_spp_model_trial"

#' North Sea 3-species model: final calibrated state
#'
#' A 3-species MizerParams object for the North Sea after full calibration.
#'
#' @format A MizerParams object with 3 species (Cod, Haddock, Sprat).
#' @name ns_3_spp_model_final
#' @docType data
#' @keywords datasets
"ns_3_spp_model_final"

#' Fishing deaths data frame
#'
#' Total weight of dead discards per gear per area, used in discard rate
#' calculations.
#'
#' @format A data frame with columns:
#' \describe{
#'   \item{species}{Species name}
#'   \item{total_weight_dead_gear_per_area}{Total weight of dead fish per
#'     gear per area (t/km^2/year)}
#'   \item{gear}{Gear category}
#' }
#' @name fishing_deaths
#' @docType data
#' @keywords datasets
"fishing_deaths"

#' Biomass discard rates
#'
#' Discard rates by gear category and species, derived from STECF data.
#'
#' @format A grouped data frame with columns:
#' \describe{
#'   \item{gear_category}{Gear category}
#'   \item{species}{Species name}
#'   \item{country}{Country}
#'   \item{target_assemblage}{Target assemblage}
#'   \item{total_landings}{Total landings (t)}
#'   \item{tot_discards_tonnes}{Total discards (t)}
#'   \item{discard_rate}{Discard rate (proportion)}
#' }
#' @name biomass_discard_rates
#' @docType data
#' @keywords datasets
"biomass_discard_rates"

#' Catch size distribution data
#'
#' Observed catch size distributions for Celtic Sea species by gear.
#'
#' @format A data frame with columns:
#' \describe{
#'   \item{species}{Species name matching MizerParams}
#'   \item{Scientific_name}{Scientific name}
#'   \item{length}{Length class (cm)}
#'   \item{dl}{Width of the length bin (cm)}
#'   \item{gear}{Gear identifier}
#'   \item{catch}{Observed catch count}
#' }
#' @name catch_distribution
#' @docType data
#' @keywords datasets
"catch_distribution"

#' Celtic Sea age-at-size data
#'
#' Age-at-size observations from DATRAS survey data for Celtic Sea species,
#' used for fitting diffusion-based growth models.
#'
#' @format A data frame with columns including Survey, Quarter,
#'   Scientific_name, Year, Month, Day, a, b, WgtClass, LngtClass, IndWgt,
#'   Maturity, Age, and CANoAtLngt.
#' @name cs_age_size
#' @docType data
#' @keywords datasets
"cs_age_size"

#' Celtic Sea calculated life history parameters
#'
#' Estimated length at maturity and age at maturity for Celtic Sea species,
#' derived from survey data.
#'
#' @format A data frame with columns:
#' \describe{
#'   \item{Scientific_name}{Scientific name}
#'   \item{l_mat}{Length at maturity (cm)}
#'   \item{age_mat}{Age at maturity (years)}
#' }
#' @name cs_calculated_histories
#' @docType data
#' @keywords datasets
"cs_calculated_histories"

#' Ecopath diet matrix (raw import)
#'
#' A raw Ecopath diet composition matrix as imported from an Ecopath export
#' file, with 67 groups and numeric column headers.
#'
#' @format A data frame with 67 rows and 69 columns. The first two columns
#'   are row numbers and group names; the remaining columns are predator
#'   group diet proportions.
#' @name diet
#' @docType data
#' @keywords datasets
"diet"

#' Celtic Sea Ecopath diet matrix
#'
#' The Ecopath diet composition matrix for the Celtic Sea model of
#' Lauria et al. (2016), including juvenile and adult stanzas.
#'
#' @format A data frame with 75 rows and 66 columns. Row names are prey
#'   groups; column names are predator groups.
#' @name diets
#' @docType data
#' @keywords datasets
"diets"

#' Life history parameters from FishBase
#'
#' Life history parameters for Celtic Sea species obtained from FishBase
#' using [fillDefaultsFromFishBase()].
#'
#' @format A data frame with columns:
#' \describe{
#'   \item{species}{Common name of the species}
#'   \item{Scientific_name}{Scientific name}
#'   \item{age_mat}{Age at maturity (years)}
#'   \item{l_mat}{Length at maturity (cm)}
#'   \item{l_max}{Maximum length (cm)}
#'   \item{a}{Length-weight coefficient}
#'   \item{b}{Length-weight exponent}
#' }
#' @name life_history_fishbase
#' @docType data
#' @keywords datasets
"life_history_fishbase"

#' Species parameters for North Sea 14-species model
#'
#' Species parameters data frame used to build a 14-species North Sea
#' mizer model.
#'
#' @format A data frame with 14 rows and columns including species,
#'   Scientific_name, age_mat, l_mat, l_max, a, b, w_mat, w_max, n, p, d,
#'   alpha, and biomass_cutoff.
#' @name sp_params
#' @docType data
#' @keywords datasets
"sp_params"

#' Stomach content data fit
#'
#' Parameters of fitted distributions for stomach content analysis, used in
#' diet matching.
#'
#' @format A data frame with 14 rows and columns including species, alpha,
#'   ll, ul, lr, ur, distribution, and min_w_pred.
#' @name stomach_data_fit
#' @docType data
#' @keywords datasets
"stomach_data_fit"

#' Survey length distribution
#'
#' Observed length distributions from survey data for Celtic Sea species.
#'
#' @format A data frame with columns:
#' \describe{
#'   \item{length}{Length class (cm)}
#'   \item{dl}{Width of the length bin (cm)}
#'   \item{gear}{Survey gear identifier}
#'   \item{catch}{Observed count}
#'   \item{Scientific_name}{Scientific name}
#'   \item{species}{Species name matching MizerParams}
#' }
#' @name survey_length_distribution
#' @docType data
#' @keywords datasets
"survey_length_distribution"

#' Ecopath diet matrix from Lauria et al. (2016)
#'
#' A group-level diet composition matrix derived from Lauria et al. (2016),
#' used as input to `reduceEcopathDiet()`. Columns and rows include juvenile and adult stanzas.
#'
#' @format A data frame with predator groups as row names and prey groups as columns.
#' @name ecopath_diet_example
#' @docType data
#' @keywords datasets
"ecopath_diet_example"

#' Example reduced diet matrix for mizer
#'
#' The output of applying `reduceEcopathDiet()` to the example diet matrix.
#' Rows and columns are matched to `species_params_example$species`.
#'
#' @format A matrix with species as both row and column names, giving the proportion of each prey
#' in each predator’s diet.
#' @name diet_matrix_example
#' @docType data
#' @keywords datasets
"diet_matrix_example"
