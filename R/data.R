#' Example species parameters for use with mizerEcopath
#'
#' A curated subset of species parameters derived from the Celtic Sea model of Lauria et al. (2016),
#' formatted for compatibility with mizer and mizerEcopath. Includes mapped Ecopath groups, life history
#' traits, and allometric scaling parameters.
#'
#' @format A data frame with 12 rows and 18 columns. Key columns include:
#' \describe{
#'   \item{species}{Common name of the species}
#'   \item{w_max}{Maximum body weight (g)}
#'   \item{w_mat}{Maturity weight (g)}
#'   \item{a, b}{Length-weight relationship parameters}
#'   \item{n, p, d}{Allometric exponents for consumption, metabolism, and mortality}
#'   \item{alpha}{Assimilation efficiency}
#'   \item{ecopath_groups}{List column mapping Ecopath stanza groups}
#' }
#' @name species_params_example
#' @docType data
#' @keywords datasets
species_params_example

#' Example Ecopath diet matrix
#'
#' A group-level diet composition matrix derived from Lauria et al. (2016),
#' used as input to `reduceEcopathDiet()`. Columns and rows include juvenile and adult stanzas.
#'
#' @format A data frame with predator groups as row names and prey groups as columns.
#' @name ecopath_diet_example
#' @docType data
#' @keywords datasets
ecopath_diet_example

#' Example reduced diet matrix for mizer
#'
#' The output of applying `reduceEcopathDiet()` to the Celtic Sea Ecopath diet matrix.
#' Rows and columns are matched to `species_params_example$species`.
#'
#' @format A matrix with species as both row and column names, giving the proportion of each prey
#' in each predator’s diet.
#' @name diet_matrix_example
#' @docType data
#' @keywords datasets
diet_matrix_example
