#' Calculate Size-Spectrum Fluxes and production-proportional Fishing Mortality
#'
#' @description
#' This function calculates the advective and diffusive fluxes through the size
#' spectrum. In the context of Ecosystem-Based Fisheries Management (EBFM),
#' these fluxes represent the "Production" (*P*) of the system—the total biomass
#' produced by an ecological group per unit time.
#'
#' By calculating the total flux (*J*), this function allows for the
#' implementation of species-size-level Balanced Harvesting (ssBH), where
#' fishing mortality (*F*) is set proportional to production (*P*) to
#' maintain ecosystem structure and stability.
#'
#' @param params A \code{MizerParams} object.
#' @param c_value A numeric proportionality constant (*c*). In ssBH theory,
#' this constant determines the overall fishing intensity while ensuring
#' mortality remains aligned with the ecosystem's natural energy flow.
#'
#' @details
#' ### Theoretical Context
#' Traditional fisheries management often targets specific high-value species or
#' large size classes, which can lead to "growth overfishing" and trophic
#' cascades. **Balanced Harvesting (BH)** seeks to distribute fishing
#' mortality in line with the ecosystem's productivity.
#'
#' This function implements the refined **ssBH** rule, where fishing mortality
#' is proportional to Production (*P_i(w)*):
#' \deqn{F_i(w) = c \cdot P_i(w)}
#'
#' Since Production is the instantaneous flux of individuals moving through
#' weight class *w* multiplied by weight, we use:
#' \deqn{P_i(w) = w \cdot J_i(w)}
#'
#' ### Mathematical Derivation
#'
#' **1. Total Flux (*J*):**
#' The total flux is derived by analogy with the Fokker-Planck equation,
#' accounting for both mean growth (advection) and individual variability
#' (diffusion):
#' \deqn{J_i(w) = J_{adv, i}(w) - J_{diff, i}(w)}
#'
#' **2. Advective Flux (*J_{adv}*):**
#' Represents mass flow due to deterministic growth *g_i(w)* and density *n_i(w)*:
#' \deqn{J_{adv, i}(w) = g_i(w) n_i(w)}
#'
#' **3. Diffusive Flux (*J_{diff}*):**
#' Accounts for stochasticity in growth, mitigating survivorship bias:
#' \deqn{J_{diff, i}(w) = \frac{1}{2} \frac{\partial}{\partial w} [D_i(w) n_i(w)]}
#' where *D_i(w)* is the diffusion coefficient calculated from the
#' \code{d_over_g} parameter.
#'
#' **4. Resulting Fishing Mortality (*f*):**
#' \deqn{f_i(w) = c \cdot w \cdot J_i(w)}
#'
#' @return A list containing:
#'  \code{J_adv}: Advective flux (biomass flow from mean growth).
#'  \code{J_diff}: Diffusive flux (stochastic growth spread).
#'  \code{J}: Total flux (Net biomass production rate).
#'  \code{f}: Fishing mortality aligned with ssBH rules.
#'  \code{c_value}: The proportionality constant used.
#'
#'
#' @references
#' Law, R., & Plank, M.J. (2022). Balanced exploitation and coexistence.
#' Zhou et al. (2019). Balanced harvest: Concept, policies and implementation.
#' @export
#'
#' @examples (cannot be run yet because the model is not yet saved)
#' results <- calculate_mizer_flux(params, c_value = 0.2)

compute_flux <- function(params, c_value = 1, n=initialN(params),
                         n_pp=params@initial_n_pp,
                         n_other=params@initial_n_other) {
    #Basic parameters
    w <- params@w
    species_names <- species_params(params)$species
    num_species <- length(species_names)

    #Growth and Population Density
    g <- getEGrowth(params, n = n, n_pp = n_pp, n_other = n_other)
    n <- n

    #Advective Flux (J_adv)
    J_adv <- g[, -ncol(g)] * n[, -ncol(n)]
    RDD <- as.matrix(getRDD(params, n = n, n_pp = n_pp, n_other = n_other))
    J_adv <- cbind(RDD, J_adv)
    colnames(J_adv) <- w

    #Diffusion Setup
    # We update a local copy of the params to calculate d
    temp_params <- params
    for (sp in species_names) {
        d_over_g <- species_params(temp_params)[sp, "d_over_g"]
        growth <- g[sp, ]
        n_exp <- species_params(temp_params)[sp, "n"]
        g_0 <- growth[1] / w[1]^n_exp
        d_0 <- d_over_g * g_0
        diffusion(temp_params)[sp, ] <- d_0 * w^(n_exp + 1)
    }
    d_mat <- temp_params@diffusion

    #Diffusive Flux (J_diff)
    J_diff <- matrix(NA, nrow = num_species, ncol = length(w))
    colnames(J_diff) <- w
    rownames(J_diff) <- species_names

    for(i in 1:num_species) {
        for(j in 2:length(w)) {
            # Central difference approximation for flux
            J_diff[i, j] <- (0.5) * ((d_mat[i, j] * n[i, j] - d_mat[i, j-1] * n[i, j-1]) / (w[j] - w[j-1]))
        }
    }

    # set the first size class to 0
    J_diff[, 1] <- 0

    #Total Flux (J) and Fishing Mortality (f)
    J <- J_adv - J_diff

    #Calculate f using the provided c_value
    F_Mort <- matrix(NA, nrow = num_species, ncol = length(w))
    colnames(F_Mort) <- w
    rownames(F_Mort) <- species_names

    for(i in 1:num_species) {
        for(j in 1:length(w)) {
            # f = c * J * w
            F_Mort[i, j] <- c_value * J[i, j] * w[j]
        }
    }

    return(list(
        J_adv = J_adv,
        J_diff = J_diff,
        J = J,
        f = F_Mort,
        c_value = c_value
    ))
}


#' Plot Mizer Flux and Balanced Harvest Mortality
#'
#' @description
#' Visualizes the components of the size-spectrum flux. This is particularly
#' useful for verifying that fishing mortality ($f$) correctly follows the
#' "Production" curve of the species, a key requirement for Balanced Harvesting.
#'
#' @param flux_list Output from \code{calculate_mizer_flux}.
#' @param variable The variable to plot:
#' \itemize{
#'   \item \code{"J_adv"}: Growth-based advection.
#'   \item \code{"J_diff"}: Growth-based diffusion.
#'   \item \code{"J"}: Net flux (Production proxy).
#'   \item \code{"f"}: The resulting Balanced Harvesting mortality.
#' }
#'
#' @details
#' In a PssBH (Partial species-size-level Balanced Harvesting) model, one
#' would expect the \code{"f"} plot to show zero mortality for weights
#' below the 10cm threshold, and a curve mirroring the total flux \code{"J"}
#' for larger individuals.
#'
#' @return A ggplot object with facets for each species, allowing for
#' comparison of biomass flow across different life histories.
#'
#' @examples (cannot be run yet because the model is not yet saved)
#'
#' # Calculate fluxes for a Celtic Sea model
#' fluxes <- calculate_mizer_flux(params, c_value = 0.05)
#'
#' # Visualize the Total Flux (Production)
#' plot_mizer_flux(fluxes, variable = "J")
#'
#' # Visualize the resulting Balanced Harvest fishing mortality
#' plot_mizer_flux(fluxes, variable = "f")
#'
#'
#' @export

plot_flux <- function(flux_list, variable = "f") {
    library(ggplot2)
    library(tidyr)
    library(dplyr)

    #Extract the specific matrix
    mat <- flux_list[[variable]]

    if (is.null(mat)) {
        stop("Variable not found in flux list. Choose 'J_adv', 'J_diff', 'J', or 'f'.")
    }

    #Title mapping for clarity
    titles <- c("J_adv" = "Advective Flux",
                "J_diff" = "Diffusive Flux",
                "J" = "Total Flux (J)",
                "f" = paste("Fishing Mortality (c =", flux_list$c_value, ")"))

    df_plot <- as.data.frame(mat)
    df_plot$Species <- rownames(mat)

    df_long <- df_plot %>%
        pivot_longer(cols = -Species,
                     names_to = "Weight",
                     values_to = "Value") %>%
        mutate(Weight = as.numeric(Weight))

    ggplot(df_long, aes(x = Weight, y = Value, color = Species)) +
        geom_line(linewidth = 1) +
        facet_wrap(~Species, scales = "free_y") +
        scale_x_log10() +
        theme_minimal() +
        labs(title = titles[variable],
             x = "Weight (g)",
             y = "Value")
}
