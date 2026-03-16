#' Create a Blended Fishing Mortality Function for Balanced Harvesting Transitions
#'
#' @description
#' This "Function Factory" returns a customized mortality function for use within
#' a \code{mizer} simulation. It creates a linear transition from
#' traditional gear-based fishing mortality to species-size-level Balanced
#' Harvesting (ssBH) based on ecosystem production flux.
#'
#' @details
#' ## How it Works
#' The function uses a blending parameter \eqn{\alpha} that is tied to the
#' simulation time (\eqn{t}).
#' \itemize{
#'   \item At \eqn{t = 0}, \eqn{\alpha = 0} and fishing is 100% determined by
#'   standard gear selectivity (\code{f_standard}).
#'   \item As \eqn{t} approaches \eqn{t_{max\_sim}}, \eqn{\alpha} moves toward
#'   \code{alpha_max}, phasing in the Balanced Harvest mortality (\code{f_ssBH}).
#'   \item The transition is calculated as: \eqn{F = (1-\alpha)f_{standard} + \alpha f_{ssBH}}.
#' }
#'
#' ## Implementation Note
#' This function must be registered using \code{setRateFunction()} before
#' running \code{project()}. Because it is a closure, it "captures" the
#' \code{t_max_blend}, \code{target_c}, and \code{alpha_max} values.
#'
#' @param t_max_blend Numeric. The total duration of the simulation (in years).
#' This defines the "speed" of the transition.
#' @param target_c Numeric. The proportionality constant for the Balanced
#' Harvest rule (\eqn{F = c \cdot Production}).
#' @param t_steady Numeric. The time it takes for the simulation to reach steady
#'  state with normal fishing employed
#' @param alpha_max Numeric (0 to 1). The maximum extent of the transition.
#' 1.0 represents a total shift to ssBH; values < 1.0 maintain a hybrid regime.
#'
#' @return A function compatible with \code{mizer}'s \code{FMort} rate slot.
#' @export
#'
#' @examples
#' # Define the transition: 20 years with a ssBH harvest intensity of 0.2
#' # Only 70% percent of fishing will be imparted by ssBH harvest, 30% will
#' # remain from current harvest practice.
#' blender <- make_blended_ssBH_FMort(t_max_blend = 20, target_c = 0.2,
#' alpha_max = 0.7)
#'
#' # Register the function
#' params <- setRateFunction(params, "FMort", blender)
#'
#' # Run the simulation
#' sim <- project(params, t_max = 20, effort = 1)
make_blended_ssBH_FMort <- function(t_max_blend, target_c, alpha_max, t_steady){
    function(params, n, n_pp, n_other, t, effort, ...) {
        # Calculate the 'Standard' Mizer Mortality (Gear-based)
        # We assume a base effort of 1 for the underlying gears for this calculation
        f_standard <- mizerFMort(params, n = n, n_pp = n_pp, n_other = n_other,
                                 t = t, effort = 1, ...)

        # Calculate the 'ssBH' Mortality
        # We use a fixed c_value (e.g., 0.2) for the target BH intensity
        flux_results <- compute_flux(params,
                                     c_value = target_c,
                                     n=n, n_pp=n_pp, n_other=n_other)
        f_ssBH <- flux_results$f

        #identify species if negative mortality for species
        neg_species <- rownames(f_ssBH)[apply(f_ssBH < 0, 1, any)]
        if(length(neg_species) > 0) {
            message("Negative fishing mortality for species: ",
                    paste(neg_species, collapse = ", "))
        }

        #set neg mortality to 0
        f_ssBH[f_ssBH < 0] <- 0

        # 3. Blend them based on the 'effort' parameter passed to project()
        # 'effort' here acts as our transition alpha (0 to 1)
        alpha <- alpha_max*((t-t_steady)/(t_max_blend-t_steady))
        alpha <- max(0, min(1, alpha)) #safety clamp

        if(t>t_steady){
            f_combined <- (1-alpha) * f_standard + (alpha) * f_ssBH
        }
        else {
            f_combined <- f_standard
        }

        return(f_combined)
    }
}
