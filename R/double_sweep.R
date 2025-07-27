#' Helper function for double sweep method (Thomas algorithm)
#'
#' This function solves a tridiagonal system of equations of the form
#' $$Ax = b$$
#' where the matrix \eqn{A} has the entries \eqn{D_1, D_2, \dots, D_N} on the
#' main diagonal, \eqn{L_2, L_3, \dots, L_N} on the lower diagonal, and
#' \eqn{U_1, U_2, \dots, U_{N-1}} on the upper diagonal.
#'
#' For more details on the double sweep method, see the vignette
#' \code{vignette("diffusion")}
#'
#' @param U Upper diagonal entries. Even though the matrix only uses
#'   \eqn{U_1\dots U_{N-1}}, you need to supply a full vector of length \eqn{N}.
#'   The last element will be ignored.
#' @param L Lower diagonal entries. Even though the matrix only uses
#'   \eqn{L_2\dots L_N}, you need to supply a full vector of length \eqn{N}.
#'   The first element will be ignored.
#' @param D Main diagonal entries
#' @param b Right-hand side vector
#' @return A vector of solutions
#' @examples
#' U <- c(1, 2, 3, NA)
#' L <- c(NA, 4, 5, 6)
#' D <- c(7, 8, 9, 10)
#' b <- c(11, 12, 13, 14)
#' solve_double_sweep(U, L, D, b)
#' @export
#'
solve_double_sweep <- function(U, L, D, b) {
    # U: upper diagonal (length N-1)
    # L: lower diagonal (length N-1)
    # D: main diagonal (length N)
    # b: right-hand side (length N)
    N <- length(D)
    if (length(U) != N || length(L) != N || length(b) != N) {
        stop("U, L, D, and b must have the same length.")
    }
    # Initialize arrays for alpha and beta coefficients
    alpha <- numeric(N)
    beta <- numeric(N)

    # Initial condition
    alpha[2] <- -U[1]/D[1]
    beta[2] <- b[1]/D[1]

    # Forward sweep - calculate alpha and beta coefficients
    for (i in 2:(N-1)) {
        denom <- alpha[i] * L[i] + D[i]
        alpha[i + 1] <- -U[i] / denom
        beta[i + 1] <- -(beta[i] * L[i] - b[i]) / denom
    }

    # Initialize n array for solution
    x <- numeric(N)
    # Set last value of x using the last beta and alpha
    x[N] <- (b[N] - beta[N] * L[N]) / (alpha[N] * L[N] + D[N])

    # Work backwards to get all x values
    for (i in (N-1):1) {
        x[i] <- alpha[i + 1] * x[i + 1] + beta[i + 1]
    }
    return(x)
}
