

# ===== Example usage =====
# Suppose:
# - Survey on June 30 -> theta ~ 0.5
# - Annulus becomes countable around May 15 -> alpha ~ 0.37
# - Spawning centered at μ = 0.1 (mid-Feb), concentration κ = 5
# - Observed K = 2 rings
# - Soft edge sigma = 0.03 (about ~11 days on the circle)
set.seed(1)
res <- age_density_given_k(
    k = 0,
    theta = 0.825,
    mu = 0.25, kappa = 1,
    alpha = 0.3,
    sigma = 0.1,
    R_max = 8
)

# A quick plot of the age density:
plot(res$density$a, res$density$dens, type = "l",
     xlab = "Age at survey (years)", ylab = "p(A | K=2)")

# Inspect posterior over completed years R
print(res$post_R)

# Example: posterior mean and 95% interval
dx <- diff(res$density$a)[1]
dens <- res$density$dens
a    <- res$density$a
post_mean <- sum(a * dens) * dx
cdf <- cumsum(dens) * dx
qlo <- approx(cdf, a, xout = 0.025, ties = "ordered")$y
qhi <- approx(cdf, a, xout = 0.975, ties = "ordered")$y
cat(sprintf("E[A|K=2]=%.3f, 95%% CI [%.3f, %.3f]\n", post_mean, qlo, qhi))
