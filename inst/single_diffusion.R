library(dplyr)
library(mizerEcopath)
library(ggplot2)
# library(deSolve)
library(bvpSolve)
source("inst/newParams.R")
params <- newParams(no_w = 500, sigma = 2, lambda = 2.2, ks = 0)
params <- setPredKernel(params, pred_kernel = getPredKernel(params))
#params <- setPredKernel(NS_params, pred_kernel = getPredKernel(NS_params))

dm <- getDiffusion(params)
#gm <- getEGrowth(params)
gm <- (getDiffusion(params, order = 1) - metab(params)) * (1 - params@psi)
mum <- getMort(params)
names(dimnames(mum)) <- c("sp", "w")
w <- w(params)

df <- melt(dm)
df$type <- "diffusion"
dg <- melt(gm)
dg$type <- "growth"
dmu <- melt(mum)
dmu$type <- "mortality"
da <- rbind(df, dg, dmu)
ggplot(da, aes(x = w, y = value, colour = type)) +
    geom_line() +
    scale_x_log10() +
    scale_y_log10() +
    facet_wrap(~sp, scales = "free")

# Start near 1g
sel <- w > 1
# Choose one species
sp <- 1
d <- dm[sp, sel]
g <- gm[sp, sel]
mu <- mum[sp, sel]
w <- w[sel]
N.guess <- params@initial_n[sp, sel]
w0 <- w[1]
wMax <- w[length(w)]

# We create spline fits for d, g, mu so we can get
# d'(w), d''(w), etc. from the spline. This might be easiest to keep
# everything smooth if your data is on a well-defined grid.

d.spline  <- splinefun(w, d,  method = "natural")
g.spline  <- splinefun(w, g,  method = "natural")
mu.spline <- splinefun(w, mu, method = "natural")
N.guess.spline <- splinefun(w, N.guess, method = "natural")

# Then we can also get d'(w) by d.spline'(w), d''(w) by d.spline''(w), etc.

# ODE function for bvpSolve
odefun <- function(w, y, pars) {
    N  <- y[1]
    Np <- y[2]  # derivative

    d   <- d.spline(w)
    dp  <- d.spline(w, deriv = 1)
    ddp <- d.spline(w, deriv = 2)

    g   <- g.spline(w)
    gp  <- g.spline(w, deriv = 1)

    mu  <- mu.spline(w)

    # The second derivative from the expanded ODE:
    Npp <- (  2*gp*N + 2*g*Np + 2*mu*N
              - ddp*N   - 2*dp*Np ) / d

    list( c(Np, Npp) )  # Return a list with derivatives
}

# define an initial guess function
guessfun <- function(x, pars) {
    m <- rbind(N.guess.spline(x), N.guess.spline(x, deriv = 1))
    m[m<0] <- 0
    return(m)
}

# Positivity constraint for N(w)
posbound <- function(x, y, pars) {
    y[1]  # Constrain the first dependent variable N(w) to be >= 0
}

N0 <- N.guess[1]
sol <- bvptwp(yini = c(N = N0, Np = NA),
              yend = c(N = 0, Np = NA),
              x = w,
              func = odefun,
              posbound = posbound,
              xguess = w,
              yguess = guessfun(w),
              islin = TRUE,
              allpoints = FALSE,
              atol = 1)

wn <- sol[, 1]
Nn <- sol[, 2]
Npn <- sol[, 3]

plot(wn, Nn * wn, type = "l", col = "blue", lwd = 2,
     xlab = "w", ylab = "B(w)")
N <- N.guess
lines(wn, N * wn, col = "red", lwd = 2)


# Using bvpcol did not work

bcfun <- function(ya, yb, pars) {
    bc1 <- ya[1] - N0      # N(w0) - N0 = 0
    bc2 <- yb[1]           # N(wMax) = 0
    return( c(bc1, bc2) )
}

sol <- bvpcol(
    x        = c(w0, wMax),   # domain
    func     = odefun,
    bound    = bcfun,
    leftbc   = 1,
    ncomp    = 2
)

# The returned object sol has the numerical solution on the mesh
# used by bvpcol. You can evaluate it or plot it at specific points:
Nout <- bvpsol(sol, x = w)[,1]  # first component = N


# juvenile exponent
n <- params@species_params[sp, "n"]
g0 <- g[[1]] / w[[1]]^n
mu0 <- mu[[1]] / w[[1]]^(n - 1)
d0 <- d[[1]] / w[[1]]^(n + 1)
a <- -n + (g0 - d0 / 2 - sqrt((g0 - d0 / 2)^2 + 2 * d0 * mu0)) / d0
a_old <- -n - mu0 / g0

N0 <- 1
Np0 <- a / w[[1]] * N0

##############################
# 2) Build approximate derivatives for d, g
##############################
n <- length(w)
dprime   <- numeric(n)
ddprime  <- numeric(n)
gprime   <- numeric(n)

# Central difference for interior
for (i in 2:(n-1)) {
    dprime[i]  <- (d[i+1] - d[i-1]) / (w[i+1] - w[i-1])
    gprime[i]  <- (g[i+1] - g[i-1]) / (w[i+1] - w[i-1])
}
# One-sided at boundaries
dprime[1] <- (d[2] - d[1]) / (w[2] - w[1])
dprime[n] <- (d[n] - d[n-1]) / (w[n] - w[n-1])

gprime[1] <- (g[2] - g[1]) / (w[2] - w[1])
gprime[n] <- (g[n] - g[n-1]) / (w[n] - w[n-1])

# Second derivative of d
for (i in 2:(n-1)) {
    ddprime[i] <- (dprime[i+1] - dprime[i-1]) / (w[i+1] - w[i-1])
}
ddprime[1] <- ddprime[2]
ddprime[n] <- ddprime[n-1]

##############################
# 3) Define a2, a1, a0 on the grid
##############################
a2 <- 0.5 * d
a1 <- dprime - g
a0 <- 0.5 * ddprime - gprime - mu
plot(w, a0,
     ylim = c(min(a0, a1, a2), max(a0, a1, a2)),
     type="l", col="blue", lwd=2,
     ylab="value", xlab="w")
lines(w, a1, col="red", lwd=2)
lines(w, a2, col="green", lwd=2)
legend("bottomright", legend = c("a0", "a1", "a2"),
       col = c("blue", "red", "green"), lwd = 2,
       cex = 0.5, bty = "o")

plot(w, a2, type="l", col="green", lwd=2,
     ylab="value", xlab="w", log = "y")

plot(w, a0 / a2, type="l", col="blue", lwd=2,
     ylab="ratios", xlab="w", log = "x")
lines(w, a1 / a2, col="red", lwd=2)
min(a0 / a2)
max(a0 / a2)
max(a1 / a2)
min(a1 / a2)

##############################
# 4) Interpolate them
##############################
a1fun <- approxfun(w, a1 / a2, rule=2)
a0fun <- approxfun(w, a0 / a2, rule=2)

##############################
# 5) ODE system function
##############################
odeSystem <- function(s, state, parms) {
    N  <- state[1]
    Np <- state[2]

    A1 <- a1fun(s)
    A0 <- a0fun(s)

    dN   <- Np
    dNp  <- -(A1 * Np + A0 * N)

    return(list(c(dN, dNp)))
}

##############################
# 6) Solve with initial conditions
##############################
y0  <- c(N = N0, Np = -1.213)       # [ N(w1), N'(w1) ]
out <- ode(y = y0, times = w, func = odeSystem, parms = NULL,
           method = "rk4")

##############################
# 7) Extract solution
##############################
N_sol  <- out[, "N"]
Np_sol <- out[, "Np"]

##############################
# 8) Plot
##############################
plot(w, N_sol, type="l", col="blue", lwd=2,
     xlab="w", ylab="Solution", main="N(w) and N'(w)",
     log = "xy")
lines(w, Np_sol, col="red", lwd=2)
legend("topright", legend=c("N(w)", "N'(w)"),
       col=c("blue","red"), lwd=2)

