library(dplyr)
library(mizerEcopath)

params <- setPredKernel(NS_params, pred_kernel = getPredKernel(NS_params))
dm <- getDiffusion(params)
gm <- getEGrowth(params)
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

# Choose one species
sp <- "Cod"
d <- dm[sp, ]
g <- gm[sp, ]
mu <- mum[sp, ]

##############################
# 1) Load library
##############################
library(deSolve)

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

##############################
# 4) Interpolate them
##############################
a2fun <- approxfun(w, a2, rule=2)
a1fun <- approxfun(w, a1, rule=2)
a0fun <- approxfun(w, a0, rule=2)

##############################
# 5) ODE system function
##############################
odeSystem <- function(s, state, parms) {
    N  <- state[1]
    Np <- state[2]

    A2 <- a2fun(s)
    A1 <- a1fun(s)
    A0 <- a0fun(s)

    dN   <- Np
    dNp  <- - (A1*Np + A0*N) / A2

    return(list(c(dN, dNp)))
}

##############################
# 6) Solve with initial conditions
##############################
y0  <- c(N0, diffN0)         # [ N(w1), N'(w1) ]
out <- ode(y = y0, times = w, func = odeSystem, parms = NULL,
           method = "rk4")

##############################
# 7) Extract solution
##############################
N_sol  <- out[, "y1"]
Np_sol <- out[, "y2"]

##############################
# 8) Plot
##############################
plot(w, N_sol, type="l", col="blue", lwd=2,
     xlab="w", ylab="Solution", main="N(w) and N'(w)")
lines(w, Np_sol, col="red", lwd=2)
legend("topright", legend=c("N(w)", "N'(w)"),
       col=c("blue","red"), lwd=2)

