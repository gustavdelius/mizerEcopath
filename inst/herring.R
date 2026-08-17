SD <- "/tmp/claude-1000/-home-gustav-Git-mizerEcopath/d6030879-de38-4344-8d0d-bc47513c7f29/scratchpad/vb"
suppressMessages({library(mizerEcopath); library(mizerExperimental); library(dplyr)})
source(file.path(SD, "tune2.R")); source(file.path(SD, "knobs.R"))
MULT <- exp(seq(log(0.15), log(2.6), length.out = 7))
yield_curve <- function(params, s, F_t, mult = MULT) {
    y <- suppressWarnings(suppressMessages(
        getYieldVsF(params, s, gear = "commercial", F_range = F_t * mult,
                    tol = 0.01, t_max = 60)))
    y[order(y$F), ]
}
gpeak <- function(params, s, F_t) peak_ratio(yield_curve(params, s, F_t), F_t)
p <- readRDS(file.path(SD, "C_final.rds")); Ft <- F_targets(p)
for (rl in c(0.02, 0.05, 0.10, 0.15)) {
    r <- gpeak(set_rl(p, "Herring", rl), "Herring", Ft[["Herring"]])
    cat(sprintf("Herring rl %.3f  peak/target %.2f %s\n", rl, r[["ratio"]],
                ifelse(r[["edge"]]==1, "(edge)", "")))
}
