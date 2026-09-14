# Maximum likelihood fit of the SIR model to the 1978 influenza outbreak
# (04_maximum_likelihood.qmd), with a Poisson and a negative binomial
# observation model.
#
# Writes:
#   figures/04_profile_likelihood.png   profile likelihoods
#   figures/04_mle_fit.png              MLE and least squares fits against the data
#   results/04_mle_estimates.txt        MLEs with standard errors and R0
#   results/04_mle_confint.txt          profile likelihood confidence intervals
#   results/04_mle_vs_ls.txt            MLE against the least squares estimates
#   results/04_model_comparison.txt     AIC of the two observation models
#
# The slides pull code from the regions between `# start snippet <name>` and
# `# end snippet <name>` markers (see filters/include-code.lua). Everything
# outside the markers is plumbing that stays off the slides, so keep the code
# inside them exactly as it should appear on the slide.
#
# Run from the project root: Rscript scripts/04_maximum_likelihood.R

# Builds on the least squares fit: the data, the model and ls_estimates
source("scripts/03_least_squares.R")

fig <- function(name, width = 10, height = 6) {
  png(file.path("figures", name), width = width, height = height,
      units = "in", res = 150, pointsize = 15)
}

# start snippet mle_implementation
# Negative log-likelihood, with the parameters on the log scale: the optimiser
# can then move freely while beta and gamma stay positive, and the two
# parameters are on a comparable scale
nll_function <- function(log_beta, log_gamma) {
  predicted <- predict_cases(beta = exp(log_beta), gamma = exp(log_gamma))
  -sum(dpois(flu$cases, lambda = predicted, log = TRUE))
}
# end snippet mle_implementation

# start snippet mle_fitting
library(bbmle)

# Start from the least squares estimates
fit_mle <- mle2(nll_function,
                start = list(log_beta = log(ls_estimates[["beta"]]),
                             log_gamma = log(ls_estimates[["gamma"]])))

mle_estimates <- exp(coef(fit_mle))              # back to the natural scale
names(mle_estimates) <- c("beta", "gamma")
mle_se <- sqrt(diag(vcov(fit_mle)))              # standard errors, log scale
# end snippet mle_fitting

sink("results/04_mle_estimates.txt", split = TRUE)
cat("beta  =", signif(mle_estimates[["beta"]], 3),
    " (log-scale SE", signif(mle_se[["log_beta"]], 2), ")\n")
cat("gamma =", signif(mle_estimates[["gamma"]], 3),
    " (log-scale SE", signif(mle_se[["log_gamma"]], 2), ")\n")
cat("R0 = beta N / gamma =", round(mle_estimates[["beta"]] * N / mle_estimates[["gamma"]], 2), "\n")
sink()

fig("04_profile_likelihood.png")
# start snippet mle_uncertainty
# Profile likelihood for each parameter
prof <- profile(fit_mle)
plot(prof, absVal = TRUE)
# end snippet mle_uncertainty
dev.off()

sink("results/04_mle_confint.txt", split = TRUE)
# start snippet mle_ci
# 95% confidence intervals from the profiles, back on the natural scale
ci <- exp(confint(fit_mle, level = 0.95))
rownames(ci) <- c("beta", "gamma")
print(ci)
# end snippet mle_ci
sink()

sink("results/04_mle_vs_ls.txt", split = TRUE)
cat("               beta    gamma    R0\n")
cat(sprintf("Least squares  %.5f  %.3f   %.2f\n", ls_estimates[["beta"]], ls_estimates[["gamma"]],
            ls_estimates[["beta"]] * N / ls_estimates[["gamma"]]))
cat(sprintf("MLE (Poisson)  %.5f  %.3f   %.2f\n", mle_estimates[["beta"]], mle_estimates[["gamma"]],
            mle_estimates[["beta"]] * N / mle_estimates[["gamma"]]))
sink()

fine_times <- seq(0, 14, by = 0.1)
mle_curve <- ode(y = init, times = fine_times, func = sir_model,
                 parms = c(beta = mle_estimates[["beta"]], gamma = mle_estimates[["gamma"]]))[, "I"]

fig("04_mle_fit.png")
par(mar = c(4.5, 4.5, 3, 1))
plot(flu$day, flu$cases, pch = 16, cex = 1.3, ylim = c(0, 320),
     xlab = "Day", ylab = "Pupils in bed", main = "Least squares and maximum likelihood fits")
lines(fine_times, ls_curve, col = "tomato", lwd = 3, lty = 2)
lines(fine_times, mle_curve, col = "steelblue", lwd = 3)
legend("topright", bty = "n", pch = c(16, NA, NA), lwd = c(NA, 3, 3), lty = c(NA, 2, 1),
       col = c("black", "tomato", "steelblue"), legend = c("data", "least squares", "MLE, Poisson"))
dev.off()

# start snippet model_comparison
# A second observation model: negative binomial counts, which allow more
# spread around the trajectory than Poisson (extra parameter phi)
nll_nb <- function(log_beta, log_gamma, log_phi) {
  predicted <- predict_cases(beta = exp(log_beta), gamma = exp(log_gamma))
  -sum(dnbinom(flu$cases, mu = predicted, size = exp(log_phi), log = TRUE))
}
# end snippet model_comparison

# start snippet model_fit
fit_nb <- mle2(nll_nb,
               start = list(log_beta = coef(fit_mle)[["log_beta"]],
                            log_gamma = coef(fit_mle)[["log_gamma"]],
                            log_phi = log(10)))
# end snippet model_fit

sink("results/04_model_comparison.txt", split = TRUE)
cat("Poisson AIC:          ", round(AIC(fit_mle), 1), "\n")
cat("Negative binomial AIC:", round(AIC(fit_nb), 1), "\n")
cat("Delta AIC:            ", round(AIC(fit_nb) - AIC(fit_mle), 1), "\n")
sink()
