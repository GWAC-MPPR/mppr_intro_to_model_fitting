# Maximum likelihood example for 04_maximum_likelihood.qmd.
#
# Fits the synthetic SIR data from the least squares example by maximum
# likelihood with a Poisson and a negative binomial observation model, and
# writes:
#   figures/04_profile_likelihood.png   profile likelihood plot
#   results/04_mle_estimates.txt        MLEs with standard errors and R0
#   results/04_mle_confint.txt          profile likelihood confidence intervals
#   results/04_mle_vs_true.txt          MLEs against the true values
#   results/04_model_comparison.txt     AIC of the two observation models
#
# The slides pull code from the regions between `# start snippet <name>` and
# `# end snippet <name>` markers (see _extensions/quarto-ext/include-code-files).
# Everything outside the markers is plumbing that stays off the slides, so keep
# the code inside them exactly as it should appear on the slide.
#
# Run from the project root: Rscript scripts/04_maximum_likelihood.R

# The MLE example fits the same data as the least squares example, so build on
# it: this defines sir_model, init_conds, true_params and data.
source("scripts/03_least_squares.R")

# start snippet mle_implementation
# Define negative log-likelihood function
nll_function <- function(beta, gamma, data) {
  # Simulate model
  out <- ode(y = init_conds, times = data$time,
             func = sir_model, parms = c(beta = beta, gamma = gamma))

  # Model predictions (scaled to cases)
  predicted <- out[,"I"] * 1000

  # Poisson negative log-likelihood
  nll <- -sum(dpois(data$observed, lambda = predicted, log = TRUE))

  return(nll)
}
# end snippet mle_implementation

# start snippet mle_fitting
# Use optimization to find MLE
library(bbmle)
fit_mle <- mle2(nll_function,
                start = list(beta = 0.2, gamma = 0.1),
                data = list(data = data),
                method = "L-BFGS-B",
                lower = c(0.01, 0.01),
                upper = c(1.0, 0.5))

# Extract results
mle_params <- coef(fit_mle)
mle_se <- sqrt(diag(vcov(fit_mle)))
# end snippet mle_fitting

# Console output between the sink() calls also goes to the results file, which
# the slides include verbatim.
sink("results/04_mle_estimates.txt", split = TRUE)
cat("Beta:", round(mle_params[1], 3), "±", round(mle_se[1], 3), "\n")
cat("Gamma:", round(mle_params[2], 3), "±", round(mle_se[2], 3), "\n")
cat("R0:", round(mle_params[1]/mle_params[2], 2), "\n")
sink()

png("figures/04_profile_likelihood.png", width = 10, height = 6, units = "in", res = 150)
# start snippet mle_uncertainty
# Profile likelihood for uncertainty
prof <- profile(fit_mle)
plot(prof, absVal = TRUE, main = "Profile Likelihood")
# end snippet mle_uncertainty
dev.off()

sink("results/04_mle_confint.txt", split = TRUE)
# start snippet mle_ci
# Confidence intervals
print(confint(fit_mle, level = 0.95))
# end snippet mle_ci
sink()

sink("results/04_mle_vs_true.txt", split = TRUE)
# Compare with true values
cat("True Beta:", true_params[1], "\n")
cat("MLE Beta:", round(mle_params[1], 3), "\n")
cat("True Gamma:", true_params[2], "\n")
cat("MLE Gamma:", round(mle_params[2], 3), "\n")
sink()

# start snippet model_comparison
# Fit different models and compare
# Model 1: SIR with Poisson
# Model 2: SIR with negative binomial

# Negative binomial likelihood
nll_nb <- function(beta, gamma, phi, data) {
  out <- ode(y = init_conds, times = data$time,
             func = sir_model, parms = c(beta = beta, gamma = gamma))

  # Extract and scale predictions to cases
  predicted <- out[,"I"] * 1000

  # Negative Binomial negative log-likelihood
  nll <- -sum(dnbinom(data$observed, mu = predicted, size = phi, log = TRUE))

  return(nll)
}
# end snippet model_comparison

# start snippet model_fit
# Fit NB model
fit_nb <- mle2(nll_nb,
               start = list(beta = 0.2, gamma = 0.1, phi = 10),
               data = list(data = data),
               method = "L-BFGS-B",
               lower = c(0.01, 0.01, 0.1),
               upper = c(1.0, 0.5, 100))
# end snippet model_fit

sink("results/04_model_comparison.txt", split = TRUE)
# Compare models
cat("Poisson AIC:", AIC(fit_mle), "\n")
cat("Negative Binomial AIC:", AIC(fit_nb), "\n")
cat("Delta AIC:", AIC(fit_nb) - AIC(fit_mle), "\n")
sink()
