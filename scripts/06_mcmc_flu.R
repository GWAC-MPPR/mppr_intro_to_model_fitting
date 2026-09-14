# MCMC fit of the SIR model to the 1978 boarding-school influenza outbreak, for
# the "Bayesian inference and MCMC" section (06_bayesian_mcmc.qmd). The same
# analysis is worked through step by step in 02_mcmc_practicals.qmd.
#
# Writes:
#   figures/06_flu_prior_predictive.png  epidemics simulated from the prior
#   figures/06_grid_1d.png               the posterior of beta on a grid of 15 values
#   figures/06_grid_2d.png               the posterior on a 15 x 15 grid of (beta, gamma)
#   figures/06_flu_grid_posterior.png    the posterior on a (beta, gamma) grid
#   results/06_grid_cost.txt             grid size and run time against number of parameters
#   figures/06_flu_traces.png            four chains, after burn-in
#   figures/06_flu_posterior.png         marginal posteriors and the joint
#   figures/06_flu_fit.png               posterior trajectories against the data
#   results/06_flu_mcmc_summary.txt      estimates, intervals, diagnostics
#
# The slides pull code from the regions between `# start snippet <name>` and
# `# end snippet <name>` markers (see filters/include-code.lua). Everything
# outside the markers is plumbing that stays off the slides, so keep the code
# inside them exactly as it should appear on the slide.
#
# Run from the project root: Rscript scripts/06_mcmc_flu.R   (about a minute)

source("scripts/00_flu_sir_model.R")   # the data, the SIR model, predict_cases()
library(MASS)     # mvrnorm, for a proposal with a covariance matrix

dir.create("figures", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)
set.seed(1978)

fig <- function(name, width = 10, height = 6) {
  png(file.path("figures", name), width = width, height = height,
      units = "in", res = 150, pointsize = 15)
}

# ---- Prior, likelihood, posterior ------------------------------------------

# start snippet flu_posterior
# Sample theta = (log beta, log gamma): the chain roams freely, the rates stay positive

# Prior: normal on the logs (lognormal on beta and gamma), centred on R0 = 2 and a
# 2-day infectious period, with sd 1 on the log scale (a factor of ~2.7 either way)
log_prior <- function(theta) {
  dnorm(theta[1], mean = log(0.0013), sd = 1, log = TRUE) +
    dnorm(theta[2], mean = log(0.5), sd = 1, log = TRUE)
}

# Likelihood: Poisson counts around the model trajectory, as in the MLE section
log_likelihood <- function(theta) {
  predicted <- predict_cases(beta = exp(theta[1]), gamma = exp(theta[2]))
  if (any(!is.finite(predicted)) || any(predicted < 0)) return(-Inf)  # solver hiccup
  sum(dpois(flu$cases, lambda = predicted, log = TRUE))
}

# Posterior up to a constant: all the Metropolis ratio needs
log_posterior <- function(theta) log_prior(theta) + log_likelihood(theta)
# end snippet flu_posterior

# ---- Prior predictive check ------------------------------------------------

fig("06_flu_prior_predictive.png")
par(mar = c(4.5, 4.5, 3, 1))
plot(flu$day, flu$cases, type = "n", ylim = c(0, 763), xlab = "Day",
     ylab = "Pupils in bed", main = "50 epidemics simulated from the prior")
for (k in 1:50) {
  b <- rnorm(1, log(0.0013), 1)
  g <- rnorm(1, log(0.5), 1)
  lines(1:14, predict_cases(exp(b), exp(g)), col = adjustcolor("steelblue", 0.4))
}
points(flu$day, flu$cases, pch = 16, cex = 1.3)
dev.off()

# ---- The posterior on a grid -----------------------------------------------

grid_posterior <- function(beta_grid, gamma_grid) {
  lp <- outer(beta_grid, gamma_grid,
              Vectorize(function(b, g) log_posterior(c(log(b), log(g)))))
  exp(lp - max(lp))   # relative to the highest point on the grid
}

wide_beta <- seq(0.0005, 0.005, length.out = 60)
wide_gamma <- seq(0.1, 1.0, length.out = 60)
zoom_beta <- seq(0.00205, 0.0024, length.out = 60)
zoom_gamma <- seq(0.42, 0.54, length.out = 60)

fig("06_flu_grid_posterior.png", height = 5)
par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3, 1))
image(wide_beta, wide_gamma, grid_posterior(wide_beta, wide_gamma),
      col = hcl.colors(50, "YlOrRd", rev = TRUE),
      xlab = expression(beta), ylab = expression(gamma), main = "Over the range of the prior")
rect(min(zoom_beta), min(zoom_gamma), max(zoom_beta), max(zoom_gamma), border = "grey30", lty = 2)
image(zoom_beta, zoom_gamma, grid_posterior(zoom_beta, zoom_gamma),
      col = hcl.colors(50, "YlOrRd", rev = TRUE),
      xlab = expression(beta), ylab = expression(gamma), main = "Zoomed in (the dashed box)")
contour(zoom_beta, zoom_gamma, grid_posterior(zoom_beta, zoom_gamma),
        levels = c(0.01, 0.1, 0.5, 0.9), add = TRUE, col = "grey30")
dev.off()

# ---- Building up the grid idea ---------------------------------------------

# One parameter: gamma fixed at its MLE, 15 values of beta
gamma_fixed <- 0.476
beta_15 <- seq(0.00205, 0.0024, length.out = 15)
lp_15 <- sapply(beta_15, function(b) log_posterior(c(log(b), log(gamma_fixed))))
beta_fine <- seq(0.00205, 0.0024, length.out = 300)
lp_fine <- sapply(beta_fine, function(b) log_posterior(c(log(b), log(gamma_fixed))))
post_15 <- exp(lp_15 - max(lp_fine)); post_fine <- exp(lp_fine - max(lp_fine))

fig("06_grid_1d.png")
par(mar = c(4.5, 4.5, 3, 1))
plot(beta_fine, post_fine, type = "l", col = "grey70", lwd = 2, xlab = expression(beta),
     ylab = "Posterior (relative to its peak)", main = expression("15 values of " * beta * ", with " * gamma * " fixed: 15 model runs"))
segments(beta_15, 0, beta_15, post_15, col = "steelblue", lwd = 3)
points(beta_15, post_15, pch = 16, cex = 1.4, col = "steelblue")
points(beta_15, rep(0, 15), pch = 4, cex = 1.2, col = "grey30")
legend("topright", bty = "n", pch = c(4, 16, NA), lwd = c(NA, 3, 2), col = c("grey30", "steelblue", "grey70"),
       legend = c("grid points", "posterior evaluated at each", "what the grid is trying to picture"))
dev.off()

# Two parameters: 15 x 15 grid
beta_g <- seq(0.00205, 0.0024, length.out = 15)
gamma_g <- seq(0.42, 0.54, length.out = 15)
grid2 <- expand.grid(beta = beta_g, gamma = gamma_g)
grid2$post <- exp(mapply(function(b, g) log_posterior(c(log(b), log(g))), grid2$beta, grid2$gamma))
grid2$post <- grid2$post / max(grid2$post)
n_useful <- sum(grid2$post > 0.01)

fig("06_grid_2d.png")
par(mar = c(4.5, 4.5, 3, 1))
plot(grid2$beta, grid2$gamma, pch = 21, cex = 1.9, lwd = 0.5,
     bg = hcl.colors(50, "YlOrRd", rev = TRUE)[pmax(1, ceiling(grid2$post * 50))],
     xlab = expression(beta), ylab = expression(gamma),
     main = "15 values of each: 15 x 15 = 225 model runs")
dev.off()

# How the cost grows: one model run timed, then multiplied up
run_time <- system.time(for (i in 1:200) predict_cases(0.0022, 0.476))[["elapsed"]] / 200
cost <- data.frame(parameters = c(1, 2, 3, 4, 6, 10))
cost$grid_points <- 15^cost$parameters
cost$time <- cost$grid_points * run_time
fmt_time <- function(s) {
  ifelse(s < 60, sprintf("%.1f seconds", s),
  ifelse(s < 3600, sprintf("%.0f minutes", s / 60),
  ifelse(s < 86400, sprintf("%.0f hours", s / 3600),
  ifelse(s < 3.2e7, sprintf("%.0f days", s / 86400), sprintf("%.0f years", s / 3.15e7)))))
}
sink("results/06_grid_cost.txt")
cat(sprintf("One model run takes about %.1f ms on this laptop\n\n", 1000 * run_time))
cat(sprintf("%-11s %18s %16s\n", "parameters", "grid points (15^p)", "time"))
for (i in seq_len(nrow(cost))) {
  cat(sprintf("%-11d %18s %16s\n", cost$parameters[i], format(cost$grid_points[i], big.mark = ","), fmt_time(cost$time[i])))
}
sink()

# ---- Metropolis sampler for several parameters -----------------------------

# start snippet metropolis
metropolis <- function(log_target, n_iter, init, proposal_cov) {
  samples <- matrix(NA, nrow = n_iter, ncol = length(init))
  theta <- init
  log_p <- log_target(theta)
  if (!is.finite(log_p)) stop("the starting value has zero posterior density")
  n_accepted <- 0
  for (i in seq_len(n_iter)) {
    proposal <- theta + MASS::mvrnorm(1, mu = rep(0, length(init)),
                                      Sigma = proposal_cov)
    log_p_proposal <- log_target(proposal)
    if (log(runif(1)) < log_p_proposal - log_p) {
      theta <- proposal
      log_p <- log_p_proposal
      n_accepted <- n_accepted + 1
    }
    samples[i, ] <- theta
  }
  list(samples = samples, acceptance_rate = n_accepted / n_iter)
}
# end snippet metropolis

# start snippet flu_pilot
# Pilot run: a rough, uncorrelated proposal to find the region of high posterior
pilot <- metropolis(log_posterior, n_iter = 3000,
                    init = c(log(0.0013), log(0.5)),
                    proposal_cov = diag(c(0.1, 0.1)^2))

# Tune: scale the covariance of the pilot samples (Roberts & Rosenthal's 2.38^2/d)
tuned_cov <- 2.38^2 / 2 * cov(pilot$samples[-(1:1000), ])
# end snippet flu_pilot

# start snippet flu_chains
# Four chains from starting points spread well beyond the posterior (a factor
# of ~1.6 either side of the pilot centre), so that agreement between them means
# something
centre <- colMeans(pilot$samples[-(1:1000), ])
starts <- list(centre + c(-0.5, -0.5), centre + c(0.5, -0.5),
               centre + c(-0.5, 0.5), centre + c(0.5, 0.5))
chains <- lapply(starts, function(s) metropolis(log_posterior, n_iter = 6000,
                                                 init = s, proposal_cov = tuned_cov))
burn_in <- 1000
posterior_samples <- do.call(rbind, lapply(chains, function(ch) ch$samples[-(1:burn_in), ]))
# end snippet flu_chains

# ---- Diagnostics -------------------------------------------------------------

rhat <- function(chains) {
  n <- length(chains[[1]])
  W <- mean(sapply(chains, var))
  B <- n * var(sapply(chains, mean))
  sqrt(((n - 1) / n * W + B / n) / W)
}

# Effective sample size: n over the integrated autocorrelation time
ess <- function(x) {
  rho <- acf(x, lag.max = min(length(x) - 1, 500), plot = FALSE)$acf[-1]
  k <- which(rho < 0)[1]
  if (!is.na(k)) rho <- rho[seq_len(k - 1)]
  length(x) / (1 + 2 * sum(rho))
}

cols <- c("tomato", "steelblue", "forestgreen", "darkorange")
par_names <- c(expression(beta), expression(gamma))

fig("06_flu_traces.png", height = 6)
par(mfrow = c(2, 1), mar = c(4, 6.5, 2.5, 1))
for (j in 1:2) {
  post_chains <- lapply(chains, function(ch) exp(ch$samples[-(1:burn_in), j]))
  plot(NA, xlim = c(0, 6000 - burn_in), ylim = range(unlist(post_chains)), yaxt = "n",
       xlab = "Iteration after burn-in", ylab = "",
       main = sprintf("R-hat = %.3f", rhat(post_chains)))
  axis(2, at = pretty(range(unlist(post_chains)), 3), las = 1)
  mtext(par_names[[j]], side = 2, line = 5)
  for (k in 1:4) lines(post_chains[[k]], col = adjustcolor(cols[k], 0.7))
}
dev.off()

# ---- Posterior summaries -----------------------------------------------------

beta_post <- exp(posterior_samples[, 1])
gamma_post <- exp(posterior_samples[, 2])
R0_post <- beta_post * N / gamma_post
period_post <- 1 / gamma_post

fig("06_flu_posterior.png", height = 6.5)
par(mfrow = c(2, 2), mar = c(4.5, 4.5, 2.5, 1))
hist(beta_post, breaks = 40, col = "grey85", border = "white", freq = FALSE,
     xlab = expression(beta), main = expression("Posterior of " * beta))
hist(gamma_post, breaks = 40, col = "grey85", border = "white", freq = FALSE,
     xlab = expression(gamma), main = expression("Posterior of " * gamma))
hist(R0_post, breaks = 40, col = "grey85", border = "white", freq = FALSE,
     xlab = expression(R[0] == beta * N / gamma), main = expression("Posterior of " * R[0]))
plot(beta_post, gamma_post, pch = 16, cex = 0.4, col = adjustcolor("steelblue", 0.3),
     xlab = expression(beta), ylab = expression(gamma),
     main = sprintf("Joint posterior (correlation %.2f)", cor(beta_post, gamma_post)))
dev.off()

# ---- Posterior predictive fit ------------------------------------------------

draws <- posterior_samples[sample(nrow(posterior_samples), 500), ]
fine_times <- seq(0, 14, by = 0.1)
trajectories <- sapply(seq_len(nrow(draws)), function(i) {
  out <- ode(y = init, times = fine_times, func = sir_model,
             parms = c(beta = exp(draws[i, 1]), gamma = exp(draws[i, 2])))
  out[, "I"]
})
traj_q <- apply(trajectories, 1, quantile, probs = c(0.025, 0.5, 0.975))
predictive <- sapply(seq_len(nrow(draws)), function(i) {
  rpois(14, predict_cases(exp(draws[i, 1]), exp(draws[i, 2])))
})
pred_q <- apply(predictive, 1, quantile, probs = c(0.025, 0.975))

fig("06_flu_fit.png")
par(mar = c(4.5, 4.5, 3, 1))
plot(flu$day, flu$cases, type = "n", ylim = c(0, 400), xlab = "Day", ylab = "Pupils in bed",
     main = "Posterior fit: 500 draws from the posterior")
polygon(c(1:14, 14:1), c(pred_q[1, ], rev(pred_q[2, ])), col = adjustcolor("tomato", 0.15), border = NA)
polygon(c(fine_times, rev(fine_times)), c(traj_q[1, ], rev(traj_q[3, ])),
        col = adjustcolor("tomato", 0.35), border = NA)
lines(fine_times, traj_q[2, ], col = "tomato", lwd = 3)
points(flu$day, flu$cases, pch = 16, cex = 1.3)
legend("topright", bty = "n", pch = c(16, NA, 15, 15), lwd = c(NA, 3, NA, NA),
       col = c("black", "tomato", adjustcolor("tomato", 0.35), adjustcolor("tomato", 0.15)),
       pt.cex = c(1.3, NA, 2, 2),
       legend = c("data", "posterior median trajectory", "95% interval: trajectory",
                  "95% interval: predicted counts"))
dev.off()

# ---- Summary table -----------------------------------------------------------

summarise <- function(x) c(mean = mean(x), median = median(x),
                           quantile(x, c(0.025, 0.975)))
summary_table <- rbind(
  beta = summarise(beta_post),
  gamma = summarise(gamma_post),
  R0 = summarise(R0_post),
  infectious_period_days = summarise(period_post)
)

sink("results/06_flu_mcmc_summary.txt", split = TRUE)
cat("Posterior summaries: 4 chains x 5,000 draws after burn-in\n\n")
print(noquote(t(apply(signif(summary_table, 3), 1, format))))
cat("\nAcceptance rates:", paste(sprintf("%.0f%%", 100 * sapply(chains, `[[`, "acceptance_rate")), collapse = ", "), "\n")
cat(sprintf("R-hat: beta %.3f, gamma %.3f\n",
            rhat(lapply(chains, function(ch) ch$samples[-(1:burn_in), 1])),
            rhat(lapply(chains, function(ch) ch$samples[-(1:burn_in), 2]))))
cat(sprintf("Effective sample size: beta %.0f, gamma %.0f of %d\n",
            ess(posterior_samples[, 1]), ess(posterior_samples[, 2]), nrow(posterior_samples)))
sink()
