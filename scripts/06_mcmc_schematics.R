# Schematics for the "Bayesian inference and MCMC" section
# (06_bayesian_mcmc.qmd). Every target here is a toy distribution we know the
# answer for, so the pictures isolate one idea each:
#
#   figures/06_bayes_update.png          prior x likelihood -> posterior
#   figures/06_monte_carlo.png           a distribution from 10, 100, 10000 draws
#   figures/06_metropolis_step.png       one Metropolis step: uphill / downhill
#   figures/06_metropolis_flowchart.png  the algorithm as a loop
#   figures/06_metropolis_normal.png     the sampler on Normal(5, 2): trace + histogram
#   figures/06_proposal_tuning.png       proposal too narrow / about right / too wide
#   figures/06_burn_in.png               chains started at the mode and far from it
#   figures/06_convergence.png           four chains that mix, and four that do not
#
# The slides pull code from the regions between `# start snippet <name>` and
# `# end snippet <name>` markers (see filters/include-code.lua). Everything
# outside the markers is plumbing that stays off the slides, so keep the code
# inside them exactly as it should appear on the slide.
#
# Run from the project root: Rscript scripts/06_mcmc_schematics.R

dir.create("figures", showWarnings = FALSE)
set.seed(2026)

fig <- function(name, width = 10, height = 6) {
  png(file.path("figures", name), width = width, height = height,
      units = "in", res = 150, pointsize = 15)
}

# ---- Bayes' rule: prior x likelihood -> posterior -------------------------

x <- seq(0, 8, by = 0.01)
prior <- dlnorm(x, meanlog = log(2.5), sdlog = 0.6)
likelihood <- dnorm(x, mean = 4.2, sd = 0.45)
posterior <- prior * likelihood
posterior <- posterior / sum(posterior * 0.01)   # normalise so it is a density

fig("06_bayes_update.png")
par(mar = c(4.5, 4.5, 3, 1))
plot(x, posterior, type = "n", ylim = c(0, 1.05 * max(posterior)),
     xlab = expression(R[0]), ylab = "Density",
     main = expression(posterior %prop% likelihood %*% prior))
lines(x, prior, col = "steelblue", lwd = 3, lty = 2)
lines(x, likelihood / sum(likelihood * 0.01), col = "darkorange", lwd = 3, lty = 3)
lines(x, posterior, col = "tomato", lwd = 4)
legend("topleft", bty = "n", lwd = c(3, 3, 4), lty = c(2, 3, 1),
       col = c("steelblue", "darkorange", "tomato"),
       legend = c("prior: what we believed before the data",
                  "likelihood: what the data say",
                  "posterior: what we believe now"))
dev.off()

# ---- Monte Carlo: learn a distribution from samples -----------------------

fig("06_monte_carlo.png", height = 4.5)
par(mfrow = c(1, 3), mar = c(4.5, 4, 3, 1))
for (n in c(10, 100, 10000)) {
  s <- rnorm(n, mean = 5, sd = 2)
  hist(s, freq = FALSE, breaks = seq(-5, 15, by = 0.5), col = "grey85",
       border = "white", xlim = c(-2, 12), ylim = c(0, 0.3),
       xlab = expression(theta), main = paste(n, "samples"))
  curve(dnorm(x, 5, 2), add = TRUE, col = "tomato", lwd = 3)
  legend("topright", bty = "n",
         legend = c(sprintf("mean = %.2f", mean(s)), sprintf("sd = %.2f", sd(s))))
}
dev.off()

# ---- One Metropolis step ---------------------------------------------------

target <- function(x) 0.6 * dnorm(x, 3, 1) + 0.4 * dnorm(x, 6.5, 0.8)
x_now <- 4.3          # where the chain is
x_up <- 3.4           # a proposal with higher density: always accepted
x_down <- 5.0         # a proposal with lower density: accepted with probability r
r_down <- target(x_down) / target(x_now)

fig("06_metropolis_step.png")
par(mar = c(4.5, 4.5, 1, 1))
curve(target(x), 0, 9, n = 500, lwd = 3, col = "grey30",
      xlab = expression(theta), ylab = expression(p(theta ~ "|" ~ data)),
      ylim = c(0, 0.3))
xx <- seq(0, 9, by = 0.01)
lines(xx, 0.1 * dnorm(xx, x_now, 0.8), lty = 2, col = "steelblue", lwd = 2)
text(5.4, 0.012, expression(proposal ~ q(theta^"*" ~ "|" ~ theta)), col = "steelblue", pos = 4)
for (p in c(x_now, x_up, x_down)) segments(p, 0, p, target(p), lty = 3, col = "grey50")
points(c(x_now, x_up, x_down), target(c(x_now, x_up, x_down)),
       pch = 21, cex = 2, lwd = 2,
       bg = c("grey30", "forestgreen", "darkorange"))
text(x_now, target(x_now), expression(theta ~ "(current)"), pos = 3, offset = 0.8)
text(x_up, target(x_up), expression(theta[1]^"*"), pos = 2, col = "forestgreen")
text(x_down, target(x_down), expression(theta[2]^"*"), pos = 4, col = "darkorange")
arrows(x_now, 0.02, x_up, 0.02, length = 0.1, lwd = 2, col = "forestgreen")
arrows(x_now, 0.04, x_down, 0.04, length = 0.1, lwd = 2, col = "darkorange")
text(0.2, 0.28, adj = 0, col = "forestgreen",
     expression("uphill: " * p(theta[1]^"*") > p(theta) * "  →  accept, r = 1"))
text(0.2, 0.255, adj = 0, col = "darkorange",
     substitute("downhill: accept with probability r = " * p(theta[2]^"*") / p(theta) == v,
                list(v = round(r_down, 2))))
dev.off()

# ---- The algorithm as a flowchart -----------------------------------------

box <- function(x, y, label, w = 3.4, h = 0.9, fill = "white", col = "black") {
  rect(x - w / 2, y - h / 2, x + w / 2, y + h / 2, col = fill, border = col, lwd = 2)
  text(x, y, label, cex = 0.95)
}

fig("06_metropolis_flowchart.png", width = 9, height = 6.5)
par(mar = c(0, 0, 0, 0))
plot.new()
plot.window(xlim = c(0, 10), ylim = c(0, 10))
box(5, 9.3, expression("start somewhere: " * theta^(0)), fill = "grey92")
box(5, 7.7, expression("propose " * theta^"*" == theta + epsilon * ",  " * epsilon %~% N(0, sigma^2)), w = 5)
box(5, 6.1, expression("r = " * p(theta^"*" ~ "|" ~ y) / p(theta ~ "|" ~ y)), w = 3.4)
polygon(c(5, 6.6, 5, 3.4), c(5.2, 4.3, 3.4, 4.3), col = "lightyellow", border = "black", lwd = 2)
text(5, 4.3, expression(u %~% U(0, 1) * ";  " * u < r * " ?"), cex = 0.95)
box(2.3, 2.4, expression("accept: " * theta %<-% theta^"*"), w = 3, fill = "#dff0d8")
box(7.7, 2.4, expression("reject: keep " * theta), w = 3, fill = "#f9e0d0")
box(5, 0.8, expression("record " * theta * " and go again"), w = 4.4, fill = "grey92")
arrows(5, 8.85, 5, 8.15, length = 0.12, lwd = 2)
arrows(5, 7.25, 5, 6.55, length = 0.12, lwd = 2)
arrows(5, 5.65, 5, 5.2, length = 0.12, lwd = 2)
arrows(3.4, 4.3, 2.3, 4.3, length = 0, lwd = 2); arrows(2.3, 4.3, 2.3, 2.85, length = 0.12, lwd = 2)
arrows(6.6, 4.3, 7.7, 4.3, length = 0, lwd = 2); arrows(7.7, 4.3, 7.7, 2.85, length = 0.12, lwd = 2)
text(2.75, 4.55, "yes", cex = 0.9); text(7.25, 4.55, "no", cex = 0.9)
arrows(2.3, 1.95, 2.3, 0.8, length = 0, lwd = 2); arrows(2.3, 0.8, 2.8, 0.8, length = 0.12, lwd = 2)
arrows(7.7, 1.95, 7.7, 0.8, length = 0, lwd = 2); arrows(7.7, 0.8, 7.2, 0.8, length = 0.12, lwd = 2)
lines(c(5, 5, 9.75, 9.75, 7.6), c(0.35, 0.2, 0.2, 7.7, 7.7), lwd = 2, lty = 2)
arrows(7.7, 7.7, 7.5, 7.7, length = 0.12, lwd = 2, lty = 2)
text(9.5, 4.2, "repeat many times", srt = 90, cex = 0.9)
dev.off()

# ---- The sampler ------------------------------------------------------------

# start snippet metropolis_1d
metropolis_1d <- function(log_target, n_iter, init, proposal_sd) {
  samples <- numeric(n_iter)
  theta <- init
  log_p <- log_target(theta)
  n_accepted <- 0
  for (i in seq_len(n_iter)) {
    # 1. propose a move: current value plus symmetric noise
    proposal <- theta + rnorm(1, mean = 0, sd = proposal_sd)
    log_p_proposal <- log_target(proposal)
    # 2. accept with probability min(1, p(proposal) / p(theta)); on the log
    #    scale that is: accept if log(u) < log_p_proposal - log_p
    if (log(runif(1)) < log_p_proposal - log_p) {
      theta <- proposal
      log_p <- log_p_proposal
      n_accepted <- n_accepted + 1
    }
    # 3. record where the chain is (the old value again, if we rejected)
    samples[i] <- theta
  }
  list(samples = samples, acceptance_rate = n_accepted / n_iter)
}
# end snippet metropolis_1d

# ---- Metropolis on a known target: Normal(5, 2) ---------------------------

# start snippet metropolis_normal
log_target <- function(theta) dnorm(theta, mean = 5, sd = 2, log = TRUE)
run <- metropolis_1d(log_target, n_iter = 5000, init = 0, proposal_sd = 5)
# end snippet metropolis_normal

fig("06_metropolis_normal.png")
layout(matrix(c(1, 2), 1), widths = c(2, 1))
par(mar = c(4.5, 4.5, 3, 1))
plot(run$samples, type = "l", col = "grey30", xlab = "Iteration", ylab = expression(theta),
     main = sprintf("Trace plot (acceptance rate %.0f%%)", 100 * run$acceptance_rate))
abline(h = 5, col = "tomato", lwd = 2)
hist(run$samples, freq = FALSE, breaks = 40, col = "grey85", border = "white",
     xlab = expression(theta), main = "Sampled distribution")
curve(dnorm(x, 5, 2), add = TRUE, col = "tomato", lwd = 3)
legend("topright", bty = "n", legend = c(sprintf("mean %.2f (true 5)", mean(run$samples)),
                                         sprintf("sd %.2f (true 2)", sd(run$samples))))
dev.off()

# ---- Tuning the proposal ---------------------------------------------------

fig("06_proposal_tuning.png", height = 7)
par(mfrow = c(3, 1), mar = c(4, 4.5, 2.5, 1))
for (sd in c(0.2, 5, 100)) {
  r <- metropolis_1d(log_target, n_iter = 1000, init = 5, proposal_sd = sd)
  plot(r$samples, type = "l", col = "grey30", ylim = c(-2, 12),
       xlab = "Iteration", ylab = expression(theta),
       main = sprintf("proposal sd = %g: acceptance rate %.0f%%", sd, 100 * r$acceptance_rate))
  abline(h = 5, col = "tomato", lwd = 2)
}
dev.off()

# ---- Burn-in ---------------------------------------------------------------

fig("06_burn_in.png", height = 4.5)
par(mfrow = c(1, 3), mar = c(4.5, 4.5, 3, 1))
for (start in c(5, 40, 100)) {
  r <- metropolis_1d(log_target, n_iter = 500, init = start, proposal_sd = 5)
  plot(r$samples, type = "l", col = "grey30", ylim = c(-3, 105),
       xlab = "Iteration", ylab = expression(theta), main = paste("start at", start))
  rect(0, -10, 100, 110, col = adjustcolor("steelblue", 0.15), border = NA)
  abline(h = 5, col = "tomato", lwd = 2)
  text(100, 100, "burn-in", col = "steelblue", pos = 4)
}
dev.off()

# ---- Convergence: several chains -------------------------------------------

# start snippet rhat
# Gelman-Rubin R-hat: variance between chains relative to variance within them.
# Chains that have found the same distribution give R-hat close to 1.
rhat <- function(chains) {
  n <- length(chains[[1]])
  W <- mean(sapply(chains, var))            # mean within-chain variance
  B <- n * var(sapply(chains, mean))        # n x variance of the chain means
  sqrt(((n - 1) / n * W + B / n) / W)
}
# end snippet rhat

log_bimodal <- function(theta) log(0.5 * dnorm(theta, -4, 1) + 0.5 * dnorm(theta, 4, 1))
starts <- c(-6, -2, 2, 6)
cols <- c("tomato", "steelblue", "forestgreen", "darkorange")

fig("06_convergence.png", height = 4.5)
par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3, 1))
for (sd in c(0.5, 4)) {
  chains <- lapply(starts, function(s) metropolis_1d(log_bimodal, 2000, s, sd)$samples)
  plot(NA, xlim = c(0, 2000), ylim = c(-8, 8), xlab = "Iteration", ylab = expression(theta),
       main = sprintf("proposal sd = %g:  R-hat = %.2f", sd, rhat(chains)))
  for (k in seq_along(chains)) lines(chains[[k]], col = adjustcolor(cols[k], 0.8))
}
dev.off()
