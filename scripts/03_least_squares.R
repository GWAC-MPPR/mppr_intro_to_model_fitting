# Least squares fit of the SIR model to the 1978 influenza outbreak
# (03_least_squares.qmd).
#
# Writes:
#   figures/03_sse_surface.png       the sum of squared errors over a grid
#   figures/03_ls_fit.png            the least squares fit against the data
#   results/03_least_squares_fit.txt the estimates
#
# The slides pull code from the regions between `# start snippet <name>` and
# `# end snippet <name>` markers (see filters/include-code.lua). Everything
# outside the markers is plumbing that stays off the slides, so keep the code
# inside them exactly as it should appear on the slide.
#
# Run from the project root: Rscript scripts/03_least_squares.R

source("scripts/00_flu_sir_model.R")
dir.create("figures", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)

fig <- function(name, width = 10, height = 6) {
  png(file.path("figures", name), width = width, height = height,
      units = "in", res = 150, pointsize = 15)
}

# start snippet ls_implementation
# The objective function: how far is the model from the data?
sse_function <- function(params) {
  beta <- params[1]
  gamma <- params[2]
  predicted <- predict_cases(beta, gamma)
  sum((flu$cases - predicted)^2)
}

# The three guesses from the introduction
sse_guesses <- c(A = sse_function(c(0.0015, 0.5)),
                 B = sse_function(c(0.0030, 0.6)),
                 C = sse_function(c(0.0022, 0.4)))
# end snippet ls_implementation

sink("results/03_sse_guesses.txt", split = TRUE)
print(round(sse_guesses))
sink()

# start snippet ls_grid
# Try a grid of parameter values
test_params <- expand.grid(
  beta = seq(0.0015, 0.0030, length.out = 20),
  gamma = seq(0.30, 0.70, length.out = 20)
)

for (i in seq_len(nrow(test_params))) {
  test_params$sse[i] <- sse_function(c(test_params$beta[i], test_params$gamma[i]))
}

best <- test_params[which.min(test_params$sse), ]
# end snippet ls_grid

fig("03_sse_surface.png")
par(mar = c(4.5, 4.5, 3, 1))
sse_matrix <- matrix(log10(test_params$sse), nrow = 20)
image(unique(test_params$beta), unique(test_params$gamma), sse_matrix,
      col = hcl.colors(50, "YlOrRd"), xlab = expression(beta), ylab = expression(gamma),
      main = "Sum of squared errors over the grid (log scale, darker = smaller)")
contour(unique(test_params$beta), unique(test_params$gamma), sse_matrix, add = TRUE, col = "grey30", nlevels = 8)
points(best$beta, best$gamma, pch = 4, cex = 2, lwd = 3)
dev.off()

# start snippet ls_optim
# Let an optimiser polish the best grid point
fit_ls <- optim(par = c(best$beta, best$gamma), fn = sse_function)
ls_estimates <- c(beta = fit_ls$par[1], gamma = fit_ls$par[2])
# end snippet ls_optim

sink("results/03_least_squares_fit.txt", split = TRUE)
# start snippet ls_results
cat("Best grid point:  beta =", signif(best$beta, 3), " gamma =", signif(best$gamma, 3),
    " SSE =", round(best$sse), "\n")
cat("After optim():    beta =", signif(ls_estimates["beta"], 3),
    " gamma =", signif(ls_estimates["gamma"], 3), " SSE =", round(fit_ls$value), "\n")
cat("R0 = beta N / gamma =", round(ls_estimates["beta"] * N / ls_estimates["gamma"], 2),
    "  infectious period =", round(1 / ls_estimates["gamma"], 2), "days\n")
# end snippet ls_results
sink()

fine_times <- seq(0, 14, by = 0.1)
ls_curve <- ode(y = init, times = fine_times, func = sir_model, parms = ls_estimates)[, "I"]

fig("03_ls_fit.png")
par(mar = c(4.5, 4.5, 3, 1))
plot(flu$day, flu$cases, pch = 16, cex = 1.3, ylim = c(0, 320),
     xlab = "Day", ylab = "Pupils in bed", main = "Least squares fit")
lines(fine_times, ls_curve, col = "tomato", lwd = 3)
legend("topright", bty = "n", pch = c(16, NA), lwd = c(NA, 3), col = c("black", "tomato"),
       legend = c("data", sprintf("SIR, beta = %.5f, gamma = %.3f", ls_estimates["beta"], ls_estimates["gamma"])))
dev.off()
