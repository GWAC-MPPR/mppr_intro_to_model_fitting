# Least squares example for 03_least_squares.qmd.
#
# Simulates synthetic SIR case data, saves the "true vs observed" plot to
# figures/03_sir_true_vs_observed.png, then fits beta and gamma by a grid
# search on the sum of squared errors and writes the best fit to
# results/03_least_squares_fit.txt.
#
# The slides pull code from the regions between `# start snippet <name>` and
# `# end snippet <name>` markers (see _extensions/quarto-ext/include-code-files).
# Everything outside the markers is plumbing that stays off the slides, so keep
# the code inside them exactly as it should appear on the slide.
#
# Run from the project root: Rscript scripts/03_least_squares.R

# Load required packages
library(deSolve)
library(ggplot2)

dir.create("figures", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)

# Define SIR model
sir_model <- function(t, y, params) {
  with(as.list(c(y, params)), {
    dS <- -beta * S * I
    dI <- beta * S * I - gamma * I
    dR <- gamma * I
    list(c(dS, dI, dR))
  })
}

# Generate synthetic data
set.seed(42)
true_params <- c(beta = 0.3, gamma = 0.1)
init_conds <- c(S = 0.99, I = 0.01, R = 0)
times <- seq(0, 365, by = 1)

# Simulate true model
out <- ode(y = init_conds, times = times, func = sir_model, parms = true_params)
true_cases <- out[,"I"] * 1000  # Scale to cases

# Add noise
observed_cases <- rpois(length(times), true_cases + 5)

# Plot
plot_data <- data.frame(
  time = times,
  observed = observed_cases,
  true_model = true_cases
)

p <- ggplot(plot_data, aes(x = time)) +
  geom_point(aes(y = observed), color = "red", alpha = 0.7) +
  geom_line(aes(y = true_model), color = "blue", linewidth = 1) +
  labs(x = "Time (days)", y = "Number of Cases",
       title = "SIR Model: True vs Observed Data") +
  theme_minimal()

ggsave("figures/03_sir_true_vs_observed.png", p,
       width = 10, height = 6, dpi = 150, bg = "white")

# start snippet ls_implementation
# Define objective function
sse_function <- function(params, data) {
  beta <- params[1]
  gamma <- params[2]
  
  # Simulate model
  out <- ode(
    y = init_conds,
    times = data$time,
    func = sir_model,
    parms = c(beta = beta, gamma = gamma)
  )
  
  #  Extract and rescale the predicted cases
  predicted <- out[, "I"] * 1000
  # Calculate sum of squared errors
  sse <- sum((data$observed - predicted)^2)
  
  return(sse)
}
# end snippet ls_implementation

# start snippet ls_fitting
# Prepare data
data <- data.frame(time = times, observed = observed_cases)

# Test different parameter values
test_params <- expand.grid(
  beta = seq(0.1, 0.5, length.out = 10),
  gamma = seq(0.05, 0.2, length.out = 10)
)
# end snippet ls_fitting

# start snippet ls_fitting_2
# Calculate SSE for each combination
for (i in 1:nrow(test_params)) {
  test_params$sse[i] <- sse_function(
    c(
      test_params[i, 1],
      test_params[i, 2]
    ),
    data
  )
}
# end snippet ls_fitting_2

# Console output between the sink() calls also goes to the results file, which
# the slides include verbatim.
sink("results/03_least_squares_fit.txt", split = TRUE)
# start snippet ls_results
# Find minimum
best_idx <- which.min(test_params$sse)
best_params <- test_params[best_idx, ]
cat("True parameters: Beta =", true_params[1], ", Gamma =", true_params[2], "\n")
cat("Best fit parameters: Beta =", round(best_params$beta, 3), ", Gamma =", round(best_params$gamma, 3), "\n")
# end snippet ls_results
sink()
