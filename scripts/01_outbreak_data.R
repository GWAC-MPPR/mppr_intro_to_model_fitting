# Synthetic outbreak data for the "Example: outbreak data" slide
# (01_introduction.qmd).
#
# Simulates an SIR epidemic, adds Poisson noise to mimic reported cases and
# saves the plot to figures/01_outbreak_data.png. None of this code is shown
# on the slides; only the figure is.
#
# Run from the project root: Rscript scripts/01_outbreak_data.R

library(deSolve)

dir.create("figures", showWarnings = FALSE)

# Simulate COVID-19-like data
set.seed(123)
days <- 1:100
true_beta <- 0.3
true_gamma <- 0.1
true_R0 <- true_beta / true_gamma

# SIR model simulation
sir_sim <- function(t, y, params) {
  with(as.list(c(y, params)), {
    dS <- -beta * S * I
    dI <- beta * S * I - gamma * I
    dR <- gamma * I
    list(c(dS, dI, dR))
  })
}

out <- ode(y = c(S = 0.99, I = 0.01, R = 0),
           times = days,
           func = sir_sim,
           parms = c(beta = true_beta, gamma = true_gamma))

# Add noise to simulate real data
observed_cases <- rpois(length(days), out[,"I"] * 1000 + 10)

png("figures/01_outbreak_data.png", width = 10, height = 6, units = "in", res = 150)
plot(days, observed_cases, pch = 16, col = "red",
     xlab = "Days", ylab = "Daily Cases",
     main = "Outbreak Data")
lines(days, out[,"I"] * 1000, col = "blue", lwd = 2)
legend("topright", c("Observed", "True Model"),
       col = c("red", "blue"), pch = c(16, NA), lty = c(NA, 1))
dev.off()
