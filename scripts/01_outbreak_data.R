# Figures for the introduction (01_introduction.qmd):
#   figures/01_flu_data.png        the outbreak data
#   figures/01_candidate_fits.png  three guesses at beta and gamma over the data
#
# Run from the project root: Rscript scripts/01_outbreak_data.R

source("scripts/00_flu_sir_model.R")
dir.create("figures", showWarnings = FALSE)

fig <- function(name, width = 10, height = 6) {
  png(file.path("figures", name), width = width, height = height,
      units = "in", res = 150, pointsize = 15)
}

fig("01_flu_data.png")
par(mar = c(4.5, 4.5, 3, 1))
plot(flu$day, flu$cases, type = "b", pch = 16, lwd = 2, col = "grey20",
     xlab = "Day", ylab = "Pupils in bed",
     main = "Influenza in a boarding school, January 1978 (763 pupils)")
dev.off()

# Three guesses at the parameters, none of them fitted: the point of the
# figure is that some guesses are clearly worse than others, and that we need
# a rule for saying which is best.
guesses <- list(A = c(beta = 0.0015, gamma = 0.5),
                B = c(beta = 0.0030, gamma = 0.6),
                C = c(beta = 0.0022, gamma = 0.4))
cols <- c(A = "steelblue", B = "darkorange", C = "forestgreen")

fig("01_candidate_fits.png")
par(mar = c(4.5, 4.5, 3, 1))
plot(flu$day, flu$cases, pch = 16, cex = 1.3, ylim = c(0, 420),
     xlab = "Day", ylab = "Pupils in bed", main = "Three guesses at beta and gamma")
fine_times <- seq(0, 14, by = 0.1)
for (g in names(guesses)) {
  curve_g <- ode(y = init, times = fine_times, func = sir_model, parms = guesses[[g]])[, "I"]
  lines(fine_times, curve_g, col = cols[g], lwd = 3)
}
legend("topright", bty = "n", lwd = 3, col = cols,
       legend = sapply(names(guesses), function(g)
         sprintf("%s: beta = %.4f, gamma = %.1f", g, guesses[[g]]["beta"], guesses[[g]]["gamma"])))
dev.off()
