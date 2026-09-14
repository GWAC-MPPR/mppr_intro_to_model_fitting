# The roadmap slide shown before the learning objectives: where this lecture
# sits in the MPPR course. Dashed boxes are modules already covered or still to
# come; the solid box is today.
#
#   figures/00_course_roadmap.png
#
# Run from the project root: Rscript scripts/00_course_roadmap.R

dir.create("figures", showWarnings = FALSE)

png("figures/00_course_roadmap.png", width = 13, height = 6.2, units = "in", res = 150, pointsize = 15)
par(mar = c(0, 0, 0, 0))
plot.new()
plot.window(xlim = c(0, 13), ylim = c(0, 6.2))

box <- function(x, y, w, h, title, lines, lty = 2, border = "grey40", fill = "white",
                title_col = border, text_col = "grey20", fonts = rep(1, length(lines))) {
  rect(x - w / 2, y - h / 2, x + w / 2, y + h / 2, col = fill, border = border, lty = lty, lwd = 3)
  text(x, y + h / 2 - 0.42, title, font = 2, cex = 1.05, col = title_col)
  # one text() call per line so that individual lines can be bold
  line_y <- y + h / 2 - 0.78 - 0.22 * (seq_along(lines) - 1)
  for (i in seq_along(lines)) text(x, line_y[i], lines[i], cex = 0.82, col = text_col, font = fonts[i])
}

# ---- The row of modules ------------------------------------------------------
y_top <- 4.3
box(1.9, y_top, 3.0, 2.7, "Programming essentials",
    c("R essentials", "data frames, plots,", "loops and functions,", "reading and writing data",
      "Version control", "git and GitHub"),
    fonts = c(2, 1, 1, 1, 2, 1))
box(5.4, y_top, 3.4, 2.7, "Infectious disease\ndynamics modelling",
    c("", "compartmental models:", "SIR, SEIR and extensions", "R0, interventions,", "heterogeneity, stochasticity"))
box(9.0, y_top, 3.4, 2.7, "Model fitting and\ncalibration",
    c("", "connecting the model to data:", "least squares,", "maximum likelihood,", "Bayesian inference and MCMC"),
    lty = 1, border = "tomato", fill = "#fff3ef")
box(12.0, y_top, 1.9, 2.7, "What comes\nnext",
    c("", "scenario modelling,", "health economic", "modelling,", "science", "communication,", "capstone projects"))

arrows(3.45, y_top, 3.65, y_top, length = 0.12, lwd = 3, col = "grey40")
arrows(7.15, y_top, 7.25, y_top, length = 0.12, lwd = 3, col = "grey40")
arrows(10.75, y_top, 11.0, y_top, length = 0.12, lwd = 3, col = "grey40")

text(1.9, y_top + 1.65, "covered", col = "grey40", cex = 0.85)
text(5.4, y_top + 1.65, "covered", col = "grey40", cex = 0.85)
text(9.0, y_top + 1.65, "today", col = "tomato", font = 2, cex = 0.95)
text(12.0, y_top + 1.65, "ahead", col = "grey40", cex = 0.85)

# ---- The calibration challenge -----------------------------------------------
y_low <- 1.3
box(6.2, y_low, 3.0, 1.5, "A model", c("SIR with unknown \u03b2 and \u03b3"))
box(10.2, y_low, 3.0, 1.5, "Outbreak data", c("cases, admissions or deaths", "counted over time"))
# the model comes out of the modelling module; both feed today's question
arrows(5.4, y_top - 1.35, 5.4, y_low + 0.75, length = 0.12, lwd = 2, col = "grey40")
arrows(7.5, y_low + 0.75, 7.5, y_top - 1.35, length = 0.12, lwd = 2, col = "tomato")
arrows(9.5, y_low + 0.75, 9.5, y_top - 1.35, length = 0.12, lwd = 2, col = "tomato")
text(8.2, 0.28, "Which values of the parameters produced these data, and how sure can we be?",
     col = "tomato", font = 2, cex = 0.95)
dev.off()
