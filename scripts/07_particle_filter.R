# Particle filter sketch for the "Example application" slide in
# 07_advanced_methods.qmd.
#
# This is an illustration of the shape of a pomp model, not a runnable example:
# it refers to data (covid_times, covid_observed) and C snippets (sir_step,
# cases_measure, ...) that are not defined here, and it is never executed.
# The slides show it as is.

# start snippet particle_filter
# Particle filter for SIR model
library(pomp)

# Define SIR model with stochasticity
sir_pomp <- pomp(
  data = data.frame(time = covid_times, cases = covid_observed),
  times = "time",
  t0 = 0,
  rprocess = euler.sim(
    step.fun = "sir_step",
    delta.t = 0.1
  ),
  rmeasure = "cases_measure",
  dmeasure = "cases_dmeasure",
  initializer = "sir_init",
  paramnames = c("beta", "gamma", "sigma"),
  statenames = c("S", "I", "R")
)

# Run particle filter
pf <- pfilter(sir_pomp, Np = 1000, params = c(beta = 0.3, gamma = 0.1, sigma = 1))
# end snippet particle_filter
