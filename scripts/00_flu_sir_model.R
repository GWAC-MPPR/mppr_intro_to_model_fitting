# The running example for the whole deck: the 1978 influenza outbreak in a
# British boarding school, and the SIR model we fit to it. Every other script
# sources this one, so the data and the model are defined in one place.
#
# The slides pull code from the regions between `# start snippet <name>` and
# `# end snippet <name>` markers (see filters/include-code.lua), so keep the
# code inside them exactly as it should appear on the slide.

# start snippet flu_data
flu <- read.csv("data/influenza_england_1978_school.csv")
flu$day <- seq_len(nrow(flu))      # days 1 to 14
flu$cases <- flu$in_bed            # pupils in bed with flu on each day
# end snippet flu_data

# start snippet sir_model
library(deSolve)

sir_model <- function(t, y, params) {
  with(as.list(c(y, params)), {
    dS <- -beta * S * I
    dI <- beta * S * I - gamma * I
    dR <- gamma * I
    list(c(dS, dI, dR))
  })
}

N <- 763                            # pupils in the school
init <- c(S = 762, I = 1, R = 0)    # one index case on day 0
times <- 0:14                       # the data are days 1 to 14

# Model prediction of the number of pupils in bed on days 1 to 14
predict_cases <- function(beta, gamma) {
  out <- ode(y = init, times = times, func = sir_model,
             parms = c(beta = beta, gamma = gamma))
  out[-1, "I"]
}
# end snippet sir_model
