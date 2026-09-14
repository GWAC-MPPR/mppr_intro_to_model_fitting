// SIR model in Stan (walkthrough)
// Functions block: derivative of [S, I, R]
functions {
  vector SIR(real t, vector y, array[] real theta) {
    real S = y[1]; real I = y[2]; real R = y[3];
    real beta = theta[1]; real gamma = theta[2];
    vector[3] dydt;
    dydt[1] = -beta * S * I;
    dydt[2] =  beta * S * I - gamma * I;
    dydt[3] =  gamma * I;
    return dydt;
  }
}

// Data block: obs counts & time grid
data {
  int<lower=1> n_obs;
  int<lower=1> n_pop;
  array[n_obs] int y;
  real t0;
  array[n_obs] real ts;
}

// Parameters: beta, gamma, S0
parameters {
  array[2] real<lower=0> theta;     // {beta, gamma}
  real<lower=0,upper=1> S0;         // initial susceptible fraction
}

// Transformed params: ODE solve + Poisson rate
transformed parameters {
  vector[3] y_init;
  array[n_obs] vector[3] y_hat;
  array[n_obs] real lambda;
  y_init[1] = S0; y_init[2] = 1 - S0; y_init[3] = 0;
  y_hat = ode_rk45(SIR, y_init, t0, ts, theta);
  for (i in 1:n_obs) lambda[i] = y_hat[i, 2] * n_pop;
}

// Model: priors + likelihood
model {
  theta ~ lognormal(0, 1);
  S0    ~ beta(1, 1);
  y     ~ poisson(lambda);
}

// GQ: derived R0
generated quantities {
  real R_0 = theta[1] / theta[2];
}
