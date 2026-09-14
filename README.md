Introduction to model fitting and calibration
================

<!-- badges: start -->
<!-- badges: end -->

This repository contains materials for the “Introduction to model fitting and calibration” part of the [Modelling for Pandemic Preparedness and Response (MPPR) modular short course](https://vclass.g-wac.org/local/gwac_landing/index.php?track=mppr&page=overview) run by the German West African Centre for Global Health and Pandemic Prevention (G-WAC) at the Kwame Nkrumah University of Science and Technology. 

The MPPR course is designed to provide an introduction to infectious disease dynamics modelling, health economics modelling, and science communication. It targets individuals with a quantitative background as well as professionals in public health, epidemiology, and related disciplines. It is interactive and hands-on and focuses on practical applications of infectious disease dynamics modelling to pandemic response.

See the course's 2026 schedule [here](https://vclass.g-wac.org/local/gwac_landing/index.php?track=mppr&page=overview#schedule)

## Resources in this repository

* [Web-based slides](https://gwac-mppr.github.io/mppr_intro_to_model_fitting/#/title-slide) and [pdf version](https://github.com/GWAC-MPPR/mppr_intro_to_model_fitting/releases/latest/download/slides.pdf)

* [R practicals](https://github.com/GWAC-MPPR/mppr_intro_to_model_fitting_r_practicals):
  least squares and maximum likelihood, then Bayesian inference and MCMC with a
  hand-written Metropolis sampler, all fitting an SIR model to the 1978
  boarding-school influenza outbreak

## Building the slides

The deck holds no executable code. The R behind the worked examples lives in
`scripts/`; each script saves its figures to `figures/` and its printed
results to `results/`, and the slides include those files as they are, so
`quarto render slides.qmd` needs Quarto but not R.

To regenerate the figures and results after editing a script, run from the
project root (the renv library provides `deSolve`, `ggplot2` and `bbmle`):

```sh
Rscript scripts/00_course_roadmap.R        # the roadmap slide
Rscript scripts/01_outbreak_data.R        # the data and three guesses at the parameters
Rscript scripts/04_maximum_likelihood.R   # sources scripts/03_least_squares.R first
Rscript scripts/06_mcmc_schematics.R      # toy-target schematics for the MCMC section
Rscript scripts/06_mcmc_flu.R             # Metropolis fit to the influenza data (~30 s)
```

Every fitting example in the deck uses the 1978 boarding-school influenza
outbreak in `data/`; `scripts/00_flu_sir_model.R` loads it and defines the SIR
model, and the other scripts source it.

The code shown on the slides is pulled from the same scripts by
`filters/include-code.lua`: a code block with `include="scripts/x.R"
snippet="name"` shows the lines between the `# start snippet name` and
`# end snippet name` markers in that script, so the slides and the scripts
cannot drift apart. `scripts/06_sir_model.stan` and
`scripts/07_particle_filter.R` are illustrations only and are not run.

The materials here were prepared and taught by [Dr. James Mba Azam](https://jamesmbaazam.github.io/jamesmbaazam/)
