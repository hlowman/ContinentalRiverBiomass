## Practice STAN Model Fit Script
## June 3, 2021
## Heili Lowman

# Load packages
library(tidyverse)
library(rstan)
library(shinystan)
library(MCMCglmm)
library(rethinking)

# Run each time you load in "rstan"
rstan_options(auto_write=TRUE)
options(mc.cores=parallel::detectCores())

# Load and examine data
dat <- mtcars
plot(dat$mpg ~ dat$cyl)

# Prep data for STAN
# Must be a list
dat_stan <- list(
  N = nrow(dat),
  x = dat$cyl,
  y = dat$mpg
)

# Linear model run
# .stan file must end in a BLANK LINE, otherwise it will spit out an error
test_run <- stan("code/practice_model.stan",
                 data = dat_stan)

# Examine the results
# Model: mpg = intercept + slope*cylinders
print(test_run) # prints distributions of parameters
plot(test_run) # plots distributions of parameters, with confidence intervals explained
pairs(test_run, pars = c("intercept", "slope", "lp__")) # parameters pairs plot

# Using the rethinking package
precis(test_run)

# First time, with (0,1) priors for intercept and slope, exponential for sigma
#           mean   sd  5.5% 94.5% n_eff Rhat4
# intercept 1.28 1.01 -0.33  2.89  2693     1
# slope     2.43 0.30  1.95  2.91  2981     1
# sigma     9.62 1.09  8.05 11.52  2952     1

# Second time, with (0,0.5) priors for intercept and slope, exponential for sigma
#            mean   sd  5.5% 94.5% n_eff Rhat4
# intercept  0.52 0.49 -0.26  1.29  2907     1
# slope      2.04 0.28  1.58  2.47  2662     1
# sigma     10.38 1.26  8.56 12.50  2669     1

# Third time, with (0,0.5) prior for sigma
#           mean   sd 5.5% 94.5% n_eff Rhat4
# intercept 1.21 0.50 0.41  2.02  2902     1
# slope     2.42 0.15 2.17  2.66  2936     1
# sigma     5.14 0.24 4.78  5.53  3166     1

# Using the shiny app
launch_shinystan(test_run)

# End of script.
