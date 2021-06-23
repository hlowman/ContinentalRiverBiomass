//
// This Stan program defines a simple model, with a
// vector of values 'y' modeled as normally distributed
// with mean 'mu' and standard deviation 'sigma'.
//
// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//
// The below example will be my first attempt at writing a basic
// linear model.
// I will be using the mtcars dataset for practice.

// The input data is two vectors 'y' and 'x' of length 'N'.
data {
  int<lower=0> N; // # of observations
  vector[N] y; // outcome
  vector[N] x; // predictor
}

// The parameters accepted by the model. Our model
// accepts three parameters 'intercept', 'slope', and 'sigma'.
parameters {
  real intercept;
  real slope;
  real<lower=0> sigma; // standard error needs to be positive
  // for state-space model, also need process (sigp) and observation (sigo) error
}

// The model to be estimated. We model the output
// 'y' to be normally distributed as a function of the slope*x + intercept
// and standard deviation 'sigma'.
model {
  // Priors
  intercept ~ normal(30,5); // intercept has a normal distribution
  slope ~ normal(0,10); // broad normal distribution
  sigma ~ normal(0,0.5); // changed this around some
  // for datasets with greater variance, may need to consider exp()
  // add in sigp and sigo priors here for state-space models
  
  // Likelihood
  y ~ normal(intercept + slope*x, sigma);
  
  // - process model describes the larger mechanism
  // - observation model describes the error between process model results
  // and actual observations (in the case of GPP, on a daily basis)
}
