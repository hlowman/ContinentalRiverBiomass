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
}

// The model to be estimated. We model the output
// 'y' to be normally distributed as a function of the slope*x + intercept
// and standard deviation 'sigma'.
model {
  // Priors
  intercept ~ normal(0,0.5); // intercept has a normal distribution
  slope ~ normal(0,0.5); // same with slope
  sigma ~ normal(0,0.5); // changed this around some
  
  // Likelihood
  y ~ normal(intercept + slope*x, sigma);
}
