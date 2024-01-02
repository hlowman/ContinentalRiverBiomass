

data {
  int Ndays; // number of days
  vector [Ndays] light; // relativized to max value
  vector [Ndays] GPP; // mean estimates from posterior probability distributions
  vector [Ndays] GPP_sd; // sd estimates from posterior probability distributions
  vector [Ndays] tQ; // standardized discharge (to 10 yr flood)
  int new_e [Ndays]; // 0/1s denoting new time sequences
}

parameters {
  // Disturbance (persistence) parameters
  real<lower=0> c; // estimate of Qcrit - critical discharge - a.k.a. discharge necessary to scour benthic algal biomass
  real<lower=0> s; // steepness of the transition of the persistence curve, from P=1 (100% day-to-day algal persistence) to P=0 (0% day-to-day algal persistence)
  
  // Logistic growth parameters  
  real B [Ndays]; // Biomass
  real r; // growth rate
  real lambda; // -r/K (growth rate/carrying capacity)
  
  // Error parameters
  real<lower=0> sig_p; // sigma process error
  real<lower=0> sig_o; // sigma observation error
}

// the persistence model goes in the transformed parameters block
// persistence being the model incorporating measured daily discharge

transformed parameters {
  real pred_GPP [Ndays];
  real P [Ndays];
  
  for(i in 1:Ndays){
    P[i]=exp(-exp(s*100*(tQ[i] - c))); // persistence as a function of discharge
    pred_GPP[i] =light[i]*exp(B[i]); // predicted GPP as a function of light and biomass
  }
  
}

// the process and observation models go in the model block
// process being the model estimating daily algal biomass (the unobserved/latent
// variable) and observation being the model incoporating estimated daily GPP
// as well as measured daily light (now entirely from StreamLight)

// the format of the biomass model is a transformed variation of the Ricker model
// incorporating information from the persistence model in the block above

model {
  
  // Initial value
  // Need to set the initial value for t = 1
  B[1] ~ normal(log(GPP[1]/light[1]), 1);
  
  // Process Model - reinitialize for every new time sequence
  // If there is a gap of over 14 days in the timeseries, it will re-initialize
  // using the same formula used to establish day 1 values at every site.
  for (j in 2:(Ndays)){
    
    if (new_e[j]==1) { // 1 = TRUE
    
    B[j] ~ normal(log(GPP[j]/light[j]),1);
    
    } 
    
    else { // 0 = FALSE
    
    B[j] ~ normal((B[(j-1)] + r + lambda*exp(B[(j-1)]))*P[j], sig_p);
    
    }
    
  }
  
  // Observation model
  for (j in 2:(Ndays)) {
    GPP[j] ~ normal(pred_GPP[j], sig_o)T[0,];
  }
  
  // Error priors
  sig_p ~ normal(0,2)T[0,];
  sig_o ~ normal(mean(GPP_sd), sd(GPP_sd))T[0,];
  
  // Param priors
  c ~ normal(0.5,0.25); // revised spring 2022
  s ~ normal(1.5,1); // revised spring 2022
  r ~ normal(0.2,0.1); // r cannot be negative
  lambda ~ normal(-0.03,0.01); // lambda cannot be positive
  
}


