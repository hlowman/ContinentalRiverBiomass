

data {
  int<lower=0> sites; // number of sites (groups)
  int<lower=1> Ndays; // number of days (observations)
  vector [Ndays] light; // relativized to max value
  vector [Ndays] GPP; // mean estimates from posterior probability distributions
  vector [Ndays] GPP_sd; // sd estimates from posterior probability distributions
  vector [Ndays] tQ; // standardized discharge
  int new_e [Ndays]; // 0/1s denoting new time sequences
}

parameters {
  // Random effect parameter
  real<lower=0> Sigmasite; //random effect standard deviation site
  vector[sites] REsite; //random effect intercept for each site
  
  // Disturbance (persistence) parameters
  real<lower=0> c; // estimate of Qcrit - critical discharge
  real<lower=0> s; // steepness of the transition from P=1 to P=0 - persistence curve
  
  // Logistic growth parameters  
  real B [Ndays]; // Biomass
  real r; // growth rate
  real lambda; // r/K - growth rate/carrying capacity
  
  // Error parameters
  real<lower=0> sig_p; // sigma processes error
  real<lower=0> sig_o; // sigma observation error
}

transformed parameters {
  
  // Loop over sites
  for(h in 1:sites){ 
    
  real pred_GPP [Ndays];
  real P [Ndays];
    
  for(i in 1:Ndays){
    
    P[i]=exp(-exp(s*100*(tQ[i] - c))); // persistence as a function of discharge
    pred_GPP[i] = light[i]*exp(B[i]); // predicted GPP as a function of light and biomass
  }
  
  }
  
}


model {
  
  // Loop over sites
  for(h in 1:sites){
    
  // Initial value
  B[1] ~ normal(log(GPP[1]/light[1]), 1);
  
  // Process Model - reinitialize for every new time sequence
  for (j in 2:(Ndays)){
    
    if (new_e[j]==1) { // 1 = TRUE
    
    B[j] ~ normal(log(GPP[j]/light[j]),1);
    
    } 
    
    else { // 0 = FALSE
    
    B[j] ~ normal(((B[(j-1)] + r + lambda*exp(B[(j-1)]))*P[j]) + REsite[h]*Sigmasite, sig_p);
    // only adding RE here to pool info re: r and lambda values??
    
    }
    
  }
  
  // Observation model
  for (j in 2:(Ndays)) {
    GPP[j] ~ normal(pred_GPP[j], sig_o)T[0,];
  }
  
  // Random effect priors
  REsite[h] ~ normal(0,1);
  
  }

  Sigmasite ~ cauchy(0,1);
  
  // Error priors
  sig_p ~ normal(0,2)T[0,];
  sig_o ~ normal(mean(GPP_sd), sd(GPP_sd))T[0,];
  
  // Param priors
  c ~ normal(0.5,0.25)T[0,]; // revised
  s ~ normal(1.5,1)T[0,]; // revised
  r ~ normal(0,1);
  lambda ~ normal(0,1)T[,0];
  
}


