

data {
  int Ndays; // number of days
  vector [Ndays] light; // relativized to max value
  vector [Ndays] GPP; // mean estimates from posterior probability distributions
  vector [Ndays] GPP_sd; // sd estimates from posterior probability distributions
  vector [Ndays] tQ; // standardized discharge
  int new_e [Ndays]; // 0/1s denoting new time sequences
  // vector of integers
}

parameters {
  
  // Logistic growth parameters  
  real B [Ndays]; // Biomass
  real r; // growth rate
  real lambda; // r/K - growth rate/carrying capacity
  
  // Error parameters
  real<lower=0> sig_p; // sigma processes error
  real<lower=0> sig_o; // sigma observation error
}

transformed parameters {
  real pred_GPP [Ndays];
  
  for(i in 1:Ndays){
    pred_GPP[i] =light[i]*exp(B[i]); // predicted GPP as a function of light and biomass
  }
  
}


model {
  
  // Initial value
  B[1] ~ normal(log(GPP[1]/light[1]), 1);
  
  // Process Model - reinitialize for every new time sequence
  for (j in 2:(Ndays)){
    
    if (new_e[j]==1) { // 1 = TRUE
    
    B[j] ~ normal(log(GPP[j]/light[j]),1);
    
    } 
    
    else { // 0 = FALSE
    
    B[j] ~ normal((B[(j-1)] + r + lambda*exp(B[(j-1)])), sig_p);
    
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
  r ~ normal(0,1);
  lambda ~ normal(0,1)T[,0];
  
}


