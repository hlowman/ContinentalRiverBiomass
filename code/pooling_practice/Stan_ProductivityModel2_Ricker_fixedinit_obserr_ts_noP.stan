

data {
  
  int Ndays; // number of days
  vector [Ndays] light; // relativized to max value
  vector [Ndays] GPP; // mean estimates from posterior probability distributions
  vector [Ndays] GPP_sd; // sd estimates from posterior probability distributions
  vector [Ndays] tQ; // standardized discharge
  int new_e [Ndays]; // 0/1s denoting new time sequences
  int<lower = 0, upper = 1> p_remove; // boolean (0/1) variable for P term filter
  // and finally, create an empty P dataset to feed the model if p_remove is TRUE
  // real<lower=0> p_data[p_remove ? 0 : 1]; // p_data is size 0 if p_remove is TRUE

}

parameters {
  // parameters will change based on persistence (P) term filter
  // dummy filter currently based on mean discharge - EDIT THIS!!!
  // if mean(tQ) >= 0.02 a.k.a. p_remove = 0 KEEP P term
  
  // using ternary operator, which reads [(condition) ? (true value) : (false value)]
  
  // Disturbance (persistence) parameters
  // c is size 0 if p_remove is TRUE
  real<lower=0> c[p_remove ? 0 : 1]; // estimate of Qcrit - critical discharge
  // s is size 0 if p_remove is TRUE
  real<lower=0> s[p_remove ? 0 : 1]; // steepness of the transition from P=1 to P=0 - persistence curve
  
  // Logistic growth parameters  
  real B [Ndays]; // Biomass
  real r; // growth rate
  real lambda; // r/K - growth rate/carrying capacity
  
  // Error parameters
  real<lower=0> sig_p; // sigma processes error
  real<lower=0> sig_o; // sigma observation error
    
}


transformed parameters {
  // transformed parameters will change based on persistence (P) term filter
  // dummy filter currently based on mean discharge - EDIT THIS!!!
  // if mean(tQ) >= 0.02 a.k.a. p_remove = 0 KEEP P term
  
  real pred_GPP [Ndays];
  // P is size 0 if p_remove is TRUE
  real P [p_remove ? 0 : Ndays]; // persistence term for use in process model below
  
  if(p_remove == 0){
    
    for(i in 1:Ndays) {
    P[i]=exp(-exp(s[1]*(tQ[i]-c[1]))); // persistence as a function of discharge
    }
    
  }
  
  // GPP is calculated the same either way
  for(i in 1:Ndays){
    pred_GPP[i] =light[i]*exp(B[i]); // predicted GPP as a function of light and biomass
  }
  
  }


model {
  
  // Initial value
  B[1] ~ normal(log(GPP[1]/light[1]), 1);
  
  // Process Model - 
  
  // first nested "if else" statement to filter for models with/without P term
  if (p_remove == 1) {
    
    // second nested "if else"" to reinitialize for every new time sequence
  for (j in 2:(Ndays)){
    
    if (new_e[j]==1) { // 1 = TRUE
    
    B[j] ~ normal(log(GPP[j]/light[j]),1);
    
    } 
    
    else { // 0 = FALSE
    
    B[j] ~ normal((B[(j-1)] + r + lambda*exp(B[(j-1)])), sig_p); // without P
    
    }
    
  }
  
  }
  
  else {
    
  // second nested "if else"" to reinitialize for every new time sequence
  for (j in 2:(Ndays)){
    
    if (new_e[j]==1) { // 1 = TRUE
    
    B[j] ~ normal(log(GPP[j]/light[j]),1);
    
    } 
    
    else { // 0 = FALSE
    
    B[j] ~ normal((B[(j-1)] + r + lambda*exp(B[(j-1)]))*P[j], sig_p); // with P
    
    }
    
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
  if (p_remove == 0) {
    c[1] ~ normal(0,1)T[0,];
    s[1] ~ normal(0,200)T[0,];
  } // c and s conditional on P term being included
  r ~ normal(0,1);
  lambda ~ normal(0,1)T[,0];
  
}


