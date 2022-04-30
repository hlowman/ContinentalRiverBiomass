

data {
  int<lower=0> sites; // number of sites
  int<lower=0> Nobs; // maximum number of observations possible at a site, to bound the matrices
  
  int Ndays [sites]; // number of days for each site
  
  // matrix [rows, columns]
  matrix [Nobs, sites] light; // relativized to max value
  matrix [Nobs, sites] GPP; // mean estimates from posterior probability distributions
  matrix [Nobs, sites] GPP_sd; // sd estimates from posterior probability distributions
  matrix [Nobs, sites] tQ; // standardized discharge
  matrix [Nobs, sites] new_e; // 0/1s denoting new time sequences
  
}


parameters {
  // Random effect parameters
  // Each parameter calculated at the site-level needs its own random effect.
  real<lower=0> csigma; // critical discharge random effect standard deviation site
  vector [sites] csite; // critical discharge random effect intercept for each site
  
  real<lower=0> ssigma; // sensitivity random effect standard deviation site
  vector [sites] ssite; // sensitivity random effect intercept for each site
  
  real<lower=0> rsigma; // max growth rate random effect standard deviation site
  vector [sites] rsite; // max growth rate random effect intercept for each site
  
  real<lower=0> lsigma; // lambda random effect standard deviation site
  vector [sites] lsite; // lambda rate random effect intercept for each site
  
  // Disturbance (persistence) parameters
  // Single values because these parameters are the population mean
  real<lower=0> c; // estimate of Qcrit - critical discharge
  real<lower=0> s; // steepness of the transition from P=1 to P=0 - persistence curve
  
  // Logistic growth parameters  
  matrix [Nobs, sites] B; // Biomass
  // Single values because these parameters are the population mean
  real r; // growth rate
  real lambda; // r/K - growth rate/carrying capacity
  
  // Error parameters
  real<lower=0> sig_p; // sigma processes error
  real<lower=0> sig_o; // sigma observation error
}


transformed parameters {
  matrix [Nobs, sites] pred_GPP;
  matrix [Nobs, sites] P;
 
  // Loop over sites
  // For each of the 207 sites...
  for(h in 1:sites) { 
    
  // and for each day within each of the sites...  
  for (j in 1:(Ndays[h])) {
    
      P[j,h] = exp(-exp(ssite[h]*100*(tQ[j,h] - csite[h]))); // persistence as a function of discharge
      
      pred_GPP[j,h] = light[j,h]*exp(B[j,h]); // predicted GPP as a function of light and biomass
      
    }
  
  }
  
}


model {
  
  // Loop over sites
  // For each of the 207 sites...
  for(h in 1:sites){
    
  // Initial value
  B[1,h] ~ normal(log(GPP[1,h]/light[1,h]), 1);
  
  // and for each day within each of the sites...
  // Process Model - reinitialize for every new time sequence
  for (j in 2:(Ndays[h])) {
    
    if (new_e[j,h]==1) { // 1 = TRUE
    
    B[j,h] ~ normal(log(GPP[j,h]/light[j,h]),1);
    
    } 
    
    else { // 0 = FALSE
    
    B[j,h] ~ normal(((B[(j-1),h] + rsite[h]*rsigma + lsite[h]*lsigma*exp(B[(j-1),h]))*P[j,h]), sig_p);
    
    }
    
  }
  
  // Observation model
  for (j in 2:(Ndays[h])) {
    
    GPP[j,h] ~ normal(pred_GPP[j,h], sig_o)T[0,];
    
  }
  
  // Random effect priors
  csite[h] ~ normal(c,csigma)T[0,];
  ssite[h] ~ normal(s,ssigma)T[0,];
  rsite[h] ~ normal(r,rsigma); // or T[0,]
  lsite[h] ~ normal(lambda,lsigma); // or T[,0]
  
  }

  // Cauchy distributions are normal distributions w/ *heavy* tails.
  csigma ~ cauchy(0,100);
  ssigma ~ cauchy(0,100);
  rsigma ~ cauchy(0,100);
  lsigma ~ cauchy(0,100);
  
  // Error priors
  sig_p ~ normal(0,2)T[0,];
  sig_o ~ normal(mean(GPP_sd), sd(GPP_sd))T[0,];
  // does the mean/sd above need to be reformatted to match incoming matrices?
  
  // Param priors
  c ~ normal(0.5,0.25)T[0,]; // revised
  s ~ normal(1.5,1)T[0,]; // revised
  r ~ normal(0,1);
  lambda ~ normal(0,1)T[,0];
  
}




