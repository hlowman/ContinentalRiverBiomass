stanmodelcode ="


data {
int<lower=0> s; // number of sites

int<lower=0> m; // number of modeled size classes 
vector[m+1] M; // Modeled size class cutoffs
vector[m] Mmid; // Modeled size class midpoints

int<lower=0> mF; // number of measured size classes 

int t;// total number of years
int ST[s]; //year of seeding--starting is 1979
int MT[s]; //year of measurment--starting is 1979
int Szdata[mF,s] ; ///ending size distributions all sites
vector[s] Sddata ; ///starting seed info for all sites
real area[s]; //Search area
matrix[s,2] Xdata[t]; ///covariates 


}




parameters {

real survb0; //survival intercept
real survbS;  //survival size
real survbSM; //survival moisture
///////
real growb0; //size transitions intercept
real growbS;  //size transitions size
real growbSM; //size transitions moisture
real<lower=0> grows; //size transitions kernal standard dev
///////
real fecdb0; //fecundity intercept
real fecdbS; //fecundity size


real<lower=0> sigmaG; //random effect standard deviation growth
real<lower=0> sigmaF; //random effect standard deviation fecund
real<lower=0> sigmaS; //random effect standard deviation survival

real<lower=0> sigmaGsz; //random effect slope standard deviation growth
real<lower=0> sigmaFsz; //random effect slope standard deviation fecund
real<lower=0> sigmaSsz; //random effect slope standard deviation survival
real<lower=0> kappa; //neg binomial dispersion


vector[s] Gs; //random effect intercept for each site growth 
vector[s] Fs; //random effect intercept for each sit fecund 
vector[s] Ss; //random effect intercept for each sit survival

vector[s] Gs1; //random effect size slope for each site growth 
vector[s] Fs1; //random effect size slope for each sit fecund 
vector[s] Ss1; //random effect size slope for each sit survival

}

transformed parameters {

///output
matrix[m,s] gamlat; /// size structure for IPM

///// vital rates
vector[m]  mug; /// vector of average size transitions by size class
vector[m] surv; /// vector of survival rates by size class
vector[m]  fecd; ///vector of fecund rates by size class

/// kernal matricies
matrix[m,m] grow ; ///desity kernal of transitions for transitions from one size to another for each year
matrix[m,m] grow1 ; ///normalized desity kernal of transitions for transitions from one size to another for each year 
matrix[m,m] kern; ///desity kernal that combines all vital rates 


vector[m]  szout;  //vector to store output of sizes

//Loop over sites
for(i in 1:s){ 



matrix[m,(MT[i]-ST[i])+1]  gam; //make matrix of m size classes over every year from seeding to measurement

  
//Fill in matrix with zeros  
for (x in 1:m){
for (tt in 1:(MT[i]-ST[i])+1){
gam[x,tt]=0;
}}


//Loop over years from seeding (ST[i]) to measurement (MT[i])
for (ii in (ST[i]):MT[i]){


mug[:] = growb0+(growbS+Gs1[i]*sigmaGsz)*Mmid+growbSM*Xdata[ii,i,1]+Gs[i]*sigmaG; /// regression of sz_t+1 as a fucntion of sz_t and covariates 
fecd[:] = exp(fecdb0+(fecdbS+Fs1[i]*sigmaFsz)*Mmid+Fs[i]*sigmaF); /// regression of fecundity_t as a fucntion of sz_t 
surv[:] = inv_logit(survb0+(survbS+Ss1[i]*sigmaSsz)*Mmid+survbSM*Xdata[ii,i,1]+Ss[i]*sigmaS); /// regression of survival_t as a fucntion of sz_t and covariates 


//loop over all combos of size class x and y
for (y in 1:m){
for (x in 1:m){

grow[y,x]= exp(normal_lpdf(Mmid[y] | mug[x], grows)); //Calculate probability of transitioning (i.e. growing) from one size class [x] to another [y](exp of normal_lpdf returns to probability scale) )
} //end size indexing loops
}

//Loop back over x (columns) to normalize
for (x in 1:m){

grow1[:,x]=grow[:,x]/sum(grow[:,x]); //normalize each column (x) to sum to 1

}

//if the year == the seeding year start the simulation 
if (ii-ST[i]==0){
for (y in 1:m){
gam[y,1]=Sddata[i]*grow1[y,1]*surv[y]; //First year recruits from seed-- new recruits grow then survive. 
}}

else{
// if it is a year other than ST[i] construct a full kernal

for (y in 1:m){;  
for (x in 1:m){; //// loops over size classes to construct kernal
kern[y,x] = grow1[y,x]*surv[x]+fecd[x]*grow1[y,1]*surv[y]; ///construct kernal


}
}

for (y in 1:m){; 
for (x in 1:m){; /// loops over size classes for population projection 

szout[x]= kern[y,x]*gam[x,(ii-ST[i])];//Project population into the future-- pull gam (size structure) from the last time step

}
gam[y,(ii-ST[i])+1]=sum(szout); //sum each row to get new size stucture -- this becomes gam for the current time step
}

}




} // end time loop

gamlat[:,i]=gam[:,(MT[i]-ST[i])+1]; //save gam from the year of measurement as gamlat

} // end site loop 



}

//Likelihood and priors
model {
  //Loop over sites (s) and measured size classes (mF)
for(i in 1:s){
for(sz in 1:mF){

Szdata[sz,i]~neg_binomial_2(gamlat[sz,i]*area[i],kappa);// Likelihood ---occurance rate (gamlat) is adjusted by search area (area) 
}
//Random effects
//site intercept random effects
Gs[i]~normal(0,1);
Fs[i]~normal(0,1);
Ss[i]~normal(0,1);
//site size random effects
Gs1[i]~normal(0,1);
Fs1[i]~normal(0,1);
Ss1[i]~normal(0,1);
}

//Priors
survb0~normal(0,5); //survival intercept
survbS~normal(0,5);  //survival size
survbSM~normal(0,5); //survival moisture
///////
growb0~normal(0,5); //Size transitions intercept
growbS~normal(0,5);  //Size transitions size
growbSM~normal(0,5); //Size transitions moisture
grows~cauchy(0,5);  // Size transitions kernal standard dev
///////
fecdb0~normal(0,5); //fecundity intercept
fecdbS~normal(0,5); //fecundity size

////// Standard deviation for site intercept random effects
sigmaG~cauchy(0,1); //Size transitions
sigmaF~cauchy(0,1); //Fecundity
sigmaS~cauchy(0,1); //Survival

////// Standard deviation for site size random effects
sigmaGsz~cauchy(0,1); //Size transitions
sigmaFsz~cauchy(0,1); //Fecundity
sigmaSsz~cauchy(0,1); //Survival
//////
kappa~cauchy(0,5); //Dispersion

}


generated quantities {
  
  //Generate predictions of data
  int preddata[mF,s];

for(i in 1:s){
for(sz in 1:mF){
 preddata[sz,i] <- neg_binomial_2_rng(gamlat[sz,i]*area[i],kappa);; //Predicts data as output using random neg binom number generator
}}
}
"
