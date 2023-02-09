## Maximum algal accrual calculation comparison

## Insert dummy parameter values
r1 <- 0.1
r2 <- 0.4
K1 <- 5 
K2 <- 15

#########
## Direct Comparisons
#########

## r1 and K1 comparison between equations
## K1 is a smaller carrying capacity river
exp(r1 - (r1/K1)*exp(0.5*K1)) # 0.866191
(0.5*K1*exp(0.5*r1)) - 0.5*K1 # 0.1281777

## r1 and K2 comparison between equations & K magnitude
## K2 is a larger carrying capacity river
exp(r1 - (r1/K2)*exp(0.5*K2)) # 6.435918e-06
(0.5*K2*exp(0.5*r1)) - 0.5*K2 # 0.3845332

## r1 and K2 comparison between equations & K magnitude
exp(r2 - (r2/K1)*exp(0.5*K1)) # 0.5629303
(0.5*K1*exp(0.5*r2)) - 0.5*K1 # 0.5535069


#########
## Range of r and K values
#########

## compare across a range of r values
r_seq <- seq(0.05, 0.9, by=0.01)

max_acc1 <- c() # Bob
max_acc2 <- c() # Heili

for(i in 1:length(r_seq)){
  max_acc1[i] = exp(r_seq[i] - (r_seq[i]/K1)*exp(0.5*K1))
  max_acc2[i] = (0.5*K1*exp(0.5*r_seq[i])) - 0.5*K1
}

plot(r_seq, max_acc1)
plot(r_seq, max_acc2) # linear increase


## compare across a range of K values
K_seq <- seq(1, 30, by=0.5)

max_acc1.K <- c()
max_acc2.K <- c()

for(i in 1:length(K_seq)){
  max_acc1.K[i] = exp(r1 - (r1/K_seq[i])*exp(0.5*K_seq[i]))
  max_acc2.K[i] = (0.5*K_seq[i]*exp(0.5*r1)) - 0.5*K_seq[i]
}

plot(K_seq, max_acc1.K)
plot(K_seq, max_acc2.K) # linear increase

# Replacing Bt = K/2 with Bt = ln(K/2).

max_acc3 <- c() # Heili 2.0

for(i in 1:length(r_seq)){
  max_acc3[i] = exp(r_seq[i] - (r_seq[i]/K1)*exp(log(K1/2)-1))
}

plot(r_seq, max_acc3) # linear increase again

max_acc3.K <- c() # Heili 2.0

for(i in 1:length(K_seq)){
  max_acc3.K[i] = exp(r1 - (r1/K_seq[i])*exp(log(K_seq[i]/2)-1))
}

plot(K_seq, max_acc3.K) # flat line?

# End of script.
