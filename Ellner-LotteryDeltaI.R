### Program to compute storage effect (Delta.I, in Chesson's notation)
### using the simulation-only method described in the talk. 
### The model is the original Chesson-Warner 1981 lottery model
### with constant mortality rate delta (equal for the two species)
### and log(per capita birth rates) for the two species drawn from
### a bivarite Gaussian distribution.  

rm(list=ls())
library(MASS)
set.seed(121212)  #change me to see the sampling variability; 

# Assign parameters---
delta <- 0.25          # death rate
mu.B <- c(0.7,0.5);    # mean of log birth rate for the two species 
sigma.B <- c(0.8,0.8); # Std Dev of log birth rates 
rho <- 0.5;            # correlation of log birth rates
totT <- 10^5;        # number of generations to simulate

# Generate E and E-sharp (in this model, E is log B)
sigma <-cbind(c(sigma.B[1]^2,rho*sigma.B[1]*sigma.B[2]),
              c(rho*sigma.B[1]*sigma.B[2],sigma.B[2]^2))
B <- exp(mvrnorm(n=totT,mu=mu.B,Sigma=sigma))
B.sharp <- exp(mvrnorm(n=totT,mu=mu.B,Sigma=sigma))

# First "simulations" -------------------------------------
## Note, in the lottery model the resident species is always
## at population size N (total number of sites). In this model  
## model C1 and C2 are given by the ratio between the total number
## of competing larvae, and the number of open sites, 
## which is (E2*N)/(delta*N)=E2/delta when species 2 is resident. 
E1 <- B[,1]; E2 <- B[,2]; C1 <- C2 <- E2/delta; 
rbar.I = mean(log(1-delta + E1/C1)) 

# Second "simulations" ------------------------------------
## Note, we can compute per-capita population growth rates
## without simulating the model, because it is unstructured  
E1.sharp <- B.sharp[,1]; E2.sharp <- B.sharp[,2]; 
rsharp.I = mean(log(1-delta + E1.sharp/C1)); 
rsharp.R = mean(log(1-delta + E2.sharp/C2)); 

Delta.I = rbar.I - rsharp.I + rsharp.R;  

# output 
cat("rbar.I = ",round(rbar.I,digits=3),"\n"); 
cat("rsharp.I = ",round(rsharp.I,digits=3),"\n"); 
cat("rsharp.R = ",round(rsharp.R,digits=3),"\n"); 
cat("Delta.I = ", round(rbar.I - rsharp.I + rsharp.R, digits=3),"\n");  
