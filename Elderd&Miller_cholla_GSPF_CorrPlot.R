###########################################################
## Author: Bret Elderd
##
## Purpose: The code uses survival, growth, probability of
##          flowering, and number of flowers to calculate
##          coefficients for use in the Integral Projection Model.
##          The code also includes temporal correlation, which
##          is calculated as a model-wide year effect based on
##          Evans et al. 2010. Ecological Monographs, vol 80(4).
##
##  NOTE: (1) The plot effect is added to the likelihood and
##        not in a hierarchical fashion.
##        (2) Data of individual size should be a
##        comma delimited file named -- "Individuals.csv".
##        The file should contain columns pertaining to
##        individual size, whether or not the individual
##        survived over the course of a sampling period,
##        whether the individual flowered/reproduced, and
##        how many seeds or offspring did the individual
##        produce. If using volume instead of a single
##        size measure such as height, size at a particular
##        time should consist of height, maximum width, and
##        the width perpendicular to the maximum.
##        (3) Plot variables must be numeric for JAGS code.
##
############################################################

## ---------------- set working directory ---------- ##
setwd(getwd())

## ---------------- load packages ----------------- ##
library(R2jags); library(mcmcplots);library(lattice)

## ----------- functions ---------------------- ##

## inverse logit function
invlogit<-function(x){exp(x)/(1+exp(x))}
## calculate the volume based on size (if needed)
vol<-function(h,w,p){(1/3)*pi*h*(((w+p)/2)/2)^2}

## ---------- load demographic data ------- ##
demog.dat<-read.csv("Individuals.csv",T)

## calculate volume (if needed)
demog.dat$sizet<-vol(h=demog.dat$Height_t,w=demog.dat$Width_t,p=demog.dat$Perp_t)
demog.dat$sizet1<-vol(h=demog.dat$Height_t1,w=demog.dat$Width_t1,p=demog.dat$Perp_t1)

## --------- growth and prob. of flowering data frame (Survivors Only) ----------- ##
grow.dat <- demog.dat[demog.dat$survival==1,]

## --------- fertility data frame (Survivors who flowered Only) ----------- ##
fert.dat <- demog.dat[demog.dat$flowered==1,]

##########################
#  Bayes Model
##########################

## ----------- The jags model --------------- ##
sink("GSPFCorr.Plot.txt")
cat("
model{

# Hyperparameters
sigma.year ~ dunif(0,1)

s.beta~dnorm(0,0.001)
s.mu.plot <- 0
s.sigma.plot ~ dunif(0,1000)

g.beta~dnorm(0,0.001)
g.mu.plot <- 0
g.sigma.plot ~ dunif(0,1000)
g.eps~dunif(0,1000)

pf.beta~dnorm(0,0.001)
pf.mu.plot <- 0
pf.sigma.plot ~ dunif(0,1000)

f.beta~dnorm(0,0.001)
f.mu.plot <- 0
f.sigma.plot~dunif(0,1000)
f.sigma ~ dunif(0,10) # fert overdispersion

## Priors for n.life stages intercept
g.alpha ~ dnorm(0,0.001)
s.alpha ~ dnorm(0,0.001)
pf.alpha ~ dnorm(0,0.001)
f.alpha ~ dnorm(0,0.001)

for (i in 1:n.life){
  sigma.yr[i] ~ dunif(0,10) 
  tau.yr[i] <- 1/(sigma.yr[i]*sigma.yr[i])
 }

for(i in 1:n.plot){
  g.gamma[i] ~ dnorm(g.mu.plot,g.tau.plot)
  f.gamma[i] ~ dnorm(f.mu.plot,f.tau.plot)
  pf.gamma[i] ~ dnorm(pf.mu.plot,pf.tau.plot)
  s.gamma[i] ~ dnorm(s.mu.plot,s.tau.plot)
}

# Priors for y.beta
for (i in 1:(n.life-1)){
  y.beta[i] <- 2.0*I.beta[i] - 1.0
  I.beta[i] ~ dbern(pi.beta[i])
  pi.beta[i] ~ dunif(0,1)
}


## Year-wide effect centered on the intercept
# Effect of YEAR positively correlated with sized-based fert
y.beta[4] <- 1.0

# Prior for YEAR-wide priority effect
for (t in 1:n.years){
  YEAR[t] ~ dnorm(0.0, tau.year) # mu=0.0 from Evans et al. 2010
 }

# Survival Loop (id = 1)
for (t in 1:n.years){
   mu.yr[t,1] <- s.alpha + y.beta[1]*YEAR[t]
   YR[t,1] ~ dnorm(mu.yr[t,1], tau.yr[1])
}

# Growth Loop (id = 2)
for (t in 1:n.years){
   mu.yr[t,2] <- g.alpha + y.beta[2]*YEAR[t]
   YR[t,2] ~ dnorm(mu.yr[t,2], tau.yr[2])
}

# Prob. Flower Loop (id = 3)
for (t in 1:n.years){
   mu.yr[t,3] <- pf.alpha + y.beta[3]*YEAR[t]
   YR[t,3] ~ dnorm(mu.yr[t,3], tau.yr[3])
}

# Flower Fertility Loop (id = 4)
for (t in 1:n.years){
   mu.yr[t,4] <- f.alpha + y.beta[4]*YEAR[t]
   YR[t,4] ~ dnorm(mu.yr[t,4], tau.yr[4])
}

## Likelihood
for(i in 1:sN){
  logit(s.p[i]) <- YR[s.years[i],1] + s.gamma[s.plot[i]] + s.beta*s.x[i] 
  s.y[i]~dbern(s.p[i])
}

for(i in 1:gN){
  g.mu[i] <- YR[g.years[i],2] + g.gamma[g.plot[i]] + g.beta*g.x[i] 
  g.y[i]~dnorm(g.mu[i],g.tau)
}

for(i in 1:pfN){
  logit(f.p[i])<-YR[pf.years[i],3] + pf.gamma[pf.plot[i]] + pf.beta*pf.x[i]
  pf.y[i]~dbern(f.p[i])
}

for(i in 1:fertN){ # Fertility
  f.lambda[i]<-exp(YR[f.years[i],4] + f.gamma[f.plot[i]] + f.beta*f.x[i] + f.eps[i])
  f.y[i]~dpois(f.lambda[i])
  f.eps[i] ~ dnorm(0,f.tau)
}
## Derived quantities
g.tau<-1/(g.eps*g.eps)
g.tau.plot<-1/(g.sigma.plot*g.sigma.plot)
s.tau.plot<-1/(s.sigma.plot*s.sigma.plot)
pf.tau.plot<-1/(pf.sigma.plot*pf.sigma.plot)
f.tau.plot<-1/(f.sigma.plot*f.sigma.plot)
f.tau <- 1/(f.sigma*f.sigma) # fert overdispersion
tau.year <- 1/(sigma.year*sigma.year)

}##end model
",fill=T)
sink()

## Need to add thing to top on getting rid of the non-survivors for growth and similar for prob flower and number of flowers.

## --------- list of data needed for JAGS ----------- ##
jag.data<-list(s.y=demog.dat$Survival.t1,s.x=log(demog.dat$size.t),
               s.years=demog.dat$years, sN=nrow(demog.dat),s.plot=demog.dat$Plot,
               g.y=log(grow.dat$size.t1),g.x=log(grow.dat$size.t),
               g.years=grow.dat$years, gN=nrow(grow.dat), g.plot=grow.dat$Plot,
               pf.y=grow.dat$Flower.t, pf.plot=grow.dat$Plot,
               pf.x=log(grow.dat$size.t), pf.years=grow.dat$years,pfN=nrow(grow.dat),
               f.y=fert.dat$fertility,f.x=log(fert.dat$size.t),f.years=fert.dat$years,
               fertN=nrow(fert.dat), f.plot=fert.dat$Plot,
               n.life=4, n.years=length(unique(demog.dat$years)),
               n.plot=length(unique(fert.dat$Plot)))

## --------- Initial Values --------- ##
inits<-function(){list(s.beta=rnorm(1,1,2),s.sigma.plot=rlnorm(1),
                       g.beta=rnorm(1,1,2), g.sigma.plot=rlnorm(1), g.eps=rlnorm(1),
                       pf.beta=rnorm(1,1,2), pf.sigma.plot=rlnorm(1),
                       f.beta=rnorm(1,1,2), f.sigma.plot=rlnorm(1), f.sigma=rlnorm(1))}

## -- Parameters to be monitored (=to estimate) -- ##
parameters<-c("s.alpha","s.beta","s.gamma","s.sigma.plot",
              "g.alpha","g.beta","g.gamma","g.sigma.plot","g.eps",
              "pf.alpha","pf.beta","pf.gamma","pf.sigma.plot",
              "f.alpha","f.beta","f.gamma","f.sigma.plot","f.sigma", # "f.eps",
              "y.beta","sigma.yr","sigma.year")

## ---------- MCMC Settings -------- ##
## ni = number of iterations; nb = number of iterations for the burn-in
## nt = number of iterations to thin; nc = number of separate chains
ni<-1.8e4; nb<-4e3; nt<-100; nc<-6

## ----------- Start Gibbs sampler in JAGS --------------- ##
out.GSPFCorr.Plot<-jags(data=jag.data,inits=inits,parameters.to.save=parameters,
                        model.file="GSPFCorr.Plot.txt",
                        n.thin=nt,n.chains=nc,n.burnin=nb,n.iter=ni,
                        DIC=T,working.directory=getwd())

###########################
#   Check for Convergence
###########################
                                    
## ----- Create mcmc object ----- ##
out.GSPFCorr.Plot.mcmc <- as.mcmc(out.GSPFCorr.Plot)

## ----- Look at ALL of the output ----- ##
str(out.GSPFCorr.Plot)
print(out.GSPFCorr.Plot, dig=3)
## ----- The median ------ ##
str(out.GSPFCorr.Plot$BUGSoutput$median)

## --- Plot output, which contains R-hat --- ##
## (R-hat, the convergence stat, should be near 1.1)
plot(out.GSPFCorr.Plot)

#######################
# Plot of convergence
#######################
#  combines density, trace, and autocorrelation plots
#  outputs the result to a browser
mcmcplot(out.GSPFCorr.Plot.mcmc, parms=c("s.alpha","s.beta","s.gamma","s.sigma.plot",
                                     "g.alpha","g.beta","g.gamma","g.sigma.plot","g.eps",
                                     "pf.alpha","pf.beta","pf.gamma","pf.sigma.plot",
                                     "f.alpha","f.beta","f.gamma","f.sigma.plot",
                                     "sigma.yr","sigma.year"))

