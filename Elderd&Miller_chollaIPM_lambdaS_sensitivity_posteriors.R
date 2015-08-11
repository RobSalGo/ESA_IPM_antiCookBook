#############################################################################
## Author: Tom Miller (tom.miller@rice.edu)
##
## Purpose: Estimate posterior distributions of lambda_S and its sensitivities to vital rate coefficients
## given the joint posterior of the vital rate coefficients
##
## Notes: (1) This code loads data from the MCMC output ("cholla.all.params.post.csv") and it requires source functions for the IPM ("chollaIPM_SOURCE.R")
##        (2) At the current settings, this script takes ~8 hrs to run on a desktop computer. Run times will depend on the number of draws taken from 
##            the joint posterior distribution ("Ndraws"), the number of iterations for stochastic population dynamics ("iter"), and the number of plots
##            from which parameter values are sampled ("Nplots"). Also, run time will depend on whether or not the sensitivity analysis is included.
##            Currently, the sensitivity analysis is commented out. If you wish to run the sensitivity analysis, a supercomputer is recommended. 
############################################################


## -------- load packages ---------------------- ##

library(R2jags); library(lattice); library(mcmcplots)

## -------- set working directory ---------------------- ##

setwd(getwd())

## -------- read in MCMC output ---------------------- ##

##This file contains random draws from the joint posterior distribution of all parameters
post.params<-read.csv("cholla.all.params.post.csv")

## Number of draws to take from the joint posterior distribution of the parameters. 
## Cannot be greater than the number of draws provided in the .csv file, which is 500.
Ndraws<-min(100,nrow(post.params))
post.params<-post.params[1:Ndraws,]

## -------- load IPM source functions ---------------------- ##

source("chollaIPM_SOURCE.R")

## -------- Set up IPM parameter vector ---------------------- ##

## 'cholla' is a matrix where rows are vital rate coefficients and columns are posterior draws
## below, we will loop over columns, sending each set of coefficients into the stochastic IPM
cholla<-matrix(NA,nrow=100,ncol=Ndraws) 

## Params 01-10: Growth
cholla[1,]<-post.params$g.alpha      	  ## growth intercept
cholla[2,]<-post.params$g.beta					## growth slope
cholla[3,]<-post.params$y.beta.2.				## growth y.beta
cholla[4,]<-post.params$g.eps					  ## SD of growth function
cholla[5,]<-post.params$sigma.yr.2.			## SD of interannual growth variance
cholla[6,]<-post.params$g.sigma.plot  	## SD of among-plot growth variance

## Params 11-20: Survival params
cholla[11,]<-post.params$s.alpha 				## survival intercept
cholla[12,]<-post.params$s.beta					## survival slope
cholla[13,]<-post.params$y.beta.1.			## survival y.beta
cholla[14,]<-post.params$sigma.yr.1.		## SD of interannual survival variance
cholla[15,]<-post.params$s.sigma.plot  	## SD of among-plot survival variance

## Params 21-30: Probability of flowering params
cholla[21,]<-post.params$pf.alpha   		## pr.flower intercept
cholla[22,]<-post.params$pf.beta				## pr.flower slope
cholla[23,]<-post.params$y.beta.3.			## pr.flower y.beta
cholla[24,]<-post.params$sigma.yr.3.		## SD of interannual pr.flower variance
cholla[25,]<-post.params$pf.sigma.plot  ## SD of among-plot pr.flower variance

## Params 31-40: Fertility params
cholla[31,]<-post.params$f.alpha     		## fertility intercept
cholla[32,]<-post.params$f.beta					## fertility slope
cholla[33,]<-post.params$y.beta.4.			## fertility y.beta
cholla[34,]<-post.params$sigma.yr.4.		## SD of interannual fertility variance
cholla[35,]<-post.params$f.sigma.plot  	## SD of among-plot fertility variance
cholla[36,]<-post.params$f.sigma        ## SD of fertility overdispersion

## Params 41-50: Seed/seeling params
cholla[41,]<-post.params$seed.num      	          ## seeds per fruit
cholla[42,]<-post.params$fruit.disp				          ## survival from seed dispersal to seed bank
cholla[43,]<-post.params$intercept.one			## pr germination from 1yo seed bank
cholla[44,]<-post.params$intercept.two		## pr germination from 2yo seed bank
cholla[45,]<-post.params$precens.surv     ## survival from germination to census
cholla[46,]<-post.params$kid.mean         ## mean recruit size
cholla[47,]<-post.params$kid.sd           ## sd of recruit size

## Params 51-60: Model-wide year effect
cholla[51,]<-post.params$sigma.year     ## sd of mwye

## Params 61-70: Misc params (bounds of continuous size domain, in units of log(cm^3)). Hard coded from Miller's data. 
cholla[61,]<- -4.78  ## minsize 
cholla[62,]<- 13.65  ## maxsize 


## -------- Stochastic IPM settings ---------------------- ##
  
## size of the approximating matrix to the IPM kernel
matsize<-200

## number of years to iterate stochastic population dynamics
iter<-1000  ## how many years to sample

## number of plots to sample in each year
Nplots<-10 ## how many plots to sample

## perturbation for sensitivity analysis
perturbation<-0.01

## -------- Initiatlize vectors to store output ---------------------- ##

## these are the indices of parameters for which we estimate sensitivity to the mean (indices correspond to rows of matrix 'cholla')
params.mu<-c(1,2,4,11,12,21,22,31,32,41,42,43,44,45,46,47)

## store lambdaS values
post.lambdaS<-vector("numeric",length=Ndraws*Nplots) ## 2 rows for correlation on and off

## sensitivities to parameter means
sensSmu.post<-matrix(NA,length(params.mu),Ndraws*Nplots)

################################################################################
#LOOP over draws from the joint posterior of all parameters -- WARNING: this loop can take a long time to run!

#Progress bar 
pb.max <- Ndraws*Nplots    
pb.tally <- 0                     
ptm <- proc.time()                
pb <- txtProgressBar(style=3, min=0, initial=0, max=pb.max) 

for(i in 1:Ndraws) {
  
  ## sample a sequence of random deviates representing plot-to-plot variance in each of the four main vital rates
  plotfx <- matrix(0,4,Nplots)
    plotfx[1,] <- rnorm(n=Nplots,mean=0,sd=cholla[6,i]) # Growth
    plotfx[2,] <- rnorm(n=Nplots,mean=0,sd=cholla[15,i]) # Survival 
    plotfx[3,] <- rnorm(n=Nplots,mean=0,sd=cholla[25,i]) # Probability of flowering 
    plotfx[4,] <- rnorm(n=Nplots,mean=0,sd=cholla[35,i]) # Fertility
  
  ## sample a sequence of random deviates representing year-to-year variance in each of the four main vital rates, individually
  yrfx <- matrix(0,4,iter)  
    yrfx[1,] <- rnorm(n=iter,mean=0,sd=cholla[5,i]) # Growth
    yrfx[2,] <- rnorm(n=iter,mean=0,sd=cholla[14,i]) # Survival 
    yrfx[3,] <- rnorm(n=iter,mean=0,sd=cholla[24,i]) # Probability of flowering 
    yrfx[4,] <- rnorm(n=iter,mean=0,sd=cholla[34,i]) # Fertility
  
  ## sample a sequence of random deviates representing year-to-year variance that will link the four main vital rates
  mwye <- rnorm(n=iter,mean=0,sd=cholla[51,i]) ##model-wide year effect
  
  ## loop over plots, calculating a stochastic growth rate for each
  for(p in 1:Nplots){
    
  post.lambdaS[p+Nplots*(i-1)]<-lambda.fun(parameters=cholla[,i],plotfx=plotfx[,p],yrfx=yrfx,mwye=mwye,iter=iter,matsize=matsize)
  
  #######################################################################################
  ###The following lines conduct a sensitivity analysis. This is commented out because it
  ###adds substantially to run times.
  ######################################################################################
  ############### collect the kernels and eigenvectors for this parameter set
  #cholla.wvK<-wvK(parameters=cholla[,i],yrfx=yrfx,plotfx=plotfx[,p],mwye=mwye,iter=iter,matsize=matsize) 
  #  
  ################# loop over the parameter vector for this plot and calculate sensitivity for each
  #for(j in 1:length(params.mu)){
  #
  ################# perturb this parameter and re-calculate kernel
  #  cholla.perturb.K<-perturb.K(parameters=cholla[,i],yrfx=yrfx,plotfx=plotfx[,p],
  #                              mwye=mwye,iter=iter,
  #                              matsize=matsize,floor.extend=1,ceiling.extend=4,extra.grid=2,
  #                              perturbation=perturbation,which.perturb=params.mu[j],Kt=cholla.wvK$Kt)
  #  
  ############## calculate sensitivity
  # sensSmu.post[j,p+Nplots*(i-1)]<-sens_S_mu(lambdaS=post.lambdaS[p+Nplots*(i-1)],vt=cholla.wvK$vt,wt=cholla.wvK$wt,
  #                            Kt=cholla.wvK$Kt,dKt.dthetai=cholla.perturb.K$dKt.dthetai,
  #                            k=round(0.05*iter),iter=iter)
  #
  #}##end parameter loop
  ############################################################################################
  
  # Update the progress bar
  pb.tally<- pb.tally + 1; setTxtProgressBar(pb, pb.tally)
  
  }##end plot loop

} ## end posterior draw loop, go to the next posterior sample

# Total processing time
print(proc.time() - ptm)



##########################################################################################
###########################################################################################
########################################################################################

#######################  FIGURE: posterior distribution of lambda_S #################################
dens<-density(post.lambdaS)
par(mar=c(5,5,2,2))
plot(dens$x,dens$y,type="n",xlim=c(0.75,1.05),xlab=expression(lambda[S]),ylab="Density",cex.lab=1.6)
polygon(dens,col="gray",border="gray")
abline(v=quantile(post.lambdaS,probs=0.025),lwd=2,lty=2)
abline(v=quantile(post.lambdaS,probs=0.975),lwd=2,lty=2)
abline(v=mean(post.lambdaS),lwd=3)


#######################  FIGURE: posterior distributions of sensitivities (d.lambda_S/d.parameter) ##################
####################### Sensitivities are split into two panels, A and B, as in the main paper ######################
####################### As in the code above, the sensitivity figure is commented out here. #########################

############# this list of symbols corresponds (in order) to the parameters list 'params.mu'. It is used to draw the axis labels.
#param.symbols<-c(expression(mu[alpha]^G),expression(beta^G),expression(sigma^G),
#                 expression(mu[alpha]^S),expression(beta^S),
#                 expression(mu[alpha]^P[Fl]),expression(beta^P[Fl]),
#                 expression(mu[alpha]^F),expression(beta^F),
#                 expression(mu[Sd]),expression(delta),"g1","g2",expression(alpha^phi))

############# sensitivities of the means
#par(mfrow=c(1,2),mar=c(5,5,2,2))
#boxplot(cbind(sensSmu.post[1,],
#        sensSmu.post[2,],
#        sensSmu.post[3,],
#        sensSmu.post[4,],
#        sensSmu.post[5,]),
#        names=param.symbols[1:5],xlab="Parameter",cex.lab=1.4,outline=F,
#        ylab=expression(paste("Sensitivity to parameter mean,  ",partialdiff, lambda[S]," / ",partialdiff,theta[i])),
#        col="gray")
#title(main="A",adj=0,font=3)

#boxplot(cbind(sensSmu.post[6,],
#              sensSmu.post[7,],
#              sensSmu.post[8,],
#              sensSmu.post[9,],
#              sensSmu.post[10,],
#              sensSmu.post[11,],
#              sensSmu.post[12,],
#              sensSmu.post[13,],
#              sensSmu.post[14,]),
#        names=param.symbols[6:14],xlab="Parameter",cex.lab=1.4,outline=F,
#        ylab=expression(paste("Sensitivity to parameter  mean,  ",partialdiff, lambda[S]," / ",partialdiff,theta[i])),
#        col="gray")
#title(main="B",adj=0,font=3)

