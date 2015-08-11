################################################################################
## Author: Tom Miller (tom.miller@rice.edu)
## Purpose: Source functions for the tree cholla IPM. These functions are required to run
##          'Elderd&Miller_chollaIPM_lambdaS_sensitivity_poseriors.R'. 
## NOTES: These functions assmble the IPM kernel, calculate the stochastic growth rate via simulation, and 
##        calculate the corresponding sensitivities. Sensitivity functions are based on Rees and Ellner 2009 (Ecol. Monographs)
################################################################################

## ----------- Miscellany...we'll need an inverse logit functions ------------- ##
invlogit<-function(x){exp(x)/(1+exp(x))}

## ----------- Vital rate functions. Parameter indices are hard-coded and must correspond to rows of MCMC matrix ------------- ##

#GROWTH FROM SIZE X TO Y
gxy<-function(x,y,params,yrfx,mwye,plotfx){
  xb=pmin(pmax(x,params[61]),params[62]) #Transforms all values below/above limits in min/max size
  return(dnorm(y,mean=params[1] + params[2]*xb + params[3]*mwye + yrfx[1] + plotfx[1],sd=params[4]))
}

#SURVIVAL AT SIZE X.
sx<-function(x,params,yrfx,mwye,plotfx){
  xb=pmin(pmax(x,params[61]),params[62])
  return(invlogit(params[11] + params[12]*xb  + params[13]*mwye + yrfx[2] + plotfx[2]))
}

#SURVIVAL*GROWTH
pxy<-function(x,y,params,yrfx,mwye,plotfx){
  xb=pmin(pmax(x,params[61]),params[62])
  sx(xb,params,yrfx,mwye,plotfx)*gxy(xb,y,params,yrfx,mwye,plotfx)
}

#PRODUCTION OF 1-YO SEEDS IN THE SEED BANK FROM X-SIZED MOMS
fx<-function(x,params,yrfx,mwye,plotfx,f.eps){
  xb=pmin(pmax(x,params[61]),params[62])
  p.flow<-invlogit(params[21] + params[22]*xb + params[23]*mwye + yrfx[3] + plotfx[3]) 
  nfruits<-exp(params[31] + params[32]*xb + params[33]*mwye + yrfx[4] + plotfx[4] + f.eps)   
  seeds.per.fruit<-params[41]
  seed.survival<-invlogit(params[42])^2  ## I measured 6-month seed survival; annual survival is its square
  return(p.flow*nfruits*seeds.per.fruit*seed.survival)  
}

#SIZE DISTRIBUTION OF RECRUITS
recruits<-function(y,params){
  yb=pmin(pmax(y,params[61]),params[62])
  dnorm(x=yb,mean=params[46],sd=params[47])
}

## ----------- Construct IPM kernel ------------- ##

bigmatrix<-function(params,yrfx,plotfx,mwye,f.eps,lower,upper,matsize){  
###################################################################################################
## returns the full IPM kernel (to be used in stochastic simulation), the F and T kernels, and meshpoints in the units of size
## params,yrfx,plotfx, and mwye get passed to the vital rate functions
## f.eps is fertility overdispersion. defaults to zero (see lambda.fun())
## lower and upper are the integration limits
## matsize is the dimension of the approximating matrix (it gets an additional 2 rows and columns for the seed banks)
###################################################################################################

  n<-matsize
  L<-lower; U<-upper
  h<-(U-L)/n                   #Bin size
  b<-L+c(0:n)*h;               #Lower boundaries of bins 
  y<-0.5*(b[1:n]+b[2:(n+1)]);  #Bin midpoints
  
  # Fertility matrix
  Fmat<-matrix(0,(n+2),(n+2))
  
  # Banked seeds go in top row
  Fmat[1,3:(n+2)]<-fx(y,params=params,yrfx=yrfx,mwye=mwye,plotfx=plotfx,f.eps=f.eps)
  
  # Growth/survival transition matrix
  Tmat<-matrix(0,(n+2),(n+2))
  
  # Graduation to 2-yo seed bank = pr(not germinating as 1-yo)
  Tmat[2,1]<-1-invlogit(params[43])
  
  # Graduation from 1-yo bank to cts size = germination * size distn * pre-census survival
  Tmat[3:(n+2),1]<-invlogit(params[43])*recruits(y,params)*h*invlogit(params[45])   
  
  # Graduation from 2-yo bank to cts size = germination * size distn * pre-census survival
  Tmat[3:(n+2),2]<-invlogit(params[44])*recruits(y,params)*h*invlogit(params[45])  
  
  # Growth/survival transitions among cts sizes
  Tmat[3:(n+2),3:(n+2)]<-t(outer(y,y,pxy,params=params,yrfx=yrfx,mwye=mwye,plotfx=plotfx))*h
  
  # Put it all together
  IPMmat<-Fmat+Tmat     
  
  return(list(IPMmat=IPMmat,Fmat=Fmat,Tmat=Tmat,meshpts=y))
}

## ----------------- Function that simulates population dynamics and returns lambdaS ------------- ############

lambda.fun<-function(parameters,yrfx,plotfx,f.eps=0,mwye,iter,
                     matsize,extra.grid=2,floor.extend=1,ceiling.extend=4,stochastic=T,corr=T){
############################################################################################
## This function returns the population growth rate for a given set of parameters
## Defaults to lambdaS (stochastic=T) but can give deterministic lambda based on vital rate means
## extra.grid adds the discrete seed stages to the cts kernel
## floor extend and ceiling.extend correct eviction (see Williams, Miller, Ellner (2012), Ecology)
## corr==T (default) includes vital rate correlations. corr==F sets mwye to zero and therefore turns correlations off
############################################################################################
  
  ## IPM bounds
  lower<- parameters[61] - floor.extend
  upper<- parameters[62] + ceiling.extend
  
  if(stochastic==F){
    yrfx <- matrix(0,1,4)
    mwye <- 0
    lambda<-Re(eigen(bigmatrix(params=parameters,yrfx=yrfx,plotfx=plotfx,mwye=mwye,f.eps=f.eps,lower=lower,upper=upper,matsize=matsize)$IPMmat)$values[1])
    return(lambda)
  }
  if(corr==F){mwye <- rep(0,iter)}
  
  ## place to store annual growth rates
  rtracker <- numeric(iter)
  ## initial population vector
  n0 <- rep(1/(matsize + extra.grid),length=matsize + extra.grid)
  
  # This loop calculates stochastic lambda for iter years
  for(g in 1:(iter-1)){
    
    K.ipm<- bigmatrix(params=parameters,yrfx=yrfx[,g],plotfx=plotfx,mwye=mwye[g],f.eps=f.eps,lower=lower,upper=upper,matsize=matsize)$IPMmat 
    n0<- K.ipm %*% n0
    N <- sum(n0)
    rtracker[g]<-log(N)
    n0<-n0/N
  }
  
  burnin <- round(length(rtracker)*0.1) ## drop the first 10% of the time series
  rtracker <- rtracker[-c(1:burnin)]
  meanR <- mean(rtracker)
  
  return(exp(meanR))
  
}

## ----------------- Functions for stochastic sensitivities ------------- ############

## function to return time series of w's and v's (eigenvectors), K's (matrices), and one-time-step lambdas
wvK<-function(parameters,yrfx,plotfx,mwye,iter,f.eps=0,
              matsize,floor.extend=1,ceiling.extend=4,extra.grid=2){
##############################
## all arguments defined above
############################
  
  ## set kernel bounds
  lower<-parameters[61]
  upper<-parameters[62]
  
  ##Generate a sequence of kernels K0...KT and corresponding eigenvalues
  Kt<-list()
  lambdas<-c()
  for(t in 1:iter){
    Kt[[t]]<-bigmatrix(params=parameters,yrfx=yrfx[,t],plotfx=plotfx,mwye=mwye[t],f.eps=f.eps,
                       lower=lower-floor.extend,upper=upper+ceiling.extend,matsize=matsize)$IPMmat
    ##lambdas[t]<-Re(eigen(Kt[[t]])$values[1])
  }
  ## calculate the stochastic growth rate from this sequence of kernels
  ##lambdaS<-exp(mean(log(lambdas)))
  
  ##Generate and store a sequence of w vectors
  wt<-list()
  wt[[1]]<-rep(1/matsize,matsize+extra.grid)
  for(t in 2:iter){
    wt[[t]]<-Kt[[t-1]]%*%wt[[t-1]]/sum(Kt[[t-1]]%*%wt[[t-1]])
  }
  
  ##Generate and store a sequence of v vectors
  vt<-list()
  vt[[iter]]<-rep(1/matsize,matsize+extra.grid)
  for(t in (iter-1):1){
    vt[[t]]<-vt[[t+1]]%*%Kt[[t]]/sum(vt[[t+1]]%*%Kt[[t]]) 
  }
  
  return(list(Kt=Kt,wt=wt,vt=vt))#,lambdas=lambdas,lambdaS=lambdaS))
}

### function to return first and second derivative of the kernel wrt coefficient i
perturb.K<-function(parameters,yrfx,plotfx,mwye,iter,f.eps=0,
                    matsize,floor.extend=1,ceiling.extend=4,extra.grid=2,
                    perturbation,which.perturb,Kt){
#################################################################
## 'perturbation' is the absolute values of the perturbation.
## 'which.perturb' specific the index of the parameter we are perturbing
## 'Kt' is the array of year-specific IPM kernels
## All other arguments defined above.
######################################################################
  lower<-parameters[61]
  upper<-parameters[62]
  
  parameters.P1<-parameters
  parameters.P2<-parameters
  parameters.P1[which.perturb]<-parameters[which.perturb]+perturbation
  parameters.P2[which.perturb]<-parameters[which.perturb]-perturbation
  
  ## calculate the first derivative of the kernel wrt the perturbed parameter
  Kt.P1<-list()
  Kt.P2<-list()
  dKt.dthetai<-list()## first derivative
  dKt2.dthetai2<-list()##second derivative
  for(t in 1:iter){
    Kt.P1[[t]]<-bigmatrix(params=parameters.P1,yrfx=yrfx[,t],plotfx=plotfx,mwye=mwye[t],f.eps=f.eps,
                          lower=lower-floor.extend,upper=upper+ceiling.extend,matsize=matsize)$IPMmat
    Kt.P2[[t]]<-bigmatrix(params=parameters.P2,yrfx=yrfx[,t],plotfx=plotfx,mwye=mwye[t],f.eps=f.eps,
                          lower=lower-floor.extend,upper=upper+ceiling.extend,matsize=matsize)$IPMmat
    dKt.dthetai[[t]]<-(Kt.P1[[t]]-Kt.P2[[t]])/(2*perturbation)
    dKt2.dthetai2[[t]]<-(Kt.P1[[t]]-2*Kt[[t]]+Kt.P2[[t]])/(perturbation^2)
  } 
  return(list(dKt.dthetai=dKt.dthetai,dKt2.dthetai2=dKt2.dthetai2))
  
}

## sensitivity of lambdaS to mean coefficient theta_i
sens_S_mu<-function(lambdaS,vt,wt,Kt,dKt.dthetai,k,iter){
#################################################################
## 'lambdaS' is the stochastic growth rate.
## 'vt' and 'wt' are arrays of year-specific eigenvectors 
## 'Kt' is the array of year-specific IPM kernels
## 'dKt.dthetai' is the array of year-specific IPM kernel derivatives wrt parameter i
## 'k' is a burn-in parameter, discarding the beginning and end of the stochastic series
## All other arguments defined above.
######################################################################
  hold<-c()
  for(t in k:(iter-k)){
    hold[t]<-(vt[[t+1]]%*%(dKt.dthetai[[t]]%*%wt[[t]]))/(vt[[t+1]]%*%(Kt[[t]]%*%wt[[t]]))
  }
  sens<-lambdaS*mean(hold,na.rm=T)
  return(sens)
}
