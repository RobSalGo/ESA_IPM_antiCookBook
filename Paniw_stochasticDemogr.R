# R script accompanying the OOS 17 presentation 
# "Population viability of a disturbance-dependent plant species in natural and human-induced environments"
# at the ESA 2015, Baltimore

# This script is a DEMO consiting of 3 parts:

# PART A: construct IPMs for two covaratiates (time since fire and grazing) used in vital-rate regressions 
# PART B: run simulations to calculate the stochastic population growth rate and elasticities
# PART C: plot some results

### PART A - COnstruct IPMs (much of the analyses and code are based on Ellner&Rees. 2006. Amer. Natur. 167: 410-428)

# Continuous transitions (survival, growth, and reproduction) are are obtained from
# regressions including a state variable (size) and two categorical variables (TSF and LS)
# Time since fire (TSF) here has five levels: 0, 1, 2, 3, and > 3 years since fire
# Grazing by large livestock (LS) has two levels: HG - high grazing equivalent to LS+ in presentation; LG - low grazing equivalent to LS_ in presentation

# Seed-bank transitions are given as constants

#load necessary packages

library(Matrix)
library(plyr)
library(reshape2)
library(ggplot2)

# parameters needed for regressions

cont_param=data.frame(a0.surv=-2.5, bc.surv= 0.54, a1.surv.one=0.50, a1.surv.two=0.77, #survival
                 a1.surv.three=-0.90, a1.surv.four=-0.37, a2.surv.HG=0.76, a2.surv.LG=-0.76,
                 bcTSF.surv.one=0.27, bcTSF.surv.two=-0.086, bcTSF.surv.three=-0.10, bcTSF.surv.four=-0.08,
                 bcD.surv.HG=-0.18, bcD.surv.LG=0.18, a1a2.surv.one.HG=-0.46, a1a2.surv.one.LG=0.46,
                 a1a2.surv.two.HG=-0.65, a1a2.surv.two.LG =0.65, a1a2.surv.three.HG=0.54,  a1a2.surv.three.LG=-0.54, 
                a1a2.surv.four.HG=0.58,  a1a2.surv.four.LG=-0.58, 
                #growth mean
                a0.gr= 2.33, bc.gr=0.62, a1.gr.one=0.72, a1.gr.two=-0.26, a1.gr.three=-0.29, 
                a1.gr.four=-0.17, a2.gr.HG=0.039, a2.gr.LG=-0.039, bcTSF.gr.one=-0.03,
                bcTSF.gr.two=-0.038, bcTSF.gr.three=0.063, bcTSF.gr.four=0.009, bcD.gr.HG=-0.058, bcD.gr.LG=0.058, 
                a1a2.gr.one.HG=-0.31, a1a2.gr.one.LG=0.31, a1a2.gr.two.HG =-0.06,
                a1a2.gr.two.LG=0.06, a1a2.gr.three.HG=0.23, a1a2.gr.three.LG=-0.23, a1a2.gr.four.HG=0.15,
                a1a2.gr.four.LG=-0.15, 
                #growth variance
                a0.gr.sd=-0.79, a1.gr.sd.one=-0.25, a1.gr.sd.two=0.38, a1.gr.sd.three=-0.12, a1.gr.sd.four=-0.01, a2.gr.sd.HG= 0.19, a2.gr.sd.LG=-0.19, 
                a1a2.gr.sd.one.HG=-0.59,   a1a2.gr.sd.one.LG=0.59,   a1a2.gr.sd.two.HG=0.7,   a1a2.gr.sd.two.LG=-0.7, 
                a1a2.gr.sd.three.HG=-0.04, a1a2.gr.sd.three.LG=0.04, a1a2.gr.sd.four.HG=-0.08, a1a2.gr.sd.four.LG=0.08, 
                #seedling size mean 
                a0.sds=3.91, a1.sds.one=-0.3, a1.sds.two=0.19, a1.sds.three=-0.16, 
                a1.sds.four=0.26, a2.sds.HG=0.08, a2.sds.LG=-0.08, a1a2.sds.one.HG=0.3,
                a1a2.sds.one.LG=-0.3, a1a2.sds.two.HG=-0.28, a1a2.sds.two.LG=0.28, a1a2.sds.three.HG=-0.28,
                a1a2.sds.three.LG=0.28, a1a2.sds.four.HG=0.26, a1a2.sds.four.LG=-0.26, 
                #seedling size variance
                a0.sds.sd=-1.13, 
                a1.sds.sd.one=0.12, a1.sds.sd.two=-0.19, a1.sds.sd.three=0.46, a1.sds.sd.four=-0.39, 
                #probability of flowering
                a0.fl=-9.1, bc.fl= 1.53, 
                a1.fl.two=-1.25, a1.fl.three=2.87, a1.fl.four=-1.6, a2.fl.HG=-0.74,
                a2.fl.LG=0.74, bcTSF.fl.two=0.057, bcTSF.fl.three=-0.38, bcTSF.fl.four=0.32,
                bcD.fl.HG=0.12, bcD.fl.LG=-0.12, a1a2.fl.two.HG=0.31, a1a2.fl.two.LG=-0.31,
                a1a2.fl.three.HG=0.58, a1a2.fl.three.LG=-0.58, a1a2.fl.four.HG=-0.89, a1a2.fl.four.LG=0.89, 
                #number of flowering stalks
                a0.fs=-2.2, bc.fs=0.42, a1.fs.two=0.32, a1.fs.three=-0.2,
                a1.fs.four=-0.12, a2.fs.HG=1.1, a2.fs.LG=-1.1, bcTSF.fs.two=-0.1,
                bcTSF.fs.three=0.07, bcTSF.fs.four=0.03, bcD.fs.HG=-0.2, bcD.fs.LG=0.2, 
                a1a2.fs.two.HG=0.12, a1a2.fs.two.LG=-0.12, a1a2.fs.three.HG=-0.09, a1a2.fs.three.LG=0.09,
                a1a2.fs.four.HG=-0.03, a1a2.fs.four.LG=0.03, 
                #number of flowers per stalk
                a0.fps=-0.007, bc.fps=0.24,
                a1.fps.two=0.22, a1.fps.three=-0.09, a1.fps.four=-0.13, a2.fps.HG=0.026, a2.fps.LG=-0.026)

attach(cont_param)
# Create density or probability distributions of the vital rates based on the parameters obtained from the Bayesian (simple) model:

# A. Continuous stages

# SURVIVAL:

S.fun <- function(z,fire,grazing) {
  
  mu.surv=exp(a0.surv+get(paste("a1.surv.",fire,sep=""))+get(paste("a2.surv.",grazing,sep=""))
              +get(paste("a1a2.surv.",fire,".",grazing,sep=""))
              +(bc.surv+get(paste("bcTSF.surv.",fire,sep=""))+get(paste("bcD.surv.",grazing,sep="")))*z)
  
  return(mu.surv/(1+mu.surv))
}

# GROWTH 

GR.fun <- function(z,zz,fire,grazing){
  
  # mean
  growth.mu = (a0.gr+get(paste("a1.gr.",fire,sep=""))+get(paste("a2.gr.",grazing,sep=""))
               +get(paste("a1a2.gr.",fire,".",grazing,sep=""))+(bc.gr+get(paste("bcTSF.gr.",fire,sep=""))+get(paste("bcD.gr.",grazing,sep="")))*z)
  
  #variance (parameters were obtained on log scale)
  se2 = exp(a0.gr.sd+get(paste("a1.gr.sd.",fire,sep=""))+get(paste("a2.gr.sd.",grazing,sep=""))
            +get(paste("a1a2.gr.sd.",fire,".",grazing,sep="")))
  
  #   se2[se2<0] = 0.1
  # Density distribution function of the normal distribution
  gr1 = sqrt(2*pi*se2)
  gr2 = ((zz-growth.mu)^2)/(2*se2)
  
  density=exp(-gr2)/gr1
  
  
  return(density)
  
}

## SEEDLING SIZES (same approach as in growth function)

SDS.fun <- function(z,zz,fire,grazing){
  
  # mean
  sds.mu=(a0.sds+get(paste("a1.sds.",fire,sep=""))+get(paste("a2.sds.",grazing,sep=""))
          +get(paste("a1a2.sds.",fire,".",grazing,sep="")))
  
  # Variance (parameters were obtained on log scale)
  
  se2.sds = exp(a0.sds.sd+get(paste("a1.sds.sd.",fire,sep="")))
  
  # Density distribution function of the normal distribution
  sds1 = sqrt(2*pi*se2.sds)
  sds2 = ((zz-sds.mu)^2)/(2*se2.sds)
  
  density.sds=exp(-sds2)/sds1
  
  return(density.sds)
  
}

# PROBABILITY OF FLOWERING 

FL.fun <- function(z,fire,grazing) {
  
  mu.fl=exp(a0.fl+get(paste("a1.fl.",fire,sep=""))+get(paste("a2.fl.",grazing,sep=""))
            +get(paste("a1a2.fl.",fire,".",grazing,sep=""))+
              (bc.fl+get(paste("bcTSF.fl.",fire,sep=""))+get(paste("bcD.fl.",grazing,sep="")))*z)
  
  return(mu.fl/(1+mu.fl))
}

# NUMBER OF FLOWERING STALKS 

FS.fun <- function(z,fire,grazing) {
  
  mu.fs=exp(a0.fs+get(paste("a1.fs.",fire,sep=""))+get(paste("a2.fs.",grazing,sep=""))+
              get(paste("a1a2.fs.",fire,".",grazing,sep="")) +
              (bc.fs+get(paste("bcTSF.fs.",fire,sep=""))+get(paste("bcD.fs.",grazing,sep="")))*z)
  
  return(mu.fs)
}

# NUMBER OF FLOWERS PER STALK

FPS.fun <- function(z,fire,grazing) {
  
  mu.fps=exp(a0.fps+get(paste("a1.fps.",fire,sep=""))+get(paste("a2.fps.",grazing,sep=""))+bc.fps*z)
  
  return(mu.fps)
}

################

# B. Discrete stages

SB_param=data.frame(TSF=rep(c("zero","one","two","three","four"),2),
                    LS=rep(c("HG","LG"),each=5),
                    corr=c(1,1,0.5,0.5,0.6,1,1,0.8,0.8,0.8),
                    goCont=c(0.01,0.01,0.05,0.05,0.065,0.01,0.01,0.002,0.0002,0.0001),
                    goSB=c(0,0,0.87,0.87,0.87,0,0,0.96,0.96,0.96),
                    staySB=c(0.1,0.05,0.6,0.6,0.6,0.1,0.05,0.6,0.85,0.867),
                    outSB=c(0.43,0.06,0.06,0.06,0.06,0.68,0.06,0.03,0.03,0.009))
seeds=9.8


# Function to make IPMs

minsize=0.9*0.01 # minimum size
maxsize=1.1*8.7 # maximum size

IPMkernel<-function(n) { # n defines the size (rows and columns) of the kernel resulting from integrating the IPM 
  
  # Function to define the IPM kernel (the midpoints for the integration)
  b <- minsize+c(0:n)*(maxsize-minsize)/n # interval that each cell of the matrix covers 
  h <- 0.5*(b[1:n]+b[2:(n+1)]) # midpoint
  
 
  # The zero (right after fire) matrix
  if(fire=="zero"){
    goCont=0.01 # immediate germination (growing season following seed production). Should be 0, but we don't want to divide by 0, see below
    goSB=0 # seed-bank ingression 
    outSB <- SB_param$outSB[SB_param$TSF==fire&SB_param$LS==grazing] # seed-bank egression
    staySB <- 0.1 # seed-bank stasis
    corr <- 1 # correction for seed mortality 
    S <-matrix(0,n,n) # survival (S), growth (G), and fecundity (FecALL) are all zero
    G <- matrix(0,n,n)
    FecALL=matrix(0,n,n)
    R <- (t(outer(h,h,SDS.fun,fire="one",grazing)))# the relevant non-0 transition is the size distribution of seedlings
    
    # TSF 1
  }else if(fire=="one"){
    goCont=0.01 
    goSB=0
    outSB <- SB_param$outSB[SB_param$TSF==fire&SB_param$LS==grazing]
    staySB <-0.05
    corr <- 1
    S <- diag(S.fun(h,fire,grazing)) # Survival Matrix 
    G <- t(outer(h,h,GR.fun,fire,grazing)) # Growth Matrix
    #Recruits distribution
    R <- (t(outer(h,h,SDS.fun,fire,grazing)))
    FecALL=matrix(0,n,n)# no seeds produced, therefore 0 fecundity
    
    # scale G and R below so columns sum to 1
    G <- G/matrix(as.vector(apply(G,2,sum)),nrow=n,ncol=n,byrow=TRUE)
    
    
  }else{
    goCont=SB_param$goCont[SB_param$TSF==fire&SB_param$LS==grazing]
    goSB=SB_param$goSB[SB_param$TSF==fire&SB_param$LS==grazing]
    outSB=SB_param$outSB[SB_param$TSF==fire&SB_param$LS==grazing]
    staySB=SB_param$staySB[SB_param$TSF==fire&SB_param$LS==grazing]
    
    corr=SB_param$corr[SB_param$TSF==fire&SB_param$LS==grazing]
    S <- diag(S.fun(h,fire,grazing)) # Survival Matrix 
    G <- t(outer(h,h,GR.fun,fire,grazing)) # Growth Matrix
    
    #Recruits distribution
    R <- (t(outer(h,h,SDS.fun,fire,grazing)))
    
    #Probability of flowering
    Fec01 = (diag(FL.fun(h,fire,grazing)))
    
    #Number of Flowering Stalks 
    Fec02 = (diag(FS.fun(h,fire,grazing)))
    
    #Number of flowers per Stalk
    Fec03= (diag(FPS.fun(h,fire,grazing)))
    
    #Number of seeds per flower that survive to become offspring 
    Fec04 = (diag(rep(9.8,n)))
    
    
    FecALL= Fec01*Fec02*Fec03*Fec04*goCont*corr # add goCont and corr to fecundity
    
    # scale R and G so columns sum to 1
    G <- G/matrix(as.vector(apply(G,2,sum)),nrow=n,ncol=n,byrow=TRUE)
    
  }
  
  R <- R/matrix(as.vector(apply(R,2,sum)),nrow=n,ncol=n,byrow=TRUE)
  
  return(list(S=S,G=G,FecALL=FecALL,R=R,meshpts=h,goCont=goCont,corr=corr,goSB=goSB,outSB=outSB,staySB=staySB))
  
}

#####################
###################### MAKE IPMS

#Loop though the LS and TSF categories
LS.name=c("HG","LG")
TSF.name=c("zero","one","two","three","four")

# choose how many bins the kernel should have (here 20 (+1 for seed bank) to speed up computation)
bins=21

#Create 2 empty arrays to store IPMs for LS+ (HG) and LS-(LG)
IPM.HG=array(0,c(bins,bins,length(TSF.name)))
IPM.LG=array(0,c(bins,bins,length(TSF.name)))

count=0
for(d in 1:length(LS.name)){
  
  grazing =  LS.name[d]
  
  for(t in 1:5){
     
    fire=TSF.name[t]

    M=IPMkernel(bins-1)
    
    ### Create P and F matrices including discrete stages
    
    Pmat.cont <- as.matrix(Matrix(M$G)%*%Matrix(M$S)) # continuous part 
    Pmat.discr = c(M$staySB,M$outSB*M$R[,2]) # discrete part
    Pmat = cbind(Pmat.discr,rbind(rep(0,length(M$meshpts)),Pmat.cont)) # combine continuous and discrete
    
    Fmat.cont <- as.matrix(Matrix(M$R)%*%Matrix(M$FecALL))# continuous part 
    Fmat.discr =rep(0,length(M$meshpts)+1)# discrete part
    Fmat=cbind(Fmat.discr,rbind(diag(M$FecALL)*M$goSB/(M$goCont),Fmat.cont))# combine
    
    mat <-Pmat+Fmat
    
    if(grazing=="HG"){IPM.HG[,,t]=mat}else{IPM.LG[,,t]=mat}
    
  }
  
}
##################################################################
#############################################################

#PART B - run stochastic simulations 
# (much of the aanalyses and code are based on Tuljapurkar et al. 2003. Amer. Natur.162:489-502
# & Horvitz et al. 2010.J of Ecol. 98:268-278)

# Vectorizing function

vecc1 <- function(x){
  d <- dim(x)
  vecc <- x[,1]
  for (i in 2:d[1]){
    vecc <- c(vecc,x[,i])}
  return(vecc)}


### Add simulations 

# Simulation function (just one of many possibilities)
# This function takes a environmental transition matrix (trans) and the number of simulation (N years)
# and returns a character vector of the sequence of environments in N years 

simula <- function(trans,N) {
  # This function tells you sample the rownames (i.e., environment at t+1), each having a probability = column name that you got at previous iteration(environment at t) 
  state.at.N <- function(char,trans) {
    sample(rownames(trans),1,prob=trans[,char])
  }
  
  sim <- character(N)
  sim[1] <- "1" # start with 0 matrix and populate 1st character with "1"
  
  # Populate characters 2-N with the row names samples obtained from (state.at.N)
  for (b in 2:N) {
    sim[b] <- state.at.N(sim[b-1],trans)
  }
  
  sim
}



### Simulate environments

# Number of simulations (here only 20 to save computational time)
M.num = 20 

# Define frequency of fire 

freqarray=c(0.04,0.03,0.02,0.01) 
# run each simulation for x years  
run=2002

# Define simulation time
tr=200 #the discard time
ts=1600 # total time of simulations
t2 <- ts + 2 * tr # time counte needed for the backwards loop for left eigenvectors

# Set some values needed to for further calculations 

n2=bins^2
n=bins

# The different fire environments 
env=TSF.name
env.n=length(env)


nstage <- bins # number of stages
nstate <- length(TSF.name) # number of environments 


# Define arrays to hold stochastic lambda and elastiticties in 
lambda.s=array(0,c(length(freqarray),M.num,length(LS.name))) # stochastic lambda (mean) for each fire frequency, simulation, and disturbance category 

elasE=array(0,c(nstage,nstage,length(freqarray),M.num,length(LS.name))) # "the" stochaStic elasticity
elasA=array(0,c(nstage,nstage,length(freqarray),M.num,length(LS.name))) # stoch elasticity to mean (across environments) of transitions
elasV=array(0,c(nstage,nstage,length(freqarray),M.num,length(LS.name))) # stoch elasticity to variance (across environments) of transitions

Pqt=array(0,c(length(freqarray),M.num,length(LS.name))) #quasiextinction after t=300 years

### THIS MAIN LOOP GOES OVER M.NUM SIMULATIONS (PICKs MATRIX FOR EACH FIRE ENVIRONMENT) 
#### +
##  RUN SIMULATIONS (GO OVER T=RUN YEARS FOR EACH M.NUM PICK) 

for(dist in 1:length(LS.name)){
  
  for(si in 1:M.num){
    
    #Choose GRAZED or UNGRAZED
    if(LS.name[dist]=="HG"){ms=IPM.HG}else{ms=IPM.LG}
    
    
    # Define arrays to hold elastiticties (T)  temporarily in 
    lambda.s.temp=rep(NA,length(freqarray)) # save lambda 
    Pqt.temp=rep(NA,length(freqarray)) # save lambda 
   
    elasE.temp <- array(0,c(nstage,nstage,length(freqarray)))
    elasA.temp <- array(0,c(nstage,nstage,length(freqarray)))  
    elasV.temp <- array(0,c(nstage,nstage,length(freqarray))) 
   
    #Loop trough fire frequency
    for(f in 1: length(freqarray)){
      
      fire=freqarray[f]; 
      
      #environmental transition matrix based on fire frequency
      P=matrix(c(fire,1-fire,rep(0,3),
                 fire,0,1-fire,rep(0,2),
                 fire,0,0,1-fire,0,
                 fire,0,0,0,1-fire,
                 fire,0,0,0,1-fire),nrow=5,ncol=5,byrow=F)
      colnames(P) <- c("1","2","3", "4","5")
      row.names(P) <- c("1","2","3", "4","5")
      
      
      ## Simulations 
      simulate=simula(P,run)

      ### average transtions weighted by dominant eigenvalue of transitison matrix P 
      cmat_array=P
      
      vec_c <- eigen(cmat_array) 
      cmax  <- as.numeric(max(abs(vec_c$values)))
      loc   <- which(vec_c$values==max(abs(vec_c$values)))
      
      # The right eingevector of the environmental matrix depicts the ralative frequencies
      # of evironmental states over time for a single or (or, equivalenty, for a multitude
      # of independetly variable sites that change according to the Markov chain)
      cstar <- as.numeric(vec_c$vec[,loc])  
      cstar <- cstar/sum(cstar)
      
      #Cacluate the average matrix 
      Abar <- array(0,nstage^2)
      
      #for each environment
      
      for (i in 1:nstate) {
        
        extractmat <-  ms[,,i]
        A <- vecc1(extractmat)  
        Abar <- Abar+cstar[i]*A # this is the mean matrix averaged over environments (ctar[i]*A) 
      }  
 
      avmat <- matrix(Abar,nstage,nstage) # Matrix of average vital rates

      ## Stochastic Elasticity
      
      states <- as.numeric(simulate)
     
      
      # Save the series of:
      growth <- array(0,ts)   # lambda
      uvecs  <- array(0,c(nstage,ts)) # right eigenvectors
      vvecs  <- array(0,c(nstage,ts)) # left eigenvectors
     
      # Initial population vectors 
      vec1=c(1000,rep(1,nstage-1))
      vec1 <- vec1/sum(vec1)
      vec1 <- t(vec1) 
      vec2 <- vec1
      vec3 <- c(1000,rep(1,nstage-1)) # to monitor actual numbers 
      
      trun <- ts+tr #overall runtime (total simulations+discard)
      
      # ITERATION TO CALCULATE LAMBDA AND V AND W FOR EACH TIME STEP
      for (i  in 1:trun){
        i2 <- states[i]
        mat1 <-  ms[,,i2]
        vec1 <- mat1%*%as.numeric(vec1)
        vec3 <- mat1%*%as.numeric(vec3)
        #save population vector at t=300
        if(i==300){Pqt.temp[f]=sum(vec3)}
        growth1 <- sum(vec1)
        vec1 <- vec1/growth1
        #Forward iteration
        if( i > tr){      
          i1 <- i - tr
          uvecs[,i1] <- vec1
          growth[i1] <- growth1 # growth rate at each iteration after discard time
        }
        #Backward iteration
        j <- (t2 - i+1)
        i2 <- (states[j+1])
        mat1 <- ms[,,i2]
        vec2 <- t(as.numeric(vec2)%*%mat1)
        vec2 <- vec2/(sum(vec2))
        if (i > tr)  {     
          vvecs[, j-tr] <- vec2
        }
      }
      
   
      lambda.s.temp[f]=sum(log(growth[1:ts])) # save cumulative growth rate
      
      ####
      
      ### CALCULATE ELASTICITIES TO ELEMENTS (elasE), MEAN (elasA), AND VARIANCE (elasV)
      
      for (i in (tr+1):(tr+ts-1)){ #start after a time lag (tr)
        
        i2 <- (states[i+1])
        mat1 <- ms[,,i2]    
        itime <- i+1-tr
        i1 <- i-tr;      
        
        mat2S <- mat1 
        mat2A <- avmat
        mat2V <- mat1 - avmat
        
        # calculate elasticties (numerator of equation)
        mat2S <- diag(vvecs[,itime])%*%mat2S%*%diag(uvecs[,i1]) 
        mat2A <- diag(vvecs[,itime])%*%mat2A%*%diag(uvecs[,i1])
        mat2V <- diag(vvecs[,itime])%*%mat2V%*%diag(uvecs[,i1])
        
        # Scalar product of left and right eigennvectors 
        scale1 <-(t(vvecs[,itime]))%*%(uvecs[,itime])
        
        #Denominator part 
        mat2S <- as.numeric((1/(scale1*growth[itime])))*mat2S                  
        elasE.temp[,,f] <- elasE.temp[,,f] + mat2S
        
        mat2A <- as.numeric((1/(scale1*growth[itime])))*mat2A
        elasA.temp[,,f] <- elasA.temp[,,f] + mat2A
        
        mat2V <- as.numeric((1/(scale1*growth[itime])))*mat2V
        elasV.temp[,,f] = elasV.temp[,,f] + mat2V
      }
      # Time-averaged elasitivities 
      elasE.temp[,,f] <- elasE.temp[,,f]/(ts-1)
      elasA.temp[,,f] <- elasA.temp[,,f]/(ts-1)
      elasV.temp[,,f] <- elasV.temp[,,f]/(ts-1)
    }
    
    
    lambda.s[,si,dist]=lambda.s.temp
    Pqt[,si,dist]=Pqt.temp
    elasE[,,,si,dist]=elasE.temp
    elasA[,,,si,dist]=elasA.temp
    elasV[,,,si,dist]=elasV.temp

 
  }
}

##########################################################
###################################################
##### PART C Plot results

############## STOCHASTIC LAMBDA

logL=adply(lambda.s,c(1,2,3))
colnames(logL)=c("firefreq","si","grazing","lambda")

logL$firefreq=factor(logL$firefreq)

levels(logL$firefreq)=c("25","33","50","100")
logL$lambda=logL$lambda/ts
logL$si=factor(logL$si)
logL$grazing=factor(logL$grazing)
levels(logL$grazing)=c( expression(LS[+phantom()]), expression(LS[-phantom()]))

# Confidence intervals around estimates
max.l=aggregate(logL$lambda, by=list(logL$firefreq,logL$grazing),max)

sig=aggregate(logL$lambda, by=list(logL$firefreq,logL$grazing),quantile,0.975)
colnames(sig)=c("firefreq","grazing","p")
sig$lab=""
sig$lab[sig$p<(0)]="*"
sig$lambda=max.l[,"x"]+0.015


ggplot(logL, aes(x=firefreq, y=lambda,fill="firefreq")) + geom_boxplot(fill="grey") +
  guides(fill=FALSE)+ 
  stat_summary(fun.y=mean, geom="point", shape=5, size=5)+
  geom_hline(aes(yintercept=0),   
             color="black", linetype="dashed", size=1)+
  facet_grid(. ~ grazing,labeller=label_parsed)+
  theme_bw()+
  xlab("Fire return (years)") +
  ylab(expression(paste("log ",lambda[s],sep="")))+
  theme(panel.grid = element_blank())+
  theme(axis.text = element_text(size=18))+
  theme(axis.title = element_text(size=24))+
  theme(axis.title.x = element_text(vjust=-0.4),
        axis.title.y = element_text(vjust=0.8))+
  theme(strip.text.x = element_text(size = 22),
        strip.background =element_blank(),
        strip.text.y = element_blank())+
  geom_text(aes(x=firefreq, y=lambda, label=lab),
            data=sig,size=8)

############## PROBABILITY OF QUASI-EXTINCTION

# If summed population vector at t=300 < 100, Pqt=1, else it is 0
ext=apply(Pqt,c(1,2,3),function(x)ifelse(x>100,x<-0,x<-1))

extSimul=adply(ext,c(1,2,3))
colnames(extSimul)=c("firefreq","si","grazing","Pqt")
extSimul$si=factor(extSimul$si)
extSimul$grazing=factor(extSimul$grazing)
levels(extSimul$grazing)=c("HG","LG")

# calculate means and standard error across simulations 

mean.ext300=aggregate(extSimul$Pqt, by=list(extSimul$firefreq,extSimul$grazing),mean)
mean.ext300$se=aggregate(extSimul$Pqt, by=list(extSimul$firefreq,extSimul$grazing),function(x)(sqrt(mean(x)*(1-mean(x)))/sqrt(length(x))))[,"x"]
colnames(mean.ext300)=c("firefreq","grazing","ext","se")
mean.ext300$firefreq=as.numeric(mean.ext300$firefreq)


ggplot(mean.ext300, aes(x=firefreq, y=ext,colour=grazing))  +
  geom_line(size=1.2) + 
  geom_point(size=4)+
  geom_errorbar(data=mean.ext300,aes(ymin=ext-se, ymax=ext+se,color=grazing),width=.4,size=1)+
  guides(fill=FALSE)+ 
  theme_bw()+
  xlab("Fire return (years)") +
  ylab(expression(paste(P[q],"(300 years)",sep="")))+
  theme(panel.grid = element_blank())+
  theme(axis.text = element_text(size=18))+
  theme(axis.title = element_text(size=24))+
  theme(axis.title.x = element_text(vjust=-0.4),
        axis.title.y = element_text(vjust=0.8))+
  scale_x_continuous(breaks=seq(1,4),labels=c("25","33","50","100"))+
  scale_color_manual(name="",values=c("#009900","red"),breaks=c("HG", "LG"),
                     labels=c(expression(LS[+phantom()]),expression(LS[-phantom()])))+
  theme(legend.justification=c(0,1), legend.position=c(0,1),
        legend.title = element_text(size=16, face="bold"),
        legend.text = element_text(size=24))

############## STOCHASTIC ELASTICITIES
#Calculate the mean across simulations (see Trotter et al. 2013. MEE, )

m.elasE=apply(elasE,c(1,2,3,5),mean) 
m.elasA=apply(elasA,c(1,2,3,5),mean) 
m.elasV=apply(elasV,c(1,2,3,5),mean) 

m.lambda.s=apply(lambda.s,c(1,3),sum) 
m.lambda.s=m.lambda.s/c(M.num*ts)

#Variance of stochastic growth rate
m.lambda.v=array(0,c(length(freqarray),length(LS.name))) 

for(dist in 1:length(LS.name)){
  
  for(f in 1:length(freqarray)){
    
    for(i in 1:M.num){
      
      m.lambda.v[f,dist] = m.lambda.v[f,dist] + (lambda.s[f,i,dist] - m.lambda.s[f,dist]*ts)^2 
    }
    m.lambda.v[f,dist]=m.lambda.v[f,dist]/(M.num*ts)
    
  }
}

### Plots (example m.elasA - but can be euqally plotted for other elasticities)

#Get midpoints from the kernel

b <- minsize+c(0:n)*(maxsize-minsize)/(bins-1) 
h <- 0.5*(b[1:n]+b[2:(n+1)]) 
y=h # for axis labels in the image plot 

names=c("25","33","50","100")
set.panel()

par(mfrow=c(2,4),mar=c(1.3,1.3,1,1),oma=c(3,3,2,4)) 

set.panel( 2,4) 

count=0
#draw all your plots using usual image command
for(dist in 1:2){
  for(f in 1:length(freqarray)){
    count=count+1
    #Plot ELASTICTIES 
    
    #function to cut upper and lower limits of elasticity values to a desired ranged for easier plotting
    # note that for elasticities of log lambda_s to mean transitions, 0 is already the lower limits - we set a limit here just for demonstration
    conv=function(x){ifelse(x>0.015,x<-0.015,ifelse(0>x,x<-(0),x))}
    # #     
    V=apply(m.elasA[,,f,dist],c(1,2),conv)
    
    if(count==1){
      image(c(-1,y),c(-0.55,y),t(V), xlab="",ylab="size t+1",cex.axis=1.5,cex.lab=2,xaxt="n",zlim=c(0,0.015),col=tim.colors(900), xpd="n")
      text(1.5,8.6, expression(LS[+phantom()]),cex=2.5,col="green")
      abline(h=0.24,col="white")
      abline(v=0.2,col="white")
      title( main=paste("Fire return: ",names[f]," years", sep=""), line =1.1,cex.main=2.5,xpd="n")
      
    } else if(count==2|count==3|count==4){
      
      image(c(-1,y),c(-0.55,y),t(V), xlab="",ylab="",xaxt="n",yaxt="n",zlim=c(0,0.015),col=tim.colors(900))
      abline(h=0.24,col="white")
      abline(v=0.2,col="white")
      title( main=paste(names[f]," years", sep=""), line =1.1,cex.main=2.5,xpd="n")
      
    }else if(count==5){
      
      image(c(-1,y),c(-0.55,y),t(V), xlab="size t",ylab="size t+1",cex.axis=1.5,cex.lab=2,zlim=c(0,0.015),col=tim.colors(900), main="",xpd="n")
      text(1.5,8.6, expression(LS[-phantom()]),cex=2.5,col="red")
      abline(h=0.24,col="white")
      abline(v=0.2,col="white")
    }else{
      image(c(-1,y),c(-0.55,y),t(V), xlab="size t",ylab="",yaxt="n",cex.axis=1.5,cex.lab=2,zlim=c(0,0.015),col=tim.colors(900), main="",xpd="n")
      abline(h=0.24,col="white")
      abline(v=0.2,col="white")
    }
     
    
  }
}


par(oma=c( 0,0,0,2.5))# reset margin to be much smaller.
image.plot(legend.only=TRUE,zlim=c(0,0.015),col=tim.colors(900)) 

set.panel()

# Summed elasticities of log lambda_s to changes in mean and variance in vital rates:

#### Tmu
Tmu=NULL
for(dist in 1:2){
  for(f in 1:length(freqarray)){
    
    cont=sum(m.elasA[2:21,2:21,f,dist]) #continuous transitions
    sb.s=sum(m.elasA[1,1,f,dist]) # seed-bank stasis
    sb.out=sum(m.elasA[2:21,1,f,dist])# seed-bank egression
    sb.in=sum(m.elasA[1,2:21,f,dist])# seed-bank ingression
    fire=f
    LS=LS.name[dist]
    
    Tmu=rbind(Tmu,data.frame(cont,sb.s,sb.out,sb.in,fire,LS))
    
  }
}

Tvar=NULL
for(dist in 1:2){
  for(f in 1:length(freqarray)){
    
    cont=sum(m.elasV[2:21,2:21,f,dist])
    sb.s=sum(m.elasV[1,1,f,dist])
    sb.out=sum(m.elasV[2:21,1,f,dist])
    sb.in=sum(m.elasV[1,2:21,f,dist])
    fire=f
    LS=LS.name[dist]
    
    Tvar=rbind(Tvar,data.frame(cont,sb.s,sb.out,sb.in,fire,LS))
    
  }
}

rowSums(Tmu[,1:4])+rowSums(Tvar[,1:4]) #must sum to 1 

Tmu; Tvar

# Plot Tmu (same fot Tvar)
Tmu=melt(Tmu,id.vars=c("fire","LS"),measure.vars=c("cont","sb.s","sb.out","sb.in"))
Tmu$fire=as.factor(Tmu$fire)
levels(Tmu$fire)=c("Fire: 25 years","33 years","50 years","100 years")
Tmu$variable=as.character(Tmu$variable)
Tmu$variable=factor(Tmu$variable,levels=c("cont","sb.in","sb.s","sb.out"))

ggplot(Tmu, aes(x=variable, y=value)) +
  geom_bar(width = 0.7, aes(fill=LS), stat="identity",position= position_dodge(width=0.7)) +
  geom_bar(width = 0.7, aes(fill=LS),colour="black", stat="identity",position = position_dodge(width=0.7),show_guide=F) +
  facet_grid(. ~ fire)+
  theme_bw()+
  ylab(expression(T^mu)) +
  scale_fill_manual(values=c("#009900","red"), 
                    labels=c(expression(LS[+phantom()]),expression(LS[-phantom()])))+
  xlab("Life cycle components")+
  theme(legend.text.align = 0)+
  theme(panel.grid = element_blank())+
  theme(axis.text = element_text(size=16),
        axis.text.x = element_text(angle=45,hjust=1))+
  theme(axis.title = element_text(size=24))+
  theme(axis.title.x = element_text(vjust=-0.4),
        axis.title.y = element_text(vjust=0.8))+
  theme(strip.text.x = element_text(size = 22),
        strip.background =element_blank(),
        strip.text.y = element_blank())+ 
  scale_x_discrete(labels=c("Continuous","Ingression","Stasis","Egression"))+
  theme(legend.key = element_rect(colour = "black",size=1))+
  theme(legend.text = element_text(size=20),
        legend.title=element_blank())
