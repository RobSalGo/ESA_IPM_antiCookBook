## KhayaIPM: Effect of heterogeneous harvesting on population dynamics
## Study system: Khaya senegalensis (Meliaceae, Tree), in Benin (West Africa)
## Orou Gaoue | ogaoue@hawaii.edu, August 3, 2015
## Last revision: August 11, 2015 at 2 a.m.

## Context: Indigenous people harvest African mahogany to feed cattle in West Africa. They climb trees and cut branches (pruning) to feed cattles undereath the tree. This is a widespread activity and questions arise about how to harvest this system in a sustainable way. I tested the effect of varying the mean and/or the variance in individual tree level harvesting (pruning) rate on population dynamics. 

rm(list=ls(all=TRUE)) ## Clear all

library(popbio) ## for eigenanalysis

## setting working directory
setwd("/Users/ogaoue/Dropbox/Khaya projects/Khaya IPM/Script R/")

### PART (1) ############################
### DATA 1: Khaya data for adults (n=580)
KhayaA<-read.csv("khaya_adults04_07.csv", header=T)

str(KhayaA) ## Need to tranform pop as factor
KhayaA$pop<-as.factor(KhayaA$pop) ## pop as factor in KhayaA
str(KhayaA) ## checking transformations

## Organizing data to build Khaya adult domain   
Tp<-4
nind<-length(KhayaA[,1])

## Pruning intensity in 2006 is 0 for all indiv because that year and data collecting took place before harvesting (late pruning that year due to seasonal changes). As a result, I am using the 2005 harvest intensity for that year. Any other option may work as well (e.g., use mean harvest intensity across years).

KhayaA$prun06<-KhayaA$prun05

## Creating matrices for tag, architecture (arc) DBH,  pruning intensity (prun), debarking intensity (dbk), bark thickness (thk), fruit production (fr), year since pruning (ysp) 

tagmat<-as.matrix(KhayaA[,4:7])
arcmat<-as.matrix(KhayaA[,8:10])
dbhmat<-as.matrix(KhayaA[,11:14])
prunmat<-as.matrix(KhayaA[,15:18])
dbkmat<-as.matrix(KhayaA[,19:22])
thkmat<-as.matrix(KhayaA[,23:26])
frmat<-as.matrix(KhayaA[,27:30])
yspmat<-as.matrix(KhayaA[,33:36])

## Building size, sizenext, fruit, pruning, debarking
## Data is available from year 2004 to year 2007. So size will be from 2004-2006 and sizenext from 2005:2007. Let's create year column 

## year vector
year<-as.factor(c(rep(2004, nind), rep(2005, nind), rep(2006, nind)))

## region vector
region<-rep(KhayaA$reg,3)

## region vector
harvest<-rep(KhayaA$harvest,3)

## tag vector
tag<-as.factor(c(tagmat[,1:3]))

## pop vector
pop<-as.factor(c(rep(KhayaA$pop, 3)))

## size and sizenext vectors log-transformed
s0<-log(c(dbhmat[,1:3]))  # size at t
s1<-log(c(dbhmat[,2:4])) # size at t+1

## pruning vector
prun<-c(prunmat[,1:3])
prunnext<-c(prunmat[,2:4])

## debarking vector
dbk<-c(dbkmat[,1:3])
dbknext<-c(dbkmat[,2:4])

## fruit vector
fr<-c(frmat[,1:3])
frnext<-c(frmat[,2:4])

## Now creating a new dataframe 
dat<-data.frame(year, region, harvest, pop, tag, s0, s1, prun, dbk, fr)

## There are obs with NA for both s0 an s1 that needs to be removed
datA<-subset(dat, !(s0 %in% NA & s1 %in% NA))
## 101 observation has NA for both s0 and s1

## Adding survival to the data
datA$survA<-as.factor(1*(!is.na(datA$s1)))

## probability of fruiting
datA$pfr<-as.factor(1*(datA$fr>0)
                    
## Creating new Adult recruited
datA$recrA<-as.factor(1*(is.na(datA$s0) & !is.na(datA$s1)))

## Creating a unique pop identifier for mixed effect models
datA$popA<-with(datA, as.factor(paste(year, region, pop, sep=".")))

## checking data structure
str(datA) 

### DATA 2: Khaya data for seedlings (n=1396) ************************
KhayaS<-read.csv("khaya_seedlings04_07.csv", header=T)

str(KhayaS) # Transform pop as factor, tag06 as integer, surv as factor
KhayaS$pop<-as.factor(KhayaS$pop) ## pop as factor in KhayaS

for (i in 6:7){
  KhayaS[,i]<-as.integer(KhayaS[,i]) ## Tag05 & 06 as integer in KhayaS
}

for (i in 9:11){
  KhayaS[,i]<-as.factor(KhayaS[,i]) ## surv as factor in KhayaS
}

str(KhayaS) ## recheck

nindS<-length(KhayaS[,1])

## Creating matrices
tgmat<-as.matrix(KhayaS[,5:8])
survmat<-as.matrix(KhayaS[,9:11])
dbmat<-as.matrix(KhayaS[,12:15])
hgtmat<-as.matrix(KhayaS[,16:19])

## Creating vectors
yearS<-as.factor(c(rep(2004, nindS), rep(2005, nindS), rep(2006, nindS))) ## year vector
regS<-rep(KhayaS$region,3) ## region vector
harvS<-rep(KhayaS$harvest,3) ## harvest vector
popS<-as.factor(c(rep(KhayaS$pop, 3))) ## pop vector
tagS<-as.factor(c(tgmat[,1:3]))  ## tag vector
survS<-c(survmat[,1:3]) ## survival vector
db0<-log(c(dbmat[,1:3])) ## seedling size at t
db1<-log(c(dbmat[,2:4])) ## seedling size at t+1
ht0<-log(c(hgtmat[,1:3])) ## seedling height at t
ht1<-log(c(hgtmat[,2:4])) ## seedling height at t+1

## Now creating a new dataframe 
datS0<-data.frame(yearS, regS, harvS, popS, tagS, db0, db1, ht0, ht1, survS)
## There are obs with NA for both s0 an s1 that needs to be removed
datS<-subset(datS0, !(db0 %in% NA & db1 %in% NA))

## Creating new recruit
datS$recr<-as.factor(1*(is.na(datS$db0) & !is.na(datS$db1)))

## Creating a unique pop identifier for mixed effect models
datS$popSS<-with(datS, as.factor(paste(yearS, popS, sep=".")))


######## (PART 2): Fiting the statistical models  ############################

####### A. Survival (logistic regression for size and harvest)
###### A.1. Adults survival models
## Survival analysis with a mixed model, with population (popA) as a random effect
library(lme4)
surv1<-glmer(survA~s0*prun+(1|popA), family=binomial, data=datA)

s2<-update(surv1,~.-s0:prun)
s3<-update(s2,~.-prun)
s4<-update(s2,~.-s0)

anova(surv1, s2, s3, s4)
## model surv1 which includes prun is one of the best supported models

summary(surv1) 
## Negative effect of pruning on survival only apparent with size. Keep in mind here that harvest is size-selective and large individual are harvested
## plot(datA$prun~datA$s0, col=as.numeric(datA$harvest))

## Saving survival coefficients  ****************************
sint<-fixef(surv1)[1]  		## Adult survival intercept
ss0<-fixef(surv1)[2]			## slope for adult size
sprun<-fixef(surv1)[3]			## slope for pruning
ss0prun<-fixef(surv1)[4]   	## slope for interaction
#############################################################

##### A.2. Seedling survival models
## The difference here is that seedling are not pruned and that simplifies the model
Ssurv<-glmer(survS~db0+(1|popSS), family=binomial, data=datS)
summary(Ssurv)

## Saving survival coefficients  ****************************
Ssint<-fixef(Ssurv)[1]  		## intercept for seedling survival
Sdb0<-fixef(Ssurv)[2]			## slope for seedling size
############################################################

###### B. GROWTH
###### B.1. Adult survival models
## Do we have size dependent variance?
with(datA, plot(s1~s0, col=as.numeric(harvest)))
## I did not see a strong size dependent variance and decided to model growth here using mixed effect model without adding a variance component. And model below confirms it.

library(nlme)
gr1<-gls(s1~s0*prun, na.action=na.omit, method="ML", data=datA)
grp<-update(gr1, weight=varPower(form=~fitted(.)))
gre<-update(gr1, weight=varExp(form=~fitted(.)))

anova(gr1, grp, gre)
## Model gr1 is by far the best model suggesting no size dependent variance. there is also an effect of size and also pruning on future size

gr1<-update(gr1, method="REML")
summary(gr1)

## Saving growth coefficients ********************
gint<-summary(gr1)$coef[1]     	## intercept for adult growth
gs0<-summary(gr1)$coef[2]   		## slope for s0
gprun<-summary(gr1)$coef[3] 		## slope for prun
gs0prun<-summary(gr1)$coef[4]   ## slope for s0 x prun
gsig<-summary(gr1)$sigma^2 		  ## overall (redisual) variance
####################################################

###### B.1. Seedling survival models
## Do we have size dependent variance?
with(datS, plot(db1~db0, col=as.numeric(harvS)))
## Yes there seems to be an exponential size dependent variance. Let's test this

gS1<-gls(db1~db0, na.action=na.omit, method="ML", data=datS)
gSp<-update(gS1, weight=varPower(form=~fitted(.)))
gSe<-update(gS1, weight=varExp(form=~fitted(.)))

anova(gS1, gSp, gSe)
## Model gSe is best suppoted confirming exponential variance component

gSe<-update(gSe, method="REML")
summary(gSe)

## Saving growth coefficients ********************
gSint<-summary(gSe)$coef[1]       ## intercept for seedling growth
gSs0<-summary(gSe)$coef[2]   		  ## slope for db0
gSsig<-summary(gSe)$sigma^2 		  ## overall (redisual) variance
varExp<-gSe$modelStruct$varStruct  ## size-dependent variance
####################################################

###### C. FERTILITY 
#### C.1. Probability of fruiting by adults
prf1<-glmer(pfr~s0*prun+(1|popA), family=binomial, data=datA)
prf2<-update(prf1,~.-s0:prun)
prf3<-update(prf2,~.-prun)
prf4<-update(prf2,~.-s0)

anova(prf1, prf2, prf3, prf4) ## Model prf1 and 2 are best supported: Size and pruning effect
summary(prf1)

## Saving probability of fruiting coefficients  *************
prfint<-fixef(prf1)[1]    ## Adult survival intercept
prfs0<-fixef(prf1)[2]			## slope for adult size
prfprun<-fixef(prf1)[3]		## slope for pruning
prfs0prun<-fixef(prf1)[4] ## slope for interaction
#############################################################


#### C.2. Number of fruits produced by adults
plot(s0, fr, col=as.numeric(harvest), log="x")

## First try with poisson error stucture
nfr<-glmer(fr~s0*prun+(1|popA), family=poisson, data=datA)
summary(nfr) ## overdispersion deviance/df = 50788.4/1518. Using negative binomial error

nfr1<-glmer.nb(fr~s0*prun+(1|popA), data=datA)
nfr2<-glmer.nb(fr~s0+prun+(1|popA), data=datA)
nfr3<-glmer.nb(fr~s0+(1|popA), data=datA)
nfr4<-glmer.nb(fr~prun+(1|popA), data=datA)

anova(nfr1, nfr2, nfr3, nfr4) ## nfr1 is one of the best supported: Size and pruning effect
summary(nfr1)

## Saving fruit production coefficients  ********************
fint<-fixef(nfr1)[1]    ## Adult survival intercept
fs0<-fixef(nfr1)[2]  		## slope for adult size
fprun<-fixef(nfr1)[3]		## slope for pruning
fs0prun<-fixef(nfr1)[4] ## slope for interaction
#############################################################


##### C.3. Fruit/new seedlings ratio "p.est". 
## We used the total number of fruit produced in year t (tfr), and the number of new seedling produced in year t+1 (nsdl) as total number of new recruit from data datS.

tfr<-sum(datA$fr, na.rm=T)
nsdl<-length(subset(datS, datS$recr=="1")[,1])

## Fruit/new seedling ration p.est ***********
p.est<-nsdl/tfr   ## p.est = 0.01796096
##########################################

###### C.4. Size distribution of new seedlings 
## Size distribution (histogram) for new seedlings.
with(datS, hist(db1[recr=="1"]))

## Mean and standard deviation of new seedling size
mu.sdl<-with(datS, mean(db1[recr=="1"], na.rm=T)) ## *********
sig.sdl<-with(datS, sqrt(var(db1[recr=="1"], na.rm=T)))


#### PART (3) ###########################
### Matrix of all the regression coefficients

p.vec<-matrix(0,26,1) ## 26 coefficients

p.vec[1]<-sint      ## intercept for adult survival 
p.vec[2]<-ss0   		## size slope for adult survivak
p.vec[3]<-sprun			## pruning slope for adult survival
p.vec[4]<-ss0prun   ## s0 x prun interaction slope for adult survival

p.vec[5]<-Ssint    	## intercept for seedling survival
p.vec[6]<-Sdb0			## size db0 slope for seedling survival

p.vec[7]<-gint       ## intercept for adult growth
p.vec[8]<-gs0 		   ## size s0 slope for adult growth
p.vec[9]<-gprun 		 ## pruning slope for adult growth
p.vec[10]<-gs0prun   ## s0 x prun interaction slope for adult growth
p.vec[11]<-gsig 		 ## overall (redisual) variance for adult growth

p.vec[12]<-gSint      ## intercept for seedling growth
p.vec[13]<-gSs0    	  ## size db0 slope for seedling growth
p.vec[14]<-gSsig 		  ## overall (residual) variance for seedling growth
p.vec[15]<-varExp     ## size-dependent variance for seedling growth

p.vec[16]<-prfint     ## intercept for probability of fruiting
p.vec[17]<-prfs0  		## size s0 slope for prob of fruiting
p.vec[18]<-prfprun		## prun slope for prob of fruiting
p.vec[19]<-prfs0prun  ## s0 x prun interaction slope for prob of fruiting

p.vec[20]<-fint     ## intercept for fruit production
p.vec[21]<-fs0  	  ## size s0 slope for fruit production
p.vec[22]<-fprun		## prun slope for fruit production
p.vec[23]<-fs0prun  ## s0 x prun interaction slope for fruit prod

p.vec[24]<-p.est    ## fruit/new seedling ratio
p.vec[25]<-mu.sdl   ## mean size for new seedling
p.vec[26]<-sig.sdl  ## stdev size for new seedling


##### PART (4) ##########################
## The IPM FUNCTIONS: s(x), g(y,x), f(y,x), c(y,x), p(y,x), K(y,x)

### A. SURVIVAL FUNCTIONS
### A1. SURVIVAL function for adult s(x)
sx<-function(x, h, pvec){
	xbeta<-pvec[1] + pvec[2]*x + pvec[3]*h + pvec[4]*h*x;
	s<-exp(xbeta)/(1+exp(xbeta))
	return(s);
}

### A2. SURVIVAL function for seedling s(x)
sx.S<-function(x, pvec){
  xbeta<-pvec[5] + pvec[6]*x;
  s<-exp(xbeta)/(1+exp(xbeta))
  return(s);
}

### B. GROWTH function g(y,x)
### B1. Growth for adults
gyx<-function(y, x, h, pvec){
	mux<-pvec[7] + pvec[8]*x + pvec[9]*h + pvec[10]*h*x;
	sigmax<-sqrt(pvec[11])
	g<-dnorm(y, mux, sigmax)
	return(g)
}

### B2. Growth for seedling
gyx.S<-function(y, x, pvec){
  mux<-pvec[12] + pvec[13]*x;
  sigmax2<-pvec[14]*exp(2*pvec[15]*mux)
  sigmax<-sqrt(sigmax2)
  g<-dnorm(y, mux, sigmax)
  return(g)
}

### C. The SURVIVAL-GROWTH function P(y, x)
### C1. Survival-growth for adults
pyx<-function(y,x, h, pvec){ 
	p<-sx(x, h, pvec)*gyx(y, x, h, pvec)
	return(p) 
	}

### C2. Survival-growth for seedlings
pyx.S<-function(y, x, pvec){ 
  p<-sx.S(x, pvec)*gyx.S(y, x, pvec)
  return(p) 
}

### D. FERTILITY function f(x,y)

fyx<-function(y, x, h, pvec){
	xbeta<-pvec[16] + pvec[17]*x + pvec[18]*h + pvec[19]*h*x;
	probF<-exp(xbeta)/(1+exp(xbeta)) ## prob of fruiting
  NFr<-exp(pvec[20] + pvec[21]*x + pvec[22]*h + pvec[23]*h*x); ## log(mu)=a + bx + ch +dxh
	scd<-dnorm(y, pvec[25], pvec[26])  ## Size class dist of new seedling
	f<-probF*NFr*pvec[24]*scd
	return(f)
}

### E. The KERNEL: K(y,x)= p(y,x) + f(y,x) + c(y,x)
### E1. Kernel for adult domain
Kyx<-function(y, x, h, pvec){
	k<-pyx(y, x, h, pvec)+fyx(y, x, h, pvec)
	return(k) 
}

### E2. Kernel for seedling domain
Kyx.S<-function(y, x, pvec){
  ks<-pyx.S(y, x, pvec)
  return(ks) 
}

##### PART (5): Numerical integration #############################
### NUMERICAL INTEGRATION 

## The big matrix function

bigmat<-function(bigM, H, pvec){
	## Set matrix size for seedling + adult 
	min.sz<-0.9*min((range(datS$db0,na.rm=T))) ## min size for seedling
	max.sz<-1.1*max((range(datA$s1,na.rm=T)))  ## max size for adult
  
	# Compute meshpoints iteration matrix KD 
	h=(max.sz-min.sz)/(bigM+1); 
	y=(seq(min.sz, max.sz, length=bigM) + seq(min.sz+h, max.sz+h, length=bigM))/2;  
  ys=y[y<max((range(datS$db1,na.rm=T)))];  
	
  ## Adult Kyx function for y and y, H, pvec (=p.vec)
	K=outer(y,y, h=H, Kyx, pvec=p.vec);
	KD.A=h*K;
  
	## Seedling Kyx.S function for ys and ys, pvec (=p.vec)
	KS=outer(ys,ys,Kyx.S, pvec=p.vec);
	KD.S=h*KS;
  
  ## bolting to get the full kernel KD
  KD<-matrix(0, length(y), length(y))
  KD[1:length(y), 1:length(y)]<-KD.A  ## Adult kernel including fxy
	KD[1:length(ys), 1:length(ys)]<-KD.S  ## Seedling kernel 
	return(KD);  		  
}


##### PART (6): Eigen analysis ##########################

## Checks: Ploting lambda for various big matrix size to see when it level off:
bigm<-seq(50, 1000, 100) ## range of big matrix dimension tested
lam<-matrix(0, length(bigm), 1) ## vector to receive associated pop g.r.

## To run the bigmat function to identify the parsimonous big matrix size, I used random normal distribution of pruning intensity of equal lenght as the big matrix dimension. 

mH<-mean(datA$prun, na.rm=T) ## mean pruning intensity: 17.68
sdH<-sqrt(var(datA$prun, na.rm=T)) ## SD pruning intensity: 36.76
#H<-rnorm(length(bigm[i]), mean=mH, sd=sdH)

## However, given that harvest is in percentage, using normal distribution may not be appropriate here and a beta distribution would be recommended to estimate mean and variance. Maybe fit a beta distribution to pruning intensity (between 0,1) to get B(a,b), and then derive the mean mu and variance var using the following functions:
  ## mu= a/(a+b) and 
  ## var= a*b/([(a+b)^2*(a+b+1)] 

## For now, I will continue with the rnorm pruning intensity and later test what difference it makes to use a beta distribution

for (i in 1:length(bigm)){
	A<-array(0, c(bigm[i], bigm[i], length(bigm)))
	#H<-seq(minH, maxH, length=bigm[i])
	H<-rnorm(length(bigm[i]), mean=mH, sd=sdH)
	A[,,i]<-bigmat(bigm[i], H=H, pvec=p.vec)
	lam[i]<-eigen.analysis(A[,,i])$lambda1
}	

## Plot to identify the big matrix dimension
plot(bigm, lam, type="l", xlab="Size of the big matrix", ylab=expression(paste("Population growth rate, ", lambda)), col="red", lwd=3)

## This generate a really wierd plot where pop g.r. is not increasing smoothly as expected and as a result there is no asymptote and it is difficult to select a big mat dim based on this. From Zuidema et al (2010) JEcol discussion and Metcalf et al (2009), two of the rare tree IPM studies, a diameter interval of ~0.65 cm can be a good candidate, even though  dbh intervalfrom 1-10cm can also yield varying results. The goal of this paper at this point is to tes the effect of the variation in pruning intensity. Therefore, I will use the 0.65 cm interval which for the maximun tree size of 139 cm recorded in my study system yields a big matrix dim of 214. I will use 200 instead here.  


###### PART (7): FIGURES for IPM Functions ###########################

## Assuming a normal dist for pruning intensity (which is obviously wrong, will switch to beta dist later), the mean pruning intensity is 17.68% with SD 36.76. I seek to answer two questions here: (1) what is the effect of increasing mean pruning intensity on pop dynamics? and (2) what is the effect of increasing pruning intensity variance on population dynamics

## I will vary mean individual tree pruning intensity from 7.68 to 27.68 (that is within +/- 10%) and the variance from 26.76 to 46.76 (that is within +/- 10%).

mprun<-seq(7.68, 37.68, length=30)
sdprun<-seq(16.76, 56.76, length=30)

bigM<-200 ## choice based on litterature
mat<-array(0, c(bigM, bigM, length(mprun)))
LamM<-LamSD<-vector("numeric")

for (i in 1:length(mprun)){
  H<-rnorm(bigM, mean=mprun[i], sd=sdH)
	mat[,,i]<-bigmat(bigM, H=H, pvec=p.vec)
	LamM[i]<-eigen.analysis(mat[,,i])$lambda1
}

for (i in 1:length(sdprun)){
  H<-rnorm(bigM, mean=mH, sd=sdprun[i])
  mat[,,i]<-bigmat(bigM, H=H, pvec=p.vec)
  LamSD[i]<-eigen.analysis(mat[,,i])$lambda1
}

require(plotrix) # package plotrix is needed for function "ablineclip""

## Plot
op<-par(mfrow=c(1,2))
plot(mprun,LamM, xlab=expression(paste("Mean pruning intensity, ", mu[p])), ylab=expression(paste("Population growth rate, ", lambda)), col="black", pch=21, bg = "grey", cex = 1.5)
mtext(a)

reg1 <- lm(LamM~mprun)
ablineclip(reg1, lwd=2, col="red") 

plot(sdprun,LamSD, xlab=expression(paste("SD pruning intensity, ", sigma[p])), ylab=expression(paste("Population growth rate, ", lambda)), col="black", pch=21, bg = "grey", cex = 1.5)
mtext(b)

reg2 <- lm(LamSD~sdprun)
ablineclip(reg2, lwd=2, col="blue") 
par(op)

### Results here show that with increasing mean individual level harvesting (pruning) rate, population growth rate decreases, which is expected. Now, what was not expected is that with increasing variance within population harvesting rate, population growth rate increases for a given mean harvesting rate (here the observed mean harvesting rate; will test later if this is true for different values of mean). This suggest that homogenizing harvesting intensity at population level (say, harvest every individual at 50%, thus driving the variance to 0) may not be a good harvesting strategy. One would rather increase harvesting variance, with some individual not harvested and other heavily harvested and some lighly harvested to offer a mosaic of trees with different harvesting burden and therefore different fitness which may compensate each other. Example for such mamagement system exist in fire ecology or grazing ecology where a spatial scale mosaic of vegetation at different grazing or fire burning level can ensure a better primary and secondary productivity. Here I am showing for the first time that this may be possible for non-timber forest product harvest.

############################# END OF #############################
