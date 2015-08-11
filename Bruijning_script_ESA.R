##########################################################################
#################### Function to perform cost-benefit analysis ###########
##########################################################################

# R-code accompanying a presentation at ESA in Baltimore, 11 Aug 2015:
# "Surviving in a co-sexual world: A cost-benefit analysis of dioecy in tropical trees"
# Marjolein Bruijning, Marco D. Visser, Helene C. Muller-Landau, S. Joseph Wright, Liza S. Comita, Stephen P. Hubbell, Hans de Kroon, Eelke Jongejans 
# (c) Marjolein Bruijning, Radboud University, Nijmegen, the Netherlands: M.Bruijning@science.ru.nl

### Required functions:
# createkernels() # function to obtain all functions and parameters to construct the IPM (based on estimated coefficients)
# buildipm() # function to contruct IPM
# predictmodel() # function to extract parameters and function based on function form (linear or logistic), trait values and coefficients

### Required objects:
# testModel... # model averaged parameters for all intercept effects ($coeffInt) and slope effects ($coeffSlope)

### Arguments:
# traits: list with trait values for reference species (default set at average hermaphrodite)
# maxdbh: boundary of the IPM (default at 158 mm diamter for an average hermaphrodite species)
# sexratio: if dioecious, proportion of females (default at 0.5)
# seedproducing, fecundity, psts, sdlgrowth, sdlsurv, adultsurv, adultgrowth, prepr, distrecruits, dbhheight: all vital rates, if replace by dioecious set TRUE
# traitsDIO: list with trait value to compare to (default equal to 'traits')
# numberclasses_sdl: number of seedling classes in IPM (default at 200)
# numberclasses_ad: number of tree classes in IPM (default at 600)

performCostBenefit=function(traits=list(wd=0,seed=0,breedingDIO=0,smax=0),maxdbh=150,sexratio=0.5,
	seedproducing=FALSE,fecundity=FALSE,psts=FALSE,sdlgrowth=FALSE,sdlsurv=FALSE,
	adultsurv=FALSE,adultgrowth=FALSE,prepr=FALSE,distrecruits=FALSE,dbhheight=FALSE,
	traitsDIO=traits,numberclasses_sdl=200,numberclasses_ad=600)
	{
		# Get all coefficients and functions for an average hermaphrodite
		kernels=createkernels(traits,sexratio=sexratio) # List containing kernels$param and kernels$functions
		
		if(seedproducing==TRUE) {
			# Effect of having males (i.e. decreased reproduction probability)
			kernels$param[[12]] <- sexratio
		}
		if (fecundity==TRUE) {
			# Effect of dioecious fecundity
			traits <- traitsDIO
			traits$breedingDIO <- 1
			coeffInt <- testModelFec$coeffInt
			coeffSlope <- testModelFec$coeffSlope
			mod <- predictmodel(fun="funlin",coeffInt=coeffInt,coeffSlope=coeffSlope,xas="seed",fact=traits)
			kernels$param[[1]] <- as.numeric(mod$fun(mod$param,x=unlist(traits[names(traits)=="seed"]))) # Fecundity if dioecious
		}
		if(prepr==TRUE) {
			# Reproduction probability
			traits <- traitsDIO
			traits$breedingDIO <- 1
			coeffInt <- testModelRepr$coeffInt
			coeffSlope <- testModelRepr$coeffSlope
			mod <- predictmodel(fun="funlog",coeffInt=coeffInt,coeffSlope=coeffSlope,xas="dbh",fact=traits)
			kernels$functions[[2]] <- mod$fun
			kernels$param[[2]] <- mod$param # Reproduction parameters if dioecious
		}
		if (psts==TRUE) {
			# Effect of decreased psts
			traits <- traitsDIO
			traits$breedingDIO <- 1
			coeffInt <- testModelPsts$coeffInt
			coeffSlope <- testModelPsts$coeffSlope
			mod <- predictmodel(fun="funlin",coeffInt=coeffInt,coeffSlope=coeffSlope,xas="seed",fact=traits)
			kernels$param[[4]] <- plogis(as.numeric(mod$fun(mod$param,x=unlist(traits[names(traits)=="seed"])))) # Psts if dioecious
		}
		if(sdlgrowth==TRUE) {
			# Effect of dioecious seedling growth
			traits <- traitsDIO		
			traits$breedingDIO=1
			coeffInt <- testModelSdlgrowth$coeffInt
			coeffSlope <- testModelSdlgrowth$coeffSlope
			mod <- predictmodel(fun="funlin",coeffInt=coeffInt,coeffSlope=coeffSlope,xas="size",fact=traits)
			kernels$param[[6]] <- mod$param # Seedling growth if dioecious
			
			# Seedling growth variation
			coeffInt <- testModelsSeedlingRes$coeffInt
			coeffSlope <- testModelsSeedlingRes$coeffSlope
			mod <- predictmodel(fun="funlin",coeffInt=coeffInt,coeffSlope=coeffSlope,
				xas="seed",fact=traits)
			sdpar <- mod$fun(mod$par,as.numeric(traits[names(traits)=='seed']))
			kernels$param[[7]] <- c(0,sdpar)
			kernels$functions[[7]] <- function(x,param) {dnorm(x,param[1],param[2])}
		}
		if(sdlsurv==TRUE) {
			# Effect of dioecious seedling survival
			traits <- traitsDIO
			traits$breedingDIO <- 1
			coeffInt <- testModelSdlsurv$coeffInt
			coeffSlope <- testModelSdlsurv$coeffSlope
			mod <- predictmodel(fun="funlog",coeffInt=coeffInt,coeffSlope=coeffSlope,xas="size",fact=traits) # seedling survival if dioecious
			kernels$param[[5]] <- mod$param
		}
		if(adultsurv==TRUE) {
			# Effect of dioecious tree survival
			traits <- traitsDIO
			traits$breedingDIO <- 1
			coeffInt <- testModelAdultsurv$coeffInt
			coeffSlope <- testModelAdultsurv$coeffSlope
			mod <- predictmodel(fun="funlog",coeffInt=coeffInt,coeffSlope=coeffSlope,xas="dbh1",fact=traits) # tree survival if dioecious
			kernels$param[[8]] <- mod$param
		}
		if(adultgrowth==TRUE) {
			# Effect of dioecious tree growth
			traits <- traitsDIO
			traits$breedingDIO <- 1
			coeffInt <- testModelAdultgrowth$coeffInt
			coeffSlope <- testModelAdultgrowth$coeffSlope
			mod <- predictmodel(fun="funlin",coeffInt=coeffInt,coeffSlope=coeffSlope,xas="dbh1",fact=traits)
			kernels$param[[9]] <- mod$param # tree growth if dioecious
			kernels$functions[[9]] <- mod$fun
						
			# Tree growth variation
			mod <- predictmodel(fun="funlin",coeffInt=testModelsAdultRes$coeffInt,coeffSlope=testModelsAdultRes$coeffSlope,
				xas="wd",fact=traits)
			sdpar <- mod$fun(mod$par,as.numeric(traits[names(traits)=='wd']))
			kernels$param[[10]] <- c(0,sdpar)
			kernels$functions[[10]] <- function(x,param) {dnorm(x,param[1],param[2])}
		}
		if(distrecruits==TRUE) {
			# Recruit size distribution
			traits <- traitsDIO
			traits$breedingDIO <- 1
			# Shape
			coeffInt <- testModelsNewRecruitsShape$coeffInt
			coeffSlope <- testModelsNewRecruitsShape$coeffSlope
			mod <- predictmodel(fun="funlin",coeffInt=coeffInt,coeffSlope=coeffSlope,
				xas="seed",fact=traits)
			shapepar <- mod$fun(mod$param,as.numeric(traits[names(traits)=='seed']))	
			# Scale
			coeffInt <- testModelsNewRecruitsScale$coeffInt
			coeffSlope <- testModelsNewRecruitsScale$coeffSlope
			mod <- predictmodel(fun="funlin",coeffInt=coeffInt,coeffSlope=coeffSlope,
				xas="seed",fact=traits)
			scalepar <- mod$fun(mod$param,as.numeric(traits[names(traits)=='seed']))	
			kernels$param[[3]] <- c(shapepar,scalepar)
			kernels$functions[[3]] <- function(x,param) {dweibull(x,param[1],param[2])}
		}
		if(dbhheight==TRUE) {
			# Seedling transition
			traits <- traitsDIO
			traits$breedingDIO <- 1	
			coeffInt <- testModelDbh$coeffInt
			coeffSlope <- testModelDbh$coeffSlope
			mod <- predictmodel(fun="funlin",coeffInt=coeffInt,coeffSlope=coeffSlope,xas="size",fact=traits)
			kernels$param[[11]] <- mod$param
			kernels$functions[[11]] <- mod$fun
		}
	
	# BUILD IPM
	model <- buildipm(param=kernels$param,functions=kernels$functions,
		numberclasses_sdl=numberclasses_sdl,numberclasses_ad=numberclasses_ad,
		maxdbh=maxdbh)
	lambda <- eigenAnalysis(model$ipm)$lambda1
	return(lambda)
}

