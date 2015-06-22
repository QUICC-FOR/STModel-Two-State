
## need to work this in to this script somewhere
# remove extremely long sampling intervals from the dataset 
# (to avoid the possibility of 2 transitions in the space we missed)
intervals = transitionData.scaled$year2 - transitionData.scaled$year1
indices = which(intervals <= 15)
transitionData.scaled = transitionData.scaled[indices,]






#!/usr/bin/Rscript
library(GenSA)
targetInterval = 5
modType = "rf"

parameters = commandArgs(trailingOnly = TRUE)
#parameters = scan("inp.txt", what=character())
spName = parameters[1]
tempVar = parameters[2]
precipVar = parameters[3]
# adding a one to the beginning for the intercept
colDesign = c(1, sapply(1:nchar(parameters[4]), function(i) as.integer(substr(parameters[4],i,i))))
extDesign = c(1, sapply(1:nchar(parameters[5]), function(i) as.integer(substr(parameters[5],i,i))))

design_err = function(x)
{
	err = length(x) != 7 | (!setequal(x, c(0,1)) & !setequal(x, 1))
}
if(design_err(colDesign) | design_err(extDesign)) 
	stop("Error in colonization or extinction design vectors: must be a vector of length 7 containing only 0 or 1")

minus_log_likelihood = function(params, dat, parlist)
{
	# parlist is a vector of 1s and zeros indicating which positions in the likelihood
	# should be included; the number of 1s should be equal to the length of params
	numPars = 14
	if(sum(parlist) != length(params)) 
		stop("Number of parameters must equal the number specified in parlist")
	if(length(parlist) < numPars) 
		stop(paste("Parlist must have", numPars, 
				"items (the number of parameters in the full model)"))
	p = rep(0, numPars)
	p[parlist == 1] = params

	
	# full model linear predictors
	logit_gamma_annual = p[1] + p[2]*dat$env1 + p[3]*dat$env2 + p[4]*dat$env1^2 + 
			p[5]*dat$env2^2 + p[6]*dat$env1^3 + p[7]*dat$env2^3 
	logit_epsilon_annual = p[8] + p[9]*dat$env1 + p[10]*dat$env2 + p[11]*dat$env1^2 + 
			p[12]*dat$env2^2 + p[13]*dat$env1^3 + p[14]*dat$env2^3 

	# model macro parameters - gamma and epsilon
	gamma_annual = plogis(logit_gamma_annual)
	epsilon_annual = plogis(logit_epsilon_annual)

	gamma_interval =  1 - (1 - gamma_annual)^(dat$interval/targetInterval)
	epsilon_interval =  1 - (1 - epsilon_annual)^(dat$interval/targetInterval)

	# likelihood
	colonizations = which(dat$state1 == 0 & dat$state2 == 1)
	extinctions = which(dat$state1 == 1 & dat$state2 == 0)
	presences = which(dat$state1 == 1 & dat$state2 == 1)
	absences = which(dat$state1 == 0 & dat$state2 == 0)
	
	liklihood = numeric(nrow(dat))

	liklihood[colonizations] = gamma_interval[colonizations] * 
			dat$expectedPresence[colonizations]
	liklihood[extinctions] = epsilon_interval[extinctions]
	liklihood[presences] = 1 - epsilon_interval[presences]
	liklihood[absences] = 1 - gamma_interval[absences] * dat$expectedPresence[absences]
	
	# guard against likelihood of zero
	liklihood[liklihood==0] = .Machine$double.xmin
	
	return(-1 * sum(log(liklihood)))
}


transitionData = readRDS(paste(spName, "_transitions_projected_subsample.rds", sep=""))

modelData = data.frame(state1 = transitionData$state1, 
		state2 = transitionData$state2,
		interval = (transitionData$year2 - transitionData$year1))
		
if(tempVar == "NA") {
	modelData$env1 = 0
} else {
	modelData$env1 = transitionData[,tempVar]
}

if(precipVar == "NA") {
	modelData$env2 = 0
} else {
	modelData$env2 = transitionData[,precipVar]
}

if(modType == "gam")
{
	modelData$expectedPresence = transitionData$expectedGAM
	prevalenceVar = "GAM"
} else if(modType == "rf")
{
	modelData$expectedPresence = transitionData$expectedRF
	prevalenceVar = "RF"
} else
{
	modelData$expectedPresence = transitionData$expectedGLM
	prevalenceVar = "GLM"
}
		

# set up initial values for the parameters
# calculate average colonization and extinction probability
# 1 - (1 - C/P)^(1/i); where C is colonization prob (i.e., gamma * Presence)
pr.c = with(modelData, 1 - (1 - (sum(state1==0 & state2==1)/length(state1))/
		mean(expectedPresence))^(1/mean(interval)))  
pr.e = with(modelData, 1 - (1 - sum(state1==1 & state2==0)/length(state1))^(1/mean(interval)))

parameters = rep(0, sum(colDesign)+sum(extDesign))
if(colDesign[1] == 1) parameters[1] = pr.c
if(extDesign[1] == 1) parameters[sum(colDesign)+1] = pr.e


print("Starting anneal with the following parameters:")
print(paste("Species: ", spName, sep=""))
print(paste("Env1: ", tempVar, sep=""))
print(paste("Env2: ", precipVar, sep=""))
print(paste("Model Design: ", paste(colDesign, collapse=""), " + ",
		paste(extDesign, collapse=""), sep=""))


controlPars =list(verbose = TRUE, max.time = (90*60), smooth=TRUE)
annealParams = GenSA(par = parameters, fn = minus_log_likelihood, 
		lower = rep(-50, length(parameters)), upper = rep(50, length(parameters)), 
		 control = controlPars, dat = modelData, parlist=c(colDesign, extDesign))
print("Anneal finished")

mll = minus_log_likelihood(annealParams$par, modelData, c(colDesign, extDesign))
annealResults = list(
	params = annealParams$par,
	AIC = 2*mll + 2*length(annealParams$par),
	BIC = 2*mll + log(nrow(modelData))*length(annealParams$par),
	env1 = tempVar,
	env2 = precipVar,
	colDesign = colDesign,
	extDesign = extDesign,
	species = spName,
	prevalenceVar = prevalenceVar,
	annealParams = annealParams)

dir.create('out/', showWarnings=FALSE, recursive=TRUE)
filename = paste(spName, annealResults$env1, annealResults$env2, 
		paste(annealResults$colDesign, collapse=""), 
		paste(annealResults$extDesign, collapse=""), sep="-")
saveRDS(annealResults, paste('out/', filename, ".rds", sep=""))
print(paste("Saved file ", paste(spName, annealResults$env1, annealResults$env2, 
		paste(annealResults$colDesign, collapse=""), 
		paste(annealResults$extDesign, collapse=""), sep="-"), sep=""))

