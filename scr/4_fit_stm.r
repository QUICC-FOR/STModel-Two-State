#!/usr/bin/Rscript
library(argparse)
parser = ArgumentParser()
parser$add_argument("-s", "--species", default="28731-ACE-SAC", help="desired species code")
parser$add_argument("-r", "--rf", action="store_true", default=FALSE, help="random forest for prevalence (default: GLM)")
parser$add_argument("-m", "--gam", action="store_true", default=FALSE, help="use GAM for prevalence (default: GLM)")
parser$add_argument("-f", "--fraction", default=0.5, type="double", help="proportion of data to use for model fitting")
parser$add_argument("-i", "--interval", default=5, type="double", help="how many years should the parameterization interval be")
parser$add_argument("-t", "--tempvar", default="annual_mean_temp", help="Temperature variable to use")
parser$add_argument("-p", "--precipvar", default="tot_annual_pp", help="Precipitation variable to use")
parser$add_argument("-c", "--coldesign", default="1111111", help="7-character design vector for colonization model")
parser$add_argument("-e", "--extdesign", default="1111111", help="7-character design vector for extinction model")


argList = parser$parse_args()
spName = argList$species
targetInterval = argList$interval
colDesign = sapply(1:nchar(argList$coldesign), function(i) as.integer(substr(argList$coldesign,i,i)))
extDesign = sapply(1:nchar(argList$extdesign), function(i) as.integer(substr(argList$extdesign,i,i)))

design_err = function(x)
{
	err = length(x) != 7 | setequal(x, c(0,1))
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


transitionData = readRDS(paste("dat/", spName, "/", spName, "_transitions_projected.rds", sep=""))

NO! this doesnt work. need to subsample externally, so that all of the models use the SAME
subsample within a given species. otherwise AIC is invalid

# subsample data, if desired
if(argList$fraction < 1)
{
	sel = sample(nrow(transitionData), as.integer(argList$fraction*nrow(transitionData)))
	transitionData.subsample = transitionData[sel,]
} else {
	transitionData.subsample = transitionData
}


modelData = data.frame(state1 = transitionData.subsample$state1, 
		state2 = transitionData.subsample$state2,
		interval = (transitionData.subsample$year2 - transitionData.subsample$year1), 
		env1 = transitionData.subsample[,argList$tempvar], 
		env2 = transitionData.subsample[,argList$precipVar])

if(argList$gam)
{
	modelData$expectedPresence = transitionData.subsample$expectedGAM
	prevalenceVar = "GAM"
} else if(argList$rf)
{
	modelData$expectedPresence = transitionData.subsample$expectedRF
	prevalenceVar = "RF"
} else
{
	modelData$expectedPresence = transitionData.subsample$expectedGLM
	prevalenceVar = "GLM"
}
		
modelData = within(modelData, 
{
	if(argList$gam)
	{
		expectedPresence = transitionData.subsample$expectedGAM
	} else if(argList$rf) {
		expectedPresence = transitionData.subsample$expectedRF
	} else {
		expectedPresence = transitionData.subsample$expectedGLM
	}
})

# set up initial values for the parameters
# calculate average colonization and extinction probability
# 1 - (1 - C/P)^(1/i); where C is colonization prob (i.e., gamma * Presence)
pr.c = with(modelData, 1 - (1 - (sum(state1==0 & state2==1)/length(state1))/
		mean(expectedPresence))^(1/mean(interval)))  
pr.e = with(modelData, 1 - (1 - sum(state1==1 & state2==0)/length(state1))^(1/mean(interval)))

parameters = rep(0, sum(colDesign)+sum(extDesign))
if(colDesign[1] == 1) parameters[1] = pr.c
if(extDesign[1] == 1) parameters[sum(colDesign)+1] = pr.e

### perhaps try setting (in control)
### nb.stop.improvement = n
### how much for n? depends on how quickly the objective function runs
### a few hundred? thousand? ten?

library(GenSA)
# controlPars =list(verbose = TRUE, max.time = (60*20), smooth=TRUE)
controlPars =list(verbose = TRUE, nb.stop.improvement = (200), smooth=TRUE)
annealParams = GenSA(par = parameters, fn = minus_log_likelihood, 
		lower = rep(-50, length(parameters)), upper = rep(50, length(parameters)), 
		 control = controlPars, dat = modelData, parlist=c(colDesign, extDesign))

mll = minus_log_likelihood(annealParams$par, modelData, c(colDesign, extDesign))
annealResults = list(
	params = annealParams$par,
	AIC = 2*mll + 2*length(annealParams$par),
	BIC = 2*mll + log(nrow(modelData))*length(annealParams$par),
	env1 = argList$tempvar,
	env2 = argList$precipVar,
	colDesign = colDesign,
	extDesign = extDesign,
	species = spName,
	prevalenceVar = prevalenceVar,
	annealParams = annealParams)

resultDir = paste("results/", spName, "/anneal", sep="")
dir.create(resultDir, showWarnings=FALSE, recursive=TRUE)
saveRDS(annealResults, paste(resultDir, "/", spName, "_", annealResults$env1, "_", 
		annealResults$env2, "_", annealResults$colDesign, "_", annealResults$extDesign, 
		".rds", sep=""))
