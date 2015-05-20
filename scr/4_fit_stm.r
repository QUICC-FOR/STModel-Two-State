#!/usr/bin/Rscript
library(argparse)
parser = ArgumentParser()
parser$add_argument("-s", "--species", default="28731-ACE-SAC", help="desired species code")
parser$add_argument("-c", "--cubic", action="store_true", default=FALSE, help="use cubic terms in model (default is square only)")
parser$add_argument("-r", "--rf", action="store_true", default=FALSE, help="random forest for prevalence (default: GLM)")
parser$add_argument("-m", "--gam", action="store_true", default=FALSE, help="use GAM for prevalence (default: GLM)")
parser$add_argument("-f", "--fraction", default=0.1, type="double", help="proportion of data to use for model fitting")
parser$add_argument("-i", "--interval", default=5, type="double", help="how many years should the parameterization interval be")
argList = parser$parse_args()
spName = argList$species
targetInterval = argList$interval

minus_log_likelihood = function(p, dat, useCubic)
{
	inv_logit = function(logitVal)
	{
		val = numeric(length(logitVal))
		val[logitVal > 0] = 1.0 / (1.0 + exp(-logitVal[logitVal > 0]))
		val[logitVal <= 0] = exp(logitVal[logitVal <= 0]) / (1.0 + exp(logitVal[logitVal <= 0]))
		
		# guard against big numbers
		val[is.infinite(exp(logitVal))] = 1
		
		return(val)
	}
	
	# full model linear predictors
	if(useCubic)
	{
		logit_gamma_annual = p[1] + p[2]*dat$env1 + p[3]*dat$env2 + p[4]*dat$env1^2 + 
				p[5]*dat$env2^2 + p[6]*dat$env1^3 + p[7]*dat$env2^3 
		logit_epsilon_annual = p[8] + p[9]*dat$env1 + p[10]*dat$env2 + p[11]*dat$env1^2 + 
				p[12]*dat$env2^2 + p[13]*dat$env1^3 + p[14]*dat$env2^3 
	} else {
		logit_gamma_annual = p[1] + p[2]*dat$env1 + p[3]*dat$env2 + p[4]*dat$env1^2 + 
				p[5]*dat$env2^2 
		logit_epsilon_annual = p[6] + p[7]*dat$env1 + p[8]*dat$env2 + p[9]*dat$env1^2 + 
				p[10]*dat$env2^2
	}

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

# set seed - drawn from sample(1:1e6, 1)
## set.seed(588533)

# subsample data, if desired
# transitionData.subsample = transitionData
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
		env1 = scale(transitionData.subsample$annual_mean_temp), 
		env2 = scale(transitionData.subsample$tot_annual_pp))
		
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

if(argList$cubic) 
{
	parameters = c(
		qlogis(pr.c),0,0,0,0,0,0,		## colonization
		qlogis(pr.e),0,0,0,0,0,0)		## extinction
} else {
	parameters = c(
		qlogis(pr.c),0,0,0,0,		## colonization
		qlogis(pr.e),0,0,0,0)		## extinction
}


library(GenSA)
annealParams = GenSA(par = parameters, fn = minus_log_likelihood, 
		lower = rep(-50, length(parameters)), upper = rep(50, length(parameters)), 
		 control = list(verbose = TRUE, max.time = (60*20), smooth=TRUE), dat = modelData,
		 useCubic = argList$cubic)
		 
saveRDS(annealParams, paste("results/", spName, "/", spName, "_anneal_parameters.rds", sep=""))
