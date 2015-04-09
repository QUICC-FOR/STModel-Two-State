setwd("/Users/mtalluto/Dropbox/work/projects/STModel-Two-State_git")

dat = readRDS("dat/28731-ACE-SAC_transitions_projected.rds")
annealPars = readRDS("results/28731-ACE-SAC_anneal_parameters_0.1.rds")

cubic = (length(annealPars$par) > 10)
if(cubic)
{
	parNames = c('g0','g1','g2','g3','g4','g5','g6','e0','e1','e2','e3','e4','e5','e6')
} else
{
	parNames = c('g0','g1','g2','g3','g4','e0','e1','e2','e3','e4')
}
mcmcDat = with(dat, data.frame(
	initial = state1,
	final = state2,
	env1 = annual_mean_temp,
	env2 = tot_annual_pp,
	interval = year2-year1,
	prevalence1 = expectedGAM))
	
parInfo = data.frame(
	name = parNames,
	initialValue = annealPars$par,
	priorMean = 0,
	priorSD = 10000,
	priorDist = "Normal")
	
# take a random sample for testing
ind = sample(nrow(mcmcDat), 0.1*nrow(mcmcDat), replace=F)
write.csv(mcmcDat[ind,], "dat/mcmc_trans_28731-ACE-SAC_0.1.txt", row.names=FALSE)
write.csv(parInfo, "dat/mcmc_pars_28731-ACE-SAC.txt", row.names=FALSE)