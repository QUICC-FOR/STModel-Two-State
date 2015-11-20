# do this for each species
## spName = '' # set the species by hand here

spPosterior = readRDS(paste0('res/posterior/',spName,'_posterior.rds'))[[1]]
rcRes = 100
info = speciesInfo[speciesInfo$spName == spName,]
calibDat = readRDS(file.path('dat', 'stm_calib', paste0(spName,'stm_calib.rds')))

# set x values for the response curves; for constants (env1 is constant in the curve 
# for env2, we lookup the values in info; if missing we use 0 (the mean)
env1 = with(calibDat, seq(min(annual_mean_temp), max(annual_mean_temp), length.out=rcRes))
env1_c = if(is.null(info$rc_tval) || is.na(info$rc_tval)) 
	{0} else {info$rc_tval}
env2 = with(calibDat, seq(min(tot_annual_pp), max(tot_annual_pp), length.out=rcRes))
env2_c = if(is.null(info$rc_pval) || is.na(info$rc_pval))
	{0} else {info$rc_pval}
	
quartz(w=18, h=6)
par(mfcol=c(2,9), mar=c(3,3,0,0), oma=c(0,0,2,0))
for(env_c in seq(-2,2,.5))
{
e1Preds_e = e1Preds_c = e2Preds_e = e2Preds_c = matrix(NA, nrow=nrow(spPosterior), ncol=rcRes)
env1_c = env2_c = env_c
for(i in 1:nrow(spPosterior))
{
	e1Preds_e[i,] = compute_e(spPosterior[i,], env1, env2_c)
	e1Preds_c[i,] = compute_c(spPosterior[i,], env1, env2_c)
	e2Preds_e[i,] = compute_e(spPosterior[i,], env1_c, env2)
	e2Preds_c[i,] = compute_c(spPosterior[i,], env1_c, env2)
}	
plot(env1, colMeans(e1Preds_e), col='red', type='l', ylim=c(0,0.5))
lines(env1, colMeans(e1Preds_c), col='blue')
mtext(env_c, side=3, cex=0.75)
plot(env2, colMeans(e2Preds_e), col='red', type='l', ylim=c(0,0.5))
lines(env2, colMeans(e2Preds_c), col='blue')
}
