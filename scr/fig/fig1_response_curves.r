library(coda)

# temporary
spName = '18032-ABI-BAL'
env1 = 'annual_mean_temp'
env2 = 'tot_annual_pp'

compute_e = function(p, env1, env2)
{
	plogis(p[8] + env1*p[9] + env2*p[10] + env1^2*p[11] + env2^2*p[12] + env1^3*p[13] + env2^3*p[14])
}


compute_c = function(p, env1, env2)
{
	plogis(p[1] + env1*p[2] + env2*p[3] + env1^2*p[4] + env2^2*p[5] + env1^3*p[6] + env2^3*p[7])
}


posterior = readRDS(file.path('species', spName, 'res', paste(spName, 'posterior.rds', sep='_')))
climGrid = readRDS('dat/climateGrid_scaled.rds')


x1 = seq(min(climGrid[,env1]), max(climGrid[,env1]), length.out=100)
x2 = seq(min(climGrid[,env2]), max(climGrid[,env2]), length.out=100)

# now compute the y values for each x value for every posterior sample
# first we have to choose a value at which to fix the second environmental variable
# do this with the mean values of the parameters
meanParams = colMeans(posterior)
meanParams = c(meanParams[1:5], 0, 0, meanParams[6:10], 0, 0)
# find the point at which env1 maximizes lambda
env1.yc = compute_c(meanParams, x1, rep(0, length(x1)))
env1.ye = compute_e(meanParams, x1, rep(0, length(x1)))
env1.lam = env1.yc - env1.ye
env1.maxx = x1[which(env1.lam == max(env1.lam))[1]]

# find the point at which env2 maximizes lambda with the max from the prev step
env2.yc = compute_c(meanParams, rep(env1.maxx, length(x2)), x2)
env2.ye = compute_e(meanParams, rep(env1.maxx, length(x2)), x2)
env2.lam = env2.yc - env2.ye
env2.maxx = x2[which(env2.lam == max(env2.lam))[1]]

# re-do the computation of x1 with the new max from x2
env1.yc = compute_c(params, x1, rep(env2.maxx, length(x1)))
env1.ye = compute_e(params, x1, rep(env2.maxx, length(x1)))
env1.lam = env1.yc - env1.ye
env1.maxx = x1[which(env1.lam == max(env1.lam))[1]]


y1.c = sapply(1:nrow(posterior), function(x)
{
	
	params = c(posterior[x,1:5], 0, 0, posterior[x,6:10], 0, 0)
	compute_c(params, x1, rep(env2.maxx, length(x1)))
})
y1.e = sapply(1:nrow(posterior), function(x)
{
	params = c(posterior[x,1:5], 0, 0, posterior[x,6:10], 0, 0)
	compute_e(params, x1, rep(env2.maxx, length(x1)))
})
y2.c = sapply(1:nrow(posterior), function(x)
{
	params = c(posterior[x,1:5], 0, 0, posterior[x,6:10], 0, 0)
	compute_c(params, rep(env1.maxx, length(x2)), x2)
})
y2.e = sapply(1:nrow(posterior), function(x)
{
	params = c(posterior[x,1:5], 0, 0, posterior[x,6:10], 0, 0)
	compute_e(params, rep(env1.maxx, length(x2)), x2)
})



env1.yc = cbind(rowMeans(y1.c), apply(y1.c, 1, quantile, 0.025), apply(y1.c, 1, quantile, 0.975))
env1.ye = cbind(rowMeans(y1.e), apply(y1.e, 1, quantile, 0.025), apply(y1.e, 1, quantile, 0.975))
env2.yc = cbind(rowMeans(y2.c), apply(y2.c, 1, quantile, 0.025), apply(y2.c, 1, quantile, 0.975))
env2.ye = cbind(rowMeans(y2.e), apply(y2.e, 1, quantile, 0.025), apply(y2.e, 1, quantile, 0.975))


# now put the x variables back on their original scales
climScale = readRDS("dat/climate_scaling.rds")
x1.us = (x1 * climScale$scale[env1]) + climScale$center[env1]
x2.us = (x2 * climScale$scale[env2]) + climScale$center[env2]

varNames = readRDS("dat/climVariableNames.rds")
e1.lab = varNames[which(varNames[,1] == env1),2]
e2.lab = varNames[which(varNames[,1] == env2),2]

pdf(w=6.5, h=8, file="img/response_curves.pdf")
par(mfrow=c(4,2), bty='n', mar=c(3,3,0.5, 0.5), mgp=c(2,0.5,0), tck=-0.03)
for(i in 1:4)
{
#par(fig=c(0, 1, 0, 1), bty='o', mar=c(0.5, 0.5, 0.5, 0.5))
#plot(0,0, type='n',xaxt='n', yaxt='n', xlab='', ylab='')
#par(fig=c(0.04, 0.5, 0.04, 0.96), bty='n', new=TRUE, mar=c(4,4, 0.1, 0.1))
plot(x1.us, env1.yc[,1], xlab=e1.lab, ylab="Probability", ylim=c(0,0.25), col='blue', type='l')
polygon(c(x1.us, rev(x1.us)), c(env1.yc[,2], rev(env1.yc[,3])), col="#0000FF22", border=NA)
lines(x1.us, env1.ye[,1], col='red')
polygon(c(x1.us, rev(x1.us)), c(env1.ye[,2], rev(env1.ye[,3])), col="#FF000022", border=NA)

#par(fig=c(0.5, 0.96, 0.04, 0.96), bty='n', new=TRUE)
plot(x2.us, env2.yc[,1], xlab=e2.lab, ylab="", ylim=c(0,0.25), col='blue', type='l')
polygon(c(x2.us, rev(x2.us)), c(env2.yc[,2], rev(env2.yc[,3])), col="#0000FF22", border=NA)
lines(x2.us, env2.ye[,1], col='red')
polygon(c(x2.us, rev(x2.us)), c(env2.ye[,2], rev(env2.ye[,3])), col="#FF000022", border=NA)
}

dev.off()