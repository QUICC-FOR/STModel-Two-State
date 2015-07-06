library(coda)

# temporary
setwd("~/Dropbox/work/projects/STModel-Two-State_git/")
spName = '19447-QUE-VEL'

p1.raw = read.csv(file.path('species', spName, 'res', 'mcmc1', 'posterior.csv'))
posterior = mcmc(p1.raw)

plot(posterior, ask=T)
summary(posterior)

pquant = summary(posterior)$quantiles
minpar = which(pquant[,1] == min(pquant[,1]))
maxpar = which(pquant[,5] == max(pquant[,5]))

minsamp = posterior[posterior[,minpar] == min(posterior[,minpar]),]
maxsamp = posterior[posterior[,maxpar] == min(posterior[,maxpar]),]

findsum = function(x, fun)
{
	if(length(dim(x)) == 2)
	{
		res = x[which(apply(x, 1, sum) == fun(apply(x, 1, sum))),]
		if(is.matrix(res)) res = res[1,]
		return(res)
	}
	if(is.vector(x))
		return(x)
	stop("Need a vector or matrix-like object")
}
	

initMax = findsum(maxsamp, max)
initMin = findsum(minsamp, min)
ylim=c(min(c(initMin, initMax)), max(c(initMin, initMax)))
plot(initMax, col='blue', ylim=ylim, pch=16)
points(initMin, col='red', pch=16)
abline(h=0)

i2path = file.path('species', spName, 'dat', 'mcmc_inits2.txt')
i2 = read.csv(i2path)
i2$initialValue = initMax
write.csv(i2, i2path, row.names=FALSE)

i2path = file.path('species', spName, 'dat', 'mcmc_inits3.txt')
i2 = read.csv(i2path)
i2$initialValue = initMin
write.csv(i2, i2path, row.names=FALSE)