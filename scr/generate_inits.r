#!/usr/bin/Rscript

parameters = commandArgs(trailingOnly = TRUE)
spName = parameters[1]

p1.raw = read.csv(file.path('species', spName, 'res', 'mcmc1', 'posterior.csv'))

no.variance = function(posterior)
{
	vars = apply(posterior, 2, var)
	return(all.equal(vars, rep(0, length(vars))))
}

# drop half
p1.half = p1.half[1:floor(nrow(p1.raw)/2),]
samp.points = floor(c(0.25, 0.5, 0.75) * nrow(p1.half))

# check that the sampler has actually moved between each sample point
st = 1
for(sp in samp.points)
{
	if(no.variance(p1.half[st:sp,])
		stop(paste("Error in selecting inits for ", spName, 
				"\nno variance in posterior between samples", nrow(p1.raw)+st, "and",
				nrow(p1.raw)+sp))
}

i0path = file.path('species', spName, 'dat', 'mcmc_inits.txt')
i0 = read.csv(i0path)

i1path = file.path('species', spName, 'dat', 'mcmc_inits1.txt')
i1 = i0
i1$initialValue = unlist(p1.half[samp.points[1],])
write.csv(i1, i1path, row.names=FALSE)

i2path = file.path('species', spName, 'dat', 'mcmc_inits2.txt')
i2 = i0
i2$initialValue = unlist(p1.half[samp.points[2],])
write.csv(i2, i2path, row.names=FALSE)

i3path = file.path('species', spName, 'dat', 'mcmc_inits3.txt')
i3 = i0
i3$initialValue = unlist(p1.half[samp.points[3],])
write.csv(i3, i3path, row.names=FALSE)

pdf(file=file.path('species', spName, 'img', 'mcmc0_trace.pdf'), w=5, h=5)
nms = c('g0','g1','g2','g3','g4','g5','g6','e0','e1','e2','e3','e4','e5','e6')
par(mfrow=c(2,ncol(p1.half)/2))
for(i in 1:ncol(p1.half))
	plot(1:nrow(p1.half), p1.half[, i], type='l', xlab="iteration", ylab=nms[i])
dev.off()
