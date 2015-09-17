#!/usr/bin/Rscript

parameters = commandArgs(trailingOnly = TRUE)
spName = parameters[1]

p1.raw = read.csv(file.path('species', spName, 'res', 'mcmc1', 'posterior.csv'))

samp.points = c(2500, 5000, 7500)

i0path = file.path('species', spName, 'dat', 'mcmc_inits.txt')
i0 = read.csv(i0path)

i1path = file.path('species', spName, 'dat', 'mcmc_inits1.txt')
i1 = i0
i1$initialValue = p1.raw[samp.points[1],]
write.csv(i1, i1path, row.names=FALSE)

i2path = file.path('species', spName, 'dat', 'mcmc_inits2.txt')
i2 = i0
i2$initialValue = p1.raw[samp.points[2],]
write.csv(i2, i2path, row.names=FALSE)

i3path = file.path('species', spName, 'dat', 'mcmc_inits3.txt')
i3 = i0
i3$initialValue = p1.raw[samp.points[3],]
write.csv(i3, i3path, row.names=FALSE)

pdf(file=file.path('species', spName, 'img', 'mcmc0_trace.pdf'), w=5, h=5)
nms = c('g0','g1','g2','g3','g4','g5','g6','e0','e1','e2','e3','e4','e5','e6')
par(mfrow=c(2,ncol(p1.raw)/2))
for(i in 1:ncol(p1.raw))
	plot(1:nrow(p1.raw), p1.raw[, i], type='l', xlab="iteration", ylab=nms[i])
dev.off()
