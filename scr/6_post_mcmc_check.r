library(coda)

setwd("~/Dropbox/work/projects/STModel-Two-State_git/")
spName = '18037-PIN-TAE'

design = '11110001111000'
design = sapply(1:nchar(design), function(i) as.integer(substr(design, i, i)))
constCols = which(design == 0)

st1 = st2 = st3 = 1
# special instructions for individual species due to poor starts
if(spName == '19408-QUE-RUB')
	st2 = 14001
if(spName == '27821-NYS-SYL')
	st2 = 17001
if(spName == '19280-QUE-NIG')
	st3 = 19001

p1 = read.csv(file.path('species', spName, 'res', 'mcmc1', 'posterior.csv'))[,-constCols]
p2 = read.csv(file.path('species', spName, 'res', 'mcmc2', 'posterior.csv'))[,-constCols]
p2 = p2[st2:nrow(p2),]
p3 = read.csv(file.path('species', spName, 'res', 'mcmc3', 'posterior.csv'))[,-constCols]
p3 = p3[st3:nrow(p3),]
posterior = mcmc.list(mcmc(p1), mcmc(p2), mcmc(p3))
lapply(posterior, nrow)

plot(posterior, ask=T)
summary(posterior)
gelman.diag(posterior)

st = floor(nrow(p1)/2) + 1
en = nrow(p1)
th.int = 10
th = seq(st, en, th.int)
posterior.burnin = mcmc.list(mcmc(p1[st:en,], start=st), mcmc(p2[st:en,], start=st), mcmc(p3[st:en,], start=st))
posterior.thinned = mcmc.list(mcmc(p1[th,], start=st, thin=th.int), mcmc(p2[th,], start=st, thin=th.int), mcmc(p3[th,], start=st, thin=th.int))
saveRDS(posterior.burnin, file.path('species', spName, 'res', paste(spName, 'posterior.rds', sep='_')))
saveRDS(posterior.thinned, file.path('species', spName, 'res', paste(spName, 'posterior_thinned.rds', sep='_')))
