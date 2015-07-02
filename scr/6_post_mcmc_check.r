library(coda)

# temporary
setwd("~/Dropbox/work/projects/STModel-Two-State_git/")
spName = '28731-ACE-SAC'
constCols = c(5,6,7,12,13,14)

p1 = read.csv(file.path('species', spName, 'res', 'mcmc1', 'posterior.csv'))[,-constCols]
p2 = read.csv(file.path('species', spName, 'res', 'mcmc2', 'posterior.csv'))[,-constCols]
p3 = read.csv(file.path('species', spName, 'res', 'mcmc3', 'posterior.csv'))[,-constCols]
posterior = mcmc.list(mcmc(p1), mcmc(p2), mcmc(p3))

plot(posterior, ask=T)
summary(posterior)
gelman.diag(posterior)

st = floor(nrow(p1)/2) + 1
en = nrow(p1)
posterior.burnin = mcmc.list(mcmc(p1[st:en,], start=st), mcmc(p2[st:en,], start=st), mcmc(p3[st:en,], start=st))

saveRDS(posterior.burnin, file.path('species', spName, 'res', paste(spName, 'posterior.rds', sep='_')))