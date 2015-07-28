library(coda)

setwd("~/Dropbox/work/projects/STModel-Two-State_git/")
spName = 'NA-CAR-ALB'

modSel = readRDS(file.path('species', spName, 'res', paste(spName, 'modelSelection.rds', sep='_')))
print(modSel[1,])

design = modSel[1,'design']

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


# throw out first half of the samples
# everyone gets thinned down to 10,000 samples for the 'full' posterior and 1000 for maps
st = floor(nrow(p1)/2) + 1
en = nrow(p1)

full.len = 10000
full.th = floor((en - st)/(full.len - 1))
full.seq = seq(st, en, by = full.th)

thin.len = 1000
thin.th = floor((en - st)/(thin.len - 1))
thin.seq = seq(st, en, by = thin.th)

posterior.burnin = list(p1[full.seq,], p2[full.seq,], p3[full.seq,])
posterior.burnin = as.mcmc.list(lapply(posterior.burnin, mcmc, start=st, thin=full.th))

posterior.thinned = list(p1[thin.seq,], p2[thin.seq,], p3[thin.seq,])
posterior.thinned = as.mcmc.list(lapply(posterior.thinned, mcmc, start=st, thin=thin.th))

saveRDS(posterior.burnin, file.path('species', spName, 'res', paste(spName, 'posterior.rds', sep='_')))
saveRDS(posterior.thinned, file.path('species', spName, 'res', paste(spName, 'posterior_thinned.rds', sep='_')))
