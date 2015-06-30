library(coda)

# temporary
setwd("~/Dropbox/work/projects/STModel-Two-State_git/")
spName = '28731-ACE-SAC'

p1.raw = read.csv(file.path('species', spName, 'res', 'mcmc1', 'posterior.csv'))
p2.raw = read.csv(file.path('species', spName, 'res', 'mcmc2', 'posterior.csv'))
p3.raw = read.csv(file.path('species', spName, 'res', 'mcmc3', 'posterior.csv'))

p1 = mcmc(p1.raw)
p2 = mcmc(p2.raw)
p3 = mcmc(p3.raw)

posterior = mcmc.list(p1, p2, p3)

plot(posterior, ask=T)