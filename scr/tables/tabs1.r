spList = readRDS('dat/speciesList.rds')


for(spName in spList) {
	trdat = rbind(readRDS(file.path('dat', 'stm_calib', paste0(spName, 'stm_calib.rds'))),
		readRDS(file.path('dat', 'stm_valid', paste0(spName, 'stm_valid.rds'))))
	prdat = readRDS(file.path('dat', 'presence', paste0(spName, '_presence.rds')))

	type = rep('pres', nrow(trdat))
	type[trdat$state1 == 0 & trdat$state2 == 0] = 'abs'
	type[trdat$state1 == 0 & trdat$state2 == 1] = 'col'
	type[trdat$state1 == 1 & trdat$state2 == 0] = 'ext'
	cat(spName, '\n')
	print(table(type))
	cat(sum(prdat[,spName]), '\n')
}

obs = sapply(spList, function(spName)
{
prdat = readRDS(file.path('dat', 'presence', paste0(spName, '_presence.rds')))
sum(prdat[,spName])
})

trdat = lapply(spList, function(spName) {
	td = rbind(readRDS(file.path('dat', 'stm_calib', paste0(spName, 'stm_calib.rds'))),
		readRDS(file.path('dat', 'stm_valid', paste0(spName, 'stm_valid.rds'))))
	type = rep('pres', nrow(td))
	type[td$state1 == 0 & td$state2 == 0] = 'abs'
	type[td$state1 == 0 & td$state2 == 1] = 'col'
	type[td$state1 == 1 & td$state2 == 0] = 'ext'
	type
})

pres = sapply(trdat, function(x) sum(x == 'pres'))
abse = sapply(trdat, function(x) sum(x == 'abs'))
colo = sapply(trdat, function(x) sum(x == 'col'))
ext = sapply(trdat, function(x) sum(x == 'ext'))
