#!/usr/bin/env Rscript
# to be run after the initial short chain to select some inits
# uses a simple distance/entropy computation to select the inits that are farthest
# from each other (starting with one random selection)


## dependencies:
## res/mcmc/spname/0[i]/*
## dat/mcmc/spname_mcmcInit[_int].txt

## creates
## dat/mcmc/spname_mcmcInit[_int][1,2,3].txt


library(foreach)
library(iterators)

spList = c('18034-PIC-RUB', '183295-PIC-GLA', '183319-PIN-BAN', '183375-PIN-RES', 
		'183385-PIN-STR', '183397-TSU-CAN', '183412-LAR-LAR', '19287-QUE-MAC', 
		'19462-FAG-GRA', '22463-POP-GRA', '32931-FRA-AME', '32945-FRA-NIG', 
		'505490-THU-OCC')

shan = function(init, orig,length.out=3)
{
	while(nrow(init) < length.out)
	{
		dm = as.matrix(dist(rbind(init, orig)))[,1:nrow(init), drop=FALSE]
		sh = apply(dm, 1, function(x) -sum(x * log(x)))
		sh[is.nan(sh)] = 0
		select = ((which(sh == min(sh))) - nrow(init))[1]
		init = rbind(init, orig[select,])
	}
	init
}

save_inits = function(samp, sp, int.only = FALSE)
{
	suffix = "_mcmcInit"
	if(int.only) suffix = "_mcmcInit_int"

	# select a random row
	inits = samp[sample(nrow(samp),1),,drop=FALSE]


	new.inits = shan(inits, samp)
	old.inits = read.csv(paste0("dat/mcmc/", sp, suffix, ".txt"))

	for(i in 1:nrow(new.inits))
	{
		res = old.inits
		for(param in old.inits$name)
			res[res$name == param,'initialValue'] = new.inits[i,param]
		write.csv(res, paste0("dat/mcmc/", sp, suffix, i, ".txt"), row.names=FALSE)
	}
}

for(sp in spList)
{
	isamp = read.csv(file.path('res','mcmc', sp, '0i', 'posterior.csv'))
	fsamp = read.csv(file.path('res','mcmc', sp, '0', 'posterior.csv'))
	save_inits(fsamp, sp)
	save_inits(isamp, sp, TRUE)
}