library(foreach)
library(xtable)
spList = readRDS('dat/speciesList.rds')
spInfo = read.csv('dat/speciesInfo.csv')
get.info = function(sp)  spInfo[spInfo$spName == sp,]

# get the 2 datasets
trdat = foreach(spName = spList, .final=function(x) {names(x) = spList; x}) %do% {
	readRDS(file.path('dat', 'stm_calib', paste0(spName, '_stm_calib.rds')))
}
prdat = foreach(spName = spList, .final=function(x) {names(x) = spList; x}) %do% {
	readRDS(file.path('dat', 'presence', paste0(spName, '_presence.rds')))
}

trans.stat = foreach(sp = spList, .final=function(x) {names(x) = spList; x}) %do%
{
	tab = trdat[[sp]]
	type = rep('pres', nrow(tab))
	type[tab$state1 == 0 & tab$state2 == 0] = 'abs'
	type[tab$state1 == 0 & tab$state2 == 1] = 'col'
	type[tab$state1 == 1 & tab$state2 == 0] = 'ext'
	table(type)
}


spTab = data.frame(
	"Scientific name"=foreach(sp = spList, .combine=c) %do% paste(get.info(sp)$genus, get.info(sp)$species),
	"English name"=foreach(sp = spList, .combine=c) %do% as.character(get.info(sp)$comname),
#	"Occurrences"=foreach(sp = spList, .combine=c) %do% sum(prdat[[sp]][,sp]),
	"Presences"=foreach(sp = spList, .combine=c) %do% trans.stat[[sp]]['pres'],
	"Absences"=foreach(sp = spList, .combine=c) %do% trans.stat[[sp]]['abs'],
	"Colonizations"=foreach(sp = spList, .combine=c) %do% trans.stat[[sp]]['col'],
	"Extinctions"=foreach(sp = spList, .combine=c) %do% trans.stat[[sp]]['ext'],
	check.names=F
)

spXTab = xtable(spTab, label="tab:species_list", align='lllllll',
	caption="List of species studied and the number of observations in the database")
	
print(spXTab, file="res/table/si_tab_s1_species.tex", caption.placement="top", booktabs=TRUE, include.rownames=FALSE)

