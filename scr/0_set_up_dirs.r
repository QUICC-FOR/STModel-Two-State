#!/usr/bin/Rscript
# this is a quick utility script that will build the full working tree of directories
# and set up a few other things that are used globally. Should be run once on each
# computer that will be running any analyses

# create the directory tree needed by the project, if it doesn't already exist
spList = c(
# temperate species that occur in the transition, ordered by occurrence
# quealb is borderline
'28728-ACE-RUB', '28731-ACE-SAC', '19290-QUE-ALB', '19408-QUE-RUB', '19049-ULM-AME',
# '32931-FRA-AME', '19462-FAG-GRA', '183385-PIN-STR', '32929-FRA-PEN', '183397-TSU-CAN',
# '505490-THU-OCC', '18034-PIC-RUB', '22463-POP-GRA', '32945-FRA-NIG',

# boreal species
'18032-ABI-BAL', '19489-BET-PAP', '195773-POP-TRE', '19481-BET-ALL', '183302-PIC-MAR', 
# '183295-PIC-GLA',

# temperate and southern species
# '18037-PIN-TAE', '19027-LIQ-STY', '18086-LIR-TUL', '19447-QUE-VEL', '27821-NYS-SYL',
# '19280-QUE-NIG', 'NA-CAR-ALB',  '19422-QUE-STE', '19231-CAR-GLA', '19277-QUE-FAL',
# '18048-JUN-VIR', '183335-PIN-ECH', '19051-ULM-ALA', 'NA-QUE-PRI', '19242-CAR-OVA'
# '19288-QUE-COC', '19254-JUG-NIG', '19050-ULM-RUB', '23690-OXY-ARB', '19511-OST-VIR'
)

subDirs = c('dat', 'img', 'res')
resDirs = c('anneal', 'mcmc1', 'mcmc2', 'mcmc3')
for(sp in spList)
{
	for(subD in subDirs)
		dir.create(file.path('species', sp, subD), recursive=TRUE)
	for(subD in resDirs)
		dir.create(file.path('species', sp, 'res', subD), recursive=TRUE)
}

# directories for global data and results
for(sDir in subDirs)
	dir.create(file.path(sDir), recursive=TRUE)
	
saveRDS(spList, file.path('dat', 'speciesList.rds'))
