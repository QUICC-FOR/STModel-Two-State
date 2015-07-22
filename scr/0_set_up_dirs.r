#!/usr/bin/Rscript
# this is a quick utility script that will build the full working tree of directories
# and set up a few other things that are used globally. Should be run once on each
# computer that will be running any analyses

# create the directory tree needed by the project, if it doesn't already exist
spList = c('18032-ABI-BAL', '18037-PIN-TAE', '18086-LIR-TUL', '19049-ULM-AME', 
		'19280-QUE-NIG', '19290-QUE-ALB', '19408-QUE-RUB', '19447-QUE-VEL', 
		'19462-FAG-GRA', '19481-BET-ALL', '19489-BET-PAP', '27821-NYS-SYL', 
		'28728-ACE-RUB', '28731-ACE-SAC', '32929-FRA-PEN', '32931-FRA-AME', 
		'183295-PIC-GLA', '183302-PIC-MAR', '183319-PIN-BAN', '183385-PIN-STR', 
		'195773-POP-TRE')
subDirs = c('dat', 'img', 'res')
for(sp in spList)
{
	for(subD in subDirs)
	{
		dir.create(file.path('species', sp, subD), recursive=TRUE)
	}
	dir.create(file.path('species', sp, "res", "anneal"), recursive=TRUE)
}

# directories for global data and results
for(sDir in subDirs)
	dir.create(file.path(sDir), recursive=TRUE)
	
saveRDS(spList, file.path('dat', 'speciesList.rds'))
