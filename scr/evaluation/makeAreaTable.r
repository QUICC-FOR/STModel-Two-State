#!/usr/bin/env Rscript

# depends:
#    res/eval/SPNAME_MOD__rdeAreas.rds

# produces:
#    res/si_tab_range_areas.tex

speciesList = readRDS('dat/speciesList.rds')
speciesInfo = read.csv('dat/speciesInfo.csv', stringsAsFactors=FALSE)

modName = '0'
spSkip = '[0.3em]'
tabLns = c(
	'\\begin{table}[tb]',
	'\\label{tab:range_size}',
	'\\caption{Median range sizes (in 1000s of km\\textsuperscript{2}).',
	'Parenthetical values are 95\\% credible intervals.}',
	'\\begin{tabular}{lccc}',
	'Species & Present Distribution & Equilibrium Distribution & Change in Range Size \\\\',
	'\\toprule'
)

areaFactor = 1000
digits = 2
for(spName in speciesList)
{
	info = speciesInfo[speciesInfo$spName == spName,]
	areas = readRDS(file.path('res','eval', paste0(spName, '_', modName, '_rdeAreas.rds')))
	curRange = round(mean(areas[,1] + areas[,3], na.rm=TRUE) / areaFactor, digits)
	eqRange = round(median(areas[,1] + areas[,2], na.rm=TRUE) / areaFactor, digits)
	eqRange.q = round(quantile(areas[,1] + areas[,2], c(0.025, 0.975), na.rm=TRUE) / areaFactor, digits)
	delta = round(median(areas[,2] - areas[,3], na.rm=TRUE) / areaFactor, digits)
	delta.q = round(quantile(areas[,2] - areas[,3], c(0.025, 0.975), na.rm=TRUE) / areaFactor, digits)


	tabLns = c(tabLns, 
		paste("\\multirow{2}{*}{\\it", info$genus, info$species, '} & \\multirow{2}{*}{', curRange, '} &', eqRange, '&', delta, '\\\\'),
		paste0(' & & (', eqRange.q[1], ', ', eqRange.q[2], ') & (', delta.q[1], ', ', delta.q[2], ') \\\\', spSkip)
	)
}

tabLns = c(tabLns, '\\bottomrule', '\\end{tabular}', '\\end{table}')
write(tabLns, file="res/si_tab_range_areas.tex")
