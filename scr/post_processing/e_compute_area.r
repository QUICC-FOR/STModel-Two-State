library(doParallel)
library(foreach)
source('scr/stm_functions.r')
load("dat/map_projections.rdata")
suppressWarnings(dir.create(file.path('res', 'areas'), recursive=TRUE))
speciesList = readRDS('dat/speciesList.rds')
models = c('0', 'i0', 'g', 'ig')

registerDoParallel(cores=7)

for(spName in speciesList)
{
	for(mod in models)
	{
		spGrid = readRDS(file.path('res','maps',paste0(spName,'_',mod,'_maps.rds')))
		postGrid = readRDS(file.path('res','posteriorGrid',paste0(spName,'_', mod, '_posteriorGrid.rds')))
		grPres = postGrid$stmPres

		areas = foreach(pres = iter(grPres, by='row'), .combine=rbind) %dopar% {
			rde = (1 * (pres & spGrid$sdm.pres)) + 
						(2 * (pres & !spGrid$sdm.pres)) + 
						(3 * (!pres & spGrid$sdm.pres))
			rde[rde == 0] = NA
			rde = rde - 1
			rdeRas = make_raster(rde, spGrid[,1:2], P4S.latlon, stmMapProjection)
			freq(rdeRas)[1:3,2] * prod(res(rdeRas)/1000)/1000
		}
		saveRDS(areas, file.path('res','areas',paste0(spName, '_', mod, '_areas.rds')))
	}	
}