#!/usr/bin/env Rscript

# depends:
#    res/eval/SPNAME_MOD__rdeAreas.rds

# produces:
#    res/si_tab_range_areas.tex
library(foreach)
library(xtable)

library(raster)
speciesList = readRDS('dat/speciesList.rds')
spInfo = read.csv('dat/speciesInfo.csv', stringsAsFactors=FALSE)
mapProj = readRDS("dat/map_projections.rds")
get.info = function(sp)  spInfo[spInfo$spName == sp,]

# get areas
digits = 2
uncertainRange = c(0.25, 0.75)
areas = sapply(speciesList, simplify=FALSE, USE.NAMES=TRUE, FUN=function(spName) {
	a = readRDS(file.path('res','areas',paste0(spName, '_0_areas.rds')))
	spGrid = readRDS(file.path('res','rangemaps',paste0(spName,'_rangemaps.rds')))
	spUncertainty = rasterFromXYZ(data.frame(x=spGrid$x, y=spGrid$y, 
			u=(spGrid$stm > uncertainRange[1] & spGrid$stm < uncertainRange[2])))
	proj4string(spUncertainty) = mapProj$projected
	spRange = rasterFromXYZ(data.frame(x=spGrid$x, y=spGrid$y, r=spGrid$stm.pres))
	u.area = sum(values(spUncertainty), na.rm=T) * prod(res(spUncertainty)/1000)/1000
	range.area = sum(values(spRange), na.rm=T) * prod(res(spRange)/1000)/1000
	proj4string(spRange) = mapProj$projected
	c(curRange = round(mean(a[,1] + a[,3], na.rm=TRUE), digits),
		eqRange = round(median(a[,1] + a[,2], na.rm=TRUE), digits),
		eqRange.q = round(quantile(a[,1] + a[,2], c(0.025, 0.975), na.rm=TRUE), digits),
		delta = round(median(a[,2] - a[,3], na.rm=TRUE), digits),
		delta.q = round(quantile(a[,2] - a[,3], c(0.025, 0.975), na.rm=TRUE), digits),
		uncertain.area = u.area,
		uncertain.percent = round(100 * u.area / range.area, digits)
	)
})

a2 = do.call(rbind, areas)
medianUncertain = median(a2[,'uncertain.percent'])
print(medianUncertain)

makerows = function(x, sp)
{	data.frame(
		"Species"=c(paste(get.info(sp)$genus, get.info(sp)$species),""),
		"Present area"=c(paste(x['curRange']), ""),
		"Equilibrium area" = c(paste(x['eqRange']), paste0("(",x['eqRange.q.2.5%'], ",", x['eqRange.q.97.5%'],")")),
		"Change in area" = c(paste(x['delta']), paste0("(",x['delta.q.2.5%'], ",", x['delta.q.97.5%'],")")),
		check.names=FALSE)
}

spTab = foreach(sp=speciesList, .combine=rbind) %do% makerows(areas[[sp]], sp)
spXTab = xtable(spTab, label="tab:dic", align='llccc',
	caption="Table S3: Median range sizes (in 1000s of km2). Parenthetical values are 90\\% credible intervals.")
	
print(spXTab, file="res/table/si_tab_s3_area.tex", caption.placement="top", booktabs=TRUE, include.rownames=FALSE)
