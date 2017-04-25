#!/usr/bin/env Rscript
if(interactive())
	setwd("/Users/mtalluto/Documents/work/projects_active/stm/STModel-Two-State")

# all species > 0.25
spList = c('32945-FRA-NIG', '19287-QUE-MAC', '183397-TSU-CAN', '183375-PIN-RES', '22463-POP-GRA', '183412-LAR-LAR')

# additionally species > 0.2
## spList = c(spList, '19462-FAG-GRA', '183385-PIN-STR', 

# spatial coverage, environmental variable coverage, sample size
cex <- c(a= 0.1, p = 0.1, c = 1, e = 1)
pch <- c(a = '.', p = '.', c = '•', e = '•')
col = c(a = '#aaaaaa', p = '#888888', c = '#0178d8', e = '#a900a9')
cex.txt=0.7



quartz(w=24, h=9)
layout(matrix(c(1,3,5,7,9,11,2,4,6,8,10,12), nrow=2, ncol=6, byrow=TRUE))
par(mar=c(4,4,1,0), oma=c(0.5,0.5,0.5,0.5))
sapply(spList, function(spName) {
	calib <- readRDS(paste0('dat/stm_calib/', spName, '_stm_calib.rds'))
	# absences are not spatially biased, so we drop them
	calib <- calib[!(calib$state1 == 0 & calib$state2 == 0),]

	# set up unique plotting colors and symbols
	calib$type <- 'p'
	calib$cex <- cex['p']
	calib$pch <- pch['p']
	calib$col <- col['p']
	calib <- within(calib, {
		type[state1 == 0 & state2 == 1] <- 'c'
		type[state1 == 1 & state2 == 0] <- 'e'
	})
	for(tp in c('c', 'e')) {
		calib$cex[calib$type == tp] <- cex[tp]
		calib$pch[calib$type == tp] <- pch[tp]
		calib$col[calib$type == tp] <- col[tp]
	}
	
	with(calib, {
		plot(annual_mean_temp, tot_annual_pp, cex=cex, pch=pch, col=col, xlab="Mean Annual Temperature",
			ylab = "Total Annual Precipitation", main=spName)
		plot(x, y, cex=cex, pch=pch, col=col, xlab="X", ylab = "Y")
	})
	
})