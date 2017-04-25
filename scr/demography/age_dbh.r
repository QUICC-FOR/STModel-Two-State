library(reshape2)

trees <- readRDS('dat/raw/dbh_trees_20151026.rds')
trees <- rbind(trees, readRDS('dat/raw/dbh_trees_20160408.rds'))
trees <- trees[!is.na(trees$dbh) & trees$is_dead == 'f',]
trees_sm <- trees[trees$dbh < 127, ]

trees_full <- trees[trees$id_spe %in% unique(trees_sm$id_spe),c(1, 2, 4, 5, 9, 10)]

plots <- dcast(trees_full, plot_id + year_measured + x + y ~ id_spe, value.var = "dbh")
spNames <- colnames(plots)[5:ncol(plots)]
spNames <- spNames[-3]
plots[,spNames] <- as.integer(plots[,spNames] > 0)


library(parallel)
transitionData <- mclapply(spNames, function(spName) {
	stateData = plots[,c(1:4, which(colnames(plots) == spName))]

	# drop plots with one sample
	freqs = table(stateData$plot_id)
	select = which(freqs > 1)
	dfselect = which(stateData$plot_id %in% as.numeric(names(freqs[select])))
	transSamples = stateData[dfselect,]
	transSamples.wide = dcast(transSamples, plot_id ~ year_measured, value.var = spName)
	trReturn = data.frame(plot = numeric(), year1 = numeric(), year2 = numeric(),
			state1 = numeric(), state2 = numeric())

	for(i in 2:ncol(transSamples.wide)) {
		subSamp = transSamples.wide[!is.na(transSamples.wide[,i]),]
		state1 = subSamp[,i]
		n = length(state1)
		year1 = rep(as.numeric(colnames(subSamp)[i]), n)
		state2 = year2 = rep(NA, n)
		j = i + 1
		remain = n
		while(nrow(subSamp) > 0 & j <= ncol(subSamp)) {
			select = which(!is.na(subSamp[,j]))
			n = length(select)
			if(n > 0) {
				trns = data.frame(
					plot = subSamp$plot_id[select],
					year1 = rep(as.numeric(colnames(subSamp)[i]), n),
					year2 = rep(as.numeric(colnames(subSamp)[j]), n),
					interval = rep(as.numeric(colnames(subSamp)[j]) - as.numeric(colnames(subSamp)[j]), n),
					state1 = subSamp[select, i],
					state2 = subSamp[select, j]
				)
				trns = trns[!is.na(trns$year2) & !is.na(trns$state2),]
				subSamp = subSamp[-select,]
				trReturn = rbind(trReturn, trns)
			}
			j = j + 1
		}
	}
	trReturn
})
names(transitionData) <- spNames


par(mfrow=c(3,3))
for(spName in spNames) {
	calibDat <- readRDS(paste0("dat/stm_calib/", spName, "_stm_calib.rds"))
	calibDat$trans = "A"
	calibDat <- within(calibDat, {
		interval <- year2 - year1
		trans[state1 == 0 & state2 == 1] <- "C"
		trans[state1 == 1 & state2 == 1] <- "P"
		trans[state1 == 1 & state2 == 0] <- "E"
	})
	allDat <- transitionData[[spName]]
	allDat$trans = "A"
	allDat <- within(allDat, {
		interval <- year2 - year1
		trans[state1 == 0 & state2 == 1] <- "C"
		trans[state1 == 1 & state2 == 1] <- "P"
		trans[state1 == 1 & state2 == 0] <- "E"
	})
	allDat <- merge(allDat, calibDat, by=c(1,2,3))
	
	mod.o.e <- glm(as.integer(!(state2)) ~ interval + annual_mean_temp + tot_annual_pp + I(annual_mean_temp^2) + 
				I(tot_annual_pp^2), data=calibDat, family=binomial)
	mod.a.e <- glm(as.integer(!(state2.x)) ~ interval.x + annual_mean_temp + tot_annual_pp + I(annual_mean_temp^2) + 
				I(tot_annual_pp^2), data=allDat, family=binomial)
	mod.o.c <- glm(state2 ~ interval + annual_mean_temp + tot_annual_pp + I(annual_mean_temp^2) + 
				I(tot_annual_pp^2), data=calibDat, family=binomial)
	mod.a.c <- glm(state2.x ~ interval.x + annual_mean_temp + tot_annual_pp + I(annual_mean_temp^2) + 
				I(tot_annual_pp^2), data=allDat, family=binomial)
				
## 	s.o.e = summary(mod.o.e)$coefficients[,1:2]
## 	s.a.e = summary(mod.a.e)$coefficients[,1:2]
## 	yl = c(min(s.o.e[,1] - s.o.e[,2], s.a.e[,1] - s.a.e[,2]), max(s.o.e[,1] + s.o.e[,2], s.a.e[,1] + s.a.e[,2]))
## 	plot(1:6, s.o.e[,1], pch=20, col='blue', main=spName, xlim=c(1,7), ylim=yl)
## 	segments(1:6, s.o.e[,1] + 2*s.o.e[,2], 1:6, s.o.e[,1] - 2*s.o.e[,2], col='black')
## 	points(1:6+0.1, s.a.e[,1], col='red', pch=20)
## 	segments(1:6+0.1, s.a.e[,1] + 2*s.a.e[,2], 1:6+0.1, s.a.e[,1] - 2*s.a.e[,2], col='black')

	xx <- data.frame(annual_mean_temp = seq(-3,3, length.out=500), tot_annual_pp = 0, interval = 5, interval.x = 5)
	yoe <- predict(mod.o.e, newdata=xx, type='response')
	yoc <- predict(mod.o.c, newdata=xx, type='response')
	yae <- predict(mod.a.e, newdata=xx, type='response')
	yac <- predict(mod.a.c, newdata=xx, type='response')
	plot(xx[,1], yoe, ylim=c(0,1), col='red', type='l')
	lines(xx[,1], yae, col='red', lty=2)
	lines(xx[,1], yoc, col='blue')
	lines(xx[,1], yac, col='blue', lty=2)

}







