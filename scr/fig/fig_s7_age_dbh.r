#!/usr/bin/env Rscript
library(xtable)

dat <- read.csv("dat/raw/age_dbh_NEE_reviewer.csv", sep=';')
drops <- which(dat$species == 'Quercus macrocarpa' | dat$species == 'Fraxinus americana' | dat$dbh < 20)
dat <- dat[-drops,]
dat$species <- factor(dat$species)

# put the species in order
spList <- c("Acer rubrum", "Acer saccharum", "Fagus grandifolia", "Fraxinus nigra", "Quercus rubra", "Tsuga canadensis",
			"Pinus resinosa", "Pinus strobus", "Abies balsamea", "Betula papyrifera", "Larix laricina", "Picea glauca",
			"Picea mariana","Picea rubens", "Pinus banksiana", "Thuja occidentalis")

yax = rep(c(TRUE, rep(FALSE, 3)), 4)
xax <- c(rep(FALSE, 12), rep(TRUE, 4))

## to do: change to PDF, save, and move on with my life
figure.width = 6.5
figure.height= 6.5
filename = file.path('img', 'figs', 'fig_s1.png')

pdf(w=figure.width, h=figure.height, file="img/figs/fig_siX_age-dbh.pdf")
par(mfrow=c(4,4), mar=c(0.2, 1, 2, 0), oma=c(3.5,3.5,0,1), cex.axis=0.9, cex.main=0.9, bty='l', tcl=-0.1, mgp=c(2,0.5,0), las=1)
fits <- t(mapply(function(spName, xa, ya) {
	dd <- dat[dat$species == spName,]
	mod <- lm(age ~ I(dbh/10), dat=dd)
	plLab = bquote(italic(.(spName)))
	plot(dd$dbh/10, dd$age, pch=20, cex=0.3, main=plLab, ylim=c(0,300), xlim=c(0, 70), xaxt='n', yaxt='n')
	lines(range(dd$dbh)/10, predict(mod, newdata=data.frame(dbh=range(dd$dbh))), col='blue', lwd=1.5)
	if(xa) axis(side=1, at=seq(0, 60, 20))
	if(ya) axis(side=2, at=seq(0,300,100))
	fit <- predict(mod, newdata=data.frame(dbh=127), se.fit=TRUE, interval='predict')
	val <- c(fit$fit, fit$se.fit)
	names(val) <- c('fit', 'lower', 'upper', 'se')
	val
}, spList, xax, yax))
mtext("DBH (cm)", side=1, outer=TRUE, cex=0.8, line=2)
mtext("Age (years)", side=2, outer=TRUE, cex=0.8, line=2, las=0)
dev.off()
rownames(fits) <- spList

dbh_xt <- xtable(fits[,c(1,4)], digits=2, align='lcc', caption="Predicted mean age (in years) with standard error for all species at 12.7 cm diameter at breast height (DBH), the minimum tree size considered in our analyses")
print(dbh_xt, file="res/table/si_tab_sX_dbh.tex", caption.placement="top", booktabs=TRUE, include.rownames=TRUE)
