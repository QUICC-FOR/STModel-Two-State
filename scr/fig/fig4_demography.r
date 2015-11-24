#!/usr/bin/env Rscript
library(rstan)
speciesInfo = read.csv('dat/speciesInfo.csv')
## deadPath = file.path('res', 'demography', 'deadStanMod.rds')
mortalityPath = file.path('res', 'demography', 'mortality_stan.rds')
recruitmentPath = file.path('res', 'demography', 'recruitment_stan.rds')

## dead = readRDS(deadPath)
mort = readRDS(mortalityPath)
recruit = readRDS(recruitmentPath)

borealSp = c('18032-ABI-BAL', '19489-BET-PAP', '183412-LAR-LAR', '183295-PIC-GLA', 
			'183302-PIC-MAR', '18034-PIC-RUB', '183319-PIN-BAN', '195773-POP-TRE', 
			'505490-THU-OCC')
tempSp = c('28728-ACE-RUB', '28731-ACE-SAC', '19462-FAG-GRA', '32931-FRA-AME',
			'32945-FRA-NIG', '19287-QUE-MAC', '19408-QUE-RUB', '183397-TSU-CAN')
transitionalSp = c('19481-BET-ALL', '183375-PIN-RES', '183385-PIN-STR', '22463-POP-GRA')

make_points = function(dat)
{
	samples = extract(dat$results)
	predictions = with(samples, list(inrange=plogis(alpha), outrange=plogis(alpha+beta)))
	res = do.call(rbind, lapply(1:ncol(predictions[['inrange']]), function(i) {
		data.frame(
			species=dat$species$name[dat$species$value == i],
			ir.mean = mean(predictions[['inrange']][,i]),
			ir.upper = quantile(predictions[['inrange']][,i], 0.975),
			ir.lower = quantile(predictions[['inrange']][,i], 0.025),
			or.mean = mean(predictions[['outrange']][,i]),
			or.upper = quantile(predictions[['outrange']][,i], 0.975),
			or.lower = quantile(predictions[['outrange']][,i], 0.025)
	)}))
	# add a null row to give some padding to the plot
	rownames(res) = res$species
	res
}


mort.points = make_points(mort)
mort.points = mort.points[c(borealSp, transitionalSp, tempSp),]
recruit.points = make_points(recruit)
recruit.points = recruit.points[c(borealSp, transitionalSp, tempSp),]


# plot settings
space = c(0,1)
cols=c('#a6cee3', '#fb9a99')
xlim=c(-0.15, 0.25)
eby.ir = seq(space[2] + 0.5, length.out = nrow(mort.points), by = 2 + space[2])
eby.or = eby.ir+1
cex.axis = 0.7
cex.annotate = 0.7
mtext.line= 0
mtext.at = c(-0.075, 0.075)
cex.splab=0.5
cex.type = 0.65
cex.axtitle=0.75
cex.legend=0.5
lwd.segment=0.75
lwd.axis=0.75


dpi = 600
figure.width = 5.5
figure.height = 6.5
filename = file.path('img', 'figs', paste0('fig4.pdf'))
fontsize=9
pdf(width=figure.width/2.54, height=figure.height/2.54, pointsize=fontsize, file=filename)
## png(width=as.integer(dpi*figure.width), height=as.integer(dpi*figure.height),
## 	file=filename, pointsize=fontsize, res=dpi)


par(mar=c(3,3.5,1,0.5), mgp=c(1.3,0.4,0), tcl=-0.2)



# draw mortality bars
with(mort.points, {
	barplot(t(cbind(ir.mean, or.mean)), beside=TRUE, horiz=TRUE, xlim=xlim, space=space, border=NA, col=cols, xaxt='n')
	segments(x0=ir.lower, y0=eby.ir, x1=ir.upper, y1=eby.ir, lwd=lwd.segment)
	segments(x0=or.lower, y0=eby.or, x1=or.upper, y1=eby.or, lwd=lwd.segment)	
})



# draw recruitment bars
par(new=TRUE)
with(recruit.points, {
	barplot(t(-1*cbind(ir.mean, or.mean)), beside=TRUE, horiz=TRUE, xlim=xlim, space=space, border=NA, col=cols, xaxt='n')
	segments(x0=-1*ir.lower, y0=eby.ir, x1=-1*ir.upper, y1=eby.ir, lwd=lwd.segment)
	segments(x0=-1*or.lower, y0=eby.or, x1=-1*or.upper, y1=eby.or, lwd=lwd.segment)
})


# axes are nice
axis(side=2, pos=0, at=c(min(eby.ir)-2, max(eby.or)+2), tcl=0, mgp=c(3,1,0), labels=rep("", 2), lwd=lwd.axis)
ax.at = pretty(c(-0.1, 0.2),4)
ax.labs=as.character(round(ax.at, 2))
axis(side=1, at=ax.at, labels=ax.labs, cex.axis=cex.axis, lwd=lwd.axis)




# draw rectangles around species groups
## rectcol="#777777"
## xpad = 0.01
## ypad = c(1, 2)
## numsp = length(eby.or)
## te.ind = c(numsp - length(tempSp) + 1, numsp)
## tr.ind = c(te.ind[1] - length(transitionalSp), te.ind[1] - 1)
## b.ind = c(1, tr.ind[1] - 1)
## ylims.boreal = eby.ir[b.ind] + c(-1, 1) *ypad
## ylims.trans = eby.ir[tr.ind] + c(-1, 1) *ypad
## ylims.temp = eby.ir[te.ind] + c(-1, 1) *ypad
## xlims.rect = xlim
## xlims.rect = c(-max(recruit.points$or.upper)-xpad, max(mort.points$or.upper)+xpad)
## polygon(rep(xlims.rect, each=2), c(ylims.boreal, rev(ylims.boreal)), border=rectcol)
## polygon(rep(xlims.rect, each=2), c(ylims.trans, rev(ylims.trans)), border=rectcol)
## polygon(rep(xlims.rect, each=2), c(ylims.temp, rev(ylims.temp)), border=rectcol)

## legend(0.11, 12, legend=c("Equilibrium Range", "Extinction Debt"), fill=cols, cex=cex.legend, box.lwd=0.65, border="#777777")


# add text annotations
mtext(c("Recruitment", "Mortality"), side=3, line=mtext.line, at=mtext.at, cex=cex.annotate)

par(xpd=TRUE)
xtext = -0.28
spNames.y = rowMeans(cbind(eby.or, eby.ir))
for(i in 1:nrow(recruit.points))
{
	sp = as.character(recruit.points$species[i])
	info = speciesInfo[speciesInfo$spName == sp,]
	plLab = bquote(italic(.(paste0(substr(as.character(info$genus),1,1),'.'))~.(as.character(info$species))))
	text(xtext, spNames.y[i], plLab, cex=cex.splab, pos=4)
}

## text(0.26, max(spNames.y[b.ind]), "Boreal", cex=cex.type, pos=2)
## text(0.26, max(spNames.y[tr.ind]), "Transitional", cex=cex.type, pos=2)
## text(0.26, max(spNames.y[te.ind]), "Temperate", cex=cex.type, pos=2)
text(0,grconvertY(0.055, from="ndc", to="user"), bquote(Rate~"("*individual^-1~year^-1*")"), adj=c(0.5, NA), cex=cex.axtitle)

dev.off()
