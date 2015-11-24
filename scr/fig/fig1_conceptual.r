#!/usr/bin/env Rscript
col.color = '#377eb8'
ext.color = '#e41a1c'
cex.axis=0.75
lwd=1.5
lwd.lims = 0.25
lab.line=1.5
cex.annote=0.7
ytop.adj=-0.03
ylett.adj=-0.01
cex.id=1
fut.xadj = 1.3
fut.yadj = -0.005
arr.yadj = 0.11
arr.xsadj = -0.3
arr.xeadj = -0.1
arrtxt.xadj = -0.7
arrtxt.yadj = -0.077

rc = readRDS("res/resp_curve/28731-ACE-SAC_0_respCurve.rds")
clim.map = readRDS('dat/clim/climateGrid_unscaled_unprojected.rds')

# fit a transformation from temperature to latitude. turns out it's pretty linear
library(gam)
latmod = gam(lat ~ s(annual_mean_temp,3), data=clim.map)
warming = 2.5
T = rc$temp
lat = predict(latmod, newdata=list(annual_mean_temp = T))
lat.warm = predict(latmod, newdata=list(annual_mean_temp = T-warming))

C = rc$col.temp
C.upp = rc$col.temp.upper
C.low = rc$col.temp.lower
E = rc$ext.temp
E.upp = rc$ext.temp.upper
E.low = rc$ext.temp.lower
xlims = c(35,51)
ylims = c(0,0.28)

presentNorthernTempLimit = 1.216216
presentSouthernTempLimit = 12.56756757
rts = predict(latmod, newdata=list(annual_mean_temp=c(presentNorthernTempLimit, presentSouthernTempLimit)))
label.col = function() mtext(expression(bold(Colonization)), side=2, line=lab.line, col=col.color)
label.ext = function() mtext(expression(bold(Extinction)), side=4, line=lab.line, col=ext.color)

figure.width = 12/2.54
figure.height = 6.5/2.54
filename = file.path('img', 'figs', 'fig1.pdf')
pdf(width=figure.width, height=figure.height, file=filename, pointsize = 9)

layout(matrix(c(1,2), nrow=2))
par(mgp=c(1,.3,0), tcl=-0.25, mar=c(0,1,0.5,0), oma=c(3,2,0,0.5), cex.axis=cex.axis)
plot(0,0, type='n', ylim=ylims, ylab="", bty='l', xlim=xlims, xlab="", xaxt='n', yaxt='n')
axis(2, cex.axis=cex.axis, at=c(0,0.1,0.2))
## axis(4, cex.axis=cex.axis)
polygon(rep(rts, each=2), c(ylims, rev(ylims)), col="#77777744", border=NA)
lines(lat, C, col=col.color, lwd=lwd)
lines(lat, C.upp, col=col.color, lwd=lwd.lims)
lines(lat, C.low, col=col.color, lwd=lwd.lims)
lines(lat, E, col=ext.color, lwd=lwd)
lines(lat, E.upp, col=ext.color, lwd=lwd.lims)
lines(lat, E.low, col=ext.color, lwd=lwd.lims)
## label.col()
## label.ext()
text(mean(rts), ylims[2]+ytop.adj, "Present range\n(at equilibrium)", cex=cex.annote)
text(min(xlims), ylims[2]+ylett.adj, "A.", cex=cex.id)

plot(0,0, type='n', ylim=ylims, ylab="", bty='l', xlim=xlims, xlab="", yaxt='n')
axis(2, cex.axis=cex.axis, at=c(0,0.1,0.2))
## axis(4, cex.axis=cex.axis)
ext.x = rep(c(rts[2], rts[2] + warming), each=2)
pres.x = rep(c(rts[2] + warming, rts[1]), each=2)
# bit of manual adjustment here to make the zones line up better
col.x = rep(c(rts[1], rts[1]+warming-0.4), each=2)
polygon(ext.x, c(ylims, rev(ylims)), col=paste0(ext.color, "22"), border=NA)		
polygon(pres.x, c(ylims, rev(ylims)), col="#55555544", border="#55555544")
polygon(col.x, c(ylims, rev(ylims)), col=paste0(col.color, "22"), border=NA)
lines(lat.warm, C, col=col.color, lwd=lwd)
lines(lat.warm, C.upp, col=col.color, lwd=lwd.lims)
lines(lat.warm, C.low, col=col.color, lwd=lwd.lims)
lines(lat.warm, E, col=ext.color, lwd=lwd)
lines(lat.warm, E.upp, col=ext.color, lwd=lwd.lims)
lines(lat.warm, E.low, col=ext.color, lwd=lwd.lims)
mtext("Latitude", side=1, cex=0.9, line=1.5)
## label.col()
## label.ext()
text(mean(pres.x), ylims[2]+ytop.adj, "Continued\npresence", cex=cex.annote)
text(mean(ext.x)+fut.xadj, ylims[2]+ytop.adj+fut.yadj, "Future\nextinction", pos=2, cex=cex.annote)
text(mean(col.x)-fut.xadj, ylims[2]+ytop.adj+fut.yadj, "Future\ncolonization", pos=4, cex=cex.annote)
text(min(xlims)+arrtxt.xadj, mean(ylims)+arrtxt.yadj, "Warming moves the\nrange poleward", pos=4, cex=cex.annote)
arrows(min(xlims) + arr.xsadj, mean(ylims)-arr.yadj, rts[2]+arr.xeadj, mean(ylims)-arr.yadj , lwd=1.5, length=0.07)
text(min(xlims), ylims[2]+ylett.adj, "B.", cex=cex.id)

mtext("Colonization/Extinction Rate", side=2, outer=TRUE, line=0.5)

dev.off()
